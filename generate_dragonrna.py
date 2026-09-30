#!/usr/bin/env python3
"""
Generate a DragonRNA with a designed 3' self-priming foldback and no other
self-structure.

Layout (5' -> 3'), each block's length set by its own flag:

    [ prefix ][ randomer (N) ][ spacer ][ arm ][ loop ][ arm' ]
        20           15            20      6      4       6     = 71 nt

arm' is the Watson-Crick reverse complement of arm, so the 3' terminus folds
back and pairs with arm.  The 3' terminal base pairs with the 5'-most base of
arm, so extension from the 3' end reads straight back through the spacer, the
randomer and the prefix.

--prefix and --spacer take either a length (random sequence) or a literal
sequence to hold fixed.  --randomer takes a length; those positions are the
degenerate (N) block of a library.

A candidate is only reported if all of these hold:

  1. Alphabet-level: no (--max-internal-foldback + 1)-mer anywhere in the
     molecule can pair with any non-overlapping window elsewhere, the designed
     stem excepted.  This bounds every contiguous duplex that could ever form,
     not just the one the energy model happens to pick.  G-U wobble counts as
     a pair unless --no-wobble; ignoring it silently permits 4+ bp helices
     that ViennaRNA does form.
  2. The ViennaRNA MFE structure of the whole molecule is exactly the designed
     stem: those base pairs and no others.
  3. The molecule minus its 3' arm folds to nothing (all dots).

Reported but deliberately not enforced: --randomer-trials resamples the
randomer block and reports the worst foldback found, because only the emitted
instance can be constrained -- the rest of the library cannot.  With a full
ACGU randomer most library members do fold beyond the designed stem; a
narrower --randomer-alphabet such as AG or AC trades library size for far
less structure.

Folding always uses RNA energy parameters; --dna (the default) only reports
T in place of U.

Usage:
    python3 generate_dragonrna.py
    python3 generate_dragonrna.py --randomer-alphabet AG
    python3 generate_dragonrna.py --prefix ACAACAACCAACCAAACCAA --spacer 20
    python3 generate_dragonrna.py --prefix 40 --randomer 0 --spacer 0 --stem-len 0
"""
import argparse
import itertools
import random
import sys
import time

import RNA

BASES = "ACGU"
WATSON_CRICK = {"A": "U", "C": "G", "G": "C", "U": "A"}
WITH_WOBBLE = {"A": "U", "C": "G", "G": "CU", "U": "AG"}


def reverse_complement(seq):
    return "".join(WATSON_CRICK[c] for c in reversed(seq))


def partner_kmers(kmer, partners):
    """Every window that could pair along the whole length of `kmer`.

    Pairing window A against window B lines a_t up with b_{k-1-t}, so B is
    drawn from the partners of `kmer` read backwards.
    """
    return ("".join(p) for p in
            itertools.product(*(partners[c] for c in reversed(kmer))))


def mfe_pairs(structure):
    """0-based (i, j) base pairs read off a dot-bracket string."""
    stack, pairs = [], set()
    for i, c in enumerate(structure):
        if c == "(":
            stack.append(i)
        elif c == ")":
            pairs.add((stack.pop(), i))
    return pairs


def longest_foldback(seq, ignore, partners):
    """Longest contiguous duplex the sequence admits, ignoring given (i, j) pairs.

    Independent of the energy model and of how the sequence was built, so it
    re-derives check 1 rather than trusting the generator.
    """
    n = len(seq)
    best = 0
    for i in range(n):
        for j in range(i + 1, n):
            run = 0
            while (i + run < j - run
                   and seq[j - run] in partners[seq[i + run]]
                   and (i + run, j - run) not in ignore):
                run += 1
            best = max(best, run)
    return best


def build(total, forced, alphabets, arm1_start, arm2_start, stem_len, k,
          max_homopolymer, max_kmer_reuse, partners, rng, budget):
    """Randomised DFS over positions, forcing arm' to the reverse complement of arm.

    A base is rejected the moment it would complete a k-mer able to pair with a
    non-overlapping window already placed, so check 1 holds by construction.
    Returns None if the step budget runs out or the space is exhausted.
    """
    seq = []
    kmer_starts = {}
    kmer_chosen = {}   # only windows this function completed itself
    steps = [0]

    def in_arm1(s):
        return bool(stem_len) and arm1_start <= s and s + k <= arm1_start + stem_len

    def in_arm2(s):
        return bool(stem_len) and arm2_start <= s and s + k <= total

    def exempt(s, j):
        """Is this pairable k-mer pair harmless?

        Either it *is* the designed stem, or the windows overlap: pairing
        seq[j:j+k] against seq[s:s+k] puts (j, s+k-1) outermost and (j+k-1, s)
        innermost, so a full k-bp duplex needs them disjoint.  Overlapping
        windows top out below k and must not be rejected; longest_foldback
        applies the same geometry.
        """
        if abs(s - j) < k:
            return True
        return (in_arm1(s) and in_arm2(j)) or (in_arm2(s) and in_arm1(j))

    def fixed(p):
        if stem_len and p >= arm2_start:
            return WATSON_CRICK[seq[arm1_start + stem_len - 1 - (p - arm2_start)]]
        return forced[p]

    def place(p):
        if p == total:
            return True
        steps[0] += 1
        if steps[0] > budget:
            return False
        f = fixed(p)
        if f:
            options = [f]
        else:
            pool = alphabets[p]
            options = rng.sample(pool, len(pool))
            if p >= k - 1:
                # Reusing a k-mer already present costs nothing, while a fresh one
                # rules out its partners for the rest of the build.  Trying reuses
                # first is what makes long molecules tractable.
                head = "".join(seq[p - k + 1:])
                options.sort(key=lambda b: head + b not in kmer_starts)
        for b in options:
            if (max_homopolymer and len(seq) >= max_homopolymer
                    and all(c == b for c in seq[-max_homopolymer:])):
                continue
            seq.append(b)
            s = p - k + 1
            added = None
            ok = True
            if s >= 0:
                kmer = "".join(seq[s:])
                if (not f and max_kmer_reuse
                        and kmer_chosen.get(kmer, 0) >= max_kmer_reuse):
                    # Without this the reuse-first ordering above happily emits a
                    # tandem repeat of one short k-mer cycle.  It governs only bases
                    # this function chose; a supplied literal block is the caller's
                    # business and must not eat the allowance.
                    ok = False
                if ok:
                    for mate in partner_kmers(kmer, partners):
                        for j in kmer_starts.get(mate, ()):
                            if not exempt(s, j):
                                ok = False
                                break
                        if not ok:
                            break
                if ok:
                    kmer_starts.setdefault(kmer, []).append(s)
                    if not f:
                        kmer_chosen[kmer] = kmer_chosen.get(kmer, 0) + 1
                    added = kmer
            if ok and place(p + 1):
                return True
            if added is not None:
                kmer_starts[added].pop()
                if not kmer_starts[added]:
                    del kmer_starts[added]
                if not f:
                    kmer_chosen[added] -= 1
            seq.pop()
        return False

    return "".join(seq) if place(0) else None


def survey_randomer(seq, start, rand_len, alphabet, want, partners, trials, rng):
    """What the rest of the library does if only the randomer varies."""
    worst_foldback = 0
    with_extra_pairs = 0
    for _ in range(trials):
        draw = "".join(rng.choices(alphabet, k=rand_len))
        cand = seq[:start] + draw + seq[start + rand_len:]
        worst_foldback = max(worst_foldback,
                             longest_foldback(cand, want, partners))
        struct, _ = RNA.fold(cand)
        if mfe_pairs(struct) != want:
            with_extra_pairs += 1
    return worst_foldback, with_extra_pairs


def preflight_block(name, seq, limit, partners):
    """Why a supplied literal block can never satisfy the checks, if it cannot."""
    problems = []
    worst = longest_foldback(seq, set(), partners)
    if worst > limit:
        problems.append(f"its longest internal foldback is {worst} nt, over the "
                        f"{limit} nt limit")
    struct, mfe = RNA.fold(seq)
    if struct.count(".") != len(struct):
        problems.append(f"it folds on its own ({struct}, {mfe:.2f} kcal/mol), so the "
                        f"designed stem could not be the only structure")
    return [f"--{name} {seq.replace('U', 'T')}: {p}" for p in problems]


def parse_block(value, name, ap):
    """A block flag is either a length or a literal sequence to hold fixed."""
    if value.isdigit():
        return int(value), None
    seq = value.strip().upper().replace("T", "U")
    bad = sorted(set(seq) - set(BASES))
    if not seq or bad:
        ap.error(f"--{name}: not a length or an ACGT/ACGU sequence "
                 f"(offending: {''.join(bad) or 'empty'})")
    return len(seq), seq


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prefix", default="20",
                    help="5' block: length, or a literal sequence (default: 20)")
    ap.add_argument("--randomer", type=int, default=12,
                    help="length of the degenerate (N) block (default: 12)")
    ap.add_argument("--spacer", default="15",
                    help="block between randomer and stem: length, or a literal "
                         "sequence (default: 15)")
    ap.add_argument("--stem-len", type=int, default=6,
                    help="designed 3' self-priming foldback, in bp; 0 for none "
                         "(default: 6)")
    ap.add_argument("--loop-len", type=int, default=4,
                    help="hairpin loop between the two stem arms (default: 4)")
    ap.add_argument("--max-internal-foldback", type=int, default=2,
                    help="longest contiguous duplex allowed anywhere outside the "
                         "designed stem (default: 2)")
    ap.add_argument("--wobble", action=argparse.BooleanOptionalAction, default=True,
                    help="count G-U as a pair when bounding foldbacks, as ViennaRNA "
                         "does (default: on)")
    ap.add_argument("--stem-min-gc", type=int, default=3,
                    help="minimum G+C in the stem; only a cheap prefilter, since check 2 "
                         "is what decides whether the stem really folds (default: 3)")
    ap.add_argument("--block-alphabet", default="ACGU",
                    help="bases the generated prefix/spacer may use. No 3-letter code "
                         "is pair-free, but AC, AG and CU are: with G-U counted, a "
                         "2 nt foldback limit over ~70 nt effectively needs one "
                         "(default: ACGU)")
    ap.add_argument("--randomer-alphabet", default="ACGU",
                    help="bases the degenerate block may use; a narrower code such as "
                         "AG or AC trades library size for far less foldback "
                         "(default: ACGU)")
    ap.add_argument("--max-kmer-reuse", type=int, default=3,
                    help="times one k-mer may repeat, which is what keeps the sequence "
                         "from collapsing into a tandem repeat; 0 for no limit "
                         "(default: 2)")
    ap.add_argument("--max-homopolymer", type=int, default=3,
                    help="longest run of one base allowed; 0 for no limit (default: 3)")
    ap.add_argument("--dna", action=argparse.BooleanOptionalAction, default=True,
                    help="report T in place of U (default: on)")
    ap.add_argument("--randomer-trials", type=int, default=200,
                    help="resample the randomer this many times and report the worst "
                         "foldback; 0 to skip (default: 200)")
    ap.add_argument("--max-attempts", type=int, default=20_000)
    ap.add_argument("--timeout", type=float, default=60.0,
                    help="give up after this many seconds (default: 60)")
    ap.add_argument("--budget", type=int, default=5_000,
                    help="DFS steps before restarting an attempt; restarting cheaply "
                         "beats searching a stuck prefix harder (default: 5000)")
    ap.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    args = ap.parse_args()

    partners = WITH_WOBBLE if args.wobble else WATSON_CRICK
    prefix_len, prefix_seq = parse_block(args.prefix, "prefix", ap)
    spacer_len, spacer_seq = parse_block(args.spacer, "spacer", ap)
    stem_len, loop_len, rand_len = args.stem_len, args.loop_len, args.randomer

    k = args.max_internal_foldback + 1
    if k < 2:
        ap.error("--max-internal-foldback must be at least 1")
    if min(prefix_len, rand_len, spacer_len, stem_len, loop_len) < 0:
        ap.error("block lengths must be non-negative")
    if stem_len:
        if loop_len < 3:
            ap.error("--loop-len must be at least 3 for a hairpin")
        if args.stem_min_gc > stem_len:
            ap.error("--stem-min-gc cannot exceed --stem-len")
    else:
        loop_len = 0
    total = prefix_len + rand_len + spacer_len + 2 * stem_len + loop_len
    if total < k:
        ap.error(f"the layout is only {total} nt, too short to check")

    rand_alphabet = args.randomer_alphabet.strip().upper().replace("T", "U")
    bad = sorted(set(rand_alphabet) - set(BASES))
    if not rand_alphabet or bad:
        ap.error(f"--randomer-alphabet: not an ACGT/ACGU subset "
                 f"(offending: {''.join(bad) or 'empty'})")

    rand_start = prefix_len
    arm1_start = prefix_len + rand_len + spacer_len
    arm2_start = total - stem_len
    want = {(arm1_start + n, total - 1 - n) for n in range(stem_len)}

    block_alphabet = args.block_alphabet.strip().upper().replace("T", "U")
    bad = sorted(set(block_alphabet) - set(BASES))
    if not block_alphabet or bad:
        ap.error(f"--block-alphabet: not an ACGT/ACGU subset "
                 f"(offending: {''.join(bad) or 'empty'})")

    alphabets = [block_alphabet] * total
    for i in range(arm1_start, total):
        alphabets[i] = BASES
    for i in range(rand_start, rand_start + rand_len):
        alphabets[i] = rand_alphabet
    forced = [None] * total
    for i, c in enumerate(prefix_seq or ""):
        forced[i] = c
    for i, c in enumerate(spacer_seq or ""):
        forced[prefix_len + rand_len + i] = c

    problems = []
    if prefix_seq:
        problems += preflight_block("prefix", prefix_seq,
                                    args.max_internal_foldback, partners)
    if spacer_seq:
        problems += preflight_block("spacer", spacer_seq,
                                    args.max_internal_foldback, partners)
    if problems:
        for line in problems:
            print(line, file=sys.stderr)
        print("No molecule can satisfy the checks with these fixed blocks. Raise "
              "--max-internal-foldback, or supply a block that does not fold.",
              file=sys.stderr)
        sys.exit(2)

    rng = random.Random(args.seed)
    started = time.monotonic()
    for attempt in range(1, args.max_attempts + 1):
        if time.monotonic() - started > args.timeout:
            print(f"Gave up after {args.timeout:g}s and {attempt - 1} attempts. Try "
                  f"raising --timeout or --max-internal-foldback, or narrowing "
                  f"--randomer-alphabet.", file=sys.stderr)
            sys.exit(1)

        seq = build(total, forced, alphabets, arm1_start, arm2_start, stem_len, k,
                    args.max_homopolymer, args.max_kmer_reuse, partners, rng,
                    args.budget)
        if seq is None:
            continue

        if stem_len:
            arm = seq[arm1_start:arm1_start + stem_len]
            if sum(c in "GC" for c in arm) < args.stem_min_gc:
                continue

        full_struct, full_mfe = RNA.fold(seq)
        if mfe_pairs(full_struct) != want:
            continue

        free = seq[:total - stem_len]
        free_struct, free_mfe = RNA.fold(free)
        if free_struct.count(".") != len(free_struct):
            continue

        worst = longest_foldback(seq, want, partners)
        if worst > args.max_internal_foldback:
            continue

        out = seq.replace("U", "T") if args.dna else seq
        blocks = [("prefix", 0, prefix_len, "P"),
                  ("randomer", rand_start, rand_len, "N"),
                  ("spacer", arm1_start - spacer_len, spacer_len, "S"),
                  ("arm", arm1_start, stem_len, "<"),
                  ("loop", arm1_start + stem_len, loop_len, "-"),
                  ("arm'", arm2_start, stem_len, ">")]

        print(f"Sequence ({total} nt): {out}")
        print(f"Blocks:            {''.join(m * n for _, _, n, m in blocks)}")
        print(f"MFE structure:     {full_struct}  ({full_mfe:.2f} kcal/mol)")
        print(f"Minus 3' arm:      {free_struct}  ({free_mfe:.2f} kcal/mol)")
        print()
        for name, start, n, mark in blocks:
            if n:
                print(f"  {name:9s} [{start + 1:3d}-{start + n:3d}] {mark * 2}  "
                      f"{out[start:start + n]}")
        print()
        if stem_len:
            print(f"3' self-priming foldback: {stem_len} bp, loop {loop_len} nt; "
                  f"3' terminal base pairs with position {arm1_start + 1}")
        print(f"Longest other foldback:   {worst} nt (limit "
              f"{args.max_internal_foldback}, G-U "
              f"{'counted' if args.wobble else 'ignored'})")
        if rand_len and args.randomer_trials:
            r_worst, r_extra = survey_randomer(seq, rand_start, rand_len,
                                               rand_alphabet, want, partners,
                                               args.randomer_trials, rng)
            print(f"Randomer survey ({args.randomer_trials} draws over "
                  f"{rand_alphabet.replace('U', 'T') if args.dna else rand_alphabet}):"
                  f" worst foldback {r_worst} nt, {r_extra} draw(s) folded beyond "
                  f"the designed stem")
        print(f"Attempts needed: {attempt}")
        return

    print(f"Failed to find a qualifying sequence in {args.max_attempts} attempts. "
          f"Try raising --max-internal-foldback or narrowing --randomer-alphabet.",
          file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
