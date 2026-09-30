#!/usr/bin/env python3
"""
Design a self-priming molecule laid out as

    [ X const ][ P degenerate (N) ][ Z const ][ B const ][ L const ][ D ]

with D the reverse complement of B, so the 3' end folds back onto B and primes
extension back through Z, P and X.

Two separate constraints do the work, and conflating them is what produces a
degenerate sequence:

  1. D must not find a partner anywhere else -- nowhere outside B may a window
     be even --min-complementarity nt complementary to D.  Because D is the
     reverse complement of B, such a duplex means that window holds a 2-mer of
     B, so the rule collapses to a small set of forbidden 2-mers derived from B
     once up front.

  2. The constant regions must not fold.  This does NOT require an alphabet
     that cannot pair -- banning two bases outright is what leaves a scaffold
     with no T in it at all.  It requires only that no helix can stack up:
     no window may be more than --max-duplex nt complementary to any
     non-overlapping window elsewhere.  All four bases stay available, and T
     appears wherever it cannot stack into a real helix.

Both are enforced during a single left-to-right walk, so a candidate is built
rather than sampled and rejected; the iteration then only retries for the
folding checks, which cannot be enforced by construction.  --min-base-fraction
holds every base above a floor so the result is a real sequence, not a
two-letter one.

G-U counts as a pair unless --no-wobble, since ViennaRNA forms G-U and ignoring
it admits helices the fold will find.

Checked on every candidate:
  a. No window outside B is >= --min-complementarity nt complementary to D.
  b. No duplex anywhere exceeds --max-duplex nt, the designed B:D stem aside.
  c. X and Z, the constant regions before B, each fold to nothing alone.
  d. The molecule minus D folds to nothing.
  e. The MFE structure of the whole molecule is exactly the designed B:D stem.
  f. G+C within [--min-gc, --max-gc], and every base above --min-base-fraction.

Reported, not enforced: P is degenerate, so only the emitted draw can be
constrained.  --p-trials resamples P and reports how much of the library breaks
rule (a).  For a 15 nt N region that is most of it -- 14 windows against 16
possible 2-mers -- but a 2 bp duplex is worth ~0 kcal/mol and cannot prime, so
watch the "fold beyond the stem" figure instead.

Usage:
    python3 design_dragonrna.py
    python3 design_dragonrna.py -x 20 -p 15 -z 20 -b 6 -l 4
    python3 design_dragonrna.py --max-duplex 2 --min-base-fraction 0.12
"""
import argparse
import itertools
import math
import random
import sys
import time

import RNA

BASES = "ACGU"
WATSON_CRICK = {"A": "U", "C": "G", "G": "C", "U": "A"}
WITH_WOBBLE = {"A": "U", "C": "G", "G": "CU", "U": "AG"}

REGION_NAMES = "XPZBLD"
DEGENERATE_REGIONS = {"P"}   # drawn per molecule, so ordered as N


def reverse_complement(seq):
    return "".join(WATSON_CRICK[c] for c in reversed(seq))


def partner_mers(window, partners):
    """Every window that could pair along the whole length of `window`.

    Pairing A against B lines a_t up with b_{n-1-t}, so B is drawn from the
    partners of A read backwards.
    """
    return ("".join(combo) for combo in
            itertools.product(*(partners[c] for c in reversed(window))))


def pairable_mers(seq, partners, n):
    """Every n-mer that could form an n-bp duplex with some window of `seq`."""
    out = set()
    for i in range(len(seq) - n + 1):
        out.update(partner_mers(seq[i:i + n], partners))
    return out


def longest_duplex(left, right, partners):
    """Longest contiguous antiparallel duplex between two sequences."""
    best = 0
    for i in range(len(left)):
        for j in range(len(right)):
            run = 0
            while (i + run < len(left) and j - run >= 0
                   and right[j - run] in partners[left[i + run]]):
                run += 1
            best = max(best, run)
    return best


def longest_self_duplex(seq, ignore, partners):
    """Longest contiguous duplex a sequence admits with itself."""
    best = 0
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            run = 0
            while (i + run < j - run
                   and seq[j - run] in partners[seq[i + run]]
                   and (i + run, j - run) not in ignore):
                run += 1
            best = max(best, run)
    return best


def register_slip(b_seq, d_seq, partners):
    """Longest B:D duplex outside the intended register.

    A contiguous B:D run pairs b_seq[u] with d_seq[v] and keeps u+v constant, so
    each register is one value of u+v and the intended one is len(b)-1.  Any
    other register lets the 3' end anneal a position or two over and prime in
    the wrong place: with B=CCCGGC the repeated CC gives D=GCCGGG a second
    2 bp register.  For Watson-Crick this reduces to "every n-mer of B is
    distinct", but the scan also covers G-U.
    """
    b = len(b_seq)
    worst = 0
    for offset in range(2 * b - 1):
        if offset == b - 1:
            continue
        run = 0
        for u in range(b):
            v = offset - u
            if 0 <= v < b and d_seq[v] in partners[b_seq[u]]:
                run += 1
                worst = max(worst, run)
            else:
                run = 0
    return worst


def offregister_duplex(seq, d_span, intended, partners):
    """Longest duplex involving D anywhere in the molecule, intended stem aside.

    Unlike the per-region check this includes B itself, so a slipped register
    against the designed partner cannot hide.
    """
    d0, d1 = d_span
    worst = 0
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            run = 0
            while (i + run < j - run
                   and seq[j - run] in partners[seq[i + run]]
                   and (i + run, j - run) not in intended):
                run += 1
            if not run:
                continue
            lo = [d0 <= i + t < d1 for t in range(run)]
            hi = [d0 <= j - t < d1 for t in range(run)]
            if all(lo) and all(hi):
                continue          # D against itself; it cannot hairpin in d nt
            if any(lo) or any(hi):
                worst = max(worst, run)
    return worst


def chimera_notation(seq, regions, dna_regions, degenerate=False):
    """Render the molecule in IDT chimera syntax: rA/rC/rG/rU for the RNA blocks,
    plain A/C/G/T for the blocks named in `dna_regions`.

    With `degenerate` the blocks in DEGENERATE_REGIONS come out as N runs, which
    is what you order for the library; without it they carry the bases actually
    drawn.  Which blocks are degenerate is a property of the design, not of the
    chemistry, so this is independent of `dna_regions` -- a degenerate RNA block
    is rN, and every constant block stays literal however it is rendered.
    """
    parts = []
    for name, start, length, _ in regions:
        if not length:
            continue
        block = seq[start:start + length]
        if degenerate and name in DEGENERATE_REGIONS:
            block = "N" * length
        if name in dna_regions:
            parts.append(block.replace("U", "T"))
        else:
            parts.append("".join("r" + base for base in block))
    return "".join(parts)


def mfe_pairs(structure):
    stack, pairs = [], set()
    for i, c in enumerate(structure):
        if c == "(":
            stack.append(i)
        elif c == ")":
            pairs.add((stack.pop(), i))
    return pairs


def unpaired(seq):
    struct, mfe = RNA.fold(seq)
    return struct.count(".") == len(struct), struct, mfe


def draw_molecule(total, alphabets, forced, b_span, d_span, forbidden2, k,
                  partners, max_homopolymer, quota, rng, budget):
    """Walk the molecule left to right, enforcing both constraints as it goes.

    Backtracks on a dead end, and prefers bases still owed quota so no base is
    squeezed out of the sequence altogether.
    """
    seq = []
    kmer_starts = {}
    used = dict.fromkeys(BASES, 0)
    steps = [0]

    def inside(i, span):
        return span[0] <= i < span[1]

    def spans_window(s, span):
        return span[0] <= s and s + k <= span[1]

    # Pairing window [s,s+k) against [j,j+k) matches s+t with j+k-1-t, so the two
    # windows sit in the intended B:D register only when s+j is this constant.
    # Exempting every B-window/D-window pair instead is what let the stem slip.
    in_register = b_span[0] + d_span[1] - k

    def designed_stem(s, j):
        if not ((spans_window(s, b_span) and spans_window(j, d_span))
                or (spans_window(s, d_span) and spans_window(j, b_span))):
            return False
        return s + j == in_register

    def place(i):
        if i == total:
            return True
        steps[0] += 1
        if steps[0] > budget:
            return False
        if forced[i]:
            options = [forced[i]]
        else:
            options = list(alphabets[i])
            rng.shuffle(options)
            options.sort(key=lambda base: used[base] - quota[base])
        for base in options:
            if (max_homopolymer and len(seq) >= max_homopolymer
                    and all(c == base for c in seq[-max_homopolymer:])):
                continue
            # Constraint 1. A 2-mer wholly inside B is the source of the
            # forbidden set, and one wholly inside D is D against itself, which
            # is not "the rest of the sequence".
            if i and not ((inside(i - 1, b_span) and inside(i, b_span))
                          or (inside(i - 1, d_span) and inside(i, d_span))):
                if seq[-1] + base in forbidden2:
                    continue
            seq.append(base)
            used[base] += 1
            ok = True
            added = None
            s = i - k + 1
            if s >= 0:
                kmer = "".join(seq[s:])
                # Constraint 2. Pairing seq[j:j+k] against seq[s:s+k] puts
                # (j, s+k-1) outermost and (j+k-1, s) innermost, so a full k-bp
                # duplex needs the windows disjoint; overlapping ones top out
                # below k and are not violations.
                for mate in partner_mers(kmer, partners):
                    for j in kmer_starts.get(mate, ()):
                        if abs(s - j) < k or designed_stem(s, j):
                            continue
                        ok = False
                        break
                    if not ok:
                        break
                if ok:
                    kmer_starts.setdefault(kmer, []).append(s)
                    added = kmer
            if ok and place(i + 1):
                return True
            if added is not None:
                kmer_starts[added].pop()
                if not kmer_starts[added]:
                    del kmer_starts[added]
            used[base] -= 1
            seq.pop()
        return False

    return "".join(seq) if place(0) else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-x", "--x", type=int, default=20,
                    help="5' constant region, nt (default: 20)")
    ap.add_argument("-p", "--p", type=int, default=15,
                    help="degenerate N region, nt (default: 15)")
    ap.add_argument("-z", "--z", type=int, default=20,
                    help="constant region after P, nt (default: 20)")
    ap.add_argument("-b", "--b", type=int, default=6,
                    help="constant region that D reverse-complements, nt (default: 6)")
    ap.add_argument("-l", "--l", type=int, default=4,
                    help="constant region between B and D, nt (default: 4)")
    ap.add_argument("-d", "--d", type=int, default=None,
                    help="3' region, reverse complement of B; defaults to B's length")
    ap.add_argument("--min-complementarity", type=int, default=2,
                    help="duplex length with D that must not occur outside B "
                         "(default: 2, i.e. not even 2 nt)")
    ap.add_argument("--max-duplex", type=int, default=3,
                    help="longest duplex allowed anywhere outside the designed B:D "
                         "stem. Lower forces a blander sequence; 2 is very tight at "
                         "70+ nt (default: 3)")
    ap.add_argument("--slip-wobble", action=argparse.BooleanOptionalAction,
                    default=False,
                    help="also count G-U when checking the stem's register. Off by "
                         "default: what makes a register ambiguous is a Watson-Crick "
                         "repeat in B, and a 2 bp off-register duplex held together "
                         "by a wobble will not prime (default: off)")
    ap.add_argument("--wobble", action=argparse.BooleanOptionalAction, default=True,
                    help="count G-U as a pair, as ViennaRNA does (default: on)")
    ap.add_argument("--b-min-gc", type=int, default=4,
                    help="minimum G+C in B, so the stem actually folds (default: 4)")
    ap.add_argument("--min-base-fraction", type=float, default=0.10,
                    help="every base must be at least this fraction of the molecule, "
                         "so the sequence does not collapse to two letters "
                         "(default: 0.10)")
    ap.add_argument("--min-gc", type=float, default=0.30,
                    help="minimum G+C fraction overall; never 0 (default: 0.30)")
    ap.add_argument("--max-gc", type=float, default=0.70,
                    help="maximum G+C fraction overall (default: 0.70)")
    ap.add_argument("--min-gc-scaffold", type=float, default=0.10,
                    help="minimum G+C fraction across X+P+Z alone, which has to fold "
                         "to nothing and so runs lower (default: 0.10)")
    ap.add_argument("--max-homopolymer", type=int, default=3,
                    help="longest run of one base allowed; 0 for no limit (default: 3)")
    ap.add_argument("--const-alphabet", default="ACGU",
                    help="bases X, Z and L may use (default: ACGU)")
    ap.add_argument("--p-alphabet", default="ACGU",
                    help="bases the degenerate region may use (default: ACGU)")
    ap.add_argument("--dna-regions", default=None,
                    help="comma-separated region names (X,P,Z,B,L,D) to render as DNA "
                         "in the chimera line; the rest come out as rA/rC/rG/rU "
                         "(default: P)")
    ap.add_argument("--rna-regions", default=None,
                    help="the inverse: name the RNA regions and the rest become DNA. "
                         "--rna-regions P is a lone ribonucleotide in an otherwise-DNA "
                         "molecule. Give this or --dna-regions, not both")
    ap.add_argument("--dna", action=argparse.BooleanOptionalAction, default=True,
                    help="report T in place of U (default: on)")
    ap.add_argument("--p-trials", type=int, default=200,
                    help="resample P this many times to survey the library; 0 to skip "
                         "(default: 200)")
    ap.add_argument("--budget", type=int, default=20_000,
                    help="DFS steps before abandoning an attempt. Restarting cheaply "
                         "beats deep backtracking here (default: 20000)")
    ap.add_argument("--max-attempts", type=int, default=100_000)
    ap.add_argument("--timeout", type=float, default=90.0, help="seconds (default: 90)")
    ap.add_argument("--seed", type=int, default=None)
    args = ap.parse_args()

    x, p, z, b, l = args.x, args.p, args.z, args.b, args.l
    d = args.d if args.d is not None else b
    if min(x, p, z, b, l, d) < 0:
        ap.error("region lengths must be non-negative")
    if b < 1:
        ap.error("-b must be at least 1 for a foldback")
    if d != b:
        ap.error(f"-d is the reverse complement of -b, so they must match "
                 f"(got d={d}, b={b})")
    if l < 3:
        ap.error("-l must be at least 3: a hairpin loop shorter than 3 nt cannot form")
    n = args.min_complementarity
    if n < 2:
        ap.error("--min-complementarity must be at least 2")
    k = args.max_duplex + 1
    if args.max_duplex < 1:
        ap.error("--max-duplex must be at least 1")
    if not 0 < args.min_gc <= args.max_gc <= 1:
        ap.error("need 0 < --min-gc <= --max-gc <= 1")
    if args.min_gc_scaffold <= 0:
        ap.error("--min-gc-scaffold must be above 0: the scaffold cannot be 0% GC")
    if not 0 <= args.min_base_fraction <= 0.25:
        ap.error("--min-base-fraction must be between 0 and 0.25")
    if args.b_min_gc > b:
        ap.error("--b-min-gc cannot exceed -b")

    def alphabet(value, name):
        seq = value.strip().upper().replace("T", "U")
        bad = sorted(set(seq) - set(BASES))
        if not seq or bad:
            ap.error(f"--{name}: not an ACGT/ACGU subset "
                     f"(offending: {''.join(bad) or 'empty'})")
        return seq

    def region_set(spec, flag):
        names = {r.strip().upper() for r in spec.split(",") if r.strip()}
        unknown = sorted(names - set(REGION_NAMES))
        if unknown:
            ap.error(f"--{flag}: unknown region(s) {','.join(unknown)}; "
                     f"choose from {','.join(REGION_NAMES)}")
        return names

    if args.dna_regions is not None and args.rna_regions is not None:
        ap.error("give --dna-regions or --rna-regions, not both; each implies the other")
    if args.rna_regions is not None:
        dna_regions = set(REGION_NAMES) - region_set(args.rna_regions, "rna-regions")
    else:
        dna_regions = region_set(
            args.dna_regions if args.dna_regions is not None else "P", "dna-regions")

    const_alphabet = alphabet(args.const_alphabet, "const-alphabet")
    p_alphabet = alphabet(args.p_alphabet, "p-alphabet")
    partners = WITH_WOBBLE if args.wobble else WATSON_CRICK
    slip_partners = WITH_WOBBLE if args.slip_wobble else WATSON_CRICK

    total = x + p + z + b + l + d
    b_start = x + p + z
    d_start = total - d
    b_span, d_span = (b_start, b_start + b), (d_start, total)
    want = {(b_start + i, total - 1 - i) for i in range(b)}
    floor = math.ceil(args.min_base_fraction * total)
    quota = dict.fromkeys(BASES, floor)

    alphabets = [const_alphabet] * total
    for i in range(x, x + p):
        alphabets[i] = p_alphabet
    for i in range(b_start, total):
        alphabets[i] = BASES

    rng = random.Random(args.seed)
    started = time.monotonic()
    rejected = dict.fromkeys(
        ["draw", "slip", "gc", "scaffold_gc", "composition", "xz", "free", "mfe",
         "dup"], 0)

    def fraction_gc(part):
        return sum(c in "GC" for c in part) / len(part) if part else 0.0

    for attempt in range(1, args.max_attempts + 1):
        if time.monotonic() - started > args.timeout:
            print(f"Gave up after {args.timeout:g}s and {attempt - 1} attempts "
                  f"({rejected}). Try --max-duplex {args.max_duplex + 1}, a lower "
                  f"--min-base-fraction, or --no-wobble.", file=sys.stderr)
            sys.exit(1)

        # Place B's G/C directly rather than sampling and filtering: at 4 of 6 a
        # uniform draw is rejected two thirds of the time.  Retry B in this inner
        # loop rather than burning a whole outer attempt, since most draws slip.
        b_seq = d_seq = None
        for _ in range(1000):
            n_gc = rng.randint(args.b_min_gc, b)
            gc_at = set(rng.sample(range(b), n_gc))
            cand = "".join(rng.choice("GC") if i in gc_at else rng.choice("AU")
                           for i in range(b))
            # A B that D can pair with in more than one register lets the 3' end
            # prime a position or two over.  Under Watson-Crick this is exactly a
            # repeated n-mer in B, which is what CCCGGC (CC twice) had.
            if register_slip(cand, reverse_complement(cand), slip_partners) < n:
                b_seq, d_seq = cand, reverse_complement(cand)
                break
        if b_seq is None:
            rejected["slip"] += 1
            continue
        forbidden2 = pairable_mers(d_seq, partners, n)

        forced = [None] * total
        for i, c in enumerate(b_seq):
            forced[b_start + i] = c
        for i, c in enumerate(d_seq):
            forced[d_start + i] = c

        seq = draw_molecule(total, alphabets, forced, b_span, d_span, forbidden2,
                            k, partners, args.max_homopolymer, quota, rng,
                            args.budget)
        if seq is None:
            rejected["draw"] += 1
            continue

        counts = {base: seq.count(base) for base in BASES}
        if min(counts.values()) < floor:
            rejected["composition"] += 1
            continue
        gc_all = fraction_gc(seq)
        if not args.min_gc <= gc_all <= args.max_gc:
            rejected["gc"] += 1
            continue
        scaffold = seq[:b_start]
        gc_scaffold = fraction_gc(scaffold)
        if gc_scaffold < args.min_gc_scaffold:
            rejected["scaffold_gc"] += 1
            continue

        x_seq, z_seq = seq[:x], seq[x + p:x + p + z]
        x_ok, x_struct, x_mfe = unpaired(x_seq) if x else (True, "", 0.0)
        z_ok, z_struct, z_mfe = unpaired(z_seq) if z else (True, "", 0.0)
        if not (x_ok and z_ok):
            rejected["xz"] += 1
            continue

        free_ok, free_struct, free_mfe = unpaired(seq[:d_start])
        if not free_ok:
            rejected["free"] += 1
            continue

        full_struct, full_mfe = RNA.fold(seq)
        if mfe_pairs(full_struct) != want:
            rejected["mfe"] += 1
            continue

        l_seq = seq[b_start + b:d_start]
        worst_d = max((longest_duplex(d_seq, part, partners), name) for name, part
                      in [("X", x_seq), ("P", seq[x:x + p]), ("Z", z_seq),
                          ("L", l_seq)] if part)
        worst_any = longest_self_duplex(seq, want, partners)
        worst_off = offregister_duplex(seq, d_span, want, slip_partners)
        if worst_d[0] >= n or worst_any > args.max_duplex or worst_off >= n:
            rejected["dup"] += 1
            continue

        placed_b, placed_d = seq[b_start:b_start + b], seq[d_start:]
        if placed_d != reverse_complement(placed_b):
            raise AssertionError(f"D is not the reverse complement of B: B={placed_b} "
                                 f"D={placed_d} expected {reverse_complement(placed_b)}")

        out = seq.replace("U", "T") if args.dna else seq
        regions = [("X", 0, x, "X"), ("P", x, p, "N"), ("Z", x + p, z, "Z"),
                   ("B", b_start, b, "<"), ("L", b_start + b, l, "-"),
                   ("D", d_start, d, ">")]
        print(f"Molecule ({total} nt = {x}+{p}N+{z}+{b}+{l}+{d}):")
        print(f"  {out}")
        print(f"  {''.join(m * ln for _, _, ln, m in regions)}")
        print(f"  {full_struct}   MFE {full_mfe:.2f} kcal/mol")
        print()
        rna_names = [nm for nm, _, ln, _ in regions if ln and nm not in dna_regions]
        print(f"  chimera, this instance (RNA: {'+'.join(rna_names) or 'none'}; "
              f"DNA: {'+'.join(sorted(dna_regions)) or 'none'}):")
        print(f"    {chimera_notation(seq, regions, dna_regions)}")
        print(f"  chimera, to order the library:")
        print(f"    {chimera_notation(seq, regions, dna_regions, degenerate=True)}")
        print()
        for name, start, ln, _ in regions:
            if ln:
                print(f"  {name} [{start + 1:3d}-{start + ln:3d}] {out[start:start + ln]}")
        print()
        print(f"  D = reverse complement of B: {out[b_start:b_start + b]} / {out[d_start:]}")
        print(f"  3' terminal base pairs with position {b_start + 1}; loop {l} nt")
        print(f"  Longest duplex, D vs outside B: {worst_d[0]} nt (in {worst_d[1]}; "
              f"must stay under {n})")
        print(f"  Longest duplex, D off-register anywhere (B included): "
              f"{worst_off} nt (must stay under {n})")
        print(f"  Longest duplex anywhere else:   {worst_any} nt (limit "
              f"{args.max_duplex}; G-U {'counted' if args.wobble else 'ignored'})")
        shown = {base.replace("U", "T") if args.dna else base: c
                 for base, c in counts.items()}
        print(f"  Composition: {shown}, each >= {floor}; "
              f"G+C {gc_all * 100:.0f}% overall, {gc_scaffold * 100:.0f}% across X+P+Z")
        print(f"  X alone: {x_struct or '-'} ({x_mfe:.2f})")
        print(f"  Z alone: {z_struct or '-'} ({z_mfe:.2f})")
        print(f"  Molecule minus D: all unpaired ({free_mfe:.2f})")
        if p and args.p_trials:
            breaks = extra = 0
            for _ in range(args.p_trials):
                draw = "".join(rng.choices(p_alphabet, k=p))
                cand = seq[:x] + draw + seq[x + p:]
                if longest_duplex(d_seq, draw, partners) >= n:
                    breaks += 1
                st, _ = RNA.fold(cand)
                if mfe_pairs(st) != want:
                    extra += 1
            print(f"  Library survey ({args.p_trials} P draws): {breaks} give D a "
                  f">={n} nt partner in P, {extra} fold beyond the B:D stem")
        print(f"\n  attempts: {attempt}  rejected: {rejected}")
        return

    print(f"No molecule found in {args.max_attempts} attempts ({rejected}).",
          file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
