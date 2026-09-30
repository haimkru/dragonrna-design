#!/usr/bin/env python3
"""
Generate a random-sequence oligo that has no self-structure anywhere in its
length, and additionally verify the 3' end has no self-structure on its own.

"No self-structure" is defined via ViennaRNA MFE folding: the minimum free
energy structure must contain zero base pairs (an all-dots dot-bracket
string). This is checked (1) for the full-length molecule, and (2) for the
last --three-prime-len nt folded in isolation.

Usage:
    python3 generate_nonfolding_oligo.py
    python3 generate_nonfolding_oligo.py --length 40 --three-prime-len 15 --dna
"""
import argparse
import random
import sys

import RNA


def random_sequence(length, alphabet, three_prime_len=0, three_prime_alphabet=None):
    """Build the sequence with an optionally restricted alphabet for the 3' end."""
    if three_prime_len <= 0 or three_prime_alphabet is None:
        return "".join(random.choices(alphabet, k=length))
    three_prime_len = min(three_prime_len, length)
    body_len = length - three_prime_len
    body = "".join(random.choices(alphabet, k=body_len))
    tail = "".join(random.choices(three_prime_alphabet, k=three_prime_len))
    return body + tail


def is_unstructured(seq):
    """Fold with ViennaRNA MFE; require an all-dots structure (no base pairs)."""
    structure, mfe = RNA.fold(seq)
    return structure.count(".") == len(structure), structure, mfe


def max_run_length(seq, base):
    longest = current = 0
    for c in seq:
        current = current + 1 if c == base else 0
        longest = max(longest, current)
    return longest


def has_homopolymer_run(seq, max_run, bases="ACGU"):
    return any(max_run_length(seq, base) > max_run for base in bases)


COMPLEMENT = {"A": "U", "U": "A", "C": "G", "G": "C"}


def reverse_complement(seq):
    return "".join(COMPLEMENT[c] for c in reversed(seq))


def three_prime_revcomp_found_elsewhere(seq, n):
    """True if the reverse complement of the last n nt occurs anywhere else in seq
    (i.e. some other n-mer could base-pair with the 3' terminus)."""
    tail_start = len(seq) - n
    target = reverse_complement(seq[tail_start:])
    for i in range(len(seq) - n + 1):
        if i == tail_start:
            continue
        if seq[i:i + n] == target:
            return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--length", type=int, default=40, help="total oligo length (default: 40)")
    ap.add_argument("--three-prime-len", type=int, default=15,
                     help="length of 3' end window checked in isolation (default: 15)")
    ap.add_argument("--dna", action="store_true",
                     help="emit T instead of U (folding still uses RNA energy params)")
    ap.add_argument("--no-t-in-3prime", action="store_true",
                     help="exclude T/U from the last --t-exclusion-len nt (3' end)")
    ap.add_argument("--t-exclusion-len", type=int, default=10,
                     help="length of the 3' end window with no T/U (default: 10)")
    ap.add_argument("--no-polya-3prime", action="store_true",
                     help="reject a poly-A run longer than --polya-max-run in the last "
                          "--polya-window-len nt (3' end)")
    ap.add_argument("--polya-max-run", type=int, default=2,
                     help="longest allowed run of consecutive A's in the 3' end window "
                          "(default: 2, i.e. AAA+ is rejected)")
    ap.add_argument("--polya-window-len", type=int, default=10,
                     help="length of the 3' end window checked for poly-A runs (default: 10)")
    ap.add_argument("--no-homopolymer-3prime", action="store_true",
                     help="reject any run (of A, C, G, or T/U) longer than --homopolymer-max-run "
                          "in the last --homopolymer-window-len nt (3' end); supersedes "
                          "--no-polya-3prime if both are set")
    ap.add_argument("--homopolymer-max-run", type=int, default=2,
                     help="longest allowed run of any single base in the 3' end window "
                          "(default: 2, i.e. any XXX+ is rejected)")
    ap.add_argument("--homopolymer-window-len", type=int, default=10,
                     help="length of the 3' end window checked for homopolymer runs (default: 10)")
    ap.add_argument("--no-3prime-revcomp", action="store_true",
                     help="reject sequences where the reverse complement of the last "
                          "--revcomp-check-len nt occurs anywhere else in the sequence "
                          "(prevents the 3' terminus from self-priming/annealing elsewhere)")
    ap.add_argument("--revcomp-check-len", type=int, default=2,
                     help="length of the 3' terminal window checked for a reverse-complement "
                          "match elsewhere (default: 2)")
    ap.add_argument("--max-attempts", type=int, default=200_000)
    ap.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    args = ap.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    alphabet = list("ACGU")
    restricted_alphabet = list("ACG") if args.no_t_in_3prime else None
    restricted_len = args.t_exclusion_len if args.no_t_in_3prime else 0

    for attempt in range(1, args.max_attempts + 1):
        seq = random_sequence(args.length, alphabet, restricted_len, restricted_alphabet)

        full_ok, full_struct, full_mfe = is_unstructured(seq)
        if not full_ok:
            continue

        three_prime_seq = seq[-args.three_prime_len:]
        tp_ok, tp_struct, tp_mfe = is_unstructured(three_prime_seq)
        if not tp_ok:
            continue

        if args.no_polya_3prime:
            polya_window = seq[-args.polya_window_len:]
            if max_run_length(polya_window, "A") > args.polya_max_run:
                continue

        if args.no_homopolymer_3prime:
            homopolymer_window = seq[-args.homopolymer_window_len:]
            if has_homopolymer_run(homopolymer_window, args.homopolymer_max_run):
                continue

        if args.no_3prime_revcomp:
            if three_prime_revcomp_found_elsewhere(seq, args.revcomp_check_len):
                continue

        out_seq = seq.replace("U", "T") if args.dna else seq
        print(f"Sequence ({len(out_seq)} nt): {out_seq}")
        print(f"Full-length structure:  {full_struct}  (MFE {full_mfe:.2f} kcal/mol)")
        print(f"3' end ({args.three_prime_len} nt) structure: {tp_struct}  (MFE {tp_mfe:.2f} kcal/mol)")
        print(f"Attempts needed: {attempt}")
        return

    print(f"Failed to find a qualifying sequence in {args.max_attempts} attempts.", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
