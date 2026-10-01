# DragonRNA oligo design

Sequence design for self-priming foldback oligos and for oligos that must not
fold at all. Candidates are drawn at random and then filtered against
[ViennaRNA](https://www.tbi.univie.ac.at/RNA/) folding plus a set of
alphabet-level combinatorial constraints, so the reported molecule is one that
provably cannot form long duplexes anywhere — not merely one whose minimum-free-
energy structure happens to look clean.

Code was written by : Haim Krupkin+Claude 

## Installation

Requires Python 3.9 or newer. The only dependency is
[ViennaRNA](https://www.tbi.univie.ac.at/RNA/), pulled in automatically by
pip.

### Install from GitHub

```bash
pip install git+https://github.com/haimkru/dragonrna-design
```

### Install from a clone

```bash
git clone https://github.com/haimkru/dragonrna-design
cd dragonrna-design
pip install .
```
The scripts have no imports beyond the standard library and ViennaRNA, so you
can clone and run them in place:

```bash
pip install ViennaRNA          # or: conda install -c bioconda viennarna
python3 design_dragonrna.py --help
```
```bash
conda install -c bioconda viennarna
```

## Example

Design a 49 nt chimera whose constant regions are DNA and whose single
degenerate position is RNA:

```bash
design-dragonrna -x 20 -p 1 -z 12 -b 6 -l 4 --rna-regions P --seed 7
```

(or `python3 design_dragonrna.py -x 20 ...` from a clone, if you skipped the
install)

```
Molecule (49 nt = 20+1N+12+6+4+6):
  TACAGATACAGAGATACATAGACATACAGACAACCGGCTAGAGAGCCGG
  XXXXXXXXXXXXXXXXXXXXNZZZZZZZZZZZZ<<<<<<---->>>>>>
  .................................((((((....))))))   MFE -9.90 kcal/mol

  chimera, this instance (RNA: P; DNA: B+D+L+X+Z):
    TACAGATACAGAGATACATArGACATACAGACAACCGGCTAGAGAGCCGG
  chimera, to order the library:
    TACAGATACAGAGATACATArNACATACAGACAACCGGCTAGAGAGCCGG
```

Drop `--seed 7` to draw a different molecule meeting the same constraints.

The last line is the sequence to order: constant positions literal, the
degenerate position as `rN`, in IDT chimera syntax.

## Layout

`design_dragonrna.py` builds a molecule of six blocks, each sized by its own
flag:

```
[ X const ][ P degenerate (N) ][ Z const ][ B const ][ L loop ][ D ]
     -x            -p               -z         -b        -l     -d
```

`D` is the reverse complement of `B`, so the 3' end folds back onto `B` and
primes extension back through `Z`, `P` and `X`. Everything else must stay
unpaired.


### Choosing the chemistry

`--dna-regions` and `--rna-regions` are inverses; pass one or the other.

```bash
--dna-regions P      # default: P is DNA, X/Z/B/L/D are RNA
--rna-regions P      # P is RNA, X/Z/B/L/D are DNA
--rna-regions ''     # every block DNA
```

## The other scripts

- **`generate_dragonrna.py`** — an earlier foldback designer. Same idea, but
  its regions are named `--prefix / --randomer / --spacer / --arm / --loop`,
  `--prefix` and `--spacer` accept a literal sequence to hold fixed, and it
  bounds internal structure with a single `--max-internal-foldback` k-mer rule
  rather than the per-region checks above.
- **`generate_nonfolding_oligo.py`** — draws an oligo with no self-structure
  at all: the MFE structure of the full molecule must be all dots, and the
  last `--three-prime-len` nt must also fold to nothing in isolation.

## License

MIT — see [LICENSE](LICENSE).

For questiosn reach out to hkrupkin@stanford.edu
