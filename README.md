# DragonRNA oligo design

Sequence design for self-priming foldback oligos and for oligos that must not
fold at all. Candidates are drawn at random and then filtered against
[ViennaRNA](https://www.tbi.univie.ac.at/RNA/) folding plus a set of
alphabet-level combinatorial constraints, so the reported molecule is one that
provably cannot form long duplexes anywhere — not merely one whose minimum-free-
energy structure happens to look clean.

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

Use `pip install -e .` instead if you want your edits to take effect without
reinstalling.

Either route puts three commands on your `PATH`, runnable from any directory:

| Command | Script |
| --- | --- |
| `design-dragonrna` | `design_dragonrna.py` |
| `generate-dragonrna` | `generate_dragonrna.py` |
| `generate-nonfolding-oligo` | `generate_nonfolding_oligo.py` |

### Without installing

The scripts have no imports beyond the standard library and ViennaRNA, so you
can clone and run them in place:

```bash
pip install ViennaRNA          # or: conda install -c bioconda viennarna
python3 design_dragonrna.py --help
```

### If ViennaRNA fails to install

pip builds it from source where no wheel is available, which needs a C
compiler. Conda ships a prebuilt binary and is the easier route on clusters:

```bash
conda install -c bioconda viennarna
```

Then install this package without re-resolving the dependency:

```bash
pip install --no-deps git+https://github.com/haimkru/dragonrna-design
```

Developed against ViennaRNA 2.7.0 and Python 3.12.

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

A candidate is emitted only if all of these hold:

1. No window outside `B` is ≥ `--min-complementarity` nt complementary to `D`.
   Because `D` is the reverse complement of `B`, such a duplex means the window
   holds a 2-mer of `B`, so the rule reduces to a small set of forbidden 2-mers.
2. No duplex anywhere exceeds `--max-duplex` nt, the designed `B:D` stem aside.
3. `X` and `Z` each fold to nothing on their own.
4. The molecule minus `D` folds to nothing.
5. The MFE structure of the whole molecule is exactly the designed `B:D` stem.

`P` is degenerate, so only the emitted draw can be constrained.
`--p-trials` resamples `P` and reports how much of the library breaks rule 1
or 5.

### Choosing the chemistry

`--dna-regions` and `--rna-regions` are inverses; pass one or the other.

```bash
--dna-regions P      # default: P is DNA, X/Z/B/L/D are RNA
--rna-regions P      # P is RNA, X/Z/B/L/D are DNA
--rna-regions ''     # every block DNA
```

These affect how the chimera line is rendered, nothing else.

> **Folding always uses RNA energy parameters.** ViennaRNA is called with its
> RNA model regardless of which blocks you mark as DNA, so for a
> DNA-dominant design the reported MFE is the wrong thermodynamic model. The
> combinatorial constraints (rules 1–2) are pure sequence logic and hold
> either way. G·U wobble is counted as a pair by default and is not a DNA
> pair — consider `--no-wobble` for a mostly-DNA molecule, and re-check the
> final candidate in a tool with DNA parameters.

Run `python3 design_dragonrna.py --help` for the full flag list.

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
