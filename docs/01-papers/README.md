---
navigation:
  order: 100
---

# Papers

## Paper I — Physics from Existence I: The Equation

*The Binary-Entropy Equation $V = -H$ and the Structure of Standard-Model Parameters*

| | |
| --- | --- |
| Author | Hidekazu Kondo, Kondo Research Institute, Tokyo, Japan |
| Current version | **v4.1.1** (2026-09-01) |
| Zenodo record | <https://zenodo.org/records/22217837> · DOI 10.5281/zenodo.22217837 |
| All versions | DOI [10.5281/zenodo.19451887](https://doi.org/10.5281/zenodo.19451887) |
| Length | 38 pages (English), 45 pages (Japanese) |

### Abstract

No theory has simultaneously derived all 26 free parameters of the Standard Model. In this paper, extending the quantum indefiniteness of existence to the existence of the universe itself, we posit three axioms about existence — existence is bivalent ($n \in \{0,1\}$, i.e. $n^2 = n$), existence is uncertain ($\sigma \neq 1$), and non-existence is forbidden ($\sigma \neq 0$) — and show that they uniquely determine the single equation $V = -H$ as an algebraic identity, and that the mathematical structure derived from its exact algebraic expansion corresponds to the structures of the Standard Model. As a result, the 26 physical constants of the Standard Model are reproduced with no adjustable parameters. Of these, the nineteen derived in this paper agree with the established measured values to 0.0006%–2.6%.

Derived from the axioms of existence, $V = -H$ is the canonical free energy of a single binary degree of freedom $n \in \{0,1\}$, determined as the Legendre transform of the Fermi–Dirac distribution. From this potential and the kinetic term of the canonical normalisation ($g = 1$, §3.1), the Schrödinger equation is derived as a transfer matrix, and we show that it has exactly three bound states (discrete solutions, $N = 3$), corresponding to the three generations of fermions. The spatial dimension $d = 3$ follows from two conditions — preservation of the well and propagating gravity — and the Lie-algebra structure unfolding on the finite geometry constructed from $N = 3$ corresponds to the structure of the gauge group $\mathrm{SU}(3) \times \mathrm{SU}(2) \times \mathrm{U}(1)$ (the coincidence of the two numbers, $N = d = 3$, is an independent consistency).

Twenty-six physical constants are derived from this structure, including $1/\alpha = 137.035999084$ (all orders, matching CODATA 2018 to all displayed digits; the leading term is 137.0368 at 0.0006%), $\alpha_s = 0.1179$ (0.07%), $\sin^2\theta_W = 0.23122$ (<0.01%), the charged-lepton mass ratio $R = 1.5293$ (0.006%), and $\sin\theta_C = 0.2245$ (0.10%). The results further predict the absolute neutrino masses ($\Sigma m_\nu = 59.37$ meV, with no oscillation input) and strictly zero neutrinoless double-beta decay, among other testable predictions. These quantities are not independent input parameters — they are all computed from the same mathematical structure fixed by $V = -H$.

All of these numerical values are reproducible with the accompanying verification code.

### Full text

The complete paper is included here as Markdown, so it can be read in place:

**[Physics from Existence I: The Equation — full text](physics_from_existence.md)**

This file is byte-identical to the manuscript sealed for v4.1.1 (SHA-256
`8caca3a8a87020c9fc6b61a53b6ddcf9876cd9eee8fafdb9d3bdbdcbdf36d73f`), and its
content parity with the published LaTeX source is checked at release. The PDF
remains the typeset version of record, and the Zenodo record is what to cite.

### Files

| File | Content |
| --- | --- |
| [physics_from_existence.md](physics_from_existence.md) | Full text in Markdown (this documentation) |
| [physics_from_existence.pdf](../../paper1/physics_from_existence.pdf) | The paper (English) |
| [physics_from_existence.tex](../../paper1/physics_from_existence.tex) | LaTeX source (English) |
| [physics_from_existence_ja.pdf](../../paper1/physics_from_existence_ja.pdf) | The paper (Japanese) |
| [physics_from_existence_ja.tex](../../paper1/physics_from_existence_ja.tex) | LaTeX source (Japanese) |
| [physics_from_existence_pure.py](../../paper1/physics_from_existence_pure.py) | Pure derivation — no experimental value appears in the code |
| [physics_from_existence_verify.py](../../paper1/physics_from_existence_verify.py) | Verification — the same derivation, compared against measured values |
| `fig_potential.pdf`, `fig_pg23_v2.pdf` | Figures referenced by both LaTeX sources |

File names are fixed across versions; each release replaces their contents and is marked with a git tag (see [Releases](../04-releases/README.md)). The Zenodo record for each version holds the English PDF and TeX and the two scripts.

### Citation

> Kondo, H. (2026). *Physics from Existence I: The Equation* (v4.1.1). Zenodo. https://doi.org/10.5281/zenodo.22217837

To cite the work independently of version, use the concept DOI, which always resolves to the latest version: https://doi.org/10.5281/zenodo.19451887

```bibtex
@article{kondo2026pfe1,
  author    = {Kondo, Hidekazu},
  title     = {Physics from Existence I: The Equation},
  year      = {2026},
  publisher = {Zenodo},
  version   = {v4.1.1},
  doi       = {10.5281/zenodo.22217837},
  url       = {https://doi.org/10.5281/zenodo.22217837}
}
```

## Paper II — Physics from Existence II: The Proof

In preparation. Several quantities shown in Paper I's tables are derived in Paper II; those rows are marked "II" in [Results](../03-results/README.md) and appear in Paper I with that attribution.
