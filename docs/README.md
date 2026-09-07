---
navigation:
  order: 0
---

# Physics from Existence — Documentation

Physics from Existence (PFE) is a research programme by Hidekazu Kondo (Kondo Research Institute, Tokyo). Its starting point is that formalizing what it means "to exist" uniquely yields a single equation, $V = -H$, with no free parameters. The mathematical structure derived from that equation is then shown to correspond to the structures of the Standard Model, and its 26 physical constants are reproduced with no adjustable parameters.

These pages accompany the published paper repository. The papers themselves live in [`paper1/`](../paper1/) as LaTeX, PDF and verification code, and every version is archived on Zenodo. The documentation explains what is in the repository, how to reproduce the numbers, and how the versions relate to each other.

## Contents

| Section | What it covers |
| --- | --- |
| [Papers](01-papers/README.md) | Paper I: abstract, files, how to cite. Paper II: status. |
| [Verification](02-verification/README.md) | Running the derivation and verification scripts, and reading their output. |
| [Results](03-results/README.md) | Comparison with the Standard Model constants, structural results, testable predictions. |
| [Releases](04-releases/README.md) | Version history with Zenodo DOIs, versioning policy, licences. |

## Where to start

- To read the paper: [Papers](01-papers/README.md), then the PDF.
- To check a number yourself: [Verification](02-verification/README.md), then compare with [Results](03-results/README.md).
- To cite: [Papers › Citation](01-papers/README.md#citation).
- Website: <https://www.kondo-lab.com/>

## About this folder

The documentation is plain Markdown in the Lunascape Docs layout: each page carries its own title and order in front matter, and Japanese pages live under [`i18n/ja/`](i18n/ja/README.md). The English pages are the canonical ones; the Japanese pages are their translations. GitHub renders everything as it is, so no build step is needed to read it.
