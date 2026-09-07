---
navigation:
  title: 論文
---

# 論文

## 論文 I — 存在からの物理学 I: The Equation

*二値エントロピー方程式 $V = -H$ と標準模型パラメータの構造*

| | |
| --- | --- |
| 著者 | 近藤秀和, Kondo Research Institute, 東京 |
| 現行版 | **v4.1.1**（2026-09-01） |
| Zenodo レコード | <https://zenodo.org/records/22217837> · DOI 10.5281/zenodo.22217837 |
| 全版共通 | DOI [10.5281/zenodo.19451887](https://doi.org/10.5281/zenodo.19451887) |
| 分量 | 英語版 38 頁、日本語版 45 頁 |

### 概要

標準模型の26個の自由パラメータを同時に導出した理論は存在しない。本論文では、量子的な存在の不確定性を宇宙の存在そのものへ拡張し、存在に関する3つの公理——存在は二値的である（$n \in \{0,1\}$、すなわち $n^2 = n$）、存在は不確実である（$\sigma \neq 1$）、非存在は禁止される（$\sigma \neq 0$）——を置くと、そこから単一の方程式 $V = -H$ が代数的恒等式として定まることを示すとともに、この方程式の厳密な代数的展開から導出される数学的構造が、標準模型の構造と対応することを示す。その結果として、標準模型の26個の物理定数が、調整可能なパラメータなしに再現される。うち本論文で導出する19個は、確定した実測値と 0.0006%–2.6% の精度で一致する。

存在の公理から導出される $V = -H$ は二値自由度 $n \in \{0,1\}$ の正準自由エネルギーであり、フェルミ・ディラック分布の Legendre 変換により定まる。このポテンシャルと正準規格化（$g = 1$、§3.1）の運動項からシュレーディンガー方程式が転送行列として導出され、この方程式がフェルミオンの3世代に対応するちょうど3つの束縛状態（離散解）を持つ（$N = 3$）ことを示す。空間次元 $d = 3$ は井戸の保存と重力伝播の2条件から導かれ、$N = 3$ から構成される有限幾何の上に展開されるリー代数の構造が、ゲージ群 $\mathrm{SU}(3) \times \mathrm{SU}(2) \times \mathrm{U}(1)$ の構造と対応する（2数の一致 $N = d = 3$ は独立な整合性である）。

この構造から、$1/\alpha = 137.035999084$（全次数、CODATA 2018 と表示桁一致；先頭項は 137.0368 で 0.0006%）、$\alpha_s = 0.1179$（0.07%）、$\sin^2\theta_W = 0.23122$（<0.01%）、荷電レプトン質量比 $R = 1.5293$（0.006%）、$\sin\theta_C = 0.2245$（0.10%）を含む26個の物理定数が導かれる。さらにニュートリノの絶対質量（振動データを用いずに $\Sigma m_\nu = 59.37$ meV）、ニュートリノなし二重ベータ崩壊の厳密なゼロを含む検証可能な予測を与える。これらの量は独立な入力パラメータではなく、すべて $V = -H$ により定まる同一の数学的構造から計算される。

これらの数値は全て、付属の検証コードによって再現可能である。

### 全文

論文の全文を Markdown で同梱しているので、そのまま読めます（英語）:

**[Physics from Existence I: The Equation — 全文](../../../01-papers/physics_from_existence.md)**

このファイルは v4.1.1 で封印された原稿とバイト単位で同一です（SHA-256
`8caca3a8a87020c9fc6b61a53b6ddcf9876cd9eee8fafdb9d3bdbdcbdf36d73f`）。
公開 LaTeX ソースとの内容パリティはリリース時に検査しています。組版の正本は
PDF であり、引用は Zenodo レコードに対して行ってください。

### ファイル

| ファイル | 内容 |
| --- | --- |
| [physics_from_existence.md](../../../01-papers/physics_from_existence.md) | 全文 Markdown（本ドキュメント、英語） |
| [physics_from_existence.pdf](../../../../paper1/physics_from_existence.pdf) | 論文（英語） |
| [physics_from_existence.tex](../../../../paper1/physics_from_existence.tex) | LaTeX ソース（英語） |
| [physics_from_existence_ja.pdf](../../../../paper1/physics_from_existence_ja.pdf) | 論文（日本語） |
| [physics_from_existence_ja.tex](../../../../paper1/physics_from_existence_ja.tex) | LaTeX ソース（日本語） |
| [physics_from_existence_pure.py](../../../../paper1/physics_from_existence_pure.py) | 純粋導出——コード中に実験値は一切現れない |
| [physics_from_existence_verify.py](../../../../paper1/physics_from_existence_verify.py) | 検証——同じ導出を実測値と比較する |
| `fig_potential.pdf`, `fig_pg23_v2.pdf` | 両 LaTeX ソースが参照する図 |

ファイル名は版によらず固定で、各リリースで中身が差し替えられ、git タグが打たれます（[リリース](../04-releases/README.md) 参照）。各版の Zenodo レコードには英語版 PDF・TeX と2本のスクリプトが収められています。

### 引用

> Kondo, H. (2026). *Physics from Existence I: The Equation* (v4.1.1). Zenodo. https://doi.org/10.5281/zenodo.22217837

版によらず引用するには、常に最新版へ解決される concept DOI を使ってください: https://doi.org/10.5281/zenodo.19451887

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

## 論文 II — 存在からの物理学 II: The Proof

準備中。論文 I の表に示される量のいくつかは論文 II で導出されるもので、[結果](../03-results/README.md) では「II」と印を付け、論文 I でもその帰属で示されています。
