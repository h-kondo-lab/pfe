---
navigation:
  title: リリース
---

# リリース

## 版管理

- 全ての版は concept DOI [10.5281/zenodo.19451887](https://doi.org/10.5281/zenodo.19451887) 配下の Zenodo レコードであり、本リポジトリに同名の git タグを持ちます。タグの `paper1/` には、その Zenodo レコードのファイルそのものに、日本語版と図を加えたものが入っています。
- ファイル名は版の間で変わりません。中身が差し替わります。
- **パッチ**版（v3.2.1、v4.1.1）は正誤のみです。表現・相互参照・表の体裁で、数値は一切変わりません。**マイナー**版は内容を改訂します。**メジャー**版は論文を再構成します。
- 公開済みの版は編集しません。公開後に必要になった訂正は、次のパッチ版として出します。

## 版履歴

| 版 | 日付 | 題名 | Zenodo | 要約 |
| --- | --- | --- | --- | --- |
| v4.1.1 | 2026-09-01 | Physics from Existence I: The Equation | [22217837](https://zenodo.org/records/22217837) | 正誤: 最終査読の所見4件（§1 の条件数、§10 の帰還項込み評価、語の置換、表 2 の $R_{\mathrm{up}}$ を2行に）。数値変更なし。 |
| v4.1.0 | 2026-08-31 | Physics from Existence I: The Equation | [22216256](https://zenodo.org/records/22216256) | $d = 3$ を2条件から導出；$\mathrm{PG}(2,\mathbb F_3)$ を $d$ を入力とせずに構成；絶対スケールを $v$ から $m_e$ へ移し、$v$ が予測量に（246.2360 GeV）；帰還項込みの評価；外部査読への対応。 |
| v4.0.0 | 2026-08-15 | Physics from Existence I: The Equation | [21941031](https://zenodo.org/records/21941031) | 第 5–8 章をパラメータ種別で再編（11 章構成）；節参照の全数監査。数値変更なし。 |
| v3.2.1 | 2026-08-12 | Physics from Existence I: The Equation | [21895908](https://zenodo.org/records/21895908) | 正誤: 相互参照の修正3件。数値変更なし。 |
| v3.2.0 | 2026-08-12 | Physics from Existence I: The Equation | [21895398](https://zenodo.org/records/21895398) | 論文 II との同期（剛性定理・一意性定理）；PDG 2026 更新；$m_W/m_Z$ を同定値へ訂正；表 2 を統合。 |
| v3.1.0 | 2026-08-05 | Physics from Existence I: The Equation | [21808864](https://zenodo.org/records/21808864) | 散文と構成の改訂；単一スケールでの較正。 |
| v3.0.0 | 2026-08-02 | Physics from Existence I: The Equation | [21757785](https://zenodo.org/records/21757785) | 証明を論文 II から移設；表を主張の種類で再編；動機の節を追加。 |
| v2.6.0 | 2026-07-29 | Physics from Existence I: The Equation | [21670655](https://zenodo.org/records/21670655) | 公開後査読への対応；strong-CP の主張を matching 点に限定；$1/\alpha$ の桁を再現可能に。 |
| v2.5.0 | 2026-07-23 | Physics from Existence I: The Equation | [21511757](https://zenodo.org/records/21511757) | 構造対応の枠組み；NuFIT 6.1・PDG 2026 更新；$N = 3$ の区間演算証明。 |
| v2.4.0 | 2026-07-02 | Physics from Existence I: The Equation | [21123893](https://zenodo.org/records/21123893) | $N \ge 3$ の厳密な変分証明；最終査読対応。 |
| v2.3.0 | 2026-06-12 | Physics from Existence I: The Equation | [20666890](https://zenodo.org/records/20666890) | Jarlskog 不変量と Dirac ニュートリノ導出の精緻化。 |
| v2.2.1 | 2026-06-11 | Physics from Existence I: The Equation | [20634925](https://zenodo.org/records/20634925) | 正誤: $\alpha_s$ の不確かさ、Weyl テンソルの証明、ビリアル比の定義。 |
| v2.2.0 | 2026-06-02 | Physics from Existence I: The Equation | [20488132](https://zenodo.org/records/20488132) | 改訂。 |
| v2.1.0 | 2026-05-31 | Physics from Existence I: The Equation | [20471724](https://zenodo.org/records/20471724) | 改訂。 |
| v2.0.0 | 2026-05-31 | Physics from Existence I: The Equation | [20470331](https://zenodo.org/records/20470331) | 大改訂。 |
| v1.2.0 | 2026-04-21 | Physics from Existence I: The Equation | [19678855](https://doi.org/10.5281/zenodo.19678855) | 改訂。 |
| v1.1.0 | 2026-04-13 | Physics from Existence I: The Equation | [19551657](https://doi.org/10.5281/zenodo.19551657) | シリーズの論文 I として改題。 |
| v1.0.2 | 2026-04-10 | Physics from Existence | [19488819](https://doi.org/10.5281/zenodo.19488819) | 改訂。 |
| v1.0.1 | 2026-04-09 | Physics from Existence | [19478719](https://doi.org/10.5281/zenodo.19478719) | 改訂。 |
| v1.0.0 | 2026-04-07 | Physics from Existence | [19451888](https://doi.org/10.5281/zenodo.19451888) | 初版公開。 |

## ライセンス

- 論文・本文・図（`.tex`、`.pdf`、図）: [CC BY 4.0](../../../../LICENSE)——出典表示のもとで自由に共有・翻案できます。
- コード（`.py`）: [MIT](../../../../LICENSE-CODE)——他のプロジェクト内での利用を含め、自由に使用・改変・再配布できます。

全ての内容は、専用ソフトウェアなしに末永く読めるよう、プレーンテキスト形式（LaTeX、Markdown、Python）で保存しています。
