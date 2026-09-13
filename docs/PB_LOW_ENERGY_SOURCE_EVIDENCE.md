# p–¹¹B 低能截面来源核查（NS/TB）

核查日期：2026-09-14。范围限于 `E_cm < 140 keV` 的 S 因子/总截面证据、Nevins–Swain（NS）来源链及 Tentori–Belloni 2023（TB23）的继承关系；不做低能外推或产物模型选择。

## 先纠正 NS 的书目信息

没有找到与“NS 1995, *Nuclear Fusion* **35**, 1813”相符的 p–¹¹B 论文。可核实的 NS 原始论文是：W. M. Nevins and R. Swain, “The thermonuclear fusion rate coefficient for p-¹¹B reactions,” *Nuclear Fusion* **40**, 865–872 (2000), DOI [10.1088/0029-5515/40/4/310](https://doi.org/10.1088/0029-5515/40/4/310)。TB23 的摘要和参考文献（pp. 1, 9）也明确写为 NS (2000), *Nucl. Fusion* 40, 865。NS 原文 PDF 本轮没有取得，因此 NS 的原始页码、原始数据表和参数协方差不能独立审计；下述 NS 参数值来自 TB23 对 NS 的逐项转录，并明确标为二手转录。

## NS 实际使用的低能数据依据

TB23 p. 1 明确写道：NS 计算反应率时使用 Becker 等的 astrophysical-factor 数据 `E < 1.1 MeV`，以及 Segel 等的 `1.1 < E < 3.5 MeV` 数据，`E` 均为质心能量。因而 NS 的低能输入不是 Sikora–Weller 的起点。

被 TB23 引用为 Becker 数据源的原始论文是 H. W. Becker, C. Rolfs and H. P. Trautvetter, “Low-energy cross sections for ¹¹B(p,3α),” *Z. Phys. A* **327**, 341–355 (1987), DOI [10.1007/BF01284459](https://doi.org/10.1007/BF01284459)。该文摘要/索引给出：

- `E_cm = 22–1100 keV`；测量了绝对截面、α 角分布和激发函数，并用经验拟合描述 `S(E)`。
- 经验拟合的零能截距为 **`S(0) = 197 ± 12 MeV b`**。这是拟合截距，不是 `E=0` 的直接测量；±12 是目前核实到的低能 S 因子不确定度。
- 22 keV 明确低于 140 keV，因此历史原始数据确实覆盖目标区间。该摘要没有给出可逐点转录的 `<140 keV` 总截面表；本轮也没有取得全文表格、机器可读文件或逐点误差。

这里的 `S(0)` 与 NS 低段的 `C₀` 数值相同，但不能自动视为含窄共振项的完整 `S₁(0)`。Becker 的正文拟合细节及 NS 是否逐点采用了 22 keV 端点，因全文缺失，保留为未核实项。

## NS 低段参数与 TB23 的继承

TB23 Eq. (3), p. 2 给出与 NS 一致的低段形式（`E`、`E_L`、`δE_L` 以 keV 代入）：

```text
S1(E) = C0 + C1 (E/1 keV) + C2 (E/1 keV)^2
        + AL / [((E−EL)/1 keV)^2 + (δEL/1 keV)^2] .
```

TB23 Table 1, p. 4 对 NS 的参数转录为：

| 项 | NS 值 | 证据边界 |
|---|---:|---|
| 低段范围 | `E ≤ 0.400 MeV` | TB23 Table 1；NS 原文页未核实 |
| `C0` | `197 MeV b` | 多项式常数；不要单独称为完整 `S1(0)` |
| `C1` | `0.240 MeV b` | TB23 Table 1 二手转录 |
| `C2` | `2.31×10⁻⁴ MeV b` | TB23 Table 1 二手转录 |
| `AL` | `1.82×10⁴ MeV b` | 148 keV 窄峰幅度，TB23 说明直接取自 NS |
| `EL` | `148 keV` | TB23 说明直接取自 NS |
| `δEL` | `2.35 keV` | TB23 说明直接取自 NS |

TB23 p. 2 明确说上述 Breit–Wigner 三项直接取自 NS，因为 SW/Spraker 数据的能量分辨率不足以包含 148 keV 窄峰；TB23 p. 4 又说明其新反应率的窄共振项仍使用这三个 NS 值，而非重新拟合。TB23 没有给这些 NS 参数附协方差或误差条；除 Becker 的 `S(0)=197±12 MeV b` 外，本轮没有核实到 NS 低能逐点不确定度。

TB23 p. 2 的“5% experimental error”是作者对 SW/BU 拟合采用的保守假设（SW 只报告统计误差，BU 数据点没有误差条），不是 NS 或 Becker 的 `<140 keV` 测量不确定度，不能移作低能数据的误差证明。

TB23 当前低段使用同一结构，但新的非共振系数是 `C0=197`、`C1=0.269`、`C2=2.54×10⁻⁴ MeV b`；这是 TB23 的新拟合值，不应回写成 NS 数据。TB23 的 SW 输入范围是 `E_p^lab=0.15–3.80 MeV`，即约 `E_cm=0.14–3.48 MeV`（TB23 p. 2）；其当前低端因此约为 140 keV，不能作为历史实验的最低能量。

## 交接结论与缺口

- `<140 keV` 的明确历史实测覆盖来自 Becker 1987，最低报告能量为 `E_cm=22 keV`；已核实的低能拟合量只有 `S(0)=197±12 MeV b`。
- NS 的低段系数及 148 keV Breit–Wigner 参数可由 TB23 Table 1/Eq. (3) 追溯，但原始 NS 2000 PDF、原始表格、参数误差/协方差和 `<140 keV` 逐点总截面尚未取得。
- SW/TB23 数据约从 `E_cm≈140 keV` 起，且 SW 漏掉 148 keV 窄峰；不能把这个起点当作所有历史 p–¹¹B 数据的起点，也不能据此声称几 keV 温度 Gamow 区已有实测覆盖。

来源：TB23 原文 [DOI](https://doi.org/10.1088/1741-4326/acda4b)，可读副本 [ResearchGate](https://www.researchgate.net/publication/371625643_Revisiting_p-_11_B_fusion_cross_section_and_reactivity_and_their_analytic_approximations)；NS 原始条目 [DOI](https://doi.org/10.1088/0029-5515/40/4/310)；Becker 原始条目 [DOI](https://doi.org/10.1007/BF01284459)。
