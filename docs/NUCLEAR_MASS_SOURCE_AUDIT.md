# pB / DD / DT mass and Q audit

审计范围：只整理可追溯的质量与反应 Q 证据，不选择产物模型或实现约定。数值中明确区分裸核质量与中性原子质量；由质量表相减得到的项目标为 `derived`，不把互不相容的极值合并成一个“推荐值”。

## 现有 C++ 状态

- `/home/cloud/research/pB11_reaction/src/fusion_constants.hpp` 已有 NIST CODATA 2022 的 `p,d,t,helion,alpha,n` 质量能（MeV）和四个两体 Q 的裸核质量组合；`q_MeV[0]` 明写为四舍五入的 `8.68`，没有 B-11 质量或电子束缚能。
- `/home/cloud/research/pB11_reaction/docs/THERMAL_NETWORK.md` 同样说明四个两体 Q 使用 NIST 裸核质量，pB 基线保留 `8.68 MeV`。这只是现有实现状态，不是本审计的独立数据源。

## 裸核质量：NIST 2022 CODATA

来源为 [NIST 2022 CODATA complete listing](https://physics.nist.gov/cuu/Constants/Table/allascii.txt)（文件首行明确为 “2022 CODATA adjustment”）；印刷版 [JPCRD2022CODATA.pdf](https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf) 的 Table XXXIII，质子 p. 46，n/d/t pp. 47–48，helion/alpha p. 48，电子 p. 45。表中“mass energy equivalent”数值按 MeV/c² 记，括号为 1σ 的末位不确定度。

| 粒子/裸核 | 质量能 (mc^2) (MeV) | 1σ (MeV) | NIST Table XXXIII |
|---|---:|---:|---|
| p (proton) | 938.27208943 | 0.00000029 | p. 46 |
| n (neutron，DD 产物所需) | 939.56542194 | 0.00000048 | p. 47 |
| D (deuteron) | 1875.61294500 | 0.00000058 | p. 47 |
| T (triton) | 2808.92113668 | 0.00000088 | p. 47 |
| ³He nucleus (helion) | 2808.39161112 | 0.00000088 | p. 48 |
| ⁴He nucleus (alpha) | 3727.3794118 | 0.0000012 | p. 48 |
| (e^-) | 0.51099895069 | 0.00000000016 | p. 45 |

相应的 NIST 质量单位值为：p `1.0072764665789(83) u`、n `1.00866491606(40) u`、D `2.013553212544(15) u`、T `3.01550071597(10) u`、helion `3.014932246932(74) u`、alpha `4.001506179129(62) u`；电子 `5.485799090441(97)×10^-4 u`。这里的 p/D/T/helion/alpha 是粒子（裸核）质量，不是带电子的中性 H/D/T/He 原子质量。

## AME2020 原子质量输入

[AME2020 官方页面](https://amdc.impcas.ac.cn/web/masseval.html) 说明 `mass_1.mas20` 是原子质量文件，且电子文件的正式引用应回到两篇原论文：[Part I, DOI 10.1088/1674-1137/abddb0](https://doi.org/10.1088/1674-1137/abddb0) 和 [Part II, DOI 10.1088/1674-1137/abddaf](https://doi.org/10.1088/1674-1137/abddaf)（[Part II PDF](https://www-nds.iaea.org/amdc/ame2020/AME2020-b.pdf)）。下表来自 [mass_1.mas20](https://amdc.impcas.ac.cn/masstables/Ame2020/mass_1.mas20) 的分析列；质量列单位是 micro-u（所以 0.013 micro-u = 1.3×10^-8 u），质量过剩列单位是 keV。

| 中性原子/核素 | atomic mass (u), 1σ | mass excess (keV), 1σ | 证据位置 |
|---|---:|---:|---|
| ¹H（用于原子质量 Q 的 p） | 1.007825031898 ± 0.000000000014 | 7288.971064 ± 0.000013 | `mass_1.mas20` row 1H |
| ²H (D) | 2.014101777844 ± 0.000000000015 | 13135.722895 ± 0.000015 | row 2H |
| ³H (T) | 3.01604928132 ± 0.00000000008 | 14949.81090 ± 0.00008 | row 3H |
| ³He | 3.01602932197 ± 0.00000000006 | 14931.21888 ± 0.00006 | row 3He |
| ⁴He | 4.00260325413 ± 0.00000000016 | 2424.91587 ± 0.00015 | row 4He |
| n | 1.00866491590 ± 0.00000000047 | 8071.31806 ± 0.00044 | row n |
| ¹¹B | **11.009305166 ± 0.000000013** | **8667.708 ± 0.012** | row 11B, plain-text line 72 |

因此 AME2020 对 ¹¹B 在此处给的是**中性原子质量**，不是裸核质量；裸核数值不能直接把该行当作核质量使用。

## ¹¹B 的五电子与总电子束缚能

AME Part II 的印刷版 `030003-1`（Eq. (1)）给出

\[
 M_N(A,Z)=M_A(A,Z)-Z m_e+B_e(Z),
\]

其中在以 u 表示质量时最后一项应按 B_e/c² 换算。Part II `030003-2`（Eq. (2)）还给出仅作近似的

\[
 B_e(Z)=14.4381 Z^{2.39}+1.55468\times10^{-6}Z^{5.35}\ {\rm eV}.
\]

AME 明确说 B_e 的计算精度没有像轻核质量那样被独立确立；该经验式和下述 NIST ASD 逐级电离能是两条来源路径，不能平均或拼接。

[NIST CODATA 2022](https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf) Sec. II.A（印刷 pp. 4–5，Eq. (1)–(2)）给出更直接的定义：中性原子质量等于核质量加电子质量再减电子束缚能；裸核对应 n=Z，且

\[
\Delta E_B(X^{n+})=\sum_{i=0}^{n-1}E_I(X^{i+}).
\]

所以 ¹¹B 需要五个电子质量 5m_e 和从 B I 到 B V 的五个逐级电离能。NIST Table XXXIII 的电子值给出

\[
5m_e=2.55499475345\ {\rm MeV}/c^2
\]

（由 0.51099895069(16) MeV/c² 乘 5；约 0.0027428995452205 u，乘法不确定度约 0.00000000080 MeV/c²）。

NIST [Atomic Spectra Database ionization output for B, all spectra](https://physics.nist.gov/cgi-bin/ASD/ie.pl?at_num_out=on&biblio=on&e_out=0&el_name_out=on&format=0&ion_charge_out=on&ion_conf_out=on&level_out=on&order=0&seq_out=on&shells_out=on&sp_name_out=on&spectra=B&submit=Retrieve+Data&unc_out=on&units=1)（ASD ver. 5.12，DOI [10.18434/T4W30F](https://doi.org/10.18434/T4W30F)）实际列出：

| 级次 | (E_I) (eV) | 1σ (eV) | ASD reference |
|---|---:|---:|---|
| B I | 8.298019 | 0.000003 | L12312 |
| B II | 25.15483 | 0.00005 | L12120 |
| B III | 37.93059 | 0.00007 | L12547 |
| B IV | (259.374379) | 0.000009 | L10054 |
| B V | (340.2260225) | 0.0000006 | L19200 |

五行的**算术和**是 B_e(B) = 670.9838405 eV（约 7.20330744×10^-7 u）。ASD 页面没有为该和另行给出协方差或合成 1σ；B IV/B V 数值在页面中带括号。因此下面的裸 ¹¹B 数值只是“采用上述 ASD 行”的条件推导，不把它冒充 AME 单独发布的裸核质量：

\[
\begin{aligned}
M_N(^{11}{\rm B})
&=11.009305166\ {\rm u}
 -5(5.485799090441\times10^{-4}\ {\rm u})\\
&\quad +670.9838405\ {\rm eV}/(uc^2)\\
&\simeq 11.0065629868\ {\rm u}.
\end{aligned}
\]

该条件计算的误差至少应保留 AME ¹¹B 原子质量的 ±0.000000013 u，并另外说明 ASD 电子项不是 ^11B 同位素专门的独立总束缚能。若不采用 ASD 五行，而采用 AME Eq. (2)，应把它标为另一套近似约定，不能与上式混用。

## 反应 Q：原始/derived 与裸核一致性

AME Part II Table III 的解释在印刷 `030003-102`（pp. 102–103），并明确这些表项由**原子质量**组合而来；`a` 表示该表项不确定度小于 5 eV，完整值在 ASCII 文件中。未四舍五入的 [rct1.mas20](https://amdc.impcas.ac.cn/masstables/Ame2020/rct1.mas20) 和 [rct2_1.mas20](https://amdc.impcas.ac.cn/masstables/Ame2020/rct2_1.mas20) 是同一 AME2020 评价的机器可读输出，正式论文页码仍以 Table III 为准。

| 反应 | AME 原子质量 Q（来源状态） | 裸核质量 Q（由 NIST 裸核质量；pB 另有条件项） | 证据/缺口 |
|---|---:|---:|---|
| p+¹¹B→3α | `Q(p,α)`(¹¹B) = 8590.0917 ± 0.0374 keV（rct2 line 73 / Table III p. 104）再加 ⁸Be `Q(α)` = 91.8399 ± 0.0354 keV（rct1 line 59 / Table III p. 105），算术和 **8681.9316 keV ≈ 8.68193 MeV** | **≈ 8.68238 MeV**，条件是采用上面 NIST ASD 五行得到的 B11 裸核质量 | AME Table III 没有单独的三 α 一行；总值是两条 AME 表项相加。裸核值未作为 AME 单行发布，且未在此重算完整协方差。 |
| D+D→T+p | **4.0326638 MeV**, `derived` from AME `mass_1.mas20` mass excess 2ME(D)−ME(T)−ME(¹H)；没有在 Table III 找到该 DD 分支单行 | **4.03266389 MeV**, 2m_D−m_T−m_p | AME 原子和 NIST 裸核结果相差小于 0.1 eV 量级；未给 derived 值另附协方差 σ。 |
| D+D→³He+n | **3.2689089 MeV**, `derived` from 2ME(D)−ME(³He)−ME(n)；没有在 Table III 找到该 DD 分支单行 | **3.26885694 MeV**, 2m_D−m_h−m_n | 约 52 eV 差异来自中性原子电子壳层束缚差；这是原子 Q 与裸核 Q 的定义差异，不应取平均。 |
| D+T→⁴He+n | **17589.2999 ± 0.0005 keV**（rct2 line 35；Table III p. 104 为 17589.30，`a`） | **17.58924794 MeV**, m_D+m_T−m_α−m_n | AME 项是原子质量 Q；NIST 项是裸核 Q，两者差约 **0.052 keV = 52 eV**。 |
| D+³He→⁴He+p | **18353.0548 ± 0.0002 keV**（rct2 line 36；Table III p. 104 为 18353.05，`a`） | **18.35305489 MeV**, m_D+m_h−m_α−m_p | 该反应两侧电子数和壳层组合基本相同，原子/裸核 Q 在列示精度内相符；NIST 项未在此传播相关不确定度。 |

上表的两个 DD 数值由 AME 质量过剩或 NIST 裸核质量分别闭合；它们不是把两条来源的极值合成一个数。AME Table III 说明 Q 误差应使用调整质量及相关矩阵；本审计没有把缺失的 DD 单行误差或 pB 三 α 总和误差自行补出来。

## 直接交付给主代理的边界

1. **8.68 MeV 是现有 pB 基线的明确四舍五入常数**。AME 原子质量组合给约 8.68193 MeV；在采用 NIST ASD B 电子项并转成裸核后约 8.68238 MeV。三者是不同精度/质量约定，不能在能量账本中静默互换。
2. B11 的 AME 数据只有中性 atomic mass；五电子质量来源是 NIST CODATA Table XXXIII，电子总束缚能来源可追溯到 NIST CODATA Eq. (1)–(2) / NIST ASD B I–V 五行，AME 还提供了无独立 σ 的多项式近似。没有一个 AME 表行直接发布“^11B 裸核质量”。
3. DT、DHe3 的 AME Q 是原子质量表项；当前 C++ 四两体 Q 是 NIST 裸核质量组合。若主代理要保持裸核能量守恒，应采用同一列裸核质量并单独处理 B11 转换；若采用 AME 原子 Q，则需把它标成原子质量定义。
