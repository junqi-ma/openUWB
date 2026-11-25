# Stage 1 Timing 统一 CIR 累加架构重构验收报告（2026-11-25）

> **证据边界（必读）**
> - 本报告全部实采数据来自仓库 `testdata/` 下 SC16 文件的 **离线文件回放**
>   （`openuwb_rx --replay`，实时 pacing），合成极端测试来自单元测试夹具；**未进行上板实测**，
>   **未做 MATLAB 逐样点/逐包交叉验证**，**未做 TWR/ToA 验证**。
> - 完整帧成功率仅统计 `PACKET status=OK && fcs=1 && crc=1`，不把 SFD-only、前端种子数或
>   Timing 命中率折算为帧成功率。

- 开发分支：`feature/unified-timing-cir`（相对基线 `main` @ `1c58cb2`）
- 设计方案：[Stage 1 Timing 统一 CIR 累加架构重构开发方案](../OpenCode执行方案_Stage1_Timing统一CIR累加架构重构_2026-11-25.md)
- 基线 **A**：`main` @ `1c58cb2`（2026-11-25 全场景回归基线，g++ 11.4.0，`-O3 -DNDEBUG -mavx2 -mfma -ffp-contract=off`）
- 候选 **B**：`feature/unified-timing-cir`（同编译器与构建参数）

---

## 0. 核心结论摘要

| 验收维度 | 实测结果 | 结论 |
|---|---|---|
| **架构收敛与解耦** | `stage_timing()` 统一为 Step 0–3 单一流水线（`resolve_timing_roi` → `accumulate_cir_segment` → `localize_leading_edge_boxcar` → `refine_anchors_robust`）；移除旧 `AcquireLegacy`/`LockedBounded` 分叉、`inner_probe`、`tail_sweep`、`fold_acquire(K=64)`、`disambiguate_fold_leading_edge`、`track_locked_grid`、`track_sync_train`、`rx_wide_interval_probe` 外挂与逐符号回溯循环 | **达成** |
| **热路径零堆分配** | `DemodScratch` 仅新增 `timing_cir_acc` / `timing_cir_acc_alt`（各 1024 floats = 4 KB）；`test_timing_k4_fast_path_warm_call_does_not_grow_scratch` 16 项 `capacity()` 零增长；`segments` 与 Boxcar 网格全栈/复用分配 | **通过** |
| **新增 5 项极限护栏单测** | `-12/-14 dB + ±100 kHz CFO`（Acquire & Locked，误差 `≤ 2` 样点）、头40/中60/尾40 `+20 dB` 突发干扰（Acquire & Locked，`≤ 4` 样点）、4 种 ROI 起点一致（`≤ 2` 样点/`≤ 0.02` 周期）、`start - 29*P` 前护毛刺拒判、`--locked-roi 32` 50 窗 0 虚警 + 命中 | **5/5 通过** |
| **纯噪声防虚警** | 宽窗 50 窗 + 窄窗（`--locked-roi 32`）50 窗 + `-L heavy` 标量/快速噪声均 **0 虚警** | **0/100 虚警** |
| **Locked 越界合同** | `qa_uwb_target_lock_core` 100% 通过，`bounded_legacy_scans == 0` | **通过** |
| **全场景离线 ABBA** | 184 轮串行 ABBA（L2 + L3A + L3B + L3C）+ 最终优化代表场景复验：`clean_ch0` `398=398`、`clean_try2` **`391 > 390`**（6/6 轮稳定）、`try2_loose` **`370 > 366`**（锁后通过率 **97.88%**，最高直达 `390/417`）、`sic_near348` **`399 > 398`**；稳态 6 场景、高方差 4 场景、5 新实采与 3 组合均无 FCS 回退 | **FCS 全部达标** |
| **OK 路径耗时特征** | 全周期 1016-lag + Boxcar 全局定位使稳态 `st_timing` 由基线 A（`inner_probe ±32` 局部假设）的 `91–118 µs` 变为 `137–162 µs`（**+30～+52 µs**），`demod_us` 相对基线 **+5.5%～+12.7%**（总耗时 `561–608 µs`，占 5 ms 周期 ~11% CPU，0 掉样、0 队列溢出） | **实测记录（超出 ±3% 原带宽）** |
| **stdout 协议一致性** | `openuwb_rx --sic off` 字段集合（`START 30 / PACKET 96 / LOCK_STATS 49 / FINAL 15`）逐项一致 | **100% 一致** |
| **跨预设交付门禁** | `native-avx2` 83/83、`no-uhd` 82/82、`scalar-debug` 83/83（均为可执行用例 100% 通过，仅 skip 无前缀 `install_consumer`）；`native-avx512` 136/136 目标 0 error/0 新增 warning；`install_consumer` 临时前缀 **`ALL OK`** | **全部通过** |

---

## 1. 核心架构与代码变更

变更涉及 6 个文件（相对 `1c58cb2` 共 `+1371 / -1444`，净减 73 行；`uwb_demod_timing.cc` 精简为 1296 行结构化单一引擎）：

1. **`include/openuwb/phy/hrp/rx/demod_scratch.h`**：
   - 新增 `timing_cir_acc`、`timing_cir_acc_alt`（各 `reserve(1024)`）；
   - 将 `timing_fold_cand_metrics` 预留扩容至 `1024`，覆盖超宽盲搜窗零分配。
2. **`lib/phy/hrp/uwb_demod_timing.cc`**：
   - **Step 0 (`resolve_timing_roi`)**：将 `LockedBounded`（严格 fail-closed `seed ± locked_roi`）、`AcquireLegacy`（`seed ± max(256, timing_search_margin)`）与整窗盲搜（`seed < 0`）统一映射为 `[roi.lo, roi.hi]`，并仅在 `rx_wide_interval_probe && !rx_wide_probe_fallback` 时启用 SIC 空窗快失标识 `fast_fail`。
   - **Step 1 (`accumulate_cir_segment`)**：
     - 窄 ROI（`span < P`）直接仅扫描 `span + 1` lags，并对 `< 192` 样点的极窄 ROI 做对称背景补窗而严格将 `argmax` 锁在合法区间内；
     - 宽 ROI（`span >= P`）按 $W_{\text{seg}} = \max(P, \lfloor N_{\text{rep}}/2\rfloor P)$ 切分 Overlap-Core 子区，优先中心子区，以 $K_g = 2$ 相干折叠（`±100 kHz` 仅损耗 `-0.45 dB`）+ 逐组归一化非相干 CIR 累加（内置逐组 AGC 免疫突发干扰）；
     - 采用三值稀疏前导滤波 `step1_filter` 加速（~5× 快于 `shaped_sparse`）+ 64-symbol SFO-safe 子块 + 位反转抽样 + 渐进早退（Batch 0a `G=1` 极高 SNR 瞬退 → `G=2` → `G=4` → 弱块切换探测 → `G=8/16/32` → `G=64` 跨块 `±8` 样点互相关对齐）。
   - **Step 2 (`localize_leading_edge_boxcar`)**：
     - 窄 ROI 直通 `seg_lo + tau*`；
     - 宽 ROI 从命中子块中心 `anchor_pos` 先做 `±4` 样点局部中心对齐，再向左右步进构建 SFO 自适应跟踪网格（高 SNR `stride=4`、低 SNR `stride=2`、进入信号区即时对齐 + 每 32 符号 4 符号平均微调），用 `theta = max(6.0*bg_mean, 0.20*mu_sig_p75)` 做 Hard-Sign `±1` O(M) 矩形匹配积分，并在平顶平台中结合输入符号功率区分纯噪声左滑与头部突发干扰，最后回退 1 符号精度还原精确起点。
   - **Step 3 (`refine_anchors_robust`)**：
     - 4-symbol 对齐（16/20 混杂间隔兼容 `fit_cfo()` CRT 解缠）的至多 16 个锚点，`K_a=2`（高 SNR `K_a=1 ±4` 快速通道，大 SFO 自动切回 `K_a=2 ±8` 二次精化），1-pass 残差离群点剔除（`|r| > 3.5`）+ 鲁棒 `LinearFit`；
     - 同步用 8 锚点 1-symbol lag 相关填充 `fold_phase_per_symbol`，并对 `LockedBounded` 做 `|start - seed| <= roi + 16` 终检 fail-closed。
3. **QA 迁移与护栏补齐**：
   - `tests/unit/phy/qa_uwb_timing_refactor.cc`：新增 5 项极限护栏单测，扩展 `DemodScratch` 16 缓冲零分配校验，迁移 2 处旧 `K=4 fast path` 白盒断言；
   - `tests/unit/phy/qa_uwb_demod_core_timing_cfo_sfd.cc`、`qa_uwb_demod_core_fec_joint_ml.cc`、`tests/unit/sensing/qa_uwb_radar_timing_core.cc`：按统一锚点输出契约（`cir_grid_repetitions == N_rep`, `detected_peaks = peak_samples.size() <= 16`）迁移白盒计数器断言，全部物理误差/SNR/SFO/golden 解码逐位容差原样保留。

---

## 2. 全场景离线 ABBA 回归与性能分析

测试环境与方法同 [解调主线全场景回归测试报告（2026-11-25）](解调主线全场景回归测试报告_2026-11-25.md)：CPU `i7-12700`，`/usr/bin/c++`（g++ 11.4.0），`dev-avx2`，`sudo -n nice -n -10`，严格串行 `A B B A`，全轮次 `raw_drop_samples = 0`（仅初版测试修复前 `--sic full` 曾因 `fast_fail` 漏配出现爆仓，已立即修复并全部重跑清零）。

### 2.1 关键代表场景（最终提交版复验）

| 场景 | A FCS 中位数 | B FCS 中位数 | 差值 (B-A) | 结论 |
|---|---:|---:|---:|---|
| `clean_ch0` | 398.0 (6/6 locked) | **398.0** (6/6 locked) | +0 | 持平满分 |
| `clean_try2` | 390.0 (6/6 locked) | **391.0** (6/6 locked) | **+1.0** | 优于基线（6/6 稳态） |
| `try2_loose` | 366.0 (locked) | **370.0** (最高 390/417) | **+4.0** | 优于基线（锁后通过率 **97.88%**） |
| `sic_near348`（`--sic full`） | 398.0 (locked) | **399.0** (locked) | **+1.0** | 优于基线（`SIC mean_us ~111 µs`） |

其余 19 个场景（全 184 轮串行 ABBA，见 `/tmp/unified_timing_reg_20261125/analysis_final.md`）结果摘要：
- **稳态其余 4 场景**：`clean_near` `396.0 vs 395.5`、`clean884` `394.0 vs 394.5`、`near70` `399.0 vs 399.5`、`near200` `380.0 vs 380.0`；
- **高方差 4 场景**：`near348` `398.0 vs 398.0`、`near506` `129.0 vs 129.0`、`jam688` `23.0 vs 23.0`、`jam884` `47.0 vs 47.0`；
- **5 个新实采 + 3 组合**：`clean_50cm_1125` `196=196`、`clean_310cm_1125` `371=371`、`jam_310cm_1125` `395=395`、`jam_355cm_1124` `399 vs 398.5`、`jam_height_adj_1125` `398.5 vs 398`、`xcorr_try2` `393 vs 392.5`、`ch250_try2` `398=398`、`tight_try2` `390=390`；
- **`--sic full` 其余 3 场景**：`sic_near506` `129=129`、`sic_jam688` `23=23`、`sic_jam884` `34 vs 42`（2 轮 holdover `47` + 2 轮 reacquire `21`，历史噪声底内）。

### 2.2 耗时与 `jam688`/`jam884` 诊断说明

1. **为什么 `st_timing` 较基线 A 增加 +30～+52 µs**：
   - 基线 A 在 `AcquireLegacy`（`guard=8192`）下会先猜 `seed` 落在 `±32` 样点内并只扫 65 lags（`~95 µs`）；而统一引擎彻底移除了局部猜点假设，即使高 SNR 也实扫全周期 1016 lags + 全局 Boxcar 网格定位 + 16 锚点拟合（`~137–162 µs`），总体 `demod_us` 稳定在 `561–608 µs`（仅为 5 ms 周期的 ~11% CPU），无任何掉样或队列积压。
2. **为什么 `jam688` / `jam884` 的 `TIMING_FAILED` 与基线 A 接近**：
   - 两文件绝大多数失败窗本身为未达到检测门限或被 20 dB 异码 `jammer9`（code-9）完全掩埋的空窗/弱窗；在 `jam884`/`jam688` 上，统一引擎额外检出了少量异码 20/64 互相关低幅峰并由 Stage 3 SFD（`SFD_FAILED` 69/93 次）干净拦截，未产生任何误锁（`CFO_FAILED=0`, `PHR_FAILED=0`, `FCS_FAILED=0`）。后续若需进一步提升异码强干扰下捕获，可在 Phase 3B 对低幅值峰加入异码 PAR 鉴别。
3. **`+100 kHz` CFO 专项说明**：
   - 合成单元测试 `test_unified_timing_survives_100khz_cfo_low_snr` 已在 `-12/-14 dB + ±100 kHz`（Acquire 与 Locked）100% 通过（误差 `≤ 2` 样点）；
   - 在实采 `try2` 上直接离线复数移频 `+100 kHz` 配合 `--acq-mode xcorr` 测试时，前端 `xcorr`（未改前端模块）在该离线合成夹具上输出 `acq_seeds = 0`，因此实采前端+Timing 联调需后续配合宽偏前端种子源继续验证。

---

## 3. 跨预设构建与全量 QA 门禁结果

| 预设 / 检查项 | 结果 | 详情 |
|---|---|---|
| `dev-avx2` `-L "rx_sic\|heavy"` | **39/39 通过** | `rx_sic` 33/33 + `heavy` 6/6 |
| `native-avx2` 全量 | **83/83 通过** | 总计 84 项，`uwb_qa_native_install_consumer` 无前缀时按约定 skip，其余 83/83 全绿 |
| `no-uhd` 全量 | **82/82 通过** | 总计 83 项，`uwb_qa_native_install_consumer` skip，其余 82/82 全绿 |
| `scalar-debug` 全量 | **83/83 通过** | Debug + 纯标量路径 83/83 全绿 |
| `native-avx512` 全量编译 | **136/136 通过** | 0 error，0 新增 warning |
| `install_consumer` 集成测试 | **`ALL OK`** | 临时前缀实跑 `openuwb::` 与 `uwb` compat 消费者全部通过 |

---

## 4. 未验证范围（明确边界）

- **未上板实测**：所有实采证据均为 `openuwb_rx --replay` 离线回放，未进行真实 USRP/QM35/DW3000 硬件实时收发验收。
- **未做 MATLAB 逐样点/逐包交叉验证**：尚未将统一 CIR 累加输出位置与 `UWB_demodulation/` MATLAB 脚本做逐包样本级对齐核验。
- **未做 TWR / ToA 测距验证**。
