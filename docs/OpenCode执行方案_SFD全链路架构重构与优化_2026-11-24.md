# OpenCode 多 Subagent 执行方案：SFD 全链路架构重构与低 SNR 灵敏度/速度优化

> 日期：2026-11-24  
> 前置基线：Timing（Stage-1）重构与收尾提交完成后的 `main`  
> 目标：在不放宽任何既有 QA 契约与纯噪声零误报前提下，对 SFD（Stage-3）全链路做模块化重构与五项针对性优化：
> 1. **稳定性（召回率）**：消除 `polarity_ok()` 6/6 一票否决对低 SNR 弱包的误杀，救回 8.84 m（尤其 `try2`）上绝大部分 `SFD_FAILED` 真包（`try2` FCS 从 380 向 390+ 逼近，`PAYLOAD_FAILED == 0`）；
> 2. **速度（尾延迟 P95 与均值）**：将 `try2` `sfd_boot` P95 从 **284 µs 压降至 ≤ 30 µs**（均值 `98.7 µs → ≤ 20 µs`），`clean` `sfd_boot` P95 从 **251 µs 压降至 ≤ 25 µs**；
> 3. **架构整洁度**：将 `stage_sfd()` 从 `uwb_demod_cfo_sfd.cc`（1,265 行）独立为 `uwb_demod_sfd.cc`，消除离散网格热路径 24 个临时 `std::vector` 堆分配与 `pipeline.cc` 失效重试死代码。

---

## 0. 问题诊断（基于代码与 8.84 m 实采日志）

### 0.1 实测瓶颈与漏包数据（Timing 重构后 `r1_B_*.log`，离线文件回放）

| 数据集 | `OK (fcs=1)` | `PAYLOAD_FAILED` | 目标窗 `SFD_FAILED` | `sfd_boot` P50 / P95 / 均值 (µs) | 落入 10 窗占比 |
|---|---:|---:|---:|---|---:|
| near 0.70 m | 399 | 0 | 0 | `6 / 8 / 5.9` | 0% |
| clean 8.84 m | 392 | 0 | 3 (~0.8%) | `9 / 251 / 25.1` | 6% (24/391) |
| try2 8.84 m (默认/宽松) | 380 | 0 | **15 (~3.8%)** | **`9 / 284 / 98.7`** | **35% (131/379)** |
| jam 8.84 m + SIC | 396 | 0 | 3 (~0.8%) | `7 / 15 / 16.9` | <2% |

### 0.2 五大代码级根因

1. **`polarity_ok()` 6/6 全符号一票否决导致误杀真包**（[uwb_demod_cfo_sfd.cc](file:///home/oi/Desktop/uwb-gnuradio/lib/phy/hrp/uwb_demod_cfo_sfd.cc#L403-L430)、[L720-L744](file:///home/oi/Desktop/uwb-gnuradio/lib/phy/hrp/uwb_demod_cfo_sfd.cc#L720-L744)）：
   - 8.84 m 几乎全部包的 `sfdmetric` 在 `0.01~0.04`（低于 `0.08`），强制经过 `polarity_ok()`；
   - `4z2` 序列有 6 个非零符号，现行代码要求 6 个符号在相干参考方向上的投影**每一个都严格 `> 0`**；
   - `try2` 失败包中大量包的 `tmetric = 0.055~0.097`、`residual = 0`、`cfo_khz` 完全准确，仅因 6 个符号中某 1 个被瞬时噪声翻号即被一票否决扔掉；
   - 且 `polarity_ok` 与 `polarity_ok_cached` 两份几乎相同的逻辑重复拷贝。
2. **`prepare_grid()` 无视 `derotated_span`，盲目重旋 32 万样点（~180 µs）**（[uwb_demod_pipeline.cc](file:///home/oi/Desktop/uwb-gnuradio/lib/phy/hrp/uwb_demod_pipeline.cc#L87-L99)）：
   - 一旦进入离散网格回退，`prepare_grid()` 不检查 Stage-2 已经旋转过的 `sync_state.derotated_span`，直接对 `[preamble_start, expected + 36*P + sfd_len]` 全量重旋（实测 `321,248` 样点 × `0.585 ns ≈ 188 µs`）；
   - 且常规 `(-36, +36)` 网格根本不需要前导前 220 个符号，只有极端晚捕获兜底 `(-255, -37)` 才需要。
3. **`scan_window` 门限（`0.02`）与 `run_discrete_grid` 门限（`0.008`）割裂，导致 35% 的 `try2` 包白扫 8 个全速率窗口（~80 µs）**（[uwb_demod_cfo_sfd.cc](file:///home/oi/Desktop/uwb-gnuradio/lib/phy/hrp/uwb_demod_cfo_sfd.cc#L489-L555)）：
   - 对 `0.008 <= sfdmetric < 0.02` 的弱包，中心窗 `scan_window(expected, ±64)` 已经算出真峰，却因 `< 0.02` 被拒；
   - 随后被迫白扫 `±1..±4` 符号共 8 个全速率 `±64` 窗口（1,032 次全 SFD 折叠相关），最后才在第 10 窗的离散网格以 `0.008` 捡回同一个 `expected`。
4. **热路径堆分配 + 标量能量滑窗 + 模板缓存未统一**：
   - `run_discrete_grid`（[L617-L619](file:///home/oi/Desktop/uwb-gnuradio/lib/phy/hrp/uwb_demod_cfo_sfd.cc#L617-L619)）每次进入在局部构造 `3 × 8 = 24` 个 `std::vector` 并 `.assign()`；
   - `scan_window`（[L451-L452](file:///home/oi/Desktop/uwb-gnuradio/lib/phy/hrp/uwb_demod_cfo_sfd.cc#L451-L452)）逐样点标量 `std::norm` 滑窗更新；
   - `sync_filter`（[L395-L396](file:///home/oi/Desktop/uwb-gnuradio/lib/phy/hrp/uwb_demod_cfo_sfd.cc#L395-L396)）未接入 Timing/CFO 已共享的 `ShapedSparseTemplate` 缓存。
5. **Pipeline 历史失效重试死代码**（[uwb_demod_pipeline.cc](file:///home/oi/Desktop/uwb-gnuradio/lib/phy/hrp/uwb_demod_pipeline.cc#L620-L641)）：
   - CFO 重构后 `run_tail_retrack_and_cfo_refit()` 已忽略 `seed_radius`（`(void)seed_radius;`），但 `final_sfd` 失败时 L620–641 仍以 `seed_radius = 0` 原封不动重跑一遍完全相同的 CFO refit + final SFD。

---

## 1. 硬约束与不可违反的红线

1. **禁止引入标量点积/标量能量循环**，禁止 `-march=native`，保持 `UwbSimd.cmake` 编译期 ISA 一致性。
2. **热路径零堆分配**：`stage_sfd()` 及其全部子函数（含 `search_sfd_discrete_grid`）预热后禁止临时 `std::vector` 堆分配，全部缓冲迁入 `DemodScratch` 并在构造函数预留。
3. **既有 QA 契约 100% 保持，不放宽任何断言或容差**：
   - `qa_uwb_demod_core_timing_cfo_sfd.cc`（L600/632/664）：中高 SNR 下 `sr.search_windows == 1u`（0 偏移）、`2u`（+1 符号）、`3u`（-1 符号），`sr.coarse_correlations > 0u`、`sr.fine_correlations > 0u`；
   - `qa_uwb_demod_core_heavy_grid.cc`（L97）：`shifted_sfd.first_threshold_backtrack_symbols == 4`；
   - `qa_uwb_demod_core_fec_joint_ml.cc`（L690–805）：保留 `Qm35825Profile::rx_seed_refine_tail` 与 `DemodResult::stage_tail_retrack_fallbacks` 公共字段及行为兼容；
   - `rx_sfd_par_only` / `rx_par_confirm` 路径行为不变；
   - 纯噪声 200 窗下 SFD 误报率必须为 **0/200**，移位 ±1..±7 符号的 SFD 旁瓣不得误锁。
4. **保持既有 `stdout` 协议、CLI 默认值与 `stage_sfd()` 旧 9 参数签名源码兼容**。
5. **所有构建与测试使用 `env -u LD_LIBRARY_PATH`**，不以全局 `install` 作为测试前置。

---

## 2. 目标架构与五项优化设计

### 2.1 模块拆分

将 `lib/phy/hrp/uwb_demod_cfo_sfd.cc`（1,265 行）拆分为两个职责单一的编译单元（同步更新 `lib/phy/hrp/CMakeLists.txt`）：

```text
lib/phy/hrp/
├── uwb_demod_cfo_sfd.cc   # (~720 行) 专注 fit_cfo / derotate_cfo_range / stage_cfo / 旋转与相位算子
└── uwb_demod_sfd.cc       # (~550 行, 新建) 专注 Stage-3 SFD 引擎：
                           #   prepare_sfd_filter()
                           #   check_sfd_polarity_soft()
                           #   scan_sfd_window()
                           #   search_sfd_discrete_grid()
                           #   stage_sfd()
```

### 2.2 五项核心优化

#### 优化 S1：相干加权软极性检验 `check_sfd_polarity_soft()`（修复漏包）

将两份重复的 `polarity_ok` / `polarity_ok_cached` 合并为统一辅助函数，接收 `const std::complex<float>* sym_corr`、非零序列 `sfd_sequence` 与长度 `n_check = min(len, 16)`：

1. 计算相干参考向量：
   $$\text{ref} = \sum_{i=0}^{n_{\text{check}}-1} s_i \cdot c_i$$
   若 $|\text{ref}|^2 \le 10^{-24}$，返回 `false`。
2. 对每个非零符号（$s_i \ne 0$，记非零符号总数为 $N_{\text{nz}}$），计算带符号投影：
   $$p_i = \big(\text{Re}(c_i)\,\text{Re}(\text{ref}) + \text{Im}(c_i)\,\text{Im}(\text{ref})\big) \cdot s_i$$
   统计：
   - `pos_cnt`：满足 $p_i > 0$ 的非零符号数；
   - `signed_sum`：$\sum p_i$；
   - `abs_sum`：$\sum |p_i|$；
   - `worst_neg`：$\max(0, \max_i(-p_i))$（最强反相符号幅值）。
3. 判定门限（严格兼顾灵敏度与抗移位旁瓣/纯噪声）：
   - 当 $N_{\text{nz}} \le 4$（如 IEEE 4 符号非零老 SFD）：保持全员同号 `pos_cnt == N_nz`；
   - 当 $N_{\text{nz}} \ge 5$（如 `4z2` 有 6 个非零符号，或 16 符号 SFD）：
     - 允许最多 1 个弱符号反相：`pos_cnt >= N_nz - 1`；
     - 且相干加权比：`signed_sum >= 0.65 * abs_sum`（若 `pos_cnt == N_nz` 则自动为 `1.0`；若恰有 1 个符号反相，则要求该反相投影绝对值不超过总投影绝对值的 `17.5%`，直接拦住一切移位旁瓣和强干扰伪峰）。
   - **验证要求**：在单测中穷举 `4z2` 所有非零移位（±1..±7 符号，含前后接 SYNC 符号）以及 200 窗纯噪声，确保移位旁瓣与纯噪声 0 误过。

#### 优化 S2：统一 `scan_sfd_window()` 与离散网格的 `0.008` 低门限（消除 35% 包白扫 8 窗）

- 在 `scan_sfd_window()`（非 `par_only` 路径）中，将低门限分支的下限从 `0.02f` 统一下调到 `kSfdGridPeakFloor`（`0.008f`）：
  - `window_best >= profile.sfd_detection_threshold`：直接通过（高 SNR，不变）；
  - `0.08f <= window_best < threshold` 且 `par >= 25.0`：直接通过（中 SNR，不变）；
  - `kSfdGridPeakFloor (0.008f) <= window_best < 0.08f` 且 `par >= 25.0` 且 `check_sfd_polarity_soft(...)` 通过：**直接通过**！
- **效果**：对 `0.008 <= sfdmetric < 0.02` 且 `residual ≈ 0` 的 8.84 m 弱包（`try2` 占 35%），第 1 个中心窗 `expected ±64` 直接命中返回（`search_windows == 1`，~9 µs），彻底跳过后续 8 个 `±1..±4` 符号全速率窗口与整个离散网格回退！
- **契约保持**：中高 SNR 偏移 ±1..±4 符号的包在中心窗无峰（`par < 25` 且 `metric < 0.008`），依然按 `+1, -1, +2, -2, +3, -3, +4, -4` 顺序命中，100% 保持 `search_windows == 2 / 3` 与 `first_threshold_backtrack_symbols == 4`。

#### 优化 S3：离散网格零堆分配 + 统一模板缓存 + 去标量滑窗能量

1. **`DemodScratch` 新增扁平化 SFD 离散网格复用缓冲**（`include/openuwb/phy/hrp/rx/demod_scratch.h`）：
   - 替换 `run_discrete_grid` 内部的 24 个局部 `std::vector`：
     - `std::vector<std::complex<float>> sfd_grid_sym_corr;`（预留 `8 * 544`）
     - `std::vector<float> sfd_grid_sym_pwr;`（预留 `8 * 544`）
     - `std::vector<uint8_t> sfd_grid_sym_valid;`（预留 `8 * 544`）
   - 每个 delta 槽 `d ∈ [0, 7]` 使用固定步距 `stride = n_sym_slots` 切分上述扁平缓冲，预热后 0 次堆分配。
2. **统一接入共享模板缓存**：在 `prepare_sfd_filter()` 中复用 `detail::ensure_shaped_sparse_template_cached()`（需注意当 `sparse_tmpl` 有效时保持现有 `sparse_tmpl` 优先级不变，仅在非 `sparse_tmpl` 时享受缓存；若 `shaped_tmpl` 带来额外提速需保证数值与 QA 一致）。
3. **消除 `scan_sfd_window` 逐样点标量能量更新**：
   - 将 L451–452 的标量 `pwr += std::norm(rx[j + sfd_len - 1]) - std::norm(rx[j - 1])` 改为预计算窗口边界功率差分或复用 `window_power_8acc` / 批量复数模平方数组（存入 `scratch.energy`），避免在热循环中逐样点串行依赖。

#### 优化 S4：Pipeline 增量补旋 `ensure_derotated_range()`（消灭 ~180 µs 全窗重复旋转）

- 在 `include/openuwb/phy/hrp/rx/demod_sfd.h` 中扩展 `stage_sfd` 支持按需区间补旋回调 `const std::function<void(size_t lo, size_t hi)>& ensure_derotated`（保留旧 9 参数 `std::function<void()>` 重载作为薄转发，100% 向后兼容）。
- 在 `search_sfd_discrete_grid(k_lo, k_hi)` 入口处，仅计算当前网格实际需要的精确样本区间 `[req_lo, req_hi)`（含 `sfd_len` 与 `±1` 精搜余量），调用 `ensure_derotated(req_lo, req_hi)`：
  - 常规 `(-36, +36)` 网格只请求 `[expected - 36*P - 4, expected + 36*P + sfd_len + 4)`，**不再旋转前导前 220 个符号**；
  - 只有 `(-36, +36)` 失败进入晚捕获兜底 `(-255, -37)` 时，才请求前导前半段区间。
- 在 `lib/phy/hrp/uwb_demod_pipeline.cc` 的 `search_sfd_with_wide_rescue()` 中实现真正的**增量补旋**：
  - 传入 `SyncBootstrapState& sync_state`；
  - 若 `sync_state.derotated_span.valid && sync_state.derotated_span.cfo_hz == res.cfo.cfo_hz`：
    - 若 `[req_lo, req_hi)` 已被 `derotated_span.range` 完全覆盖 → **直接返回，0 样点旋转**；
    - 否则仅对左缺口 `[req_lo, cur.lo)` 和/或右缺口 `[cur.hi, req_hi)` 调用 `rotate_window_range()`，并将 `derotated_span.range` 更新为合并后的并集区间 `[min(req_lo, cur.lo), max(req_hi, cur.hi))`！

#### 优化 S5：清理 `anchor_and_verify_final_sfd()` 中已失效的重复重试分支

- 核对 `tests/unit/phy/qa_uwb_demod_core_fec_joint_ml.cc`（L690–805）对 `rx_seed_refine_tail` 与 `stage_tail_retrack_fallbacks` 的全部断言；
- 确认无测试依赖失败时的重复空转后，移除 `uwb_demod_pipeline.cc` L620–641 中因为 `seed_radius` 已被忽略而无意义重跑完全相同 CFO refit + final SFD 的死代码（保留 `Qm35825Profile::rx_seed_refine_tail` 与 `DemodResult::stage_tail_retrack_fallbacks` 字段定义不动，保持 API 与日志字段兼容）。

---

## 3. 多 Subagent 并行实施计划

### Phase 0（主 Agent / 顺序，锁接口与基线）

1. 确认当前 `HEAD` 已包含 Timing 重构与收尾提交且 `dev-avx2` 编译通过；
2. 核对 `qa_uwb_demod_core_timing_cfo_sfd.cc`（L590–680）、`qa_uwb_demod_core_heavy_grid.cc`（L80–115）、`qa_uwb_demod_core_fec_joint_ml.cc`（L690–810）全部 SFD 契约；
3. 更新 `include/openuwb/phy/hrp/rx/demod_scratch.h`（新增 3 个 SFD 扁平网格缓冲并在构造函数 `reserve()`）；
4. 更新 `include/openuwb/phy/hrp/rx/demod_sfd.h`（新增带 `std::function<void(size_t, size_t)>` 的 `stage_sfd` 重载，旧签名保持兼容）。

### Phase 1（3 个 Subagent 并行开发，文件所有权严格互斥）

| Subagent | 负责文件（严格互斥） | 具体任务 |
|---|---|---|
| **Agent A（TU 剥离 + CMake 注册）** | `lib/phy/hrp/uwb_demod_cfo_sfd.cc`、`lib/phy/hrp/CMakeLists.txt` | 1. 从 `uwb_demod_cfo_sfd.cc` 移除 `stage_sfd()`（原 L300–848），保持 `fit_cfo`、`derotate_cfo_range`、`stage_cfo` 及旋转/相位算子逐位不变；<br>2. 在 `lib/phy/hrp/CMakeLists.txt` 的 `UWB_ENABLE_RX` 列表注册 `uwb_demod_sfd.cc`。 |
| **Agent B（SFD 引擎全面重构）** | `lib/phy/hrp/uwb_demod_sfd.cc`（新建） | 1. 拆分 `prepare_sfd_filter`、`check_sfd_polarity_soft`、`scan_sfd_window`、`search_sfd_discrete_grid`、`stage_sfd`；<br>2. 落地优化 S1（软极性检验）、优化 S2（`0.008` 门限统一）、优化 S3（`DemodScratch` 零堆分配 + 去标量滑窗 + 按需 `ensure_derotated(lo, hi)` 调用）；<br>3. 100% 保留 `search_windows`、`coarse/fine_correlations`、`first_threshold_backtrack_symbols`、`rx_sfd_par_only`、`rx_par_confirm` 既有契约。 |
| **Agent C（Pipeline 增量补旋与清理）** | `lib/phy/hrp/uwb_demod_pipeline.cc` | 1. 落地优化 S4：`search_sfd_with_wide_rescue()` 接入 `SyncBootstrapState&`，实现 `ensure_derotated_range(lo, hi)` 仅补旋未覆盖左右缺口并更新 `derotated_span`；<br>2. 落地优化 S5：在确认 `qa_uwb_demod_core_fec_joint_ml` 契约无冲突前提下，清理 L620–641 失效的重复 CFO refit 重试分支。 |

### Phase 2（集成联调 + 针对性单测）

1. `dev-avx2` 构建并运行 `-L rx_sic`、`-L sensing`、`-L core` 全绿；
2. 新建 `tests/unit/phy/qa_uwb_sfd_refactor.cc`（注册至 `tests/CMakeLists.txt`，打 `core;rx_sic` 标签），至少包含 6 项单测：
   - `test_sfd_soft_polarity_rescues_single_flipped_symbol`：-12/-14 dB 下人为令 6 个非零符号中 1 个弱反相，验证 SFD 正确通过且起点误差 ≤ 1；
   - `test_sfd_soft_polarity_rejects_shifted_sidelobes_and_noise`：穷举 `4z2` 移位 ±1..±7 符号旁瓣与 200 窗纯噪声，验证 **0 误报**；
   - `test_sfd_low_floor_hits_center_window_without_offset_scan`：`0.008 <= metric < 0.02` 的弱信号在 `expected` 直接命中，断言 `search_windows == 1`；
   - `test_sfd_incremental_derotate_matches_full_rotation`：验证增量补旋在已覆盖区间 0 额外旋转，扩窗时仅旋转缺口，且结果与整段旋转逐样点一致（最大误差 ≤ 1e-6）；
   - `test_sfd_late_acquisition_rescue_still_works`：Timing 晚 > 36 符号时，仍能经 `(-255, -37)` 离散网格正确找回 SFD；
   - `test_sfd_discrete_grid_warm_call_zero_allocation`：强制触发离散网格路径，二次热调用前后 `DemodScratch` 全部缓冲 `capacity()` 完全不变。

### Phase 3（离线回放 A/B 验收 + 跨预设门禁）

#### 3.1 串行 ABBA 离线回放（Timing 重构后 HEAD vs SFD 重构后，同机串行，`dev-avx2`）

必测数据集（与 Timing 验收报告同口径，必须带 `--interference off`）：
1. 8.84 m 三组（默认锁定）：`clean 8.84 m`、`try2 8.84 m`、`jammer 8.84 m + SIC`；
2. `try2 8.84 m` 宽松锁定；
3. 近距回归组：`0.70 m`、`3.48 m`、`5.06 m`。

#### 3.2 验收指标门槛

| 指标 | 基线（Timing 重构后） | 验收目标 |
|---|---:|---|
| `try2 8.84 m` 默认/宽松锁定 FCS | 380 | **≥ 388**（目标 390+，`PAYLOAD_FAILED == 0`） |
| `clean / jam 8.84 m` FCS | 392 / 396 | **持平或更高，`PAYLOAD_FAILED == 0`** |
| `try2 8.84 m` `sfd_boot` P95 / 均值 | 284 µs / 98.7 µs | **P95 ≤ 30 µs，均值 ≤ 20 µs** |
| `clean 8.84 m` `sfd_boot` P95 | 251 µs | **≤ 25 µs** |
| 近距 `0.70/3.48/5.06 m` `sfd_boot` P50 | 6–9 µs | **≤ 9 µs，FCS 零退化** |
| 纯噪声 200 窗误报 / 移位旁瓣误锁 | 0 / 0 | **0 / 0** |
| `dropped_queue` / `raw_drop_samples` | 0 / 0 | **0 / 0** |

#### 3.3 跨预设门禁（全部 `env -u LD_LIBRARY_PATH`）

因改动了公共头 `demod_scratch.h`、`demod_sfd.h` 与 `lib/phy/hrp/CMakeLists.txt`，交付前须通过：
1. `native-avx2` 全量构建 + 全量 `ctest`；
2. `no-uhd` 全量构建 + 全量 `ctest`；
3. `scalar-debug` 全量构建 + 全量 `ctest`；
4. `native-avx512` 全量编译 + 核验无重复强符号；
5. 临时前缀实跑 `uwb_qa_native_install_consumer`。

---

## 4. 交付物清单

1. 源码与头/CMake：`lib/phy/hrp/uwb_demod_sfd.cc`（新建）、`lib/phy/hrp/uwb_demod_cfo_sfd.cc`、`lib/phy/hrp/uwb_demod_pipeline.cc`、`include/openuwb/phy/hrp/rx/demod_sfd.h`、`include/openuwb/phy/hrp/rx/demod_scratch.h`、`lib/phy/hrp/CMakeLists.txt`。
2. 单测：`tests/unit/phy/qa_uwb_sfd_refactor.cc`（新建，6 例）与 `tests/CMakeLists.txt`。
3. 验收报告：`docs/performance/SFD全链路架构重构与优化验收报告_2026-11-24.md`（明确标注离线文件回放证据边界、完整帧 FCS、未覆盖项）并更新 `开发状态.md`。
4. 规范本地提交（拆分代码/测试 commit 与文档 commit，严禁提交 `testdata/`，不 push 远端）。
