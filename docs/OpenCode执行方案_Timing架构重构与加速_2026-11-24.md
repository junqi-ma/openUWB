# OpenCode 多 Subagent 执行方案：Timing（Stage-1）架构重构与低 SNR 加速

> 日期：2026-11-24  
> 基线提交：`fc5e0a7`（`main`，已含 CFO ML 架构重构与 `12fcaee` 优化）  
> 目标：在不损失任何检测灵敏度、锁定稳定性与既有 QA 契约的前提下，将 8.84 m 低 SNR 场景下 `st_timing` 从 **~720 µs 压降至 ~120–180 µs**（整包 `demod_us` 从 ~1100 µs 降至 **~550–620 µs**），同时将 2,445 行 / 1,405 行巨型函数的 `uwb_demod_timing.cc` 拆解为职责单一、热路径零堆分配的模块化实现。

---

## 0. 背景与问题诊断（基于代码与 `20261123` 实采日志）

### 0.1 实测瓶颈数据

对比不同 SNR 下的 `st_timing` P50（`12fcaee`，`dev-avx2`）：

| 场景 | `tmetric` P50 | 策略 / 折叠 | `st_timing` P50 (µs) | `coarse` + `fine` (µs) |
|---|---:|---|---:|---:|
| 极近 `ch0`（高 SNR） | 0.343 | `locked_bounded`, `fold=0` | **79** | `0 + 79` |
| `0.70–5.06 m`（中高 SNR） | 0.106–0.149 | `acquire_legacy`, `fold=0` | **89–99** | `0 + 89–99` |
| `8.84 m` clean / try2 / jam（低 SNR） | 0.074 | `acquire_legacy`, **`fold=1`** | **719–722** | **`508 + 212`** |

### 0.2 五大代码级根因

1. **低 SNR 下 `AcquireLegacy` 白烧 ~508 µs 无效盲搜**（[uwb_demod_timing.cc](file:///home/oi/Desktop/uwb-gnuradio/lib/phy/hrp/uwb_demod_timing.cc#L571-L899)）：
   - 8.84 m 锁定帧预测残差仅 ~2 样点，但因默认 `--sync acquire` 且单符号 `tmetric ≈ 0.074` 过不了第 1/2 步窄搜门限，代码被迫先跑完整个 `±131072` 样点的 `coarse_scan`（~508 µs）和全 SYNC 网格探针，全部失败后才掉到第 7 步 `try_fold_acquire`（K=4 折叠）一把救回。
2. **`fold_anchor_fit()` 做了“双倍 64 组重折叠”+ 标量能量循环（~190–212 µs）**（[uwb_demod_timing.cc](file:///home/oi/Desktop/uwb-gnuradio/lib/phy/hrp/uwb_demod_timing.cc#L1643-L1780)）：
   - 第一遍把 256 符号切成 64 组做 K=4 折叠 + 17 lag 相关做 OLS；
   - 第二遍又把 64 组全部重新折叠 + 相关 + **标量 `std::norm` 能量循环**（L1765），只为了给旧 `fit_cfo()` 填 64 个 `peak_corr`。
   - 现在 Stage-2/3 已经采用独立 ML CFO 估计器，`fold_anchor_fit` 只需在完整 256 符号跨度上取 **均匀稀疏 16 组**（且第二遍复用/按需最多 16 点），耗时可降 **~4×**（~212 µs → ~40–50 µs），且同时满足 Stage-2 legacy `fit_cfo` 交叉校验与 `qa_uwb_radar_timing_core` 契约。
3. **`stage_timing()` 每次调用重复对 1016 点模板做 FNV-1a 全量哈希**（[uwb_demod_timing.cc](file:///home/oi/Desktop/uwb-gnuradio/lib/phy/hrp/uwb_demod_timing.cc#L62-L65)）：
   - 未复用 `uwb_demod_cfo_sampler.cc` 中已引入的 `DemodScratch` 8 点探针缓存键。
4. **热路径存在隐式堆分配**：
   - `uwb_demod_timing.cc` L264 `std::vector<float> cand_metric(ncand)`；
   - L1200 `std::vector<size_t> peak_indices`；
   - L1543 `std::vector<float> lag_metrics(lags)`。
5. **代码严重臃肿，维护困难**：
   - `stage_timing()` 单个函数 1,366 行，内嵌 6 个大 lambda 与 7 层嵌套 fallback；
   - 文件尾部 L1786–2445（约 660 行）混放了通用 SIMD 向量算子、旋转斜坡、稀疏/成形模板构造与 CIR 相关，与 Timing 逻辑无关。

---

## 1. 硬约束与不可违反的红线

1. **禁止引入标量点积/标量能量循环**：全部相关、折叠、能量计算必须复用已有 AVX2/AVX-512 原语（`SyncMatchedFilter`、`scan_normalized_sync`、`scan_raw_sync`、`fold_accumulate`、`window_power_8acc`、`vector_accumulate_complex`），禁止 `-march=native`。
2. **热路径零堆分配**：`stage_timing`、`fold_acquire`、`fold_anchor_fit`、`disambiguate_fold_leading_edge` 预热后禁止新增 `std::vector` 堆分配，全部临时缓冲放入 `DemodScratch`。
3. **既有 QA 契约 100% 保持，不放宽任何断言或容差**：
   - `LockedBounded` 契约：`bounded_legacy_scans == 0`、`locked_tail_probe == false` 时必须 hard-fail（不走 rescue/fold）、`grid_single_points == preamble_repetitions`、`detected_peaks == preamble_repetitions`（高 SNR 场景）；
   - `AcquireLegacy` 契约：高 SNR 逐符号跟踪 `detected_peaks`、回溯 `backtrack_symbols`、`rx_wide_interval_probe` 快失败行为（SIC L0_Bypass）全部不变；
   - `Stage-2` 交叉校验契约：折叠路径输出的 `peak_samples / peak_corr` 仍需跨越完整前导、供 `fit_cfo()` 给出稳定的交叉校验频率，确保 `try2` 宽松锁定不回归（`≥379/433`，`PAYLOAD_FAILED == 0`）。
4. **保持既有 `stdout` 协议、CLI 参数默认值（`--sync acquire`、`--cir full`、`--fused on`、`--taps eq16`、`--sic off`）不变**。
5. **所有构建与测试必须使用 `env -u LD_LIBRARY_PATH`**，不以全局 `install` 作为测试前置。

---

## 2. 目标架构与模块拆分

将原 2,445 行的 `lib/phy/hrp/uwb_demod_timing.cc` 拆分为 3 个职责单一的编译单元（并更新 `lib/phy/hrp/CMakeLists.txt`）：

```text
lib/phy/hrp/
├── uwb_demod_simd_utils.cc      # (~660 行) 原 L1786-2445：rxopt SIMD 向量算子、模板构造、CIR 相关
├── uwb_demod_timing_fold.cc     # (~400 行) 折叠引擎：probe_fold_phase / fold_acquire / fold_probe_acquire /
│                                #           disambiguate_fold_leading_edge / fold_anchor_fit (稀疏 16 组优化)
└── uwb_demod_timing.cc          # (~550 行) 顶层编排 stage_timing() + 4 个独立子函数：
                                 #           prepare_timing_filter()
                                 #           acquire_locked_bounded()
                                 #           acquire_unbounded_seed()   <-- 低 SNR 局部 K=4 折叠前置快路径
                                 #           track_locked_grid() / track_sync_train()
```

### 2.1 核心提速设计（3 项零损优化）

#### 优化 A：`acquire_unbounded_seed()` 低 SNR 局部 K=4 折叠前置（消灭 ~508 µs 无效粗扫）

在 `AcquireLegacy` 路径中，保留第 1/2 步单符号窄搜（高/中 SNR 毫秒级命中，零变化）：

```text
[1] seed 窄搜 (narrow) + legacy 窄搜 (±256)
     ├─ 命中 -> 走原有高 SNR 跟踪（89-99 µs，行为 100% 不变）
     └─ 未命中：
         ├─ (若开启 rx_wide_interval_probe，保持原有 SIC 宽窗间隔探针与快失败逻辑不变)
         ├─ [2] 新增局部 K=4 折叠快路径（仅当 seed_start >= 0 且 !rx_wide_interval_probe）：
         │      - 在 seed ± min(margin, 8192) 上先只跑 K=4 相干折叠（fold_acquire 内部已钳制到 ±P/2，
         │        仅 1016 lag，成本 ~15-20 µs；注意此处先不跑耗时的 probe_fold_phase + K=64！）
         │      - 若 K=4 折叠命中（metric >= kFoldThreshold 或 PAR >= 12.5 dB）
         │        且前沿消歧 + 稀疏 16 组 fold_anchor_fit 成功：
         │        -> 直接返回成功！彻底跳过 TM-1、±131072 coarse_scan（508 µs）与全网格探针！
         └─ [3] 若 K=4 局部折叠未命中（例如极低 SNR 需要 K=64，或无锁盲搜 seed 偏出 ±8192）：
                -> 无缝落入原有 TM-1 -> coarse_scan -> 全网格探针 -> 完整 try_fold_acquire (含 K=64) 兜底！
```

**关键安全性检查（必须实现）**：
- **防假峰误锁（纯噪声/偏移 seed 场景）**：局部 K=4 折叠前置命中后，必须满足：
  1. `fold_acquire(K=4)` 通过既有 `kFoldThreshold` 或 `PAR >= 12.5 dB` 门限；
  2. 前沿消歧 `disambiguate_fold_leading_edge` 找到有效前沿（`max_here >= 0.003f`）；
  3. `fold_anchor_fit`（稀疏 16 组）成功拟合出合法周期（`|p - P| <= 0.005 * P`）且有效锚点组数 `>= 4`、中位相关度合理。
  - 若任一条件不满足，则**丢弃局部折叠结果并回退到 `[3]` 原有盲搜路径**，保证无锁大偏移或极低 SNR 场景的召回率 100% 不下降。
- **单符号 vs 折叠门限衔接**：当 `seed_start < 0`（无 seed 的纯盲搜）时跳过 `[2]` 直接走 `[3]`。

#### 优化 B：`fold_anchor_fit()` 稀疏 16 组 + 去标量循环（~212 µs → ~40–50 µs）

1. **均匀稀疏 16 组覆盖完整前导**：
   - 将 `try_fold_acquire` 传入 `fold_anchor_fit` 的目标组数改为 `kFoldAnchorMaxGroups = 16`（原为 64），在 `[0, total_groups - 1]` 之间**均匀采样 16 个组索引**（含首尾组，跨越完整 256 符号基线，而不是只取前 16 组！）。
   - 计算公式：`g_idx = round(i * (total_groups - 1) / (n_anchors - 1))`，其中 `total_groups = preamble_repetitions / kFoldAnchorGroupSymbols`（256 SYNC 时 `total_groups = 64`，覆盖 `[0, 63]` 全跨度，时间基线与 64 组完全相同）。
2. **第二遍按拟合网格采样同样只做这 16 组（4× 提速），并消除标量能量循环**：
   - 将 L1764–1766 的标量循环：
     ```cpp
     for (size_t k = 0; k < L; ++k) rx_energy += std::norm(sup[k]);
     ```
     替换为已有 AVX2 向量化友好的 `window_power_8acc(sup.data(), L)`；
   - 保持拟合网格采样逻辑（`lock_integer_grid` 漂移 `< 1.5` 样点锁整数周期），使输出的 16 个 `peak_samples / peak_corr` 无多径跳峰相位突变，同时精确满足 Stage-2 `fit_cfo()` 交叉校验与 `qa_uwb_radar_timing_core`（`detected_peaks <= 16`）契约。
3. **注意**：检查 `fold_anchor_fit` 是否有独立公开声明或 QA 直接调用其 `groups` 参数；若外部测试直接传 `groups`，保留按 `groups` 参数采样、当 `groups <= total_groups` 时均匀分布或兼容原签名（可新增 `stride_groups` 或内部按 `total_groups` 均匀分布）。

#### 优化 C：共享模板探针缓存 + 热路径零堆分配

1. 在 `include/openuwb/phy/hrp/rx/demod_scratch.h` 中增加 Timing 热路径复用缓冲：
   - `std::vector<float> timing_fold_cand_metrics;`（替代 L264 局部 `cand_metric`，预留 64）
   - `std::vector<float> timing_fold_lag_metrics;`（替代 L1543 局部 `lag_metrics`，预留 2048）
   - `std::vector<size_t> timing_peak_indices;`（替代 L1200 局部 `peak_indices`，预留 512）
2. 将 `uwb_demod_cfo_sampler.cc` 中的 `template_sample_probe()` + 缓存命中判断抽取为共享 inline/helper（例如放在 `include/openuwb/phy/hrp/rx/demod_sync.h` 或 `demod_scratch.h` 的 `detail::ensure_shaped_sparse_template_cached()`），供 `stage_timing()` 与 `sample_sync_symbols_simd()` 共同调用，彻底消除每包 1016 点 FNV-1a 全量哈希。

---

## 3. 多 Subagent 并行实施计划

### Phase 0（主 Agent / 顺序，锁接口与基线，~10 分钟）

- **职责**：
  1. 确认 `lib/phy/hrp/CMakeLists.txt` 源码列表与 PCH 配置；
  2. 在 `include/openuwb/phy/hrp/rx/demod_scratch.h` 中新增 3 个 scratch 缓冲字段并在构造函数/预热中 `reserve`；
  3. 抽出共享模板缓存函数 `ensure_shaped_sparse_template_cached()`（供 Timing 与 CFO 采样器复用）；
  4. 核对 `tests/unit/sensing/qa_uwb_radar_timing_core.cc`、`tests/unit/phy/qa_uwb_demod_core_timing_cfo_sfd.cc`、`tests/unit/runtime/qa_uwb_target_lock_core.cc` 中所有直接调用 `fold_anchor_fit` / `fold_acquire` 的签名。
- **交付门禁**：提交 Phase 0 骨架，`dev-avx2` 增量编译通过。

### Phase 1（3 个 Subagent 并行开发，文件所有权互斥）

| Subagent | 负责文件（严格互斥） | 具体任务 |
|---|---|---|
| **Agent A（SIMD 工具剥离 + 模板缓存接入）** | `lib/phy/hrp/uwb_demod_simd_utils.cc`（新建）、`lib/phy/hrp/uwb_demod_cfo_sampler.cc`、`lib/phy/hrp/CMakeLists.txt` | 1. 将原 `uwb_demod_timing.cc` L1786–2445 完整迁移至 `uwb_demod_simd_utils.cc`，保持逐位一致、所有 `target("avx2")` / AVX-512 分发宏与编译期 ISA 一致性不变；<br>2. 让 `uwb_demod_cfo_sampler.cc` 改调共享 `ensure_shaped_sparse_template_cached()`；<br>3. 更新 `lib/phy/hrp/CMakeLists.txt` 注册新 `.cc` 文件。 |
| **Agent B（折叠引擎独立 + 稀疏 16 组优化）** | `lib/phy/hrp/uwb_demod_timing_fold.cc`（新建）、`include/openuwb/phy/hrp/rx/demod_sync.h` | 1. 迁移 `probe_fold_phase`、`fold_acquire`、`fold_probe_acquire`、`fold_anchor_fit` 至 `uwb_demod_timing_fold.cc`；<br>2. 抽离 `disambiguate_fold_leading_edge()`，使用 `scratch.timing_fold_cand_metrics` 消除堆分配；<br>3. `fold_acquire()` 使用 `scratch.timing_fold_lag_metrics` 消除堆分配；<br>4. 实现优化 B：`fold_anchor_fit` 均匀稀疏 16 组覆盖完整前导基线 + 用 `window_power_8acc` 替换标量能量循环。 |
| **Agent C（`stage_timing` 模块化重构 + 局部折叠快路径）** | `lib/phy/hrp/uwb_demod_timing.cc` | 1. 将 1,366 行 `stage_timing()` 拆为 `prepare_timing_filter`、`acquire_locked_bounded`、`acquire_unbounded_seed`、`track_locked_grid`、`track_sync_train`、`stage_timing`；<br>2. 接入共享 `ensure_shaped_sparse_template_cached()` 与 `scratch.timing_peak_indices`；<br>3. 在 `acquire_unbounded_seed()` 实现优化 A（低 SNR 局部 K=4 折叠前置 + 严格门限校验 + 未命中无缝回退）；<br>4. 100% 保持 `LockedBounded`、`tail_sweep`、`backtrack`、`rx_wide_interval_probe` 所有既有契约。 |

### Phase 2（集成联调 + 单测与边界验证）

1. 构建 `dev-avx2` 并运行 `ctest --test-dir build/dev-avx2 -L rx_sic --output-on-failure -j` 及 `-L sensing`（重点确认 `qa_uwb_demod_core_timing_cfo_sfd`、`qa_uwb_radar_timing_core`、`qa_uwb_target_lock_core`、`qa_uwb_cfo_ml`、`qa_uwb_live_packet_processor_sic` 全绿）。
2. **新增针对性单测**（可加在 `qa_uwb_demod_core_timing_cfo_sfd.cc` 或新建 `qa_uwb_timing_refactor.cc`）：
   - 低 SNR（-8 dB / -12 dB）有准确 seed 时命中局部 K=4 折叠快路径，验证 `coarse_search_us == 0`、`sync_peaks_folded == true`、`detected_peaks <= 16` 且起点误差在容差内；
   - 低 SNR 大 seed 偏移（超出 ±8192）时正确回退至宽窗/兜底路径，召回率不下降；
   - 纯噪声窗口下局部 K=4 折叠不误报；
   - 大 CFO（±40 kHz）与大 SFO（±27 ppm）下局部 K=4 折叠未命中时无缝回退至 K=64 + `probe_fold_phase`，仍正确捕获；
   - 热路径零堆分配检查（二次调用 `stage_timing` 不增长 `scratch` 容量）。

### Phase 3（实采离线回放 A/B 验收 + 全预设门禁）

#### 3.1 串行实采回放 A/B（`fc5e0a7` vs 重构后，同机串行，`dev-avx2`）

必测数据集（与 [CFO架构重构验收报告_2026-11-24.md](file:///home/oi/Desktop/uwb-gnuradio/docs/performance/CFO%E6%9E%B6%E6%9E%84%E9%87%8D%E6%9E%84%E9%AA%8C%E6%94%B6%E6%8A%A5%E5%91%8A_2026-11-24.md) 完全同口径）：
1. **8.84 m 三组（默认锁定 `--sync acquire`）**：
   - `clean 8.84 m`
   - `try2 8.84 m`
   - `jammer 8.84 m + SIC`
2. **`try2 8.84 m` 宽松锁定**（`--lock-observations 2 --holdover 500 --reacquire 1000 --lock-fcs off`）：验证无回归（`FCS ≥ 379/433`，`PAYLOAD_FAILED == 0`）。
3. **近距离回归组（默认锁定）**：`0.70 m`、`3.48 m`、`5.06 m`（确认高 SNR `st_timing` 保持 ≤99 µs 且 FCS 不降）。
4. **附加对照**：8.84 m `clean` 跑一轮 `--sync locked` 对照，记录 `st_timing` 与 FCS 差异。

#### 3.2 验收指标门槛

| 指标 | 基线 (`fc5e0a7`) | 验收目标 |
|---|---:|---|
| 8.84 m 默认锁定 `st_timing` P50 | 719–722 µs | **≤ 200 µs**（目标 120–180 µs） |
| 8.84 m 默认锁定 `timing_coarse` P50 | 508 µs | **0 µs**（命中局部折叠快路径） |
| 8.84 m 默认锁定 `demod_us` P50 | 1078–1123 µs | **≤ 650 µs**（降幅 ≥ 40%） |
| 8.84 m 三组 FCS / `PAYLOAD_FAILED` | 391 / 380 / 396，`PF = 0` | **持平或更高，`PF == 0`** |
| `try2` 宽松锁定 FCS / `PAYLOAD_FAILED` | 379/433，`PF = 0` | **≥ 375，`PF == 0`** |
| 近距离 `0.70/3.48/5.06 m` `st_timing` P50 | 89–99 µs | **≤ 99 µs**（无退化） |
| `dropped_queue` / `raw_drop_samples` | 0 / 0 | **0 / 0** |

#### 3.3 交付前跨预设门禁（含 CFO + Timing 合并收口）

由于改动了公共头 `demod_scratch.h`、`demod_sync.h` 及 `lib/phy/hrp/CMakeLists.txt`，须执行完整门禁（全部使用 `env -u LD_LIBRARY_PATH`）：
1. `native-avx2` 全量构建 + 全量 `ctest`；
2. `no-uhd` 全量构建 + 全量 `ctest`；
3. `scalar-debug` 全量构建 + 全量 `ctest`；
4. `native-avx512` 全量编译 + 核验强符号（无 AVX-512 机器仅编译）；
5. 临时前缀实跑 `uwb_qa_native_install_consumer`（`tests/integration/install_consumer/run_native_install_consumer.sh`）。

---

## 4. 交付物清单

1. 拆分并重构后的源码：`uwb_demod_timing.cc`、`uwb_demod_timing_fold.cc`、`uwb_demod_simd_utils.cc`、`demod_scratch.h`、`demod_sync.h`、`uwb_demod_cfo_sampler.cc`、`lib/phy/hrp/CMakeLists.txt`。
2. 单元测试更新与新增用例。
3. 验收报告：`docs/performance/Timing架构重构与加速验收报告_2026-11-24.md`（包含 `fc5e0a7` vs 新提交的串行对照表、各数据集 `st_timing` / `timing_coarse` / `timing_fine` / `demod_us` / FCS / CFO 精度统计、以及未覆盖项清单）。
