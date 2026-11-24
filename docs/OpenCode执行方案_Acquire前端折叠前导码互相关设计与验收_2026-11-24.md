# OpenCode 执行方案：Acquire 前端折叠前导码互相关（±100 kHz CFO 鲁棒）设计与验收

日期：2026-11-24  
适用仓库：`/home/oi/Desktop/uwb-gnuradio`（原生 UHD/C++ 架构）  
目标角色：由 `opencode` 统筹多 Subagent 分工实施、测试并完成 7 场景离线回放 ABBA 验收。

---

## 1. 背景与问题根因

在 2.0 s 实测 DW3000 离线回放（400 包，5 ms 周期）中，当前基于延迟自相关（lag = 750 样点）的流式前端 `openuwb::live::Acquisition`（[live_acquisition.h](file:///home/oi/Desktop/uwb-gnuradio/include/openuwb/runtime/live_acquisition.h)）暴露出系统性虚警问题：

1. **自相关对非白底噪与窄带残余极度敏感**：`peak_thr = 3.0e-4` 假设输入为纯白噪声（理论底噪 `1/192000 ≈ 5.2e-6`）。实测即便经过 `-232.96 MHz` 单音消除，静默段仍存在 `-232.96 MHz` 残余（~1.1%）、`-124.1/-124.4 MHz` 窄带杂散（最高达 7% 噪声功率）以及带限有色底噪，使静默自相关值长期维持在 `3.5e-4 ~ 8.2e-2`，永久越过 `peak_thr`。
2. **机械周期性假种子挤占 `acquire` 窗口**：一旦常驻越限，`Acquisition::scan` 每隔 `holdoff_len (589,824) + fall_after (65,536) ≈ 0.89 ms` 机械发射假种子（2.0 s 产生 **2,237–2,284** 个 `acq_seeds`，真实包仅 400 个），每个假种子占用唯一 `pending_seed` 槽位 `~0.51–0.61 ms` 并施加 `0.80 ms` 盲区，导致：
   - `clean` 阶段出现 45 次 `TIMING_FAILED`、首锁延迟 24.6 ms（加一次 `reacquire` 至 44.6 ms）；
   - `try2` 出现 25 次 `TIMING_FAILED` + 2 次 `SFD_FAILED`、首锁延迟 23.6 ms；
   - `near506`（70 cm 强 Code 9 干扰 + 506 cm 目标）出现 34 次 `TIMING_FAILED`、首锁延迟 28.4 ms；
   - 前端扫描耗时 `frontend_scan_us ≈ 326–359 ms / 2.0 s`。

参考实际 DW3000 接收机的 preamble hunting 架构，本方案引入**原生 `737.28 MS/s` 域混合相干+非相干折叠前导码互相关前端（`PreambleXcorrAcquisition`）**，通过已知目标 HRP 前导码（Code 9–12）做匹配滤波与 PAR（Peak-to-Average Power Ratio）判决，从根源上免疫单音残留、`-124 MHz` 杂散、有色底噪、Wi-Fi 6E 带外突发以及异码（Code 9）强干扰。

---

## 2. 物理参数设计与 ±100 kHz CFO 实测标定

### 2.1 为什么相干叠加数必须取 `K_coh = 4`

HRP SYNC 符号时长 $T_{\text{sym}} = \frac{1016}{998.4\text{ MHz}} = 1.017628\ \mu\text{s}$，在原生采样率 `737.28 MS/s` 下对应精确周期：
$$P_{\text{nat}} = 1016 \times \frac{48}{65} = 750.276923076923\ \text{samples}$$

连续相干叠加 $K_{\text{coh}}$ 个符号时，载波频偏 $\Delta f$ 导致的功率衰减为 $\text{sinc}^2(\Delta f \cdot K_{\text{coh}} \cdot T_{\text{sym}})$：
- $K_{\text{coh}} = 16$：第一零点 **61.4 kHz**（50 kHz 以上已严重失效）；
- $K_{\text{coh}} = 8$：第一零点 **122.8 kHz**（100 kHz 衰减 **-13.3 dB**，失效）；
- **`K_coh = 4`**：第一零点 **245.7 kHz**，在 **±50 kHz 仅衰减 -0.6 dB**、**±100 kHz 仅衰减 -2.5 dB**、**±120 kHz 仅衰减 -3.5 dB**。

### 2.2 合成 AWGN（±100 kHz CFO）与真实 DW3000 采集标定结果

采用 64-tap 稀疏三元 Code 10 模板（无乘法，`stride = 1` 全 750 lag），每组 `K_coh = 4` 符号相干折叠、多组非相干能量累加（`G_full = 16` 组 = 64 符号 ≈ 65.1 µs），并引入 **两级早停**（Stage 1 `G_fast = 6` 组 → Stage 2 补算剩余 10 组）：

| 测试场景 | 配置 | 信号 PAR（中位数 / 最小） | 静默/纯噪最大 PAR | 备注 |
| :--- | :--- | :--- | :--- | :--- |
| AWGN `-14 dB` + `100 kHz` CFO（100 次随机相位） | `G_fast=6` / `G_full=16` | `18.0 / 12.2` → `17.8 / 13.9` | `3.80` → `2.53` | Stage 1 (`>=4.5`) 信号通过 100%、噪声通过 0% |
| AWGN `-16 dB` + `100 kHz` CFO（100 次随机相位） | `G_fast=6` / `G_full=16` | `11.7 / 7.6` → `12.4 / 8.7` | `3.62` → `2.26` | Stage 1 (`>=4.5`) 信号通过 100%、噪声通过 0% |
| 实测 `clean_ch0`（前 60 ms，12 包） | `K=4, G=16` | `177.0 ~ 227.9` | `1.83` | 12/12 检出，首包 1.39 ms，0 虚警 |
| 实测 `try2_884cm`（经 `-232.96 MHz` TC，前 60 ms，12 包） | `K=4, G=16` | `34.3 ~ 54.4` | `3.27` | 12/12 检出，首包 3.60 ms，0 虚警 |
| 实测 `near506_jam`（70 cm 强 Code 9 干扰 + 506 cm 目标） | `K=4, G=16` | `205.1 ~ 224.7` | `2.93` | 12/12 检出，首包 3.37 ms，Code 9 干扰完全抑制 |
| 实测 `near70_jam`（前 60 ms，12 包） | `K=4, G=16` | `206.7 ~ 232.8` | `2.85` | 12/12 检出，首包 1.20 ms，0 虚警 |

> **关键设计约束（已实测证实）**：
> 1. **禁止 `lag_stride = 2` 抽样粗扫**：稀疏三元码未做脉冲成形滤波，奇数 lag 处主瓣窄，实测 `stride = 2` 会使 `-14 dB / 100 kHz` 最小 PAR 从 `14.6` 暴跌至 `5.8`。必须保留 **`stride = 1` 全 750 lag**，改用**占空比步进（`hop_symbols`）+ 两级早停（`G_fast = 6 → G_full = 16`）**控制算力。
> 2. **必须保持精确浮点起址重对齐**：每个符号的折叠起点取 `st = llround(base + (g * K_coh + k) * P_nat)`，保证跨 16 组 64 符号（累计小数漂移 17.72 样点）时，所有组在 750 点循环坐标系内的峰值位置误差始终 `<= ±0.5` 样点，不散焦。

---

## 3. 范围边界与硬约束（遵守 `AGENTS.md`）

1. **严格限定范围**：仅针对 `openuwb_rx` 连续模式（`continuous`）的 acquire/扫描前端。**不改动** TWR、ToA、Radar/Sensing、SIC 内核、重采样器系数与后端解调算法。
2. **默认行为完全兼容（opt-in）**：
   - 新增 CLI 开关 `--acq-mode autocorr|xcorr`，**默认值保持 `autocorr`**（`--acq-mode autocorr` 下行为、stdout 协议、既有 QA 100% 保持不变；待 7 场景 ABBA 验收通过后再评估是否切换默认）。
   - 不删除或修改现有 `openuwb::live::Acquisition` 及其 QA 容差。
3. **单向依赖与头文件分层**：
   - 新模块归属 `openuwb_runtime`（依赖链 `support → dsp → phy_hrp → runtime`，`openuwb_runtime` 已链接 `openuwb_phy_hrp` 且已应用 `uwb_simd_apply(openuwb_runtime)`，可直接复用 `openuwb::demod::kPreambleCode9..12`）。
   - 公共头 `include/openuwb/runtime/preamble_xcorr_acquisition.h` 仅放声明与配置结构体；重型折叠、AVX2/标量稀疏相关与前沿回溯实现放在 `lib/runtime/uwb_preamble_xcorr_acquisition.cc`。
4. **编译期 SIMD 一致性**：
   - 禁止 `-march=native`。使用 `#if defined(__AVX2__) && defined(__FMA__)` 提供 AVX2 路径，同时提供逻辑等价的标量路径（供 `scalar-debug` 与非 AVX2 平台）。
   - 新增 QA 目标若调用 `uwb_simd_apply()`，同步按 `tests/CMakeLists.txt` 约定处理（复用 `uwb_register_qa` 即可自动对齐 PCH 与 SIMD flags）。
5. **构建、测试与数据规范**：
   - 所有命令统一前缀 `env -u LD_LIBRARY_PATH`。
   - 严禁提交大体积 `testdata/dw3000_*` 文件或未请求地 `git push`。
   - 报告中严格区分离线回放（replay）与硬件上板证据，不把前端种子数或 SFD-only 写成完整帧 FCS 成功率。

---

## 4. 核心架构与算法设计

### 4.1 类与接口：`openuwb::live::PreambleXcorrAcquisition`

- 头文件：`include/openuwb/runtime/preamble_xcorr_acquisition.h`
- 源文件：`lib/runtime/uwb_preamble_xcorr_acquisition.cc`（加入 `lib/runtime/CMakeLists.txt` 的 `_uwb_runtime_sources`）

```cpp
namespace openuwb::live {

struct PreambleXcorrConfig {
    size_t code_index = 10;          // 9..12 (IEEE 802.15.4a/z Ipatov 127)
    size_t preamble_reps = 256;      // 目标前导码重复数（用于推导默认 hop 与回溯上限）
    size_t k_coh = 4;                // 相干折叠符号数（默认 4，±100 kHz CFO 衰减 -2.5 dB）
    size_t g_fast = 6;               // Stage-1 早停组数（默认 6 组 = 24 符号）
    size_t g_full = 16;              // Stage-2 确认总组数（默认 16 组 = 64 符号）
    size_t hop_symbols = 128;        // 主扫描步进符号数（0 = 自动: clamp(preamble_reps/2, k_coh, 128)）
    double par_fast_thr = 4.5;       // Stage-1 早停门限
    double par_thr = 6.5;            // Stage-2 判决门限
    int64_t holdoff_len = 589824;    // 命中后盲区（~0.80 ms，与 Acquisition 一致）
};

struct PreambleXcorrStats {
    uint64_t windows_scanned = 0;    // Stage-1 执行次数
    uint64_t stage2_promoted = 0;    // 通过 Stage-1 进入 Stage-2 次数
    uint64_t seeds_emitted = 0;      // 最终提交 seed 次数
};

class PreambleXcorrAcquisition {
public:
    static constexpr size_t kSymLen = 750;
    static constexpr double kExactPeriod = 1016.0 * 48.0 / 65.0; // 750.276923076923

    PreambleXcorrAcquisition();
    void configure(const PreambleXcorrConfig& cfg);
    void reset(int64_t start_sample);

    // 与 Acquisition::scan 一致的流式推进接口：
    // 要求调用方保证 ring 中 [probe, probe + required_window_samples()) 均在 [valid_begin, end) 内。
    // 返回 >= 0 表示命中的前导码起始 seed（绝对样点坐标），-1 表示未命中。
    int64_t required_window_samples() const noexcept;
    int64_t probe() const noexcept { return probe_; }
    const PreambleXcorrStats& stats() const noexcept { return stats_; }

    int64_t scan(const std::vector<std::complex<int16_t>>& ring,
                 int64_t valid_begin,
                 double* out_metric = nullptr);
};

} // namespace openuwb::live
```

### 4.2 关键内核实现细节

1. **三元码表预编译（`configure` 时一次性完成，热路径零分配）**：
   - 从 `demod::kPreambleCode9..12` 提取 64 个非零 chip 索引 $c_i \in [0, 126]$ 与符号 $s_i \in \{-1, +1\}$。
   - 计算原生 750 点循环偏移：`off = llround(8.0 * c_i * 48.0 / 65.0) % 750`。
   - 拆分为 `pos_offsets_`（`+1` 偏移表）与 `neg_offsets_`（`-1` 偏移表），固定栈/成员变量数组（最大 64 项）。

2. **单组相干折叠（`fold_group`）**：
   - 给定组号 $g$ 和窗口起点 `win_start`，对 $k = 0 \dots K_{\text{coh}}-1$ 计算精确起点 `st = win_start + llround((g * K_coh + k) * kExactPeriod)`。
   - 将 `ring[(st + t) & mask]` 的 `int16_t` I/Q 转为 `float` 累加到 `fold_re_[0..749]`、`fold_im_[0..749]`（非跨 ring 边界时直接线性遍历，跨边界时拆两段 memcpy/累加）。
   - **线性双倍展开**：将 `fold_re_[0..749]` 与 `fold_im_[0..749]` 各自复制一份到 `[750..1499]`（`memcpy` 仅 3 KB × 2），彻底消除后续循环相关中的 `% 750` 取模与条件分支。

3. **64-tap 零乘法稀疏循环相关（`correlate_and_accumulate_power`）**：
   - 对 `t = 0 .. 749`（AVX2 8-wide 步进，处理 `t..t+7`，尾部标量补齐至 750）：
     - 初始化 `acc_re = 0, acc_im = 0`；
     - 遍历 `pos_offsets_`：`acc_re += fold_re[t + off]`, `acc_im += fold_im[t + off]`；
     - 遍历 `neg_offsets_`：`acc_re -= fold_re[t + off]`, `acc_im -= fold_im[t + off]`；
     - 非相干累加功率：`acc_power_[t] += acc_re * acc_re + acc_im * acc_im`。
   - **纯加减法，零乘法，L1 Cache 常驻**（`fold_re/im` 共 12 KB，`acc_power_` 仅 3 KB）。

4. **两级早停与 PAR 计算**：
   - 先执行 `g = 0 .. g_fast - 1`，计算 `peak_val`, `peak_tau`, `mean_val = sum / 750.0`、`par_fast = peak_val / max(mean_val, 1e-12)`。
   - 若 `par_fast < par_fast_thr`：立刻返回 `-1`，推进 `probe_ += hop_samples_`（95%+ 静默窗直接在此退出！）。
   - 否则 `++stats_.stage2_promoted`，继续在**同一 `acc_power_` 上累加** `g = g_fast .. g_full - 1`，重算 `peak_val, peak_tau, par_full`。
   - 若 `par_full < par_thr`：返回 `-1`，推进 `probe_ += hop_samples_`。

5. **轻量级前沿回溯（Leading-Edge Alignment）——确保 `seed` 落在第 1 个 SYNC 附近**：
   - 仅当 Stage 2 确认命中时触发（每秒 ~200 次，总耗时 < 1 ms/2 s）。
   - 已知主扫描窗起点 `win_start = probe_` 和精确循环峰值相位 `peak_tau`。
   - 因为 `hop_symbols = 128`，`win_start` 可能落在 256-SYNC 前导码的任意位置（如第 100 个符号）。若直接返回 `win_start + peak_tau`，后端切窗 `seed - pre_guard (36864)` 会裁掉前导码头部。
   - **回溯算法**（每步回退 `back_step_sym = 16` 符号，最多回退 `max_back_sym = preamble_reps` 符号，且起点不早于 `valid_begin` 与 `holdoff_`）：
     - 在候选起点 `cand_pos = win_start - m * back_step_sym * kExactPeriod`，用 `K_coh = 4, G_probe = 4`（16 符号）仅评估 `peak_tau` 及其邻域 `{(peak_tau - 1 + 750) % 750, peak_tau, (peak_tau + 1) % 750}` 的峰值功率 `p_sig`，并与均匀抽样的 32 个非峰 lag（避开 `peak_tau ± 4`）的平均功率 `p_bg` 求比值 `probe_par = p_sig / max(p_bg, 1e-12)`。
     - 只要 `probe_par >= 3.5`（说明仍在前导码内部），更新 `edge_pos = cand_pos` 并继续往前退一步；一旦跌破 `3.5`（已退到前导码之前的静默底噪），停止回溯。
     - 最终提交 `seed = max<int64_t>(0, edge_pos + peak_tau)`。
     - 设置 `holdoff_ = seed + cfg_.holdoff_len`，并将 `probe_` 推进至 `max(probe_ + hop_samples_, holdoff_)`，防止同一包后半段再次触发。
   - **几何校验**：回溯后的 `seed` 与真实第 1 个 SYNC 符号起点的偏差严格控制在 `[-16, +16]` 符号（约 `±12,000` 样点）以内，远小于 `pre_guard = 36,864` 与后端 `timing_search_margin = 131,072`，既保证完整包住第 1 个 SYNC，又保证完整包住 Payload 末尾。

---

## 5. CLI 与 `openuwb_rx` 会话集成设计

### 5.1 新增 CLI 参数（`apps/live_rx_config.{h,cc}`）

| 参数 | 类型 / 默认值 | 说明与 Fail-Closed 校验 |
| :--- | :--- | :--- |
| `--acq-mode` | `autocorr \| xcorr`（默认 `autocorr`） | 选择流式扫描前端；非法值启动即抛异常 |
| `--acq-xcorr-k-coh` | `size_t`（默认 `4`） | 相干叠加符号数，校验 `1 <= k_coh <= 32` |
| `--acq-xcorr-g-fast` | `size_t`（默认 `6`） | Stage-1 早停组数，校验 `1 <= g_fast <= g_full` |
| `--acq-xcorr-g-full` | `size_t`（默认 `16`） | Stage-2 总组数，校验 `g_full <= 64` |
| `--acq-xcorr-hop` | `size_t`（默认 `128`，`0` 表示自动） | 步进符号数，校验 `hop <= 256` |
| `--acq-xcorr-par-fast` | `double`（默认 `4.5`） | Stage-1 PAR 门限，校验 `1.0 < par_fast <= par` |
| `--acq-xcorr-par` | `double`（默认 `6.5`） | Stage-2 PAR 门限，校验 `par > 1.0` 且有限 |

### 5.2 会话接入（`apps/live_rx_session.cc`）

1. 在 `run_live_rx_session` 中，根据 `o.acq_mode` 初始化 `Acquisition acq` 或 `PreambleXcorrAcquisition xcorr_acq`（`code_index = lock_cfg.code_index` 即 `10`，`preamble_reps = lock_cfg.preamble_repetitions` 即 `256`）。
2. 在 gap reset 与 `autocorr_gate`（`--autocorr-when-locked suspend|auto|scan`）处统一调用当前激活前端的 `reset(...)` 与 `scan(...)`。
   - 注意：`autocorr_lookback` 默认值 `-1` 在 `xcorr` 模式下同样生效（`required_window_samples() + window`）。
3. 扫描循环：
   - `autocorr` 模式：保持现有 `while (acq.probe + acq.step - 1 + Acquisition::lag < end)` 零改动。
   - `xcorr` 模式：`while (xcorr_acq.probe() + xcorr_acq.required_window_samples() <= end)`，命中 `seed >= 0` 后调用同一个 `submissions.handle_autocorr_seed(seed, metric, valid_begin, cur_pre_guard, window, o.interference, id, lock, report)`。
4. 日志输出：
   - `START` 行追加 `acq_mode=<autocorr|xcorr>`（当 `acq_mode == "xcorr"` 时附带 `acq_xcorr_k_coh=... acq_xcorr_g_fast=... acq_xcorr_g_full=... acq_xcorr_hop=... acq_xcorr_par=...`，当 `acq_mode == "autocorr"` 时尽量保持既有字段不破坏现有日志解析脚本）。
   - `COST_TOTAL` 行在 `xcorr` 模式下追加 `xcorr_windows=... xcorr_stage2=...` 以便观测早停效率。

---

## 6. 多 Subagent 分工与交付清单

建议 `opencode` 按以下 4 个 Subagent 分工推进（Subagent-1 与 Subagent-2 可在确定头文件契约后并行，Subagent-3 依赖 Subagent-1，Subagent-4 依赖全部合入）：

### Subagent-1：核心算法与 SIMD 实现（`openuwb_runtime`）
- **产物**：
  - 新增 `include/openuwb/runtime/preamble_xcorr_acquisition.h`
  - 新增 `lib/runtime/uwb_preamble_xcorr_acquisition.cc`
  - 修改 `lib/runtime/CMakeLists.txt` 注册新 `.cc`
- **要求**：
  - 完整实现 `PreambleXcorrAcquisition`（支持 Code 9/10/11/12，精确 `kExactPeriod` 折叠、线性双倍缓冲、AVX2+FMA 8-wide 零乘法稀疏相关 + 标量后备、两级早停、前沿回溯、holdoff、ring wrap 边界安全、gap reset）。
  - 热路径零堆内存分配（所有工作缓冲在 `configure()` 中一次性分配）。
  - 核查 `dev-avx2`、`scalar-debug` 均无编译警告。

### Subagent-2：CLI 与 `openuwb_rx` 会话接入
- **产物**：
  - 修改 `apps/live_rx_config.h`、`apps/live_rx_config.cc`（新增 `--acq-mode` 与 `--acq-xcorr-*` 参数解析、help 文本、fail-closed 校验）
  - 修改 `apps/live_rx_session.cc`（接入 `PreambleXcorrAcquisition`、`AutocorrGate`、gap reset、`START`/`COST_TOTAL` 诊断字段）
- **要求**：
  - `--acq-mode autocorr` 默认行为与 stdout 协议保持完全兼容。
  - CLI 校验覆盖所有非法值并 fail-closed。

### Subagent-3：单元测试 QA（`tests/unit/runtime/`）
- **产物**：
  - 新增 `tests/unit/runtime/qa_uwb_preamble_xcorr_acquisition.cc`
  - 修改 `tests/CMakeLists.txt` 注册 `qa_uwb_preamble_xcorr_acquisition`（纳入 `rx_sic` 标签）
- **必测用例**：
  1. **纯白噪声与强单音/多音虚警测试**：注入纯高斯噪声 + 强 CW 单音（含 `-232.96 MHz`、`-124.4 MHz`）扫描至少 50 ms 等效样点，断言 `seeds == 0`。
  2. **大 CFO 灵敏度测试（核心！）**：使用 `testdata/reference_preamble_code10_998p4.cf32` 重采样到 `737.28 MS/s` 合成 256-SYNC 前导码，在 **CFO = 0 / ±50 kHz / ±100 kHz**、**SNR = -12 dB 与 -14 dB** 下，断言 100% 检出且无重复误报。
  3. **前沿回溯定位精度测试**：将真实前导码起点分别放在扫描周期 `hop_symbols = 128` 的 4 个不同相位（0%、25%、50%、75% 偏移），断言提交的 `seed` 与真实 `preamble_start` 误差 `|seed - true_start| <= 16 * 751` 样点（远小于 `pre_guard = 36864`）。
  4. **异码拒绝测试**：注入高 SNR Code 9 前导码，用 Code 10 配置扫描，断言 `seeds == 0`。
  5. **边界与重置测试**：覆盖 ring buffer 回绕边界（`probe` 跨越 `ring.size()`）、`holdoff` 抑制二次触发、`reset()` 清空状态。
  6. **AVX2 与标量一致性测试**：对比 AVX2 内核与标量参考实现在相同输入下的 `acc_power_` 与 `PAR`，最大相对误差 `< 1e-5`。

### Subagent-4：7 场景离线回放 ABBA 验收与跨预设门禁
- **编译与 QA 门禁**（全部 `env -u LD_LIBRARY_PATH`）：
  1. `dev-avx2` 构建 + `ctest -L rx_sic --output-on-failure -j`
  2. `native-avx2`、`no-uhd` 全量构建 + 全量 `ctest`
  3. `scalar-debug` 全量构建 + 运行 `uwb_qa_uwb_preamble_xcorr_acquisition` 与 `rx_sic`
  4. `native-avx512` 全量编译检查（无 AVX-512 机器仅编译并核验无错误）
- **7 场景 2.0 s 离线回放对照（串行、同环境、ABBA 或明确 A/B 配对）**：
  - 使用同一套 `openuwb_rx` 二进制，对比：
    - **基线 A**：`--taps ch250_eq16 --acq-mode autocorr`（以及默认 `--taps eq16 --acq-mode autocorr`）
    - **新方案 B**：`--taps ch250_eq16 --acq-mode xcorr`
  - 7 个标准离线回放场景（位于 `testdata/dw3000_*_20261123/raw_2s.sc16`，如目录存在全部跑完）：
    1. `clean_ch0`
    2. `clean_near_ch0`
    3. `clean_comm884cm_att9p3_try2`
    4. `jammer9_70cm_comm506cm_att9p3`
    5. `jammer9_70cm_comm348cm_att9p3`
    6. `jammer9_70cm_comm70cm_att9p3`
    7. `jammer9_...`（仓库内现有第 7 个 20261123 场景）
  - **额外高 CFO 离线注入验证**：用 `tools/` 或临时脚本给 `clean_ch0` / `try2` 乘以 `+100 kHz` 频偏因子生成临时文件到 `/tmp/opencode/xcorr_ab/`，验证 `--acq-mode xcorr` 在实测信号 +100 kHz CFO 下仍能稳定首包进锁。
- **核心验收指标**：
  1. **`acq_seeds` 断崖式下降**：从基线 `2,237–2,284` 降至接近真实包数（`~400` 左右，或 `--autocorr-when-locked suspend` 时仅首锁所需的 `1–3` 个）。
  2. **`acquire` 阶段虚警清零、首锁大幅提速**：7 场景 `acquire` 阶段 `TIMING_FAILED` 趋近于 `0`，首锁延迟从 `8–45 ms` 降至 `< 10 ms`（多数在第一个 5 ms 周期内完成）。
  3. **总 FCS 成功率提升且零回退**：所有 7 场景 `fcs_ok` 均 `>=` 基线 A，且 `dropped_queue == 0`、`reacquire_entries == 0`。
  4. **前端 CPU 耗时显著下降**：`frontend_scan_us` 从基线 `~330–360 ms / 2.0 s` 降至 **`< 200 ms / 2.0 s`**（目标 `~100–150 ms`）。
- **交付报告**：
  - 将完整验收数据与日志摘要写入 `docs/performance/Acquire折叠互相关前端验收报告_2026-11-24.md`（记录 commit、CPU/ISA、构建预设、7 场景 A/B 对照表、`acq_seeds`、首锁时间、`acquire` 失败数、`fcs_ok`、`frontend_scan_us`、以及 `+100 kHz` CFO 注入实测结果），所有临时日志放 `/tmp/opencode/xcorr_ab/`。

---

## 7. 给 `opencode` 的直接执行提示词

```text
请严格按照 docs/OpenCode执行方案_Acquire前端折叠前导码互相关设计与验收_2026-11-24.md 的设计与硬约束，组织多 Subagent 完成以下工作：
1. 实现 openuwb::live::PreambleXcorrAcquisition（K_coh=4, G_fast=6, G_full=16, hop=128，±100 kHz CFO 鲁棒，64-tap 稀疏三元零乘法 AVX2+Scalar 内核，带轻量级前沿回溯与 holdoff），放入 include/openuwb/runtime/preamble_xcorr_acquisition.h 与 lib/runtime/uwb_preamble_xcorr_acquisition.cc。
2. 在 apps/live_rx_config.{h,cc} 与 apps/live_rx_session.cc 中接入 --acq-mode autocorr|xcorr（默认保持 autocorr 兼容不变）及 --acq-xcorr-* 参数与日志统计。
3. 编写 tests/unit/runtime/qa_uwb_preamble_xcorr_acquisition.cc 覆盖纯噪/强单音零虚警、0/±50/±100 kHz CFO @ -14 dB SNR 灵敏度、前沿回溯定位精度、异码拒绝、ring 边界与 AVX2/Scalar 一致性。
4. 完成 dev-avx2 / native-avx2 / no-uhd / scalar-debug / native-avx512 门禁检查，以及 7 场景 2.0 s 离线回放 A/B 验收与 +100 kHz CFO 实测注入验证，输出验收报告到 docs/performance/Acquire折叠互相关前端验收报告_2026-11-24.md。

注意：遵守 AGENTS.md（始终 env -u LD_LIBRARY_PATH，禁止 -march=native，不触碰 TWR/ToA/Radar，不提交 testdata/dw3000_* 大文件，不 push 远端）。
```
