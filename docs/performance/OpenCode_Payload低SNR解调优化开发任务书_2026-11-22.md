# OpenCode 开发任务书：Payload 低 SNR 解调精度优化（多 Subagent 协作版）

- **日期**：2026-11-22
- **适用分支**：`feature/openuwb-rx-sic`
- **目标模块**：`openuwb_phy_hrp`（Payload / PHR / FEC / Fused CMF）
- **前置背景**：
  - 仓库已具备 BPM-BPSK 的 M1–M4 增强能力（`rx_bpm_known_hop`、`rx_bpm_joint_ml`、`rx_bpm_terminated_viterbi`、`rx_phr_phase_refine`、`rx_fused_sfo_track`，见 `include/openuwb/phy/hrp/common/profile.h` 与 `apps/bpm_ablation_env.h`）。
  - 本次优化针对**低 SNR 下 Payload 解调精度（FCS 通过率）**，在现有 Fused CMF + Bi-orthogonal ML + RS 擦除框架上补齐 4 个高性价比短板，同时严格保证高 SNR 路径零额外热路径开销、不破坏既有逐位一致性与 fail-open 语义。

---

## 一、全局开发约束与测试边界（必读）

### 1.1 代码与架构约束
1. **头文件分层**：重型非模板实现必须写入 `lib/phy/hrp/*.cc`，不得把大段实现写回 `include/openuwb/phy/hrp/rx/*.h`。
2. **禁止 `-march=native`**：保持 `UwbSimd.cmake` 编译期 ISA 一致性。
3. **不破坏高 SNR 热路径**：所有新增救援/重试逻辑仅在**初次译码/FCS 失败时**触发；高 SNR 一次通过的包不增加耗时。
4. **Fail-Open / 保底原则**：任何新增优化（如 CIR 去噪、PHR 相位解耦）必须保证不降低无噪/高 SNR 黄金参考包的解调结果；出现异常（如非有限值、抽头全零）时自动回退原始输入。
5. **不越界修改**：不触碰 `stage_timing`、`stage_cfo`、`stage_sfd`、SIC、TWR、Sensing、UHD/Radio 代码。

### 1.2 测试边界：**禁止过度测试！**
> **硬性要求**：每个 Subagent 以及最终集成验收，**只构建并运行与 Payload/PHR/FEC 最相关的 5 个单元测试**，**严禁**执行以下重型或无关操作：
> - ❌ **禁止**跑全量 `ctest`（无需跑 SIC、TWR、Sensing、Resampler、UHD、Console 等测试）。
> - ❌ **禁止**跑 `-L heavy` 重型网格扫描测试。
> - ❌ **禁止**切多个 CMake preset（`native-avx512`、`no-uhd`、`scalar-debug`、`install_consumer` 等）做全套门禁。
> - ❌ **禁止**做大规模数百种子慢速蒙特卡洛扫描（单测内低 SNR 验证控制在 20–30 个固定种子、单测总耗时 < 1 秒）。

**统一快速验证命令（耗时仅数秒）**：
```sh
env -u LD_LIBRARY_PATH cmake --build --preset dev-avx2 -j
env -u LD_LIBRARY_PATH ctest --test-dir build/dev-avx2 \
  -R "uwb_qa_uwb_demod_core_fec_joint_ml|uwb_qa_uwb_rx_fused_cmf|uwb_qa_uwb_rx_fused_phr|uwb_qa_uwb_rx_phase_sfo|uwb_qa_uwb_demod_core_cir_pipeline" \
  --output-on-failure -j
```
（若尚未配置 `build/dev-avx2`，仅需先执行一次 `env -u LD_LIBRARY_PATH cmake --preset dev-avx2`。）

---

## 二、多 Subagent 工作包拆解与文件所有权

为了支持 OpenCode 并行派发多个 Subagent 且最小化合并冲突，拆分为 **3 个核心开发包（WP1–WP3）** 与 **1 个集成/消融验证包（WP4）**：

| 工作包 | 核心内容 | 主要修改文件（所有权） | 并行关系 |
| :--- | :--- | :--- | :--- |
| **WP1** | FEC 深度救援（`v=8` 擦除 + Trellis-Chase 联合擦除）+ PHR/Payload 独立裁幅 | `lib/phy/hrp/uwb_demod_fec.cc`<br>`include/openuwb/phy/hrp/rx/demod_fec.h`<br>`tests/unit/phy/qa_uwb_demod_core_fec_joint_ml.cc` | **完全独立**，可与 WP2/WP3 并行 |
| **WP2** | CIR 抽头自适应噪声底收缩（Denoise），提升匹配滤波 SNR | `include/openuwb/phy/hrp/common/profile.h`（仅加 1 个开关）<br>`include/openuwb/phy/hrp/rx/rx_fused_cmf.h`<br>`lib/phy/hrp/uwb_demod_softchips_payload.cc`（仅 `denoise_cir_taps` 辅助函数与模板构建入口）<br>`tests/unit/phy/qa_uwb_rx_fused_cmf.cc` | 与 WP1 独立；与 WP3 错开函数区段 |
| **WP3** | Fused Payload 零成本 `-pol` 重试 + PHR 导频相位/频偏精修门限解耦 | `lib/phy/hrp/uwb_demod_softchips_payload.cc`（仅 `stage_phr_fused_cmf` 与 `stage_payload_fused_cmf` 内部）<br>`tests/unit/phy/qa_uwb_rx_fused_phr.cc`<br>`tests/unit/phy/qa_uwb_rx_phase_sfo.cc` | 与 WP1 独立；建议在 WP2 之后合入同一 `.cc` 文件 |
| **WP4** | 环境变量接入 + 轻量级低 SNR Payload A/B 对比验收 | `apps/bpm_ablation_env.h`<br>`tests/unit/phy/qa_uwb_demod_core_fec_joint_ml.cc`（追加 1 个轻量 A/B 单测） | 在 WP1–WP3 完成后执行 |

---

## 三、各工作包详细实施方案

### WP1：FEC 深度救援与 PHR/Payload 独立裁幅

**负责文件**：
- `lib/phy/hrp/uwb_demod_fec.cc`
- `include/openuwb/phy/hrp/rx/demod_fec.h`
- `tests/unit/phy/qa_uwb_demod_core_fec_joint_ml.cc`

**具体改动点**：
1. **RS(63,55) 补全满擦除 `v = 8`**：
   - 在 `rescue_payload_biorth()`（`lib/phy/hrp/uwb_demod_fec.cc` 约 L1190）以及 `rescue_payload_soft()`（约 L709）的 Phase 1 擦除循环中：
     - 将 `for (int v : { 2, 4, 6 })` 扩展为 `for (int v : { 2, 4, 6, 8 })`；
     - 对应栈数组 `int ep[6]` 扩大为 `int ep[8]`。
   - **依据**：RS(63,55) 最小距离 $d_{\min} = 9$，满足 $2t + v \le 8$，支持最多 8 个纯擦除（0 额外错误）。

2. **Trellis-Chase 局部翻转 + RS 擦除联合救援**：
   - 在 `rescue_payload_biorth()` Phase 2（约 L1288–L1319）的 1-node 和 2-node 局部网格翻转循环中：
     - 保持先尝试 `rs_decode_block(blk, d330.data())`（快路径）；
     - 若 `rs_decode_block` 失败且 `bc[b].cands.size() < 4`，复用 Phase 1 已经排序好的弱符号表 `srel`（或跳过被翻转节点覆盖的符号后取前 4 个最弱符号），追加一次 `rs_decode_block_erasures(blk, ep, 4, d330.data())`（也可按 `{2, 4}` 顺序，命中即停）。
   - **性能保护**：仅在 Phase 1 失败且翻转后的普通 RS 失败时运行，对高 SNR 包完全不可见。

3. **PHR 与 Payload 分开裁幅（修复 8x 幅度尺度混叠）**：
   - 在 `decode_payload_biorth_soft()`（`lib/phy/hrp/uwb_demod_fec.cc` 约 L1360）中：
     - 当前代码对拼接后的 `m0j[0..total-1]`（前 21 个为 64-chip PHR，后 `payload_symbols` 个为 8-chip Payload）统一调用 `clip_biorth_metrics(m0j.data(), m1j.data(), total)`。
     - 改为**分别裁幅**：
       ```cpp
       if (enable_clip) {
           clip_biorth_metrics(m0j.data(), m1j.data(), 21);
           clip_biorth_metrics(m0j.data() + 21, m1j.data() + 21, payload_symbols);
       }
       ```
     - 同时保证 PHR 与 Payload 各自归一化到同一噪声标准差尺度（或仅独立裁幅不改变两者信噪比权重，先做独立裁幅并观察单测结果）。

**WP1 单元测试要求**（加在 `qa_uwb_demod_core_fec_joint_ml.cc`）：
- 构造 7~8 个弱 6-bit RS 符号被破坏的码字，验证 `v = 8` 满擦除能成功恢复并通过 FCS。
- 构造需要“网格翻转 + RS 擦除”共同作用才能恢复的码字，验证 Phase 2 联合救援有效。
- 构造 PHR 幅值比 Payload 大 8 倍（模拟真实 64-chip vs 8-chip）的输入，验证独立裁幅不会误裁 PHR 或抬偏 Payload 裁幅上限。
- 仅运行：`uwb_qa_uwb_demod_core_fec_joint_ml.cc`。

---

### WP2：CIR 抽头自适应噪声底收缩（Denoise）

**负责文件**：
- `include/openuwb/phy/hrp/common/profile.h`
- `include/openuwb/phy/hrp/rx/rx_fused_cmf.h`
- `lib/phy/hrp/uwb_demod_softchips_payload.cc`
- `tests/unit/phy/qa_uwb_rx_fused_cmf.cc`

**具体改动点**：
1. **Profile 开关**（`include/openuwb/phy/hrp/common/profile.h`）：
   - 在 `Qm35825Profile` 的 BPM 算法开关区新增：
     ```cpp
     bool rx_cir_tap_denoise = false; // adaptive noise-floor shrinkage for fused CIR taps
     ```
   - 保持默认 `false`，确保既有逐位一致性 QA 零变化。

2. **实现自适应抽头去噪辅助函数**（声明放 `rx_fused_cmf.h`，实现放 `lib/phy/hrp/uwb_demod_softchips_payload.cc`）：
   - 接口建议：
     ```cpp
     // Copies `src[0..tap_count)` into `dst[0..tap_count)` and zeroes noise-only
     // taps below an adaptive floor estimated from the weakest half of taps,
     // always retaining at least `min_keep` strongest taps.
     void denoise_cir_taps(const std::complex<float>* src,
                           size_t tap_count,
                           std::complex<float>* dst,
                           float sigma_mult = 6.0f,
                           size_t min_keep = 4);
     ```
   - **算法细节**（全栈上计算，零堆分配，`tap_count <= 64`）：
     1. 计算每个抽头功率 $p_q = |h_q|^2$；
     2. 将 $p_q$ 拷贝到栈数组 `float sorted_p[64]` 并排序；
     3. 用能量最低的 **50% 抽头**（`sorted_p[0 .. tap_count/2 - 1]`）的中位数或均值作为噪声功率底估计 $\hat{\sigma}_n^2$（对于 38 taps 即取最低 19 个抽头的均值，避开真实多径主峰与反射径）；
     4. 设定功率门限 $T = \text{sigma\_mult} \cdot \hat{\sigma}_n^2$（默认 `sigma_mult = 6.0f`，约高于噪声底 7.8 dB）；
     5. 对 $p_q < T$ 的抽头置零，但强制保留能量最高的至少 `min_keep = 4` 个抽头（防止极高噪或单径场景下误删主径）；
     6. 若输入含非有限值或全零，直接原样拷贝（fail-open）。

3. **接入 Fused Payload / PHR 路径**：
   - 在 `stage_payload_fused_cmf()`（及 `stage_phr_fused_cmf()`）中，当 `profile.rx_cir_tap_denoise` 为 `true` 时，用栈数组 `std::complex<float> denoised_taps[rxfusion::kMaxTapCount]` 接收 `denoise_cir_taps(...)` 结果，再传给 `build_fractional_cir_bank` / `correlate_payload_known_hop_fused_cfo` / `correlate_payload_fused_cfo`。
   - **关键**：不修改 `cir.cir_filter_taps` 原值（不影响诊断日志与 legacy 回退），且 `tap_count` 仍保持 38，完全复用 `dot2_shared_fixed52` AVX2 快路径。

**WP2 单元测试要求**（加在 `qa_uwb_rx_fused_cmf.cc`）：
- 构造 3-tap 真实信道 + 35-tap 高斯白噪声抽头的 `cir_filter_taps`，验证 `denoise_cir_taps` 精准保留 3 个真实径、清零绝大部分噪声抽头，且相关性输出 SNR 优于未去噪的 38-tap 模板。
- 验证无噪单径与全零边界输入均 fail-open 不崩溃。
- 仅运行：`uwb_qa_uwb_rx_fused_cmf.cc`。

---

### WP3：Fused Payload 零成本 `-pol` 重试 + PHR 相位/频偏精修解耦

**负责文件**：
- `lib/phy/hrp/uwb_demod_softchips_payload.cc`
- `tests/unit/phy/qa_uwb_rx_fused_phr.cc`
- `tests/unit/phy/qa_uwb_rx_phase_sfo.cc`

**具体改动点**：
1. **Fused Payload 零成本负极性（`-pol`）重试**：
   - 在 `stage_payload_fused_cmf()`（`lib/phy/hrp/uwb_demod_softchips_payload.cc` 约 L2116–L2150）中：
     - 当 Attempt 1（硬判决）与 Attempt 2（正极性 `decode_payload_biorth_soft`）均未通过 FCS 时，在回退到 4-候选/整帧重算之前，**不重做任何复数相关**，直接将已计算的 `bm0[s], bm1[s]`（以及 PHR 指标 `m0p, m1p`）取反到已有 scratch 缓存（如 `scratch.fused.biorth_m0` / `biorth_m1`），再调用一次 `decode_payload_biorth_soft(..., terminated, /*enable_clip=*/true)`。
     - 若通过 FCS 则直接 `return true`。
   - **收益**：在低 SNR 下如果 `ns_sfd.polarity` 偶发翻号，耗时增加仅 ~10 µs（而原本会触发 ~15 ms 的全流重算且丢失 M2 增益）。

2. **PHR 导频相位精修与 CFO 频偏外推解耦**：
   - 位置：`stage_phr_fused_cmf()`（`lib/phy/hrp/uwb_demod_softchips_payload.cc` 约 L1566–L1639）。
   - **现状问题**：当前 `prod_gate = 0.997` 同时控制**常数相偏修正（`rot_at_payload_start`）**与**频偏斜率外推（`delta_cfo_w_rad`）**。低 SNR 下 `prod` 落到 `0.90 ~ 0.996` 时两者一起被禁用。
   - **改动**：
     - 保持 `cfo_prod_gate = 0.997`（支持 `UWB_BPM_REFINE_PROD` 覆盖）用于保护风险较高的**频偏斜率外推**；
     - 新增较宽松的相位门限 `phase_prod_gate = 0.90`（支持 `UWB_BPM_REFINE_PHASE_PROD` 覆盖）：
       - 当 `prod >= cfo_prod_gate`：保持原逻辑，同时应用相位 + 频偏修正；
       - 当 `phase_prod_gate <= prod < cfo_prod_gate`：清零频偏外推（`scratch.phr_refine.delta_cfo_w_rad = 0.0`），仅用 21 个 PHR 符号的**均值相偏**更新 `rot_at_payload_start`（`applied = true`），消除常数残余相偏而不放大频偏噪声；
       - 当 `prod < phase_prod_gate`：全部不应用（`applied = false`）。

**WP3 单元测试要求**（加在 `qa_uwb_rx_fused_phr.cc` / `qa_uwb_rx_phase_sfo.cc`）：
- 构造传入错误 `ns_sfd.polarity = -true_pol` 的 Fused Payload 场景，验证 `stage_payload_fused_cmf` 通过零成本取反直接解出正确 FCS（不触发 fallback）。
- 构造中等噪声下 `0.90 <= prod < 0.997` 且存在 $25^\circ$ 常数相偏的 21-symbol PHR，验证常数相位修正生效且 `delta_cfo_w_rad == 0.0`。
- 仅运行：`uwb_qa_uwb_rx_fused_phr.cc` 与 `uwb_qa_uwb_rx_phase_sfo.cc`。

---

### WP4：环境变量接入与轻量级集成验收

**负责文件**：
- `apps/bpm_ablation_env.h`
- `tests/unit/phy/qa_uwb_demod_core_fec_joint_ml.cc`

**具体改动点**：
1. 在 `apps/bpm_ablation_env.h` 中接入 `UWB_BPM_CIR_DENOISE` 环境变量（`p.rx_cir_tap_denoise = env_flag("UWB_BPM_CIR_DENOISE", level >= 2);` 或默认随 `level >= 3` 开启），方便后续快速做消融对照。
2. 在 `qa_uwb_demod_core_fec_joint_ml.cc` 末尾增加 1 个轻量级低 SNR 对比单测（固定 20 个种子，单测耗时 < 0.5 s），验证开启本次优化后 FCS 通过率不低于基线（并在低 SNR 点有净提升或零退化）。
3. 跑完第 1.2 节指定的 **5 个 Payload 核心测试**，全部通过即完成交付。

---

## 四、OpenCode 调度建议

- **第一波（并行）**：同时派发 **WP1** 与 **WP2**（两者几乎无文件重叠）。
- **第二波**：WP2 完成后派发 **WP3**（避免与 WP2 在 `uwb_demod_softchips_payload.cc` 中产生合并冲突）。
- **第三波**：WP1 + WP2 + WP3 全部合入后，派发 **WP4** 做环境变量接线与 5 个核心测试的最终验证。
- **验收标准**：上述 5 个 Payload 相关 CTest 100% 通过，不新增编译警告，不触发任何重型或无关测试套件。
