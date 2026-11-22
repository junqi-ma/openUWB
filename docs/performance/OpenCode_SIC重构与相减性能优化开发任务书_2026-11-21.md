# OpenCode 开发任务书：SIC 重构与相减阶段性能优化及问题修复

- **日期**：2026-11-21
- **适用分支**：`feature/openuwb-rx-sic`（基于当前工作区，**保留** `include/openuwb/sic/sic_unified.h` 中已修改的 `lag_lo = -4`）
- **目标范围**：仅聚焦 SIC 后半段（重构 `reconstruct`、相减 `subtract`、ROI 能量统计 `compute_roi_energies` 及 `PacketProcessor` 数据流衔接），不改动硬件、TWR 或前端定时/解调算法逻辑。
- **架构边界**（严格遵守 `AGENTS.md`）：
  - 单向依赖：`support → dsp → phy_hrp → {sic, sensing, twr_core}`，`runtime → phy_hrp + sic + sensing + radio`。
  - 重型非模板实现一律放 `lib/<domain>/*.cc`，不得写回公共头文件。
  - 严禁 `-march=native`，严格遵守 `UwbSimd.cmake` 编译期 ISA 分层；若新增 `uwb_simd_apply()` QA 目标须同步加入 `tests/CMakeLists.txt` 的 `UWB_QA_PCH_SIMD_TARGETS`。

---

## 1. 背景与性能基线

在主机 `i7-12700`（AVX2+FMA，无 AVX-512）`dev-avx2` 构建下，主路径 `Ipatov + SparsePulse` 实测基线如下：

| 场景 | `core_total` (p50) | `tx_remod` | `reconstruct` | `subtract` | `roi_energies` | 后半段合计及占比 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| Live Hit（`n=506k`, `span=198.6k`, L2, `sparse=on`） | 799.8 us | 10.6 us | **176.2 us** | **150.0 us** | **68.1 us** | **394.3 us (49.3%)** |
| 短窗（`n=135k`, `span=133.6k`, L2） | 568.0 us | 9.6 us | **121.1 us** | **59.5 us** | **45.9 us** | **226.5 us (39.9%)** |
| 短窗（`n=135k`, `span=73.2k`, L1 PreambleOnly） | 454.2 us | 0.0 us | **78.1 us** | **50.9 us** | **25.2 us** | **154.2 us (33.9%)** |

### 核心瓶颈与现存问题清单

1. **内存带宽浪费（最大瓶颈）**：506k 场景下后半段对大数组进行了 5 遍独立读写（~16 MB 内存流量，远超单核 L2）：`recon` 1.59 MB 清零 → 脉冲累加 → `out != y` 导致 307k 非 ROI 样点（4.92 MB）拷贝 → `vector_sub3_complex` 4.77 MB 读写 → `compute_roi_energies` 3.18 MB 二次重读。
2. **`vector_sub3_complex` 缺少 AVX2 实现**：`lib/phy/hrp/uwb_demod_timing.cc:1828-1843` 只有 AVX-512 分支，缺少 AVX2 实现与 `__restrict`。
3. **热循环微架构低效 & 20-tap 硬编码**：`accumulate_pulses_cfo_20tap_avx2`（`lib/sic/uwb_discrete_cir_reconstruct.cc:191-264`）内存在 `std::complex<float>::operator*`（GCC 降级为 `__mulsc3` 调用）、每脉冲 `sign_mask` 向量乘法、有符号 `/ 1024` 与 `% 1024`；且第 809 行硬编码 `model.tap_count == 20`，其他长度会静默退化为标量 `noinline` `fused_axpy`；同时未接入 AVX-512 生产核。
4. **精度缺陷（Tap 内 CFO 相位缺失）**：稀疏脉冲路径只乘了脉冲起点相位 $e^{j\omega m_i}$，漏掉了 tap 内部偏移 $l = \text{lag\_lo} + d$ 对应的 $e^{j\omega l}$，在 80–100 kHz 大频偏下会将抑制上限压低至 ~40 dB。
5. **Fail-Closed 与拷贝契约缺口**：
   - 稀疏路径跳过了输出有限性检查（`uwb_discrete_cir_reconstruct.cc:854-857`）；
   - 底层 `fail()` 无视 `copy_out_on_bypass = false`，失败时仍无条件拷贝 4 MB（若上层为 `true` 则双拷贝）；
   - L1 稀疏路径通过 `apply_discrete_cir_pulses` 传入 `transaction_prepared = false`，每帧重复校验已缓存的 `cached_preamble_events` 并标量全窗重扫 `y` 的 NaN/Inf。
6. **边界截断刚性**：`lo < 0 || hi > n` 直接整包返回 `context_truncated`，滑窗边缘干扰包（尤其是 `lag_lo = -4` 时 `q ∈ [0, 3]` 或包尾超窗）既不裁剪也不自动降级到 L1。
7. **计时与配置缺陷**：
   - `uwb_sic_unified.cc:738-741` L2 分支漏填 `recon_scale_us` / `recon_conv_us` / `recon_cfo_rot_us`，且 `zero_us` 未暴露；
   - `fullpacket_cancel_discrete`（`uwb_discrete_cir_reconstruct.cc:937`）对 Ipatov 路径仍硬校验 `profile.pulse.shape == Legacy`；
   - `SicUnifiedConfig::to_sic_config()`（`uwb_sic_unified.cc:174-197`）漏设 `validate_threshold`，且 `PathFractional` L1 计时错位（`uwb_sic_unified.cc:608-609`）。

---

## 2. 量化优化目标

在 `i7-12700`（AVX2+FMA）单核同环境 ABBA 对比下：

| 指标 | 优化前基线 | 目标值 |
| :--- | ---: | :--- |
| Live Hit 506k 后半段合计（`recon + sub + roi_e`） | ~394 us | **≤ 150 us**（降幅 ≥ 60%） |
| Live Hit 506k `core_total` (p50, `sparse=on`) | ~800 us | **≤ 560 us** |
| 短窗 135k L2 后半段合计 | ~226 us | **≤ 100 us** |
| 100 kHz CFO 下 `Ipatov + SparsePulse` 全包抑制比 | ~40 dB（受 tap 相位近似限制） | **与 NonSparse 一致（≥ 52 dB @ 20-tap）** |
| 24-tap / 32-tap 稀疏重构耗时 | 退化标量（3–5× 慢） | **与 20-tap 同量级线性缩放（≤ 1.6×）** |
| 所有既有 `rx_sic` / 全量 `ctest` 回归 | 100% pass | **100% pass，零内存分配，严格 fail-closed** |

---

## 3. 工作包拆解与严格文件所有权（面向多 Subagent 并行）

为避免并行编辑冲突，将开发拆分为 **Wave 0（接口冻结） → Wave 1（3 个独立 Subagent 并行） → Wave 2（集成联调与基准验证）**。

```text
Wave 0 (Coordinator 串行，~5 min)
  └─ 冻结公共头文件结构体扩展与函数签名，确保所有 Wave 1 并行编译契约一致

Wave 1 (3 Subagents 严格按文件所有权并行)
  ├─ Subagent A (WP1): DSP/PHY 向量相减 + 融合相减能量统计核（AVX2 & AVX-512）
  ├─ Subagent B (WP2): SIC 稀疏/非稀疏重构核心 + L1D 分块流水 + 精确 CFO 预旋 + 任意 tap_count
  └─ Subagent C (WP3): Unified SIC 编排层 + PacketProcessor 原地相减 + 计时/配置/降级修复

Wave 2 (Subagent D / Coordinator)
  └─ 跨模块集成修补、基准工具与 ABBA 性能实测、全预设门禁验收（含 install_consumer）
```

> **重要**：若 WP3 依赖 WP1/WP2 的新接口声明，则必须在 **Wave 0** 先一次性提交公共头文件的最小声明改动（或先完成 WP1+WP2 头文件再并行 WP1/WP2/WP3 的 `.cc` 与单测实现），彻底消除写冲突。

---

## 4. 各工作包详细技术规格

### Wave 0：公共接口冻结（仅头文件声明，由 Coordinator 先完成）

涉及文件：
- `include/openuwb/phy/hrp/rx/rx_sync_opt.h`
- `include/openuwb/dsp/demod_sic_avx512.h`
- `include/openuwb/sic/discrete_cir_sic.h`
- `include/openuwb/sic/sic_unified.h`

#### 0.1 新增向量相减与融合能量统计声明
1. 在 `include/openuwb/dsp/demod_sic_avx512.h` 的 `openuwb::core::simd::avx512` 命名空间新增：
   ```cpp
   // Computes dst[i] = a[i] - b[i] (supports dst == a in-place) while
   // accumulating double-precision energies e_a = sum |a[i]|^2 and
   // e_dst = sum |dst[i]|^2.  Writes *e_a_out and *e_dst_out (additive or
   // overwrite: overwrite on entry).
   void sub3_energy_cf32(gr_complex* dst,
                         const gr_complex* a,
                         const gr_complex* b,
                         size_t count,
                         double* e_a_out,
                         double* e_dst_out);

   // 20-tap sparse pulse accumulation (zero-CFO & CFO-folded).
   void accumulate_pulses_zero_cfo_20tap_cf32(
       gr_complex* recon_base,
       const gr_complex* taps_ptr,
       const uint32_t* sample_offsets,
       const int8_t* signs,
       size_t count);
   ```
   *(注：为避免 `dsp` 反向依赖 `phy_hrp` 的 `mod::PulseEvent`，AVX-512 脉冲核可选择直接放在 `lib/sic/uwb_discrete_cir_reconstruct.cc` 内部用 `#if UWB_SIMD_HAVE_AVX512` + `__attribute__((target("avx512f,avx512dq,avx512vl")))` 实现，无需跨库暴露 `PulseEvent`；而 `sub3_energy_cf32` 只用 `gr_complex`，放在 `openuwb_dsp` 非常干净，符合单向依赖！)*

2. 在 `include/openuwb/phy/hrp/rx/rx_sync_opt.h` 的 `openuwb::demod::rxopt` 命名空间新增：
   ```cpp
   // Element-wise dst[i] = a[i] - b[i] (supports dst == a in-place) and
   // computes double-precision energies sum |a[i]|^2 and sum |dst[i]|^2
   // in a single pass.
   UWB_API void vector_sub3_energy_complex(std::complex<float>* dst,
                                           const std::complex<float>* a,
                                           const std::complex<float>* b,
                                           size_t count,
                                           double& e_a,
                                           double& e_dst);
   ```

#### 0.2 扩展 `DiscreteCirResult` 与 `DiscreteCirApplyOptions`
在 `include/openuwb/sic/discrete_cir_sic.h` 中：
1. 在 `DiscreteCirResult` 末尾追加（保持默认初始化，向后兼容）：
   ```cpp
   bool roi_energy_valid = false;
   double roi_e_y = 0.0;
   double roi_e_out = 0.0;
   ```
2. 新增轻量选项结构体与带选项的重载（**不破坏**任何旧函数签名）：
   ```cpp
   struct DiscreteCirApplyOptions {
       bool transaction_prepared = false; ///< Skip full event-order/y-finite rescan (caller guarantees)
       bool copy_out_on_fail = true;      ///< Copy y->out on failure when out != y
       bool allow_partial_clip = false;   ///< Clip [lo, hi) to [0, n) instead of failing context_truncated
       bool compute_roi_energy = false;   ///< Populate result.roi_e_y / roi_e_out during fused subtract
   };

   UWB_API bool apply_discrete_cir_pulses_ex(
       const std::complex<float>* y,
       size_t n,
       const mod::PulseList& pulses,
       const std::complex<float>* reference_sync,
       size_t sync_n,
       const DiscreteCirModel& model,
       int64_t roi_lo,
       int64_t roi_hi,
       const DiscreteCirApplyOptions& opts,
       DiscreteCirScratch& scratch,
       std::complex<float>* out,
       std::complex<float>* ihat,
       DiscreteCirResult* result = nullptr);

   UWB_API DiscreteCirResult fullpacket_cancel_discrete_ex(
       const std::complex<float>* y,
       size_t n,
       const std::vector<uint8_t>& phr_bits,
       const std::vector<uint8_t>& wire_psdu,
       bool fcs_pass,
       const mod::HrpModConfig& profile,
       const DiscreteCirModel& model,
       const std::complex<float>* reference_sync,
       size_t sync_n,
       int64_t roi_lo,
       int64_t roi_hi,
       const DiscreteCirApplyOptions& opts,
       DiscreteCirScratch& scratch,
       std::complex<float>* out,
       std::complex<float>* ihat);
   ```
   原有 `apply_discrete_cir_pulses(...)` 与 `fullpacket_cancel_discrete(...)` 保持完全相同的签名与默认行为（委托给 `*_ex`）。

#### 0.3 扩展 `SicUnifiedConfig` 与 `SicUnifiedTimings`
在 `include/openuwb/sic/sic_unified.h` 中：
1. `SicUnifiedConfig` 追加（均带向后兼容默认值）：
   ```cpp
   bool allow_partial_window_clip = false; ///< Clip edge interferers to [0, n) instead of context_truncated
   double validate_threshold = 0.25;       ///< Forwarded to PathFractional SicConfig::validate_threshold
   ```
2. `SicUnifiedTimings` 追加：
   ```cpp
   double recon_zero_us = 0.0;             ///< Buffer zeroing / tile init time (child of reconstruct_us)
   ```

---

### WP1（Subagent A）：DSP / PHY 向量相减与融合能量统计内核

- **文件所有权**：
  - `lib/dsp/uwb_demod_sic_avx512.cc`
  - `lib/phy/hrp/uwb_demod_timing.cc`
  - `tests/unit/dsp/qa_uwb_demod_sic_avx512.cc`
- **禁止修改**：`lib/sic/*`、`lib/runtime/*`

#### 详细任务
1. **补齐 `vector_sub3_complex` 的 AVX2 分支**（`lib/phy/hrp/uwb_demod_timing.cc:1828`）：
   - 添加 `__attribute__((target("avx2")))` 的 `vector_sub3_complex_avx2(dst, a, b, count)`，使用 4×256-bit（每次处理 16 个复数 = 32 floats）循环展开 + 尾部清理；
   - 正确支持 `dst == a` 原地操作（先 load `a` 和 `b`，再 store `dst`）；
   - 在 `vector_sub3_complex` 中按 `AVX-512 → AVX2 → scalar` 分发。
2. **实现 `vector_sub3_energy_complex`（标量 + AVX2 + AVX-512）**：
   - 在 `lib/dsp/uwb_demod_sic_avx512.cc` 实现 `sub3_energy_cf32`：
     - 每次处理 8 个复数（16 floats = 1 个 `__m512`），计算 `vd = _mm512_sub_ps(va, vb)` 并写回 `dst`；
     - 在同一循环内将 `va*va` 和 `vd*vd` 拆成低/高 8 个 float 转成 `__m512d`（`_mm512_cvtps_pd`）累加到 4 组独立 `__m512d` 累加器（隐藏加法延迟）；
     - 尾部用标量 `double` 累加；
     - 保持 `dst[i]` 与 `a[i] - b[i]` **逐位一致**（无 FMA 重排），且支持 `dst == a` 原地操作。
   - 在 `lib/phy/hrp/uwb_demod_timing.cc` 实现 `vector_sub3_energy_complex_avx2`（参考 `uwb_sic_unified.cc:778-848` 的双精度累加结构，每次处理 8 个复数 = 2 个 `__m256`，4 组 `__m256d` 累加器，支持 `dst == a` 原地操作）及标量回退。
   - **重要**：不得丢失 `NaN/Inf` 传播（如果 `va` 或 `vd` 含 `NaN/Inf`，最终 `e_a` 或 `e_dst` 必须为非有限值，以便调用方仅靠 `std::isfinite(e_dst)` 即可发现溢出）。
3. **单元测试**（`tests/unit/dsp/qa_uwb_demod_sic_avx512.cc`）：
   - 覆盖 `count ∈ {0, 1, 3, 4, 7, 8, 15, 16, 17, 1023, 4096, 10000}`；
   - 覆盖 `dst != a` 与 `dst == a`（in-place）两种模式；
   - 验证 `dst` 逐位等于 `a - b`；
   - 验证 `e_a` / `e_dst` 与参考 `double` 累加相对误差 `< 1e-12`；
   - 验证注入单个 `NaN` 或 `Inf` 时，`std::isfinite(e_dst)` 为 `false`。

---

### WP2（Subagent B）：SIC 重构内核、精确 CFO 预旋、L1D 分块流水与鲁棒性修复

- **文件所有权**：
  - `lib/sic/uwb_discrete_cir_reconstruct.cc`
  - `tests/unit/sic/qa_uwb_discrete_cir_reconstruct.cc`
- **禁止修改**：`lib/phy/*`、`lib/dsp/*`、`lib/sic/uwb_sic_unified.cc`、`lib/runtime/*`

#### 详细任务

##### 2.1 修复 Tap 内部 CFO 相位误差（精确物理模型）
- 在 `apply_discrete_cir_pulses_impl` 中，当 `model.cfo_hz != 0.0` 时，在脉冲循环前构造局部数组：
  ```cpp
  alignas(32) std::array<cf, kDiscreteCirMaxTaps> taps_eff;
  const double omega = 2.0 * kPi * model.cfo_hz / kWorkRateHz;
  for (size_t d = 0; d < model.tap_count; ++d) {
      const int64_t lag = model.lag_lo + static_cast<int64_t>(d);
      const cd rot = std::exp(cd(0.0, omega * static_cast<double>(lag)));
      const cd tv = cd(model.taps[d]) * rot;
      taps_eff[d] = cf(static_cast<float>(tv.real()), static_cast<float>(tv.imag()));
  }
  ```
- 当 `model.cfo_hz == 0.0` 时直接使用 `model.taps.data()`（逐位不变）。
- **所有路径（AVX-512、AVX2 20-tap、AVX2 generic、标量回退）统一消费 `taps_eff`**，确保彼此等价且与精确公式 $e^{j\omega (j - q)} \sum_l h[l] x[j - q - l]$ 数学一致。
- 同步更新 `include/openuwb/sic/discrete_cir_sic.h`（若 Wave 0 未改注释则补充）以及 `tests/unit/sic/qa_uwb_discrete_cir_reconstruct.cc` 的 `sparse_pulse_apply_matches_approximate_cfo_reference`（升级为精确 CFO 参考模型），并新增 `cfo = 100000.0 Hz` 高频偏对照 `apply_discrete_cir` 的回归用例（要求相对残差 `< 1e-4`）。

##### 2.2 重构脉冲累加微架构优化 + 任意 `tap_count` SIMD + AVX-512 支持
1. **优化 `accumulate_pulses_cfo_20tap_avx2`**（以及新增泛化版 `accumulate_pulses_cfo_generic_avx2`、`accumulate_pulses_zero_cfo_generic_avx2`）：
   - **消除 `__mulsc3`**：手写 `pr = br * rr - bi * ri; pi = br * ri + bi * rr;` 及按 `e.sign` 取符号（`if (e.sign < 0) { pr = -pr; pi = -pi; }`，连浮点乘法都省了）；
   - **预乘 `sign_mask` 到 `hs0..hs4`**：循环外计算 `hs_k = _mm256_mul_ps(_mm256_permute_ps(h_k, 0xB1), sign_mask)`，循环内省去每脉冲一次 `_mm256_mul_ps(v_ci, sign_mask)` 和 1 个 YMM 寄存器；
   - **无符号位运算**：用 `uint32_t off = e.sample_offset; uint32_t b = off >> 10; uint32_t k = off & 1023u;` 替代有符号除法和取模；
   - **块相位递推 + 周期 anchor**：预计算 `cd block_step = std::exp(cd(0.0, omega * 1024.0))`，相邻块递推 `block_base_d *= block_step`，每 16 块（或跳块时）用 `std::exp` 重新 anchor，避免重复超越函数调用且零相位漂移；
   - **泛化任意 `tap_count`**：
     - `tap_count == 20` 保持 5×YMM 寄存器常驻全展开最快路径；
     - `tap_count != 20`（`1..256`）走 AVX2+FMA 泛化核（每脉冲内按 4-complex / YMM 展开 + 标量尾部，且若 `tap_count <= 32` 可按 8-YMM 常驻寄存器），彻底消除非 20-tap 静默退化到标量 `noinline` `fused_axpy` 的陷阱。
2. **接入 AVX-512 稀疏脉冲累加核**（受 `#if UWB_SIMD_HAVE_AVX512` 和 `openuwb::core::simd::use_avx512(KernelGroup::Sync)` 或 runtime 探测门控，参考 `bench_demod_sic_kernels.cc:372` 与 `UwbSimd.cmake` 约定）：
   - 20 taps = 40 floats = 2×`__m512`（16+16 floats）+ 1×`__m256`（8 floats），实现 zero-CFO 与 CFO 两个版本；
   - 遵循 `set_b3_simd_disabled(1)` 时强制走标量回退，确保 QA parity 测试继续有效。

##### 2.3 L1D 分块重构 + 融合相减流水（核心性能突破）
- **触发条件**：`ihat == nullptr` 且未触发极端溢出保底条件（见下方 Fail-Closed 协议）。当 `ihat != nullptr` 时，仍写入 `ihat` 完整支持区，保持 `ihat` 诊断语义 100% 不变（但相减阶段仍调用 `vector_sub3_energy_complex` 融合能量统计）。
- **Tile 几何**：
  - 常量 `constexpr size_t kTileSamples = 4096;`
  - 栈上 64 字节对齐缓冲区：`alignas(64) cf tile_buf[kTileSamples + kDiscreteCirMaxTaps];`（约 34 KB，完全常驻单核 48 KB L1D cache）。
- **执行流程**（覆盖全窗口有效区 `[clip_lo, clip_hi)`，其中默认不裁剪时 `clip_lo = lo, clip_hi = hi`）：
  1. 先清零 `tile_buf[0 .. kTileSamples + tap_count)`；
  2. 维护事件游标 `ev_idx`（利用 `pulses.events` 严格递增特性，跳过所有 `lo + offset + tap_count <= clip_lo` 的事件）；
  3. 按步长 `kTileSamples` 遍历 `tile_lo = clip_lo; tile_lo < clip_hi; tile_lo += kTileSamples`：
     - `tile_len = min(kTileSamples, clip_hi - tile_lo)`；
     - 找出所有起始位置 `pos = lo + e.sample_offset` 落在 `[tile_lo, tile_lo + tile_len)` 内的脉冲区间 `[ev_begin, ev_end)`，将它们的 `taps_eff` 累加到 `tile_buf + (pos - tile_lo)`（注：首块若 `clip_lo > lo`，即左边界裁剪时，需先将跨过 `clip_lo` 的最多 `tap_count` 个前置脉冲的右半截安全累加进 `tile_buf[0..tap_count)`）；
     - 此时 `tile_buf[0 .. tile_len)` 已经包含该块内所有脉冲的完整叠加结果（与全量累加逐位完全等价！）；
     - 计算该 tile 与 `[sub_lo, sub_hi)` 的交集 `[s0, s1)`：
       - 若 `s1 > s0`：调用 `demod::rxopt::vector_sub3_energy_complex(out + s0, y + s0, tile_buf + (s0 - tile_lo), s1 - s0, tile_ey, tile_eout)`（若 `!opts.compute_roi_energy`，可调 `vector_sub3_complex` 或统一调融合核并忽略能量以顺便免费做有限性检查！），累加 `total_ey += tile_ey; total_eout += tile_eout;`；
     - **滑动 Overlap**：将 `tile_buf[tile_len .. tile_len + tap_count - 1]` 拷贝到 `tile_buf[0 .. tap_count - 1]`，并清零 `tile_buf[tap_count .. tile_len + tap_count)`（只需 1 次 160 B 拷贝 + 32 KB L1D 清零！）。
  4. 若 `out != y`，拷贝非重叠两端 `[0, sub_lo)` 与 `[sub_hi, n)`（与原逻辑一致；若 `out == y`，零拷贝！）。
- **严格 Fail-Closed 溢出防护（尤其针对原地 `out == y` 模式）**：
  - 进入 tile 循环前，用 $O(\text{tap\_count})$ 计算 `float max_tap_abs`（20 个 tap 的最大幅度）；
  - 安全阈值判定：`const bool safe_for_inplace_tile = (max_tap_abs <= 1e20f);`（对于正常物理信号，`taps` 与 `y` 均在 `1e-6 ~ 1e4` 量级，单样点最多叠加 `tap_count` 个脉冲，绝不可能溢出 `FLT_MAX ≈ 3.4e38f`；若调用方在 `out == y` 模式下传 `transaction_prepared = false`，前面全窗扫描顺便确认 `max(|y|) <= 1e20f`，或当 `max_tap_abs > 1e20f` / 极端测试输入时，自动退回二阶段 `scratch.accumulator` 路径，先完整验证全部 `finite(recon[j]) && finite(y[j] - recon[j])`，确认无溢出后再写 `out`！）；
  - 正常 tile 路径结束后检查 `if (!std::isfinite(total_eout))`，若异常在非原地模式下直接 `fail("nonfinite_reconstruction")`（原地模式下由前置上界保证不可能发生，彻底杜绝半写污染）。
- **计时口径**：
  - tile 模式下将初始化与清零归入 `zero_us`，脉冲累加归入 `convolution_us`，tile 内相减+能量归入 `subtract_us`（避免在每 4096 块内频繁调用 `Clock::now()` 带来额外开销：只在整个 tile 流水前后各打一次时间戳，或记录 `channel_us` / `subtract_us` 合并统计时在注释中明确标注，注意确认既有 QA 是否断言 `convolution_us > 0` 或 `subtract_us >= 0`——经查单测主要校验 `cfo_us == 0.0` 和 `committed`）。

##### 2.4 边界裁剪（Opt-in）与 `copy_out_on_fail`、Legacy 脉冲限制修复
1. `fail()` lambda 改为仅在 `opts.copy_out_on_fail && out && y && n > 0 && out != y` 时执行 `std::copy(y, y + n, out)`。
2. 当 `opts.allow_partial_clip == false`（默认）：保持 `if (lo < 0 || hi > int64_t(n)) return fail("context_truncated");`，100% 兼容既有单测。
   当 `opts.allow_partial_clip == true`：将有效支持区间钳制为 `clip_lo = max<int64_t>(0, lo)`, `clip_hi = min<int64_t>(int64_t(n), hi)`；若 `clip_hi <= clip_lo` 返回 `no_overlap`；相减区间 `sub_lo = max(clip_lo, roi_lo)`, `sub_hi = min(clip_hi, roi_hi)`，若 `sub_hi <= sub_lo` 返回 `no_overlap`。脉冲累加安全处理跨越 `clip_lo` 和 `clip_hi` 的脉冲（仅累加落在 `[clip_lo, clip_hi)` 内的 tap 片段，绝不越界写内存）。
3. `fullpacket_cancel_discrete_ex`（第 937 行）：仅在 `!model.impulse_template`（即需要用 `profile.pulse` 调制成形波形的 NonSparse DiscreteCir 路径）时检查 `profile.pulse.shape != mod::PulseShape::Legacy`；Ipatov 路径本身使用单位冲激 `scratch.impulse_profile`，不再拒收非 Legacy 脉冲配置。
4. `apply_discrete_cir_impl`（NonSparse 路径顺带优化）：
   - scale 阶段（#L544-L549）：当 `scale_real == 1.0` 且 `waveform != scaled` 时向量化有限性检查并拷贝，或至少采用 AVX2 有限性位掩码检查 + 乘常数；
   - 相减阶段（#L664）：同步切换为 `vector_sub3_energy_complex`（当需要能量时）。

##### 2.5 单元测试（`tests/unit/sic/qa_uwb_discrete_cir_reconstruct.cc`）
- 更新 `sparse_pulse_apply_matches_approximate_cfo_reference` 为精确 CFO 参考并新增 `cfo = 100 kHz` 对比 `apply_discrete_cir` 精度用例；
- 扩展 `b2_b3_sparse_pulse_avx2_fma_20tap_and_fallback_parity` 覆盖 `tap_count ∈ {1, 7, 16, 19, 20, 24, 32, 48}`；
- 新增 `tile_path_matches_full_ihat_path_bit_exact`：对比 `ihat != nullptr`（全量缓冲）与 `ihat == nullptr`（L1D tile 路径，含 `out != y` 与 `out == y` 原地）结果逐位一致，且 `roi_e_y` / `roi_e_out` 与双精度参考一致；
- 新增 `partial_window_clip_left_and_right_edges`：测试 `q = 2, lag_lo = -4`（左边界 `lo = -2`）和 `hi = n + 500`（右边界超窗）在 `allow_partial_clip = false` 时返回 `context_truncated`，在 `allow_partial_clip = true` 时成功消除窗口内干扰且无越界（配合 `scalar-debug` / ASan 检查）；
- 新增 `copy_out_on_fail_false_skips_buffer_write`：验证 `copy_out_on_fail = false` 时失败绝不写 `out`；
- 新增 `inplace_overflow_injection_is_fail_closed`：故意注入超大 tap/信号触发溢出保护，验证原地模式 `y` 未被破坏。

---

### WP3（Subagent C）：Unified SIC 编排层修复与 `PacketProcessor` 原地相减

- **文件所有权**：
  - `lib/sic/uwb_sic_unified.cc`
  - `include/openuwb/runtime/live_packet_processor.h`
  - `lib/runtime/uwb_live_packet_processor.cc`
  - `tests/unit/sic/qa_uwb_sic_unified.cc`
  - `tests/unit/runtime/qa_uwb_live_packet_processor_sic.cc`
- **禁止修改**：`lib/dsp/*`、`lib/phy/*`、`lib/sic/uwb_discrete_cir_reconstruct.cc`

#### 详细任务

1. **L1 稀疏路径跳过重复校验与标量全窗扫描**（`lib/sic/uwb_sic_unified.cc:668`）：
   - 改为调用 `apply_discrete_cir_pulses_ex`，构造 `DiscreteCirApplyOptions opts`：
     - `opts.transaction_prepared = scratch.preamble_cache_valid && !scratch.cached_preamble_events.empty();`（仅当 `estimate_discrete_cir` 已校验过 `y` 且使用预缓存前导脉冲列表时开启）；
     - `opts.copy_out_on_fail = result.copy_out_on_bypass;`
     - `opts.allow_partial_clip = cfg.allow_partial_window_clip;`
     - `opts.compute_roi_energy = true;`
2. **L2 路径贯通选项 + 补齐细分计时**（`lib/sic/uwb_sic_unified.cc:734-741`）：
   - 改为调用 `fullpacket_cancel_discrete_ex`，传入 `opts`（`copy_out_on_fail = result.copy_out_on_bypass`、`allow_partial_clip = cfg.allow_partial_window_clip`、`compute_roi_energy = true`）；
   - 补齐漏填的 4 个细分耗时字段：
     ```cpp
     result.timings.recon_zero_us = dres.zero_us;
     result.timings.recon_scale_us = dres.scale_us;
     result.timings.recon_conv_us = dres.convolution_us;
     result.timings.recon_cfo_rot_us = dres.phase_us + dres.cfo_us;
     ```
   - L1 分支（#L673, #L714）同步补充 `result.timings.recon_zero_us = dres.zero_us;`。
3. **避免失败路径重复拷贝**：
   - `reconstruct_and_subtract_unified` 的 `fail_closed`（#L580）在底层已经根据 `copy_out_on_fail` 处理或未写 `out` 时避免重复拷贝（只需在底层失败且未做拷贝时拷贝一次，最简单：透传 `copy_out_on_fail = false` 给底层 `*_ex`，统一由 `reconstruct_and_subtract_unified` 的 `fail_closed` 根据 `result.copy_out_on_bypass` 最多拷贝一次！）。
4. **消费融合计算的 ROI 能量，消除 68 us 二次重读**（`lib/sic/uwb_sic_unified.cc:1194-1202`）：
   - 若 `reconstruct_and_subtract_unified` 已通过 `dres.roi_energy_valid` 直接算出 `suppression_db`（或存到局部/结果中），`sic_process_unified` 直接使用它，**不再调用 `compute_roi_energies(y, result.out, ...)`**；仅对未填充融合能量的后端（如 `PathFractional` / NonSparse）保留 `compute_roi_energies` 回退，且保持 `compute_roi_energies` 公开 API 和单测不变。
   - *(注意：当 `result.out == y` 原地模式时，必须在相减前/相减中算出 `e_y`，相减后再读 `y` 会拿到已修改的值！因此消费融合能量是支持原地相减的必要前提！)*
5. **`AutoFallback` 下 L2 截断自动回退 L1**：
   - 在 `sic_process_unified`（#L1183-L1191）：若 `cfg.policy == SicPolicy::AutoFallback && target_level == SicCancellationLevel::L2_FullPacket && !ok && result.reason && std::strcmp(result.reason, "context_truncated") == 0`，自动设置 `result.fallback_reason = "l2_context_truncated"` 并回退尝试 `L1_Preamble`（前导更短，常能完整落在窗口内）。
6. **`PathFractional` 统一接口小修复**：
   - `SicUnifiedConfig::to_sic_config()`（#L174-L197）补充 `c.validate_threshold = validate_threshold;`；
   - 修正 `reconstruct_and_subtract_unified` 中 `PathFractional` L1 分支（#L608-L609）的计时归属与 `copy_out_on_bypass` 一致性。
7. **`PacketProcessor` 启用原地相减（in-place）**：
   - 修改 `include/openuwb/runtime/live_packet_processor.h`（移除或按需保留空 `sic_out_`——移除 `sic_out_` 可为每个 worker 节省 ~4 MB 内存）与 `lib/runtime/uwb_live_packet_processor.cc:66, 113-132`：
     - `sic_result.out = work_.data();`（原地相减！）；
     - `sic_result.copy_out_on_bypass = false;`
     - `sic_result.i_hat = nullptr;`
     - 若需启用滑窗边缘裁剪，确认是否由 `openuwb_rx` 配置传入 `allow_partial_window_clip`（不在 `PacketProcessor` 强行覆盖调用方配置，保留默认 `false` 并在单测中显式验证开启与关闭两种模式）。
     - 无论 `sic_committed` 为 true 还是 false，`demod_in` 均为 `work_.data()`（当失败或 L0 时，`work_` 保证逐字节未被触碰，严格满足 fail-open 契约！）。
8. **单元测试**：
   - `tests/unit/sic/qa_uwb_sic_unified.cc`：
     - 验证 L2 `recon_conv_us > 0`（不再恒为 `0.0`）；
     - 验证原地模式 `result.out = y_copy.data()` 下 `suppression_db` 与异地模式完全一致（差异 `< 1e-4 dB`），且 L0/失败路径 `y_copy` 逐字节不变；
     - 验证 `AutoFallback` 在干扰包尾部超窗时自动从 L2 回退到 L1 或通过 `allow_partial_window_clip = true` 直接 L2 部分消除；
     - 验证 `PathFractional` 通过 `sic_process_unified` 不再因 `uncalibrated_validate_threshold` 失败。
   - `tests/unit/runtime/qa_uwb_live_packet_processor_sic.cc`：
     - 验证原地 SIC 在无干扰（L0 fail-open）、强干扰碰撞（L2 committed）及重复窗口下零状态残留、PRR/FCS 与此前一致。

---

### WP4（Subagent D）：基准复现工具、ABBA 性能对比与跨预设交付验收

- **文件所有权**：
  - `tools/bench_sic_reconstruct_opt.cc`（新复现工具放 `tools/`，或独立编译脚本，不污染核心库）
  - 验收报告：`docs/performance/SIC重构与相减性能优化验收报告_2026-11-21.md`

#### 详细任务
1. 编写可复现基准工具，固定随机种子，串行、同环境、ABBA 配对对比优化前后：
   - 场景 1：Live Hit `n = 506k`、`code 11`、`128 SYNC`、`decawave`、`cfo = 20 kHz`、`sparse=on` 与 `shaped_sparse=on`（含原地 vs 异地对比）；
   - 场景 2：短窗 `n = 135k`、`code 10`、`64 SYNC`、`4z2`、`L1` 与 `L2`；
   - 场景 3：高频偏精度验证 `cfo ∈ {0, 20 kHz, 50 kHz, 100 kHz}`（对比抑制比 dB）；
   - 场景 4：多 tap 宽度性能验证 `tap_count ∈ {16, 20, 24, 32}`；
   - 场景 5：边缘截断验证（`q = 2` 左边缘与 `hi = n + 2000` 右边缘）。
2. 执行 `AGENTS.md` 规定的完整门禁（**所有命令必须前缀 `env -u LD_LIBRARY_PATH`**；若本机 ccache 不可用追加 `-DUWB_USE_CCACHE=OFF`）：
   - 快闭环：`dev-avx2` 构建 + `ctest -L rx_sic --output-on-failure -j`；
   - 全量交付门禁：
     1. `native-avx2` 全量构建 + 全量 `ctest`；
     2. `no-uhd` 全量构建 + 全量 `ctest`；
     3. `native-avx512` 全量编译 + `nm` 核验新增 AVX-512 强符号存在（无 AVX-512 机器不执行）；
     4. `scalar-debug` 全量构建 + 全量 `ctest`（核验标量路径与边界安全）；
     5. 临时前缀安装消费验证：`tests/integration/install_consumer/run_native_install_consumer.sh`（因修改了公共头文件，必跑）。
3. 汇总提交哈希、CPU/ISA、构建参数、种子、实测数据表与全部门禁结果写入验收报告。

---

## 5. 关键风险与避坑指南

1. **禁止破坏既有单测契约**：
   - `qa_uwb_discrete_cir_reconstruct.cc` 显式断言默认模式下越界返回 `"context_truncated"`，所以边界裁剪必须做成默认 `false` 的 opt-in 开关 `allow_partial_clip`；
   - `set_b3_simd_disabled(1)` 必须继续屏蔽所有 AVX2/AVX-512 稀疏脉冲核，确保 SIMD vs 标量 parity 单测真实覆盖；
   - `zero-CFO` 下 AVX2 / AVX-512 / 标量以及 L1D tile 与全量缓冲必须**逐位一致**（`std::memcmp == 0`）。
2. **注意 `std::complex` 的性能陷阱**：
   - 在无 `-ffast-math` 编译选项下，`std::complex<float>` 乘法会调用 `__mulsc3`，任何热循环内都必须拆成实/虚部 `float` 运算。
3. **注意 `out == y` 原地模式下的两处陷阱**：
   - 陷阱 1：先相减再调用 `compute_roi_energies(y, out, ...)` 会把 `e_y` 算成 `e_out`，导致 `suppression_db` 恒为 `0 dB`！必须用融合相减核在写回 `out[i]` 前读取旧的 `y[i]`。
   - 陷阱 2：中途失败绝不能留下一半已修改的 `y`，因此必须执行 2.3 节规定的前置上界检查，高风险幅度自动回退二阶段提交。
4. **注意 `UWB_USE_CCACHE` 沙箱问题**：
   - 如遇到 `ccache: not found`，配置预设时加 `-DUWB_USE_CCACHE=OFF`。
5. **不要覆盖未提交的 `lag_lo = -4` 修改**：
   - `include/openuwb/sic/sic_unified.h:248` 已将 `lag_lo` 默认值改为 `-4`，开发过程中须保持该默认值。
