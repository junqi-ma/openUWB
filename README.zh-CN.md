# OpenUWB, 802.15.4a/z 兼容的 UWB 软基带

[English](README.md) | 简体中文

> [!NOTE]
> **施工告示（Work in Progress）**：本项目目前仍处于积极开发阶段，模块接口、命令行参数与功能仍在持续调整中，暂不承诺 API 稳定性，也不在此处列出性能指标。

## 为什么做 OpenUWB

现有的 UWB 研究大多基于 DW1000、DW3000 等商用专用芯片。这类芯片虽然集成度高，但底层基带固化在硬件与闭源固件中，无法提供足够的修改自由度，限制了研究者对物理层波形、同步与信道估计算法、干扰消除以及协议时序的深入探索。

我们的目标是利用现代大带宽软件无线电（例如 USRP X410、X420）来实现开放的 UWB 软基带（UWB Radio），并在物理层与协议层保持对 IEEE 802.15.4a/z 标准的兼容，使其能够直接与商用标准 UWB 模组（如 NXP、Qorvo、Apple 等）进行互操作。

我们相信，将「大带宽 SDR 的完全可编程性」与「商用标准生态的互操作性」结合起来，能够为 UWB 通信、无线感知（Sensing）、测距（TWR）以及物理层信号处理带来大量新的研究机会。

## 核心进展时间线

- **2026-08 — 实时解调、商用设备互操作与按窗重采样架构**
  - 确立原生采样率（737.28 MS/s）按窗截包 + PDU 级 `65/48` 多相重采样（998.4 MS/s 工作采样率）的主架构，打通未知起点盲捕获、周期锁定以及真实 QM35825 / DW1000 信号解调。
  - 完成 HRP BPRF 解调（前导码、SFD 搜索、CFO 补偿等）与 MATLAB 参考实现对齐。
- **2026-09 — X410 定时回波 CIR 感知、双 TX 干扰与 DPDK**
  - 在 USRP X410（CG400 / CG600）上打通定时回波与逐 repetition CIR 实时输出、运行时扫频/跳频及 UHD DPDK 高速传输支持。
  - 实现 C++ 双通道定时发射（感知信号 + 同步干扰信号）与连续 Buffer 调度架构。
- **2026-09 ~ 2026-10 — SS-TWR / DS-TWR 协议核心与离线原生 PHY 闭环**
  - 交付独立于硬件的 TWR 配置契约、ToF 数学库（已与 MATLAB oracle 对齐）、SS/DS-TWR 离线协议状态机（FSM）与原生 PHY 闭环仿真（硬件测距尚未实现）。
- **2026-11（中旬）— 低 SNR 解调、SIMD 稀疏相关与 SIC 管线统一**
  - 完成低 SNR 捕获与解调优化（折叠同步、SFD 离散符号网格救援、BPM 软判决等），并新增 AVX2 / AVX-512 稀疏三元相关加速内核。
  - 将单包连续干扰消除（SIC：解码 → CIR 估计 → 波形重构 → 相减）与标准解调管线统一（待进一步接入实时调度运行时）。
- **2026-11-20 — 去 GNU Radio 化与 OpenUWB 架构重构（Phase 0–4）**
  - 剥离 GNU Radio 运行时依赖（历史适配器归档至 `archive/gnuradio-adapter`），将原生 C++17 / UHD 工程提升至仓库根目录。
  - 项目正式更名为 **OpenUWB**，确立 `openuwb` 命名空间与单向依赖的模块化域库（`support → dsp → phy_hrp → {sic, sensing, twr_core}`，由 `openuwb_radio` 隔离 UHD、`openuwb_runtime` 统一组装），并切换为 `openuwb_*` 应用入口。
- **2026-11-21 — `openuwb_rx` 接入 SIC 与构建体系加速**
  - 将前导/Ipatov SIC 接入 `openuwb_rx`（`--sic`，默认关闭零分配、开启 fail-open），补齐经验稳态 L1 向量化与截断 SFD 救援（离线回放验证）。
  - 全预设切换为 Ninja，引入 PCH、`ccache`/`mold` 自动探测、重型头文件声明/实现解耦及按域裁剪开发预设（`dev-rx` / `dev-tx` / `dev-twr` / `dev-radar`）。
- **2026-11-24 — RX 前端与解调管线全面重构与加速（离线回放）**
  - 重构 HRP 解调核心阶段：Stage-1 Timing（局部 `K=4` 折叠快路径与稀疏锚点拟合）、Stage-2 CFO（解耦 SIMD 符号采样与复数域 ML 估计器）、Stage-3 SFD（软极性、低门限搜索与增量补旋）。
  - 新增 250 MHz 信道化 `65/48` 重采样滤波器（`--taps ch250_eq16/24/32/48`）、折叠前导码互相关捕获前端（`--acq-mode xcorr`）以及整段 AVX2 + 锁后静默期按需固定单音消除（`skip_samples` + `--tone-when-locked window`，配合 `--autocorr-when-locked suspend`）；以上均为离线 DW3000 实采回放验证（尚未上板）。

## 当前仓库导航

如果你希望了解当前的构建方式或完整开发细节，可以参考以下文档：

- 构建与运行说明：[原生 UHD/C++ 构建使用说明](docs/原生UHD与C++构建使用说明.md)
- 完整开发状态：[开发状态](开发状态.md)
- 更多文档索引：[docs/README.md](docs/README.md)
