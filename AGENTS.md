# 项目开发约定

`main` 已合入 Phase 0–4 原生架构重构（[重构验收报告](docs/performance/openUWB架构重构验收报告_2026-11-20.md)）。
当前工作分支 `feature/openuwb-rx-sic` 已完成 SIC 接入 `openuwb_rx`（纯软件/回放，未上板）及
Ninja + PCH + 重型头文件声明/实现分离重构（[构建加速验收报告](docs/performance/Ninja与头文件解耦重构验收报告_2026-11-21.md)）。
当前能力、入口和验证边界见[开发状态](开发状态.md)；构建预设与参数见
[原生 UHD/C++ 使用说明](docs/原生UHD与C++构建使用说明.md)。历史低 SNR 优化在
`feature/sic-research`；相关结果为软件/文件回放，未上板。

## 工作范围与边界

- 按用户本次任务确定范围，不自动推进硬件、TWR 或 SIC 里程碑。
- 后续开发统一使用原生 UHD/C++，从仓库根目录构建。
- 目标与目录：`openuwb_support`、`openuwb_dsp`、`openuwb_phy_hrp`、
  `openuwb_sic`、`openuwb_sensing`、`openuwb_twr_core/phy`、`openuwb_radio`、
  `openuwb_runtime`；头文件 `include/openuwb/<domain>/`，源码 `lib/<domain>/`。
  依赖严格单向：`support → dsp → phy_hrp → {sic,sensing,twr_core}`，
  `runtime → phy_hrp + sic + sensing + radio`，只有 `openuwb_radio` 可链接 UHD。
- 历史 GNU Radio 实现在 `archive/gnuradio-adapter`（`d5dfd30`）。过渡兼容层
  `include/gnuradio/uwb/`、`gr::uwb` 别名、`uwb::` 包与 `uwb_*` 应用软链全部
  标注 `TODO(openuwb-compat-removal)`，清除顺序见重构验收报告 §4.5。
- timed 无 replay 必须 fail-closed。
- `openuwb_rx` 默认 `--fused on --cir full --taps eq16 --sic off`（`--fused off` 已废弃）。
  `--sic off` 保持零分配与既有 `stdout` 协议，开启后 L0/失败路径保持 fail-open。
  Python wrapper 为兼容入口，默认 `--cir auto`；其他 C++ 参数使用 `--` 透传。
- 不把 SFD-only 进锁、前端种子数写成完整帧 FCS 成功率；区分仿真、回放和上板证据。

## 开发与验证

- 完整包消费者使用 `include/openuwb/phy/hrp/rx/demod_pipeline.h` 并链接
  `openuwb::phy_hrp`（或 `openuwb::runtime`）；单阶段 QA 可用 `demod_core.h`
  聚合入口。TX/RX 共用编码原语放 `include/openuwb/phy/hrp/common/codec.h`。
- 头文件与编译分层：重型非模板实现放在 `lib/<domain>/*.cc` 编译进静态库，不把
  大段实现体写回公共头文件。重采样内核使用 `extern template` + `.cc` 显式实例化
  （`65_48`、`65_32`、`48_65`、`32_65`），消费者须链接对应库目标。保持
  `UwbSimd.cmake` 编译期 ISA 一致性，禁止 `-march=native`；新增 `uwb_simd_apply()`
  QA 目标时同步加入 `tests/CMakeLists.txt` 的 `UWB_QA_PCH_SIMD_TARGETS`；保留
  `hrp_mod_core.h` 中 `append_ieee_fcs` 内联（`qa_uwb_twr_frame.cc` 扫描契约）。
- 修改前检查相关实现、测试和报告。保持逐位一致性约定，不放宽既有 QA 容差。
- UWB 检测/解调算法参考 `UWB_demodulation/` MATLAB 实现；检测位置最终须逐样点/
  逐包交叉验证。未做的 MATLAB 验证明确列出。
- 修改代码后构建并运行测试（全预设已切 **Ninja**，QA/lib/apps 默认启用 PCH，ccache 自动探测且
  PCH+ccache 时自动追加 `-fpch-preprocess`；
  始终使用 `env -u LD_LIBRARY_PATH` 防止命中旧安装库，不把全局 install 作为普通测试前置）：
  - 日常迭代快闭环（关闭实验目标；`rx_sic` 已排除 `heavy` 标签的重型网格扫描）：
    ```sh
    env -u LD_LIBRARY_PATH cmake --preset dev-avx2
    env -u LD_LIBRARY_PATH cmake --build --preset dev-avx2 -j
    env -u LD_LIBRARY_PATH ctest --test-dir build/dev-avx2 -L rx_sic --output-on-failure -j
    ```
    （仅跑受影响域：改 RX/SIC 只跑 `-L rx_sic`，TWR/Sensing/DSP 修改按需换用
    `-L twr` / `-L sensing` / `-L dsp_extra`；日常不连带跑 `legacy_diag`、
    `heavy` 或无关域，无真机/沙箱环境跳过 `uwb_qa_uwb_uhd_backend_device.cc`。
    重型低 SNR 网格扫描仅在改相关算法时用 `-L heavy` 显式运行。）
  - 单功能域开发可用 `UWB_DEV_SCOPE=rx|tx|twr|radar` 对应预设
    `dev-rx` / `dev-tx` / `dev-twr` / `dev-radar`（独立 `build/dev-*` 目录，
    均默认开启 `UWB_DEV_SHARED_LIBS=ON`，改库实现只重链 `.so` 而不重链全部 QA 可执行文件）；
    非 `all` 作用域不生成安装规则，也不注册 install_consumer。
    `dev-avx2` 如需共享库可追加 `-DUWB_DEV_SHARED_LIBS=ON`。
    `UWB_LINKER=AUTO|lld|mold|gold|ld` 选择链接器（先探活再启用）。
  - 交付验收与跨预设门禁（**仅在最终合入/交付或用户要求时运行，日常修改不重复跑**）：
    运行 `native-avx2`、`no-uhd` 全量构建与全量 `ctest`（无硬件/沙箱可排除
    `uwb_qa_uwb_uhd_backend_device.cc`），`native-avx512` 全量编译（无 AVX-512 机器仅编译并核验强符号），
    以及 `scalar-debug`；改动公共头、安装规则或 CMake 包导出时，须以临时前缀
    （`-DUWB_NATIVE_INSTALL_PREFIX=...` 或 `tests/integration/install_consumer/run_native_install_consumer.sh`）
    实跑 `uwb_qa_native_install_consumer`。
- 仅在有明确性能对比需求或用户要求时再做串行 ABBA 回归（避免每次小改都跑全场景大矩阵）。
  性能测试串行、同环境对照；有调度反馈的 locked 回放采用 ABBA 配对。
  记录提交、CPU、ISA、构建参数、worker/queue、夹具种子/哈希和统计口径。
- 测试信号通常在 `testdata/`；实验摘要在 `docs/performance/data_*`。
  新复现工具放入 `tools/`，生成数据放 build 或临时目录并记录来源。
- TWR 相关修改还须遵守历史快照中“验证环境”和“TWR 开发必须遵守”两节的约束；
  未获得 ToA/硬件证据不得声明测距有效。

## 历史记录

原有阶段报告、完整历史规则与 TWR 契约保存在 `feature/sic-research` 分支的
`docs/archive/AGENTS_2026-11-17历史快照.md`。2026-10-06 及更早阶段开发状态归档在
[开发状态历史记录](docs/archive/开发状态_历史记录_截至2026-10-06.md)。
其中分支目标、默认参数、测试数量与“本轮”仅属于对应阶段，不作为当前任务目标。
