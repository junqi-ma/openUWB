# OpenCode 执行提示词 — Phase 2（按域拆库、预留 HRP 4a/4z/4ab、拆巨型文件）

> 适用分支：`refactor/openuwb-arch`（前置条件：Phase 1 全绿）  
> 参考总方案：`docs/OpenCode执行方案_openUWB全仓库架构重构_2026-11-20.md`

---

Phase 1 已全绿。现在请在当前分支 `refactor/openuwb-arch` 上，严格按照 `AGENTS.md` 与 `docs/OpenCode执行方案_openUWB全仓库架构重构_2026-11-20.md` 第 2、3、5 节，**只执行 Phase 2**。

请主动使用**多 Subagent 分波次并行加速**，严格遵守写隔离规则（`.git` 操作与共享 `CMakeLists.txt` 由主控独占）。

---

### 多 Subagent 执行编排

#### Wave 1：并行只读调研（同时启动 4 个只读 Subagent）
- **Subagent-R1（全量头文件依赖与 DAG 审查）**：
  - 按总方案第 3 节映射表，扫描 `gr-uwb/include/gnuradio/uwb/` 全部头文件之间的 `#include` 关系，核查是否存在违反第 2.2 节单向依赖 DAG 的反向依赖（例如 `dsp` 引用 `phy`、`phy` 引用 `sic`/`runtime`、非 `radio` 层引用 `UHD`）。
  - 如有反向依赖，给出最小零行为变更的解耦建议。
- **Subagent-R2（`uwb_twr_config.h` 6818 行拆分设计）**：
  - 分析 `uwb_twr_config.h` 的内部结构，设计清晰的子头拆分方案（例如 profile/defaults/validation/parser 等），保留原聚合入口，保证外部引用与所有 TWR QA 零改动通过。
- **Subagent-R3（`uwb_uhd_sc16_pipeline.cc` 2387 行拆分设计）**：
  - 梳理 `UhdDeviceSc16Source`（属于 `radio/`）与 `UhdSc16Pipeline`（属于 `runtime/`）的精确切分边界，保证零性能退化、零并发语义改变、无 UHD 时 `dry-run` 正常编译。
- **Subagent-R4（HRP 802.15.4a / 4z / 4ab 预留扩展点设计）**：
  - 检查 `uwb_phy_profile.h` 与 `uwb_hrp_codec.h`，设计轻量级 `include/openuwb/phy/hrp/common/standard_rev.h`（定义 `HrpStandardRev::Rev4a` / `Rev4z` / `Rev4ab` 及轻量元数据字段，不引入空壳逻辑，不影响现有 DW1000/DW3000/QM35825 行为）。

#### Wave 2：主控串行搬迁与建立 DAG（Commit 2A + Commit 2B + Commit 2C）
- **Commit 2A（纯 `git mv`，不改内容）**：
  - 在 `gr-uwb/include/openuwb/` 下按总方案第 3 节映射表建立 `{support,dsp,phy/hrp/{common,tx,rx},sic,sensing,twr,radio,runtime}/` 目录（包含 Phase 1 新增的 `uwb_log.h` → `include/openuwb/support/log.h` 及提升的 `uwb_live_*.h` → `include/openuwb/runtime/`）。
  - 在 `gr-uwb/lib/` 下建立镜像子目录，用 `git mv` 搬迁全部 `.h` 与 `.cc`（纯移动，不混入代码修改，保留 Git rename 历史）。
- **Commit 2B（生成过渡兼容转发头）**：
  - 在 `gr-uwb/include/gnuradio/uwb/` 下为每个旧头文件生成薄转发头（`#pragma once` + `// TODO(openuwb-compat-removal): remove legacy forwarding header` + `#include <openuwb/<domain>/<file>>`）。
  - 更新 `gr-uwb/include/openuwb/` 内部头文件之间的 `#include` 为新路径 `<openuwb/...>`。
- **Commit 2C（独立 CMake target 单向 DAG）**：
  - 为 `lib/{support,dsp,phy/hrp,sic,sensing,twr,radio,runtime}` 建立独立静态/接口库 target（`openuwb_support`, `openuwb_dsp`, `openuwb_phy_hrp`, `openuwb_sic`, `openuwb_sensing`, `openuwb_twr_core`, `openuwb_twr_fake_link`, `openuwb_twr_phy`, `openuwb_radio`, `openuwb_runtime`），严格执行第 2.2 节单向依赖。
  - 导出 `uwb_core`、`uwb_native` 兼容聚合目标并标注 `TODO(openuwb-compat-removal)`。
  - 快速增量构建验证 DAG 编译通过。

#### Wave 3：无冲突并行重构（同时启动 3 个写 Subagent，禁止改根 `CMakeLists.txt` 与执行 `git`）
- **Subagent-W1（`twr/` 目录专属）**：按 Subagent-R2 方案拆分 `uwb_twr_config.h` 为多个职责清晰的子头并保留聚合头。
- **Subagent-W2（`radio/` + `runtime/` 专属）**：按 Subagent-R3 方案拆分 `uwb_uhd_sc16_pipeline.cc` 为 `radio/uhd_device_source.cc` 与 `runtime/scheduled_pipeline.cc`；在 `include/openuwb/runtime/config.h` 预留 `RxConfig` / `TxConfig` / `TwrConfig` / `RadarConfig` 统一配置结构体约定。
- **Subagent-W3（`phy/hrp/common/` 专属）**：按 Subagent-R4 方案新增 `include/openuwb/phy/hrp/common/standard_rev.h` 并补充极简单元测试 `qa_hrp_standard_rev.cc`。

#### Wave 4：主控收口提交（Commit 2D）与全量 Gate 验收
- 主控更新相关子目录 `CMakeLists.txt`、检查无遗留警告，提交 **Commit 2D**。

---

### 强制约束
- **只做 Phase 2**：暂不移除 `gr-uwb/` 外壳、暂不全局切换 `namespace gr::uwb` → `namespace openuwb`（留给 Phase 3）、严禁把 SIC 接入 `live_rx`。
- **保护热路径性能**：严禁将 `uwb_demod_*.h`、`uwb_rx_*.h`、`uwb_rational_resampler_*.h`、`uwb_cir_fir_simd.h` 中的内联 SIMD 热路径核盲目搬进 `.cc`。
- **零算法改动**：不修改任何数值逻辑，不放宽任何 QA 容差，不破坏 `Reporter` 的 `stdout` 输出契约。
- 所有构建和测试必须使用 `env -u LD_LIBRARY_PATH`。

### Phase 2 Gate 验收（由主控统一执行并如实汇报）
1. **三预设编译与全量单测**（逐个执行，全部通过）：
   - `cmake --preset native-avx2 && cmake --build --preset native-avx2 -j && env -u LD_LIBRARY_PATH ctest --test-dir build/native-avx2 --output-on-failure`
   - `cmake --preset native-avx512 && cmake --build --preset native-avx512 -j && env -u LD_LIBRARY_PATH ctest --test-dir build/native-avx512 --output-on-failure`
   - `cmake --preset no-uhd && cmake --build --preset no-uhd -j && env -u LD_LIBRARY_PATH ctest --test-dir build/no-uhd --output-on-failure`
2. **纯基带边界静态检查**：
   - 确认 `openuwb_support`、`openuwb_dsp`、`openuwb_phy_hrp`、`openuwb_sic` 四个库源码与链接项中 **零 `UHD` 依赖、零 `std::cout`/`std::cerr` 调用**。
3. **控制台回放回归**：
   - `env -u LD_LIBRARY_PATH python3 tools/rx/qa_live_dw3000_console.py`（指向 `build/native-avx2`）通过。
4. 汇报：Commit 2A~2D 的 hash、各 Subagent 产出摘要、三个预设的 CTest 精确通过/失败/跳过数、纯基带静态检查结果及最终 `git status --short`。
