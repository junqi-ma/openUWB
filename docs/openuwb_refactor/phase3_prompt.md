# OpenCode 执行提示词 — Phase 3（顶层去 `gr-uwb/` 外壳、切换 `namespace openuwb`、`tools/` 批量迁移）

> 适用分支：`refactor/openuwb-arch`（前置条件：Phase 2 全绿）  
> 参考总方案：`docs/OpenCode执行方案_openUWB全仓库架构重构_2026-11-20.md`

---

Phase 2 已全绿。现在请在当前分支 `refactor/openuwb-arch` 上，严格按照 `AGENTS.md` 与 `docs/OpenCode执行方案_openUWB全仓库架构重构_2026-11-20.md` 第 0、2、5 节，**只执行 Phase 3**。

请主动使用**多 Subagent 分波次并行加速**，尤其在 `tools/` 42 个脚本迁移与多预设验证环节充分利用并行能力；`.git` 搬迁与核心 `CMakeLists.txt` 由主控独占。

---

### 多 Subagent 执行编排

#### Wave 1：并行只读调研（同时启动 3 个只读 Subagent）
- **Subagent-R1（顶层提升与 install/package 契约盘点）**：
  - 检查 `gr-uwb/CMakeLists.txt`、`gr-uwb/cmake/`、`gr-uwb/tests/integration/install_consumer/`、`gr-uwb/python/` 中所有相对路径（如 `CMAKE_SOURCE_DIR`、`gr-uwb/` 引用），输出提升至仓库根目录后的精确修改清单。
- **Subagent-R2（`namespace gr::uwb` → `namespace openuwb` 迁移审查）**：
  - 扫描 `include/openuwb/`、`lib/`、`apps/`、`tests/` 中的 `namespace gr`、`gr::uwb` 用法，确认迁移为 `namespace openuwb` 后，如何在兼容层提供 `namespace gr { namespace uwb { using namespace ::openuwb; } }`（避免 ADL/重定义冲突，全部标注 `TODO(openuwb-compat-removal)`）。
- **Subagent-R3（`tools/` 42 个脚本硬编码路径盘点）**：
  - 按 `tools/rx/`、`tools/sic/`、`tools/perf/`+`tools/performance/`、`tools/twr/` 四组列出所有引用 `gr-uwb`、`build-*`、`uwb_dw3000_live`、`uwb_uhd_sc16_pipeline_runner` 的位置，设计统一的 `tools/_common/paths.py`（优先环境变量 `OPENUWB_BUILD_DIR` → 根目录 `build/native-avx2` → 兼容回退旧路径 + `TODO(openuwb-compat-removal)`）。

#### Wave 2：主控串行完成顶层提升与命名空间切换（Commit 3A + Commit 3B + Commit 3C）
- **Commit 3A（纯 `git mv` 提升顶层，不改内容）**：
  - 将 `gr-uwb/CMakeLists.txt`、`gr-uwb/cmake/`、`gr-uwb/include/`、`gr-uwb/lib/`、`gr-uwb/apps/`、`gr-uwb/tests/`、`gr-uwb/python/uwb`（→ `python/openuwb`）用 `git mv` 提升到仓库根目录，清理空 `gr-uwb/` 目录。
- **Commit 3B（根 CMake + `CMakePresets.json` + 命名空间切换）**：
  - 根 `CMakeLists.txt` 切换 `project(OpenUWB LANGUAGES CXX)`，导出 `OpenUWB` 包与 `openuwb::*` target，同时保留 `uwb_core`/`uwb_native` 兼容导出（标注 `TODO(openuwb-compat-removal)`）。
  - 更新 `CMakePresets.json` 使 `sourceDir` 直接指向仓库根目录。
  - 库与测试代码切换至 `namespace openuwb`，兼容头与 `namespace gr::uwb` 别名保持向后兼容并标注 `TODO(openuwb-compat-removal)`。
- **Commit 3C（应用目标统一命名 `openuwb_*` + 兼容别名）**：
  - 注册正式应用目标：`openuwb_rx`、`openuwb_pipeline_runner`、`openuwb_capture`、`openuwb_twr_fake_demo`、`openuwb_twr_native_phy_demo`（并在 `apps/openuwb_tx/` 预留清晰目录/README 占位，不引入伪造逻辑）。
  - 通过 CMake post-build symlink 或兼容别名保留 `uwb_dw3000_live`、`uwb_uhd_sc16_pipeline_runner`（标注 `TODO(openuwb-compat-removal)`）。
  - 新增 `tools/_common/paths.py` 并先完成一次 `native-avx2` 增量构建确认全绿。

#### Wave 3：无冲突并行迁移 `tools/`（同时启动 4 个写 Subagent，严格按子目录隔离，禁止执行 `git`）
- **Subagent-W1（仅修改 `tools/rx/`）**：将 `tools/rx/` 脚本接入 `tools/_common/paths.py` 并优先查找新应用名（兼容旧名 + TODO），保证 `qa_live_dw3000_console.py`、`run_distance_replay.py`、`autocorr_suspend_ab.py` 开箱即用。
- **Subagent-W2（仅修改 `tools/sic/`）**：迁移 `tools/sic/` 下脚本与 `.cc` 工具的路径/头文件引用。
- **Subagent-W3（仅修改 `tools/perf/` 与 `tools/performance/`）**：迁移性能/ABBA 基准脚本路径。
- **Subagent-W4（仅修改 `tools/twr/` 与 `tests/integration/install_consumer/`）**：迁移 TWR 验证脚本；扩展 `install_consumer` 同时验证 `find_package(OpenUWB)` 与兼容目标。

#### Wave 4：主控收口提交（Commit 3D）与全量 Gate 验收
- 主控检查全部 `tools/` 语法（`python3 -m py_compile`）并提交 **Commit 3D**。

---

### 强制约束
- **只做 Phase 3**：严禁把 SIC 接入 `openuwb_rx` / `live_rx`；严禁提前删除带 `TODO(openuwb-compat-removal)` 的兼容层。
- **零行为变更**：不改算法逻辑、不放宽 QA 容差、不改变 `Reporter` `stdout` 文本协议。
- 所有构建和测试必须使用 `env -u LD_LIBRARY_PATH`。

### Phase 3 Gate 验收（由主控统一执行并如实汇报）
1. 根目录三预设全量构建与 CTest（`native-avx2`、`native-avx512`、`no-uhd`），全部 100% 通过：
   - `cmake --preset native-avx2 && cmake --build --preset native-avx2 -j && env -u LD_LIBRARY_PATH ctest --preset native-avx2 --output-on-failure`（其余两预设同上）。
2. 运行 `tests/integration/install_consumer`（同时验证 `OpenUWB` 新包与兼容别名）。
3. 运行控制台回归：`env -u LD_LIBRARY_PATH python3 tools/rx/qa_live_dw3000_console.py`（不传旧 `gr-uwb/` 路径，验证默认自动发现）。
4. 如有 `testdata/` 基线夹具，运行一次 `tools/rx/run_distance_replay.py` 关键样本，核对解调包数/FCS 与 Phase 0 基线逐包一致（如无本地夹具则如实说明）。
5. 汇报：Commit 3A~3D hash、各 Subagent 迁移文件数、三个预设 CTest 精确计数、`install_consumer` 与回放结果及最终 `git status --short`。
