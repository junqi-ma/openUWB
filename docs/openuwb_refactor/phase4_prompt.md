# OpenCode 执行提示词 — Phase 4（串行 ABBA 性能验收、TODO 汇总与文档收尾）

> 适用分支：`refactor/openuwb-arch`（前置条件：Phase 3 全绿）  
> 参考总方案：`docs/OpenCode执行方案_openUWB全仓库架构重构_2026-11-20.md`

---

Phase 3 已全绿。现在请在当前分支 `refactor/openuwb-arch` 上，严格按照 `AGENTS.md` 与 `docs/OpenCode执行方案_openUWB全仓库架构重构_2026-11-20.md` 第 5 节，**执行最后的 Phase 4 收尾与性能验收**。

⚠️ **特别注意多 Subagent 边界**：
- **ABBA 性能测试必须严格串行、独占 CPU 执行，严禁与任何编译/测试/Subagent 并行**（否则 CPU 与缓存争抢会导致性能数据作废）。
- 文档更新与 TODO 扫描可采用多 Subagent 并行，但必须在 ABBA 基准全部跑完之后（或之前完全空闲时）进行。

---

### 多 Subagent 执行编排

#### Wave 1：严格串行 ABBA 性能回归对比（仅主控单线程独占 CPU 执行）
- 对比对象：基线快照（Phase 0 的 Commit 0A）vs 当前 Phase 3 HEAD（`refactor/openuwb-arch`）。
  - 可在临时目录或独立只读 worktree 构建基线二进制，在同一台机器、同一构建参数、同一 ISA（`native-avx2` / `native-avx512`）下按 **ABBA** 顺序串行执行。
- 必测项（均使用 `env -u LD_LIBRARY_PATH`）：
  1. 核心解调/SIC 微基准（`bench_demod_sic_kernels` 等）
  2. 稳态管线基准（`benchmark_sc16_pipeline_steady`）
  3. 实时接收回放吞吐与单包耗时（`openuwb_rx` / `uwb_dw3000_live` 对比基线夹具，记录锁包率、FCS 成功率、单包平均/P95 耗时）
- 必须按 `AGENTS.md` 规范记录：双方 commit hash、CPU 型号、ISA、CMake 预设/编译参数、worker/queue 配置、夹具路径与 SHA256 哈希、原始 ABBA 数据与统计口径。
- 如发现超过噪声水平的明显性能退化（例如 >3%~5%），**立即定位原因并汇报，不得强行跳过**。

#### Wave 2：并行扫描与文档收尾（ABBA 完成后，同时启动 3 个 Subagent）
- **Subagent-D1（只读扫描全部 `TODO(openuwb-compat-removal)`）**：
  - 全仓 `grep` 统计所有 `TODO(openuwb-compat-removal)` 标记（转发头文件、命名空间别名、CMake ALIAS、旧二进制软链、`tools/` 回退路径），生成按类别分组的清除清单表。
- **Subagent-D2（更新构建与架构文档，仅编辑 `docs/`）**：
  - 更新 `docs/原生UHD与C++构建使用说明.md`、`docs/架构_解调与窗口执行.md`、`docs/README.md`，同步新的根目录结构、`CMakePresets.json`、`openuwb_*` 入口、`openuwb/support/log.h` 用法与单向依赖 DAG。
  - 新增 `docs/performance/openUWB架构重构验收报告_2026-11-20.md`，完整写入 Wave 1 的 ABBA 性能对比数据、三预设 CTest 结果与兼容层清除清单。
- **Subagent-D3（更新根入口文档，仅编辑根目录 `README.md`、`AGENTS.md`、`开发状态.md`）**：
  - 更新 `README.md` 为 `openUWB` 定位（含新目录结构与快速构建命令，如 `docs/assets/` 已有 Logo 可一并挂入）。
  - 更新 `AGENTS.md` 中的新目录边界（`openuwb_support`/`dsp`/`phy_hrp`/`sic`/`radio`/`runtime`）与默认验证命令。
  - 更新 `开发状态.md`：记录 Phase 0~4 重构完成状态、集中登记 `TODO(openuwb-compat-removal)` 清单，并明确标注“下一步在独立分支基于 `openuwb_runtime::PacketProcessor` + `openuwb_sic` 合入 SIC、以及扩展 `openuwb_tx` / 4z / 4ab”。

#### Wave 3：主控收口提交（Commit 4A）与最终总验收
- 检查所有文档链接有效性（无断链、无残留 `gr-uwb/` 死链接）。
- 提交 **Commit 4A**（文档与验收报告）。
- 最终跑一次 `env -u LD_LIBRARY_PATH ctest --preset native-avx2 --output-on-failure` 与 `git status --short` 确认仓库洁净。

---

### Phase 4 Gate 最终交付汇报
1. ABBA 性能回归对比表（含 commit、CPU、ISA、夹具哈希、ABBA 原始值与差异百分比）。
2. 全量 CTest（`native-avx2` / `native-avx512` / `no-uhd`）+ `install_consumer` + `qa_live_dw3000_console.py` 最终结果。
3. 全仓 `TODO(openuwb-compat-removal)` 总数与分类摘要。
4. 确认分支 `refactor/openuwb-arch` 已达到合并主分支标准，并给出后续仓库改名 `openUWB` 与另开分支合入 SIC 的下一步指令。
