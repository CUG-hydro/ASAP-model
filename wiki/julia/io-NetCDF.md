# io/NetCDF — NetCDF 只读 I/O

> 源文件：`src/io/NetCDF.jl`
> Fortran 来源：`fortran/module_io.f90`（子集：`READINITIAL` L10、`READWTDNC` L309）
> 测试：`test/test_io_netcdf.jl`
> 状态：已摄取（仅 P0 子集；P1/P2 子程序留待后续 PR）

## 1. 功能概述

`module NetCDFIO` 是 ASAP-model 的 NetCDF **读取**入口。本 PR 只翻译 Fortran
`module_io.f90` 中最高优先级（**P0**）的两个子程序：

- `READINITIAL` → `read_initial(path)`：读取土壤质地、地形、根系深度因子与
  陆地掩码四个静态场；
- `READWTDNC`   → `read_wtdnc(path)`：读取地下水位初始场并施加
  `min(-raw, 0.)` 截断。

P1/P2 子程序（`READLATLON`、`READVEG`/`READHVEG`、`READFLOWDIRECTION`、
`READRIVERPARAMETERS`、`READHISTORYNC` 等）**不在本 PR 范围**，待后续 PR 增量
接入。Module 同样不实现任何写入路径。

## 2. 函数签名

```julia
module NetCDFIO
using NCDatasets
export read_initial, read_wtdnc

function read_initial(path::String) -> NamedTuple
function read_wtdnc(path::String)    -> Matrix{Float64}
end
```

## 3. 算法 / 公式

### 3.1 静态场（`read_initial`）

四个 NetCDF 变量名沿用 Fortran `nf90_inq_varid` 调用约定：

| Fortran 变量 | NetCDF 变量名 | Julia 类型 |
|---|---|---|
| 土壤质地 (STXT) | `"STXT"` | `Array{Int}` |
| 地形 (topo) | `"topo"` | `Array{Float64}` |
| 根系深度因子 (F) | `"F"` | `Array{Float64}` |
| 陆地掩码 | （派生） | `Array{Int}` |

> **注意**：原始 Fortran 代码中 `READINITIAL` 通过 *二进制直读*（`access='direct'`）
> 加载地形，NetCDF 变量名 `"topo"` 是为 Julia 端口所选定的命名约定。如生产数据
> 集使用不同变量名，需要相应调整 `src/io/NetCDF.jl` 中 `ds["topo"][:, :]`
> 这一行。

后处理步骤完全镜像 Fortran 中的标量化代码：

```fortran
landmask = 1
where (topo .lt. -1.e5) landmask = 0
where (topo .lt. -1.e5) topo = 0.

where (fdepth .lt. 1.e-6) fdepth = 100.
```

用 Julia 的 broadcast 表示：

```julia
TOPO_MISSING_THRESHOLD = -1.0e5
FDEPTH_MISSING_THRESHOLD = 1.0e-6
FDEPTH_MISSING_REPLACE   = 100.0

landmask = ifelse.(topo .< TOPO_MISSING_THRESHOLD, 0, 1)
topo     = ifelse.(topo .< TOPO_MISSING_THRESHOLD, 0.0, topo)
fdepth   = ifelse.(fdepth_raw .< FDEPTH_MISSING_THRESHOLD,
                   FDEPTH_MISSING_REPLACE, fdepth_raw)
```

### 3.2 地下水位（`read_wtdnc`）

直接对应 Fortran `module_io.f90:309`：

```fortran
varread(1:n2, 1:n3) = min(-varreadbig(nw:ne, ns:nn), 0.)
```

```julia
raw = Array{Float64}(ds["WTD"][:, :])
wtd = min.(-raw, 0.0)   # 负值=低于地表；正值（陆地表面以上）被夹到 0
```

返回的 `wtd` 矩阵保证 `≤ 0`，与 ASAP 主算法中"地下水位深度"约定一致
（参见 `wiki/conventions.md §1`）。

## 4. 关键变量与单位

| 符号 | 含义 | 单位 |
|---|---|---|
| `STXT` | USDA 土壤质地类型 1..13 | —（枚举） |
| `topo` | 地表高程（海洋/缺测格点已置 0） | m |
| `F` | 根系深度因子（`fdepth`，< 1e-6 替换为 100） | m |
| `landmask` | 0 = 海洋，1 = 陆地 | —（掩码） |
| `WTD` | 地下水位深度（已施加 `min(-raw, 0.)`） | m（≤ 0） |

阈值常量（导出为模块内部 `const`）：

| 常量 | 值 | 含义 |
|---|---|---|
| `TOPO_MISSING_THRESHOLD` | `-1.0e5` | 地形 < 此值视为海洋/缺测 |
| `FDEPTH_MISSING_THRESHOLD` | `1.0e-6` | 根深因子 < 此值视为缺测 |
| `FDEPTH_MISSING_REPLACE` | `100.0` | 缺测根深因子的替换值（m） |

## 5. 与 Fortran 对应

| Fortran 子程序 | Julia 函数 | 差异 |
|---|---|---|
| `READINITIAL` (`fortran/module_io.f90:10`) | `read_initial` | Fortran 串行 + MPI 分发；Julia 单进程 + NetCDF 一次性读取 |
| `READWTDNC` (`fortran/module_io.f90:309`) | `read_wtdnc` | 同上；并省略 MPI 广播与子区域裁剪 |
| `READLATLON` (L212) | — | **P1，未实现** |
| `READVEG` / `READHVEG` (L378 / L447) | — | **P1，未实现** |
| `READSMOIEQ` (L517) | — | **P2，未实现** |
| `READFLOWDIRECTION` (L591) | — | **P2，未实现** |
| `READRIVERPARAMETERS` (L693) | — | **P2，未实现** |
| `READHISTORYNC` 等 | — | **P2（restart），未实现** |
| `WRITEOUTPUT*` | — | **整段写出**未实现（与本 PR 范围无关） |

详细对照见 `wiki/mapping/julia-fortran-对照.md` §3 NetCDF 章节（待补）。

## 6. 引用

- 行号：
  - `src/io/NetCDF.jl:L41-L92` `read_initial`
  - `src/io/NetCDF.jl:L96-L113` `read_wtdnc`
  - `src/ASAP.jl:L31-L34` include + export
- 测试断言：
  - `test/test_io_netcdf.jl:L11-L25`  `read_wtdnc` 圆环测试
  - `test/test_io_netcdf.jl:L27-L46`  `read_wtdnc` 边界（混合正负）
  - `test/test_io_netcdf.jl:L48-L75`  `read_initial` 圆环测试
  - `test/test_io_netcdf.jl:L77-L97`  `read_initial` 根深因子夹断
  - `test/test_io_netcdf.jl:L99-L120` `read_initial` 地形/掩码
- 文档：`docs/*.typ` 暂未涉及（I/O 章节若需可后续补 Typst 推导）。

## 7. 已知问题与备注

- **NCDatasets API 选择**：本 PR 选用 `NCDatasets = "0.12, 0.13, 0.14"`，
  `ds[name][:, :]` 而非 `NCDatasets.readvar`（后者在 v0.14 已不再导出）。
  v0.14 中 `v[:]` 会展平成 1D 向量，必须用 `v[:, :]` 保留二维形状。
- **地形变量名约定**：Fortran 通过二进制直读加载地形，NetCDF 端无变量名
  约定；Julia 端口暂以 `"topo"` 为默认名，若与生产数据集不匹配需调整
  `src/io/NetCDF.jl:55`。
- **不导出读子程序内部 helper**：模块内只 `export read_initial, read_wtdnc`；
  `TOPO_MISSING_THRESHOLD` 等阈值常量以模块内 `const` 形式存在，**不 export**，
  调用方如需访问应通过 `ASAP.NetCDFIO.TOPO_MISSING_THRESHOLD`。
- **MPI 行为省略**：原 Fortran 子程序带 MPI 分发与子区域裁剪
  （`nw:ne, ns:nn`）；Julia 端口假设单进程/整域读取。生产部署如需多进程 I/O
  需重新设计。
- **未实现写入**：所有 `WRITEOUTPUT*` / `WRITEHISTORYNC*` 子程序均**未**
  翻译，与本 PR 范围分离。
- **`src/ASAP.jl` include 顺序**：NetCDF 子模块在 `RootDepth.jl` 之后引入，
  不依赖任何业务子模块；位置安全，未来可前移。
