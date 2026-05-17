# 关联几何与二聚体模型 — 射影拼图框架中的圆问题研究

基于 arXiv:2505.02229（*Incidences, tilings, and fields*）和 arXiv:2512.15499（*The dimer model and dynamical incidence geometry*），研究如何将射影几何的"拼图证明"框架扩展到涉及圆锥曲线（特别是圆）的关联定理。

## 文件说明

| 文件 | 内容 |
|------|------|
| `2026-05-17-incidences-tilings-fields.md` | arXiv:2505.02229 全文分析——三版主定理、证明强度层级、有限域反例法 |
| `2026-05-17-dimer-model-analysis.md` | arXiv:2512.15499 框架原理——回路条件、相干性、局部变换、圆问题三方案 |
| `2026-05-17-pascal-brianchon-mobius.md` | Möbius 升维分析——Pascal 与 Brianchon 定理在二聚体框架中的编码 |
| `research/pascal_brianchon_mobius.py` | Python 数值验证（蒙特卡洛 100 次，通过率 99%+） |

## 核心结果

1. **Möbius 升维是桥梁**：圆 → 射影 3-空间中的二次曲面 + 平面交线，圆定理变为线性关联定理
2. **圆条件 = 相干性条件，非回路条件**：ℙ³ 中回路最多 4 点，6 共圆点不构成回路
3. **Pascal 定理 = 3 面球面拼图**：与 Desargues 定理的 4 面四面体拼图同属一个框架
4. **Brianchon 通过对偶自动跟随**

## 术语体系

| 英文 | 中文 |
|------|------|
| circuit | 回路 |
| coherence | 相干性 |
| dimer model | 二聚体模型 |
| urban renewal | 城区改造 |
| tiling | 拼图 |
| multi-ratio | 多重比 |
| excision | 切除 |
| simplicial complex | 单纯复形 |
| spectral curve | 谱曲线 |
| cross-ratio | 交比 |
| skew field | 斜域 |

## 作者

klarisa · 独立研究 · 2026-05-17
