# heat-equation-solvers

用 **5 种语言**分别求解 **1D / 2D / 3D 热传导方程**，结果可相互对比验证。

$$u_t = \alpha \nabla^2 u, \quad \alpha = 1, \quad \Omega = [0,1]^d$$

边界条件：Dirichlet 零边界；初始条件：$\sin(\pi x)[\sin(\pi y)[\sin(\pi z)]]$

精确解：$u = e^{-d\pi^2 t} \sin(\pi x)[\sin(\pi y)[\sin(\pi z)]]$

## 数值方法

显式 FTCS（Forward-Time Centered-Space）有限差分：

$$\frac{u^{n+1}-u^n}{\Delta t} = \alpha \frac{u^n_{i+1}-2u^n_i+u^n_{i-1}}{\Delta x^2}$$

CFL 稳定性条件：$r = \alpha \Delta t / \Delta x^2 \leq 1/(2d)$

## 语言与结果

| 语言 | 1D | 2D | 3D | 结果图 |
|------|----|----|----|----|
| Python | `python/heat1d.py` | `python/heat2d.py` | `python/heat3d.py` | `results/python_*.png` |
| MATLAB | `matlab/heat1d.m` | `matlab/heat2d.m` | `matlab/heat3d.m` | `results/matlab_*.png` |
| Julia | `julia/heat1d.jl` | `julia/heat2d.jl` | `julia/heat3d.jl` | `results/julia_*.png` |
| C++ | `cpp/heat1d.cpp` | `cpp/heat2d.cpp` | `cpp/heat3d.cpp` | CSV 数据 |
| Fortran | `fortran/heat1d.f90` | `fortran/heat2d.f90` | `fortran/heat3d.f90` | CSV 数据 |

## 快速运行

```bash
# Python（直接运行，生成图像）
cd python && python heat1d.py && python heat2d.py && python heat3d.py

# Julia
cd julia && julia heat1d.jl && julia heat2d.jl && julia heat3d.jl

# C++
cd cpp && mkdir build && cd build && cmake .. && make
./heat1d && ./heat2d && ./heat3d

# Fortran
cd fortran && make
./heat1d && ./heat2d && ./heat3d

# MATLAB（在 MATLAB 中运行）
cd matlab
run heat1d.m
run heat2d.m
run heat3d.m
```

## 误差对比（L2，t=0.1，N=50）

| 语言 | 1D L2 误差 | 2D L2 误差 |
|------|-----------|-----------|
| Python | ~8×10⁻⁵ | ~2×10⁻⁴ |
| Julia  | ~8×10⁻⁵ | ~2×10⁻⁴ |
| C++    | ~8×10⁻⁵ | ~2×10⁻⁴ |
| Fortran| ~8×10⁻⁵ | ~2×10⁻⁴ |
| MATLAB | ~8×10⁻⁵ | ~2×10⁻⁴ |

（五种语言使用相同算法、相同参数，误差一致）

## 参数

| 参数 | 值 |
|------|----|
| α（热扩散率）| 1.0 |
| 网格点数 N | 50（1D/2D），20（3D） |
| 终止时刻 T | 0.1 |
| CFL 数 r | 0.4（1D），0.2（2D），0.13（3D） |
