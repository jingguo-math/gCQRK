# Runge--Kutta generalized Convolution Quadrature for sectorial problems

MATLAB implementation accompanying the paper:

> **Runge–Kutta Generalized Convolution Quadrature for Sectorial Problems**  
> Jing Guo, Maria Lopez-Fernandez  
> arXiv:2506.21242, 2025  
> https://arxiv.org/abs/2506.21242

---
## Abstract

We study the application of the generalized convolution quadrature (gCQ) based on Runge--Kutta methods to approximate the solution of an important class of sectorial problems. The gCQ  generalizes Lubich's original convolution quadrature (CQ) to variable steps. High-order versions of the gCQ  have been developed in the last decade, relying on certain Runge--Kutta methods. The Runge--Kutta based gCQ  has been studied so far in a rather general setting, which includes applications to boundary integral formulations of wave problems. The available stability and convergence results for these new methods are suboptimal compared to those known for the uniform-step CQ, both in terms of convergence order and regularity requirements of the data. Here we focus on a special class of sectorial problems and prove that in these important applications it is possible to achieve the same order of convergence as for the original CQ, under the same regularity hypotheses on the data, and for very general time meshes. In the particular case of data with some known algebraic type of singularity, we also show how to choose an optimally graded time mesh to achieve convergence with maximal order, overcoming the well-known order reduction of the original CQ in these situations. An important advantage of the gCQ  method is that it allows for a fast and memory-efficient implementation. We describe how the fast and oblivious Runge--Kutta based gCQ  can be implemented and illustrate our theoretical results with several numerical experiments.

## What this code does

This repository implements the **generalized Convolution Quadrature (gCQ)** based on Runge–Kutta (RK) methods for numerically approximating Volterra-type convolution integrals of the form

$$u(t) = \int_0^t k(t-s) f(s) ds, \quad t > 0,$$

on **non-uniform time meshes**, where $k(t)$ can be scalar-valued or operator-valued, depending on the problem setting.

---

## Repository structure

```
gCQRK/
├── quadratures/                          % Quadrature weights and nodes used by the solver
│
├── examples/
│   ├── Frac_int/
│   │   └── test_cqrk_varn0.m            % Example 1 — Fractional  integral
│   │
│   ├── ODE_heaviside/
│   │   └── test_cqrk_varn0ODERV.m       % Example 3 — Fractional ordinary diffusion equation with Heaviside forcing
│   │
│   ├── conv_genKa/
│   │   └── test_cqrk_varn0_genk_exact.m % Example 2 — Convolution with kernel k_a(t)
│   │
│   ├── conv_genKb/
│   │   └── test_cqrk_varn0_genkb.m      % Example 2 — Convolution with kernel k_b(t)
│   │
│   ├── PDE/
│   │   └── FracDiff_2D/
│   │       └── test_gcqrk_PDE2d.m       % Example 4 — 2D fractional diffusion PDE
│   │
│   └── Westervelt/
│       └── test_gCQRK_WestveltFPI.m     % Example 5 — Nonlinear Westervelt equation
│
├── LICENSE
└── README.md
```

---


## Paper examples and corresponding scripts

| Paper example | Problem | Script |
|---|---|---|
| **Example 1** | Fractional convolution integral, $K(z) = z^{-\alpha}$ | [`examples/Frac_int/test_cqrk_varn0.m`](examples/Frac_int/test_cqrk_varn0.m) |
| **Example 2** ($k_a$) | Convolution integral with kernel $k_a(t)$ | [`examples/conv_genKa/test_cqrk_varn0_genk_exact.m`](examples/conv_genKa/test_cqrk_varn0_genk_exact.m) |
| **Example 2** ($k_b$) | Convolution integral with kernel $k_b(t)$ | [`examples/conv_genKb/test_cqrk_varn0_genkb.m`](examples/conv_genKb/test_cqrk_varn0_genkb.m) |
| **Example 3** | Fractional ordinary diffusion equation with Heaviside forcing | [`examples/ODE_heaviside/test_cqrk_varn0ODERV.m`](examples/ODE_heaviside/test_cqrk_varn0ODERV.m) |
| **Example 4** | 2D fractional diffusion PDE | [`examples/PDE/FracDiff_2D/test_gcqrk_PDE2d.m`](examples/PDE/FracDiff_2D/test_gcqrk_PDE2d.m) |
| **Example 5** | Nonlinear damped Westervelt equation | [`examples/Westervelt/test_gCQRK_WestveltFPI.m`](examples/Westervelt/test_gCQRK_WestveltFPI.m) |

---

## Example descriptions

### Example 1 — Fractional convolution integral (`examples/Frac_int/test_cqrk_varn0.m`)

Approximates

$$u(t) = \int_0^t \frac{(t-s)^{\alpha-1}}{\Gamma(\alpha)} f(s) ds, \quad \alpha \in (0,1),$$

where the transfer operator is $K(z) = z^{-\alpha}$ and the spectral density is

$$G(x) = \frac{\sin(\pi\alpha)}{\pi} x^{-\alpha}.$$

---

### Example 2 — Convolution with kernels $k_a$ and $k_b$ (`conv_genKa/test_cqrk_varn0_genk_exact.m`, `conv_genKb/test_cqrk_varn0_genkb.m`)

Computes the gCQ-RK approximation of the convolution integral with the following two kernels:

$$k_a(t) = -\frac{d}{dt} E_{\alpha,1}(-t^\alpha), \qquad k_b(t) = \frac{t^{\alpha-1}}{\Gamma(\alpha)} e^{-t},$$

acting on $f(t) = t^\beta$.

---

### Example 3 — Scalar ODE with Heaviside forcing (`examples/ODE_heaviside/test_cqrk_varn0ODERV.m`)

Solves the fractional ordinary diffusion equation

$$D_t^\alpha u(t) + u(t) = f(t), \quad u(0) = 1, \quad 0 < t \leq 1,$$

with exact solution

$$u(t) = 1 + t^{\beta_1} + H(t - \sigma)(t - \sigma)^{\beta_2},$$

where $\beta_1 > -1$, $\beta_2 > \alpha$, and $H(\cdot)$ denotes the Heaviside function. The source term is

$$f(t) = u(t) + \frac{\Gamma(\beta_1 + 1)}{\Gamma(\beta_1 - \alpha + 1)} t^{\beta_1 - \alpha} + \frac{\Gamma(\beta_2 + 1)}{\Gamma(\beta_2 - \alpha + 1)} H(t - \sigma)(t - \sigma)^{\beta_2 - \alpha},$$

exhibiting piecewise smooth behavior.

---

### Example 4 — 2D fractional diffusion PDE (`examples/PDE/FracDiff_2D/test_gcqrk_PDE2d.m`)

Solves the 2D subdiffusion equation

$$\partial_t^\alpha u - u_{xx} - u_{yy} = f, \quad (x,y) \in \Omega,$$

with initial condition $u(x,y,0) = 0$ and homogeneous Dirichlet boundary conditions. The exact solution is

 $$u(x,y,t) = t^\alpha \cos\left(\frac{\pi}{2}x\right) \cos\left(\frac{\pi}{2}y\right),$$


and the source term is

$$f(x,y,t) = \left(\Gamma(\alpha+1) + \frac{\pi^2}{2} t^\alpha\right) \cos\!\left(\frac{\pi}{2}x\right) \cos\!\left(\frac{\pi}{2}y\right).$$


---

### Example 5 — Nonlinear Westervelt equation (`examples/Westervelt/test_gCQRK_WestveltFPI.m`)

Solves the nonlinear Westervelt equation

$$(1-2\kappa u)\partial_t^2 u - u_{xx} - (k \ast \partial_t u_{xx}) = 2\kappa(\partial_t u)^2 + f, \quad x \in (-8,8),\quad t \in (0,2],$$

with boundary conditions $u(\pm 8, t) = 0$ and initial conditions $u(x,0) = e^{-x^2}/2$, $\partial_t u(x,0) = 0$, where

$$k(t) = \frac{t^{\alpha-1}}{\Gamma(\alpha)} e^{-t}, \qquad f(x,t) = \bigl(1 + \log(t)\bigr) e^{-x^2}.$$

---

## Citation

If you use this code, please cite:

```bibtex
@misc{GuoLo25,
  title         = {Runge--Kutta Generalized Convolution Quadrature for Sectorial Problems},
  author        = {Jing Guo and Maria Lopez-Fernandez},
  year          = {2025},
  eprint        = {2506.21242},
  archivePrefix = {arXiv},
  primaryClass  = {math.NA},
  url           = {https://arxiv.org/abs/2506.21242}
}
```

---

## License

GPL-3.0 — see [LICENSE](LICENSE) for details.

---

## Authors

- **Jing Guo** — School of Mathematics and Statistics, Guangdong University of Technology (jingguo@gdut.edu.cn)
- **Maria Lopez-Fernandez** — University of Malaga (maria.lopezf@uma.es)
