# LinearReservoirOptimization

# Optimization Theory

The optimization problem is formulated as:

$$
\begin{aligned}
\min_{\lambda_1, \ldots, \lambda_N, \boldsymbol{\kappa}} \quad & \tilde{\boldsymbol{\epsilon}}_q^\top \tilde{\boldsymbol{\epsilon}}_q + \beta_1 \tilde{\boldsymbol{\kappa}}^\top \tilde{\boldsymbol{\kappa}} + \beta_2 \left( \frac{1}{H} \right) \\
\text{subject to} \quad
& \tilde{\boldsymbol{\epsilon}}_q = \tilde{\Omega} \tilde{\boldsymbol{\kappa}} - \tilde{y}, \\
& \tilde{y} = \begin{bmatrix}
b_1 \cos(\phi_1), & b_1 \sin(\phi_1), & b_2 \cos(\phi_2), & b_2 \sin(\phi_2), & \ldots, & b_K \cos(\phi_K), & b_K \sin(\phi_K)
\end{bmatrix}^\top, \\
& \tilde{\Omega} = \tilde{\Omega}_q, \\
& M_{ik} = \left| \frac{\gamma c_i}{j \omega_k - \gamma(\lambda_i - 1)} \right|, \\
& \theta_{ik} = \angle \frac{\gamma c_i}{j \omega_k - \gamma(\lambda_i - 1)}, \\
& \lambda_i \leq 0, \\
& H = \frac{N}{\sum_{\substack{j,z=1 \\ j \neq z}}^{N} \frac{1}{|\lambda_j - \lambda_z|}}, \\
& \omega_{\max} + \gamma(\lambda_i - 1) \leq 0.
\end{aligned}
$$
