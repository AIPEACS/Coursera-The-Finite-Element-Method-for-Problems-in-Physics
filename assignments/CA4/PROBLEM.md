# CA4: Steady and transient heat conduction

Model heat conduction in a slab with dimensions `1 m x 1 m x 0.1 m`. The volumetric heat capacity is `rho = 3.8151e6 J/(m^3 K)`, and the isotropic conductivity is `kappa = 385 W/(m K)`. The faces at `x = 0` and `x = 1 m` are held at 300 K and 310 K respectively; all other faces are insulated.

Solve the steady problem and the transient problem. For the transient initial state, set the temperature to 300 K for `x < 0.5 m`, then use `300 + 20(x - 0.5) K` for `x >= 0.5 m`.

Use linear basis functions on a 20 by 20 by 1 mesh. Advance to `t = 3000 s` with the v-method for `alpha = 0` (Forward Euler), `alpha = 0.5` (Mid-Point), and `alpha = 1` (Backward Euler). For each value, save the transient field and the L2 difference from the steady solution at the final time.
