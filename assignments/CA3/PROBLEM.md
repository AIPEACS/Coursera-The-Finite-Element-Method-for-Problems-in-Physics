# CA3: Three-dimensional elasticity

Find the displacement field in a unit cube using linear finite elements. The material is isotropic, with `E = 2e11 Pa` and Poisson's ratio `nu = 0.3`.

Clamp the entire bottom face (`z = 0`) in all three displacement components. On the top face (`z = 1 m`), apply only a vertical traction whose third component is `h_3 = 1e9 x_1 Pa`; the other traction components are zero. All other faces carry zero traction.

Use a 10 by 10 by 10 hexahedral mesh for the submitted result.
