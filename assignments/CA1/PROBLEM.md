# CA1: One-dimensional elastostatics

Solve the axial displacement of a uniform bar on `0 <= x <= 0.1 m`. The material has Young's modulus `E = 1e11 Pa` and cross-sectional area `A = 1e-4 m^2`; the body-force density varies as `f(x) = 1e11 x N/m^3`.

Implement both boundary-value cases:

1. Fix `u(0) = 0` and prescribe `u(0.1) = 0.001 m`.
2. Fix `u(0) = 0` and apply axial traction `EA u'(0.1) = 1e6 N` at the right end.

Use 100 elements for the submitted solutions. The implementation must support Lagrange orders 1, 2, and 3, integrate the element terms and L2 error using Gaussian quadrature, and produce one HDF5 result per order/case pair.
