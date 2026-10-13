# CA2: Steady heat conduction

Solve the steady conduction equation for an isotropic material with conductivity `kappa = 385 W/(m K)`. Use linear basis functions. All unspecified boundaries are insulated (zero heat flux).

## Two-dimensional case

Use a rectangular domain `0 <= x <= 0.03 m`, `0 <= y <= 0.08 m`, discretized with 15 by 40 quadrilateral elements. Prescribe

- bottom edge: `u(x, 0) = 300(1 + 31x) K`;
- top edge: `u(x, 0.08) = 310(1 + 8x^2) K`.

## Three-dimensional case

Use `0 <= x <= 0.04 m`, `0 <= y <= 0.08 m`, and `0 <= z <= 0.02 m`, with 8 by 16 by 4 hexahedral elements. Set

- left face: `u(0, y, z) = 300(1 + 31(y + z)) K`;
- right face: `u(0.04, y, z) = 310(1 + 31(y + z)) K`.

The values above follow the assignment PDF. Its two-dimensional bottom boundary formula does not match the supplied debug VTK values; this repository intentionally follows the PDF specification.
