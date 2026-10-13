# The Finite Element Method for Problems in Physics

This repository contains the C++/deal.II implementations for the four Coursera coding assignments. Each assignment keeps the supplied class and function interfaces, build file, and solution writer. The supplied debugging VTK files are under each assignment's `reference/` directory.

## Problem summaries

These independently written summaries describe the four tasks without redistributing the course PDFs:

- [CA1: One-dimensional elastostatics](assignments/CA1/PROBLEM.md)
- [CA2: Steady heat conduction](assignments/CA2/PROBLEM.md)
- [CA3: Three-dimensional elasticity](assignments/CA3/PROBLEM.md)
- [CA4: Steady and transient heat conduction](assignments/CA4/PROBLEM.md)

## Required Coursera zip contents

Create one zip per assignment and include only the listed files:

| Zip | Files to include |
|---|---|
| `CA1.zip` | `FEM1.h` and six `.h5` files, one per basis order (1, 2, 3) and problem (1, 2) |
| `CA2.zip` | `FEM2a.h`, `FEM2b.h`, `CA2a.h5`, and `CA2b.h5` |
| `CA3.zip` | `FEM3.h` and `CA3.h5` |
| `CA4.zip` | `FEM4.h` and three `.h5` files for α = 0, 0.5, and 1 |

The `.h5` names are set by the supplied `main*.cc` files. Run each assignment from its own build directory so its generated `.h5` and VTK outputs do not overwrite another assignment's outputs. Do not include the VTK files, executable, CMake files, or source `main*.cc` files in the Coursera zips.

## Submission mesh inputs

- CA1: 100 one-dimensional elements. Set the `order` and `problem` inputs in `main1.cc` to generate each of the six results.
- CA2a: 15 × 40 quadrilateral elements.
- CA2b: 8 × 16 × 4 hexahedral elements.
- CA3: 10 × 10 × 10 hexahedral elements.
- CA4: 20 × 20 × 1 hexahedral elements. Set `alpha` in `main4.cc` to 0, 0.5, and 1 and run once for each value.

## Build and verification

Build each assignment independently with a deal.II installation discoverable by CMake:

```sh
cmake -S assignments/CA1 -B build/CA1 -DDEAL_II_DIR=/path/to/deal.II
cmake --build build/CA1
```

Use the same command pattern for `assignments/CA2/CA2a`, `assignments/CA2/CA2b`, `assignments/CA3`, and `assignments/CA4`, with a separate build directory for each. The `main*.cc` inputs are left at the required submission meshes.

For debug checks, temporarily use the mesh sizes documented in the assignment PDFs and compare the generated VTK fields with `reference/`. CA1's 10-element L2 error checks are `1.66468e-07` (linear), `1.81848e-09` (quadratic), and near machine precision (cubic), for both problems. CA4's debug reference uses a 10 × 10 × 1 mesh; its listed L2 values are in the assignment PDF.

`writeSolutions.h` and `CMakeLists.txt` are kept as supplied. The source headers contain the finite-element assembly, Gaussian integration, boundary conditions, and time stepping.

## License and attribution

The MIT license in `LICENSE` applies to original repository contributions. The supplied course template code remains attributed under the CC BY 4.0 notice included with the template downloads; see the [Creative Commons Attribution 4.0 license](https://creativecommons.org/licenses/by/4.0/). The course PDFs are not included in this repository.
