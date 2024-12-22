# VectorFields.jl

Finite difference operators for vector fields. Field representation is generic and can be vector of arrays, named tuples of arrays, or dictionaries of arrays. Supports gradient, divergence, curl, and laplacian as well as elementwise arithmetic. Supports orthogonal but possibly nonuniform grids. Grid spacing, boundary condition padding and staggered vs central difference are implicit in operator construction.

- Operators:
    - `Del`
    - `Lap`
- Operations:
    - `+`, `-` overloaded to be recursively elementwise
    - `⊙` for recursively elementwise * instead of matrix multiplication
    - `⊘` for recursively elementwise / instead of matrix division


## Contributors
Paul Shen
pxshen@alumni.stanford.edu
Luminescent AI
