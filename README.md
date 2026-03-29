# matrix_algebra.h

A single-header C++ library for basic linear algebra operations — matrix arithmetic, linear systems, determinants, rank, and matrix inversion.

No dependencies beyond the C++ standard library.

---

## Usage

Just copy `matrix_algebra.h` into your project and include it:

```cpp
#include "matrix_algebra.h"
```

Requires **C++11** or later.

---

## Features

### Integer matrix operations
| Function | Description |
|---|---|
| `somma(A, B)` | Element-wise sum A + B |
| `sottrazione(A, B)` | Element-wise subtraction A - B |
| `moltiplicazione_per_scalare(A, k)` | Scalar multiplication k * A |
| `prodotto_matriciale(A, B)` | Matrix product A * B |
| `prodotto_scalare(a, b)` | Dot product of two integer vectors |

### Linear systems (double)
| Function | Description |
|---|---|
| `risolvi_sistema(A, b)` | Solves Ax = b via Gaussian elimination |
| `risolvi_equazione({a, b})` | Solves ax = b, handles singular cases |

### Matrix analysis
| Function | Description |
|---|---|
| `riduzione_gauss(A)` | Row reduces A to upper triangular form |
| `calcola_determinante(A)` | Determinant from a reduced matrix |
| `calcola_rango(A)` | Rank from a reduced matrix |
| `inversa(A)` | Matrix inverse via Gauss-Jordan elimination |

---

## Quick example

```cpp
#include <iostream>
#include "matrix_algebra.h"

int main() {
    // Solve the system:
    //  2x + y  = 5
    //  x  + 3y = 10
    MatD A = {{2, 1},
              {1, 3}};
    VecD b = {5, 10};

    VecD sol = risolvi_sistema(A, b);
    std::cout << "x = " << sol[0] << ", y = " << sol[1] << "\n";
    // x = 1, y = 3

    // Matrix product (integers)
    MatI M1 = {{1, 2}, {3, 4}};
    MatI M2 = {{5, 6}, {7, 8}};
    MatI P  = prodotto_matriciale(M1, M2);

    // Inverse
    MatD Ainv = inversa(A);

    return 0;
}
```

---

## Notes

- Integer matrix functions use `MatI` (`vector<vector<int>>`).
- Floating-point functions use `MatD` / `VecD` (`vector<vector<double>>`).
- `risolvi_sistema` returns an empty vector if the system is singular.
- `inversa` returns `{{}}` if the matrix is singular or non-square.
- Zero threshold is `EPS = 1e-9`.

---

## License

Do whatever you want with it.
