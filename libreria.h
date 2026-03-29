/**
 * @file matrix_algebra.h
 * @brief Libreria di algebra lineare: operazioni su matrici e sistemi lineari.
 *
 * Fornisce funzioni per:
 *  - Operazioni tra matrici intere (somma, sottrazione, prodotto, scalare)
 *  - Risoluzione di sistemi lineari con eliminazione gaussiana
 *  - Calcolo del determinante e del rango
 *  - Calcolo della matrice inversa (Gauss-Jordan)
 *
 * @note Le matrici intere sono rappresentate come vector<vector<int>>,
 *       quelle reali come vector<vector<double>>.
 *
 * @author  [Il tuo nome]
 * @version 1.1
 */

#pragma once

#include <vector>
#include <cmath>
#include <utility>   // pair
#include <limits>    // INFINITY

// ─────────────────────────────────────────────────────────────────────────────
//  Alias di comodità
// ─────────────────────────────────────────────────────────────────────────────

using MatI = std::vector<std::vector<int>>;
using MatD = std::vector<std::vector<double>>;
using VecD = std::vector<double>;

// ─────────────────────────────────────────────────────────────────────────────
//  Costanti interne
// ─────────────────────────────────────────────────────────────────────────────

static constexpr double EPS = 1e-9;   ///< Soglia per confronti con zero

// ─────────────────────────────────────────────────────────────────────────────
//  1. OPERAZIONI SU MATRICI INTERE
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Somma elemento per elemento due matrici intere.
 * @param A Prima matrice (m x n).
 * @param B Seconda matrice (deve avere le stesse dimensioni di A).
 * @return Matrice C = A + B, oppure matrice vuota se le dimensioni non corrispondono.
 */
inline MatI somma(const MatI& A, const MatI& B) {
    if (A.size() != B.size()) return {};

    MatI C(A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        if (A[i].size() != B[i].size()) return {};
        C[i].resize(A[i].size());
        for (size_t j = 0; j < A[i].size(); ++j)
            C[i][j] = A[i][j] + B[i][j];
    }
    return C;
}

/**
 * @brief Moltiplica ogni elemento di una matrice intera per uno scalare.
 * @param A Matrice di input.
 * @param k Scalare intero.
 * @return Matrice k * A.
 */
inline MatI moltiplicazione_per_scalare(MatI A, int k) {
    for (auto& riga : A)
        for (auto& elem : riga)
            elem *= k;
    return A;
}

/**
 * @brief Sottrae elemento per elemento la matrice B dalla matrice A.
 * @param A Matrice a sinistra.
 * @param B Matrice a destra (stesse dimensioni di A).
 * @return Matrice C = A - B.
 */
inline MatI sottrazione(const MatI& A, const MatI& B) {
    MatI neg = moltiplicazione_per_scalare(B, -1);
    return somma(A, neg);
}

/**
 * @brief Prodotto scalare tra due vettori interi (dot product).
 * @param A Primo vettore.
 * @param B Secondo vettore (stessa lunghezza di A).
 * @return Valore intero del prodotto scalare.
 */
inline int prodotto_scalare(const std::vector<int>& A, const std::vector<int>& B) {
    int risultato = 0;
    for (size_t i = 0; i < A.size(); ++i)
        risultato += A[i] * B[i];
    return risultato;
}

/**
 * @brief Prodotto righe-per-colonne (prodotto matriciale) tra due matrici intere.
 *
 * Richiede che il numero di colonne di A sia uguale al numero di righe di B.
 *
 * @param A Matrice (m x p).
 * @param B Matrice (p x n).
 * @return Matrice C = A * B di dimensione (m x n),
 *         oppure {{0}} se le dimensioni sono incompatibili.
 */
inline MatI prodotto_matriciale(const MatI& A, const MatI& B) {
    if (A.empty() || B.empty() || A[0].size() != B.size()) return {{0}};

    size_t m = A.size(), p = B.size(), n = B[0].size();
    MatI C(m, std::vector<int>(n, 0));

    // Estrai le colonne di B per accesso cache-friendly
    MatI colonneB(n, std::vector<int>(p));
    for (size_t j = 0; j < n; ++j)
        for (size_t k = 0; k < p; ++k)
            colonneB[j][k] = B[k][j];

    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            C[i][j] = prodotto_scalare(A[i], colonneB[j]);

    return C;
}

// ─────────────────────────────────────────────────────────────────────────────
//  2. UTILITY PER SISTEMI LINEARI (DOUBLE)
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Risolve un'equazione lineare ax = b e restituisce x = b/a.
 *
 * Casi particolari:
 *  - a ≈ 0 e b ≈ 0  → NAN  (equazione identica, infiniti valori)
 *  - a ≈ 0 e b ≠ 0  → INFINITY (equazione impossibile)
 *
 * @param eq Coppia (a, b) che rappresenta l'equazione ax = b.
 * @return Soluzione x, NAN oppure INFINITY.
 */
inline double risolvi_equazione(std::pair<double, double> eq) {
    if (std::fabs(eq.first) < EPS && std::fabs(eq.second) < EPS) return NAN;
    if (std::fabs(eq.first) < EPS)                                return INFINITY;
    return eq.second / eq.first;
}

/**
 * @brief Somma elemento per elemento due vettori double (usata internamente).
 * @param v1 Primo vettore.
 * @param v2 Secondo vettore (stessa lunghezza).
 * @return Vettore somma v1 + v2.
 */
inline VecD somma_vettori(const VecD& v1, const VecD& v2) {
    VecD v3(v1.size());
    for (size_t i = 0; i < v1.size(); ++i)
        v3[i] = v1[i] + v2[i];
    return v3;
}

// ─────────────────────────────────────────────────────────────────────────────
//  3. SISTEMI LINEARI
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Risolve un sistema lineare Ax = b con eliminazione gaussiana.
 *
 * Il sistema è definito dalla matrice dei coefficienti @p incognite (n x n)
 * e dal vettore dei termini noti @p soluzioni (n).
 *
 * @param incognite Matrice quadrata dei coefficienti.
 * @param soluzioni Vettore dei termini noti.
 * @return Vettore delle soluzioni x, oppure vettore vuoto se il sistema
 *         è singolare o privo di soluzione.
 */
inline VecD risolvi_sistema(MatD incognite, VecD soluzioni) {
    int n    = static_cast<int>(incognite.size());
    int cols = static_cast<int>(incognite[0].size());
    VecD sol(n);

    for (int i = 0; i < n; ++i) {

        // Pivoting parziale: cerca la riga con il valore assoluto massimo
        if (std::fabs(incognite[i][i]) < EPS) {
            for (int k = i + 1; k < n; ++k) {
                if (std::fabs(incognite[k][i]) > EPS) {
                    std::swap(incognite[i], incognite[k]);
                    std::swap(soluzioni[i], soluzioni[k]);
                    break;
                }
            }
        }

        // Matrice singolare: sistema impossibile o indeterminato
        if (std::fabs(incognite[i][i]) < EPS) return {};

        // Elimina i coefficienti sotto il pivot
        for (int M = i + 1; M < n; ++M) {
            if (std::fabs(incognite[M][i]) < EPS) continue;

            double m = -incognite[M][i] / incognite[i][i];

            VecD riga_pivot(cols + 1), riga_M(cols + 1);
            for (int c = 0; c < cols; ++c) {
                riga_pivot[c] = incognite[i][c] * m;
                riga_M[c]     = incognite[M][c];
            }
            riga_pivot[cols] = soluzioni[i] * m;
            riga_M[cols]     = soluzioni[M];

            VecD somma = somma_vettori(riga_pivot, riga_M);

            for (int c = 0; c < cols; ++c)
                incognite[M][c] = somma[c];
            soluzioni[M] = somma[cols];
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; --i) {
        double s = soluzioni[i];
        for (int j = i + 1; j < n; ++j)
            s -= incognite[i][j] * sol[j];
        sol[i] = s / incognite[i][i];
    }

    return sol;
}

// ─────────────────────────────────────────────────────────────────────────────
//  4. DETERMINANTE E RANGO
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Riduce una matrice a forma triangolare superiore (eliminazione di Gauss)
 *        e restituisce la matrice ridotta insieme al numero di scambi di righe.
 *
 * Usata internamente per calcolare determinante e rango.
 *
 * @param incognite Matrice quadrata di input.
 * @return Coppia {matrice triangolare, numero di scambi},
 *         oppure {{{}}, 0} se la matrice è singolare.
 */
inline std::pair<MatD, int> riduzione_gauss(MatD incognite) {
    int n    = static_cast<int>(incognite.size());
    int cols = static_cast<int>(incognite[0].size());
    int scambi = 0;

    for (int i = 0; i < n; ++i) {

        if (std::fabs(incognite[i][i]) < EPS) {
            for (int k = i + 1; k < n; ++k) {
                if (std::fabs(incognite[k][i]) > EPS) {
                    std::swap(incognite[i], incognite[k]);
                    ++scambi;
                    break;
                }
            }
        }

        if (std::fabs(incognite[i][i]) < EPS) return {{{}}, 0};

        for (int M = i + 1; M < n; ++M) {
            if (std::fabs(incognite[M][i]) < EPS) continue;

            double m = -incognite[M][i] / incognite[i][i];
            VecD rigaPivot(cols), rigaM(cols);

            for (int c = 0; c < cols; ++c) {
                rigaPivot[c] = incognite[i][c] * m;
                rigaM[c]     = incognite[M][c];
            }

            VecD s = somma_vettori(rigaPivot, rigaM);
            for (int c = 0; c < cols; ++c)
                incognite[M][c] = s[c];
        }
    }

    return {incognite, scambi};
}

/**
 * @brief Calcola il determinante di una matrice triangolare superiore.
 *
 * Il determinante è il prodotto degli elementi sulla diagonale principale.
 *
 * @param incognite Matrice triangolare superiore (già ridotta con riduzione_gauss).
 * @return Valore del determinante.
 *
 * @note Tenere conto del segno: se riduzione_gauss ha eseguito un numero dispari
 *       di scambi, il determinante va negato.
 */
inline double calcola_determinante(const MatD& incognite) {
    double det = 1.0;
    for (size_t i = 0; i < incognite.size(); ++i)
        det *= incognite[i][i];
    return det;
}

/**
 * @brief Verifica se un vettore double è il vettore zero.
 * @param v Vettore da controllare.
 * @return true se tutti gli elementi sono esattamente 0.0, false altrimenti.
 */
inline bool e_zero(const VecD& v) {
    for (double x : v)
        if (x != 0.0) return false;
    return true;
}

/**
 * @brief Calcola il rango di una matrice già ridotta a forma triangolare.
 *
 * Il rango è il numero di righe non nulle.
 *
 * @param incognite Matrice ridotta (output di riduzione_gauss).
 * @return Rango della matrice.
 */
inline double calcola_rango(const MatD& incognite) {
    int righe_nulle = 0;
    for (const auto& riga : incognite)
        if (e_zero(riga)) ++righe_nulle;
    return static_cast<double>(incognite.size()) - righe_nulle;
}

// ─────────────────────────────────────────────────────────────────────────────
//  5. MATRICE INVERSA (GAUSS-JORDAN)
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Crea la matrice identità n x n.
 * @param n Dimensione.
 * @return Matrice identità.
 */
inline MatD matrice_identita(int n) {
    MatD I(n, VecD(n, 0.0));
    for (int i = 0; i < n; ++i) I[i][i] = 1.0;
    return I;
}

/**
 * @brief Calcola la matrice inversa tramite il metodo di Gauss-Jordan.
 *
 * Applica le stesse operazioni elementari di riga sia alla matrice A
 * sia alla matrice identità affiancata, finché A diventa I.
 * Il risultato è l'inversa di A.
 *
 * @param A Matrice quadrata da invertire.
 * @return Matrice inversa di A, oppure {{}} se A è singolare o non quadrata.
 */
inline MatD inversa(MatD A) {
    int n    = static_cast<int>(A.size());
    int cols = static_cast<int>(A[0].size());

    if (n != cols) return {{}};   // Non quadrata

    MatD I = matrice_identita(n);

    // — Fase 1: eliminazione verso il basso (forma triangolare superiore) —
    for (int i = 0; i < n; ++i) {

        // Pivoting parziale
        if (std::fabs(A[i][i]) < 1e-12) {
            for (int k = i + 1; k < n; ++k) {
                if (std::fabs(A[k][i]) > 1e-12) {
                    std::swap(A[i], A[k]);
                    std::swap(I[i], I[k]);
                    break;
                }
            }
        }

        if (std::fabs(A[i][i]) < 1e-12) return {{}};   // Matrice singolare

        // Normalizza la riga del pivot (A[i][i] → 1)
        double pivot = A[i][i];
        for (int c = 0; c < n; ++c) {
            A[i][c] = A[i][c] / pivot;
            I[i][c] = I[i][c] / pivot;
        }

        // Elimina i coefficienti sotto il pivot
        for (int M = i + 1; M < n; ++M) {
            if (std::fabs(A[M][i]) < 1e-12) continue;
            double m = -A[M][i];
            for (int c = 0; c < n; ++c) {
                A[M][c] += m * A[i][c];
                I[M][c] += m * I[i][c];
            }
        }
    }

    // — Fase 2: eliminazione verso l'alto (forma identità) —
    for (int i = n - 1; i >= 0; --i) {
        for (int M = i - 1; M >= 0; --M) {
            if (std::fabs(A[M][i]) < 1e-12) continue;
            double m = -A[M][i];
            for (int c = 0; c < n; ++c) {
                A[M][c] += m * A[i][c];
                I[M][c] += m * I[i][c];
            }
        }
    }

    return I;
}





