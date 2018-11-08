#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "matrix.h"

Matrix make_identity_homography() {
    Matrix H = make_matrix(3, 3);
    H.data[0][0] = 1;
    H.data[1][1] = 1;
    H.data[2][2] = 1;
    return H;
}

Matrix make_translation_homography(float dx, float dy) {
    Matrix H = make_identity_homography();
    H.data[0][2] = dx;
    H.data[1][2] = dy;
    return H;
}

void free_matrix(Matrix m) {
    int i;
    for (i = 0; i < m.rows; ++i) free(m.data[i]);
    free(m.data);
}

Matrix make_matrix(int rows, int cols) {
    Matrix m;
    m.rows = rows;
    m.cols = cols;
    m.data = calloc(m.rows, sizeof(double *));
    int i;
    for (i = 0; i < m.rows; ++i) m.data[i] = calloc(m.cols, sizeof(double));
    return m;
}

Matrix copy_matrix(Matrix m) {
    int i, j;
    Matrix c = make_matrix(m.rows, m.cols);
    for (i = 0; i < m.rows; ++i) {
        for (j = 0; j < m.cols; ++j) {
            c.data[i][j] = m.data[i][j];
        }
    }
    return c;
}

Matrix augment_matrix(Matrix m) {
    int i, j;
    Matrix c = make_matrix(m.rows, m.cols * 2);
    for (i = 0; i < m.rows; ++i) {
        for (j = 0; j < m.cols; ++j) {
            c.data[i][j] = m.data[i][j];
        }
    }
    for (j = 0; j < m.rows; ++j) {
        c.data[j][j + m.cols] = 1;
    }
    return c;
}

Matrix make_identity(int rows, int cols) {
    int i;
    Matrix m = make_matrix(rows, cols);
    for (i = 0; i < rows && i < cols; ++i) {
        m.data[i][i] = 1;
    }
    return m;
}

Matrix matrix_mult_matrix(Matrix a, Matrix b) {
    assert(a.cols == b.rows);
    int i, j, k;
    Matrix p = make_matrix(a.rows, b.cols);
    for (i = 0; i < p.rows; ++i) {
        for (j = 0; j < p.cols; ++j) {
            for (k = 0; k < a.cols; ++k) {
                p.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
    return p;
}

Matrix matrix_elmult_matrix(Matrix a, Matrix b)
{
    assert(a.cols == b.cols);
    assert(a.rows == b.rows);
    int i, j;
    Matrix p = make_matrix(a.rows, a.cols);
    for(i = 0; i < p.rows; ++i){
        for(j = 0; j < p.cols; ++j){
            p.data[i][j] = a.data[i][j] * b.data[i][j];
        }
    }
    return p;
}

Matrix matrix_sub_matrix(Matrix a, Matrix b) {
    assert(a.cols == b.cols);
    assert(a.rows == b.rows);
    int i, j;
    Matrix p = make_matrix(a.rows, a.cols);
    for (i = 0; i < p.rows; ++i) {
        for (j = 0; j < p.cols; ++j) {
            p.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
    return p;
}

Matrix axpy_matrix(double a, Matrix x, Matrix y)
{
    assert(x.cols == y.cols);
    assert(x.rows == y.rows);
    int i, j;
    Matrix p = make_matrix(x.rows, x.cols);
    for(i = 0; i < x.rows; ++i){
        for(j = 0; j < x.cols; ++j){
            p.data[i][j] = a*x.data[i][j] + y.data[i][j];
        }
    }
    return p;
}

Matrix transpose_matrix(Matrix m) {
    Matrix t;
    t.rows = m.cols;
    t.cols = m.rows;
    t.data = calloc(t.rows, sizeof(double *));
    int i, j;
    for (i = 0; i < t.rows; ++i) {
        t.data[i] = calloc(t.cols, sizeof(double));
        for (j = 0; j < t.cols; ++j) {
            t.data[i][j] = m.data[j][i];
        }
    }
    return t;
}

void scale_matrix(Matrix m, double s) {
    int i, j;
    for (i = 0; i < m.rows; ++i) {
        for (j = 0; j < m.cols; ++j) {
            m.data[i][j] *= s;
        }
    }
}

double *matrix_mult_vector(Matrix m, double *v) {
    double *p = calloc(m.rows, sizeof(double));
    int i, j;
    for (i = 0; i < m.rows; ++i) {
        for (j = 0; j < m.cols; ++j) {
            p[i] += m.data[i][j] * v[j];
        }
    }
    return p;
}

void print_matrix(Matrix m) {
    int i, j;
    printf(" __");
    for (j = 0; j < 16 * m.cols - 1; ++j) printf(" ");
    printf("__ \n");

    printf("|  ");
    for (j = 0; j < 16 * m.cols - 1; ++j) printf(" ");
    printf("  |\n");

    for (i = 0; i < m.rows; ++i) {
        printf("|  ");
        for (j = 0; j < m.cols; ++j) {
            printf("%15.7f ", m.data[i][j]);
        }
        printf(" |\n");
    }
    printf("|__");
    for (j = 0; j < 16 * m.cols - 1; ++j) printf(" ");
    printf("__|\n");
}

double *LUP_solve(Matrix L, Matrix U, int *p, double *b) {
    int i, j;
    double *c = calloc(L.rows, sizeof(double));
    for (i = 0; i < L.rows; ++i) {
        int pi = p[i];
        c[i] = b[pi];
        for (j = 0; j < i; ++j) {
            c[i] -= L.data[i][j] * c[j];
        }
    }
    for (i = U.rows - 1; i >= 0; --i) {
        for (j = i + 1; j < U.cols; ++j) {
            c[i] -= U.data[i][j] * c[j];
        }
        c[i] /= U.data[i][i];
    }
    return c;
}

Matrix matrix_invert(Matrix m) {
    //print_matrix(m);
    Matrix none = {0};
    if (m.rows != m.cols) {
        fprintf(stderr, "Matrix not square\n");
        return none;
    }
    Matrix c = augment_matrix(m);
    //print_matrix(c);

    int i, j, k;
    for (k = 0; k < c.rows; ++k) {
        double p = 0.;
        int index = -1;
        for (i = k; i < c.rows; ++i) {
            double val = fabs(c.data[i][k]);
            if (val > p) {
                p = val;
                index = i;
            }
        }
        if (index == -1) {
            fprintf(stderr, "Can't do it, sorry!\n");
            free_matrix(c);
            return none;
        }

        double *swap = c.data[index];
        c.data[index] = c.data[k];
        c.data[k] = swap;

        double val = c.data[k][k];
        c.data[k][k] = 1;
        for (j = k + 1; j < c.cols; ++j) {
            c.data[k][j] /= val;
        }
        for (i = k + 1; i < c.rows; ++i) {
            float s = -c.data[i][k];
            c.data[i][k] = 0;
            for (j = k + 1; j < c.cols; ++j) {
                c.data[i][j] += s * c.data[k][j];
            }
        }
    }
    for (k = c.rows - 1; k > 0; --k) {
        for (i = 0; i < k; ++i) {
            double s = -c.data[i][k];
            c.data[i][k] = 0;
            for (j = k + 1; j < c.cols; ++j) {
                c.data[i][j] += s * c.data[k][j];
            }
        }
    }
    //print_matrix(c);
    Matrix inv = make_matrix(m.rows, m.cols);
    for (i = 0; i < m.rows; ++i) {
        for (j = 0; j < m.cols; ++j) {
            inv.data[i][j] = c.data[i][j + m.cols];
        }
    }
    free_matrix(c);
    //print_matrix(inv);
    return inv;
}

int *in_place_LUP(Matrix m) {
    int *pivot = calloc(m.rows, sizeof(int));
    if (m.rows != m.cols) {
        fprintf(stderr, "Matrix not square\n");
        return 0;
    }

    int i, j, k;
    for (k = 0; k < m.rows; ++k) pivot[k] = k;
    for (k = 0; k < m.rows; ++k) {
        double p = 0.;
        int index = -1;
        for (i = k; i < m.rows; ++i) {
            double val = fabs(m.data[i][k]);
            if (val > p) {
                p = val;
                index = i;
            }
        }
        if (index == -1) {
            fprintf(stderr, "Matrix is singular\n");
            return 0;
        }

        int swapi = pivot[k];
        pivot[k] = pivot[index];
        pivot[index] = swapi;

        double *swap = m.data[index];
        m.data[index] = m.data[k];
        m.data[k] = swap;

        for (i = k + 1; i < m.rows; ++i) {
            m.data[i][k] = m.data[i][k] / m.data[k][k];
            for (j = k + 1; j < m.cols; ++j) {
                m.data[i][j] -= m.data[i][k] * m.data[k][j];
            }
        }
    }
    return pivot;
}

Matrix random_matrix(int rows, int cols, double s)
{
    Matrix m = make_matrix(rows, cols);
    int i, j;
    for(i = 0; i < rows; ++i){
        for(j = 0; j < cols; ++j){
            m.data[i][j] = 2*s*(rand()%1000/1000.0) - s;
        }
    }
    return m;
}

double mag_matrix(Matrix m)
{
    int i, j;
    double sum = 0;
    for(i = 0; i < m.rows; ++i){
        for(j = 0; j < m.cols; ++j){
            sum += m.data[i][j]*m.data[i][j];
        }
    }
    return sqrt(sum);
}

double *sle_solve(Matrix A, double *b) {
    int *p = in_place_LUP(A);
    return LUP_solve(A, A, p, b);
}

Matrix solve_system(Matrix M, Matrix b) {
    Matrix none = {0};
    Matrix Mt = transpose_matrix(M);
    Matrix MtM = matrix_mult_matrix(Mt, M);
    Matrix MtMinv = matrix_invert(MtM);
    if (!MtMinv.data) return none;
    Matrix Mdag = matrix_mult_matrix(MtMinv, Mt);
    Matrix a = matrix_mult_matrix(Mdag, b);
    free_matrix(Mt);
    free_matrix(MtM);
    free_matrix(MtMinv);
    free_matrix(Mdag);
    return a;
}

void test_matrix() {
    int i;
    for (i = 0; i < 100; ++i) {
        int s = rand() % 4 + 3;
        Matrix m = random_matrix(s, s,10);
        Matrix inv = matrix_invert(m);
        Matrix res = matrix_mult_matrix(m, inv);
        print_matrix(res);
    }
}



