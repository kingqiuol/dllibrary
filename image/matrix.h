//
// Created by jinqiu on 2018/10/19.
//

#ifndef MATRIX_H
#define MATRIX_H

typedef struct {
    int rows, cols;
    double **data;
    int shallow;
} Matrix;

typedef struct LUP {
    Matrix *L;
    Matrix *U;
    int *P;
    int n;
} LUP;

Matrix make_identity_homography();

Matrix make_translation_homography(float dx, float dy);

void free_matrix(Matrix m);

double mag_matrix(Matrix m);

Matrix make_matrix(int rows, int cols);

Matrix copy_matrix(Matrix m);

double *sle_solve(Matrix A, double *b);

Matrix matrix_mult_matrix(Matrix a, Matrix b);

double *matrix_mult_vector(Matrix m, double *v);

Matrix matrix_elmult_matrix(Matrix a, Matrix b);

Matrix matrix_sub_matrix(Matrix a, Matrix b);

Matrix random_matrix(int rows, int cols, double s);

void print_matrix(Matrix m);

double **n_principal_components(Matrix m, int n);

void test_matrix();

Matrix solve_system(Matrix M, Matrix b);

Matrix matrix_invert(Matrix m);

Matrix transpose_matrix(Matrix m);

Matrix axpy_matrix(double a, Matrix x, Matrix y);

#endif //MATRIX_H
