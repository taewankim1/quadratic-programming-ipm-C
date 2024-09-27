#ifndef MATRIX_H
#define MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef int8_t s8;
typedef int16_t s16;
typedef int32_t s32;
typedef int64_t s64;
typedef double f64;

// typedef int32_t b32;


typedef struct {
    u32 rows;
    u32 cols;
    u8 is_square;
    double* data;
} Matrix;

typedef struct LUP_S {
    Matrix* L;
    Matrix* U;
    Matrix* P;
    u32 num_permutations;
} LUP;

typedef struct {
    Matrix* Q;
    Matrix* R;
} QR;

Matrix* create_matrix(const u32 rows, const u32 cols);
void free_matrix(Matrix* mat);

Matrix* create_random_matrix(unsigned int num_rows, unsigned int num_cols, double min, double max);
Matrix* create_identity_matrix(unsigned int size);
Matrix* create_matrix_fromfilef(FILE* f);
Matrix* copy_matrix(const Matrix* a);
int copy_matrix_inplace(const Matrix* a, const Matrix* r);
Matrix* transpose(const Matrix* a);

bool check_equality(Matrix* m1, Matrix* m2, double tolerance);

void print_matrix(const Matrix* mat);

int add_matrices_r(const Matrix* a, const Matrix* b);
int add_matrices_with_multiplier_r(const Matrix* a, const Matrix* b, f64 c);
int add_matrices_inplace(const Matrix* a, const Matrix* b, Matrix* r);
Matrix* add_matrices(const Matrix* a, const Matrix* b);
bool subtract_matrices_r(const Matrix* a, const Matrix* b);
Matrix* subtract_matrices(const Matrix* a, const Matrix* b);

double get_value(const Matrix* mat, const unsigned int row, const unsigned int col);
void set_value(Matrix* mat, const unsigned int row, const unsigned int col, double value);

Matrix* multiply_matrices(const Matrix* a, const Matrix* b);
Matrix* multiply_matrices0(const Matrix* a, const Matrix* b);
Matrix* multiply_matrices1(const Matrix* a, const Matrix* b);
Matrix* multiply_matrices2(const Matrix* a, const Matrix* b);
Matrix* multiply_matrices3(const Matrix* a, const Matrix* b);

int multiply_matrices_inplace(const Matrix* a, const Matrix* b, Matrix* r);

double f_norm(const Matrix* mat); // frobenius norm

Matrix* get_col(Matrix* mat, u32 col);
Matrix* get_row(Matrix* mat, u32 row);

void set_all(Matrix* mat, double value);
void set_diag(Matrix* mat, double value);

bool multiply_row_with_scalar_r(Matrix* a, u32 row, double num);
Matrix* multiply_row_with_scalar(const Matrix* a, u32 row, double num);
bool multiply_col_with_scalar_r(Matrix* a, u32 col, double num);
Matrix* multiply_col_with_scalar(const Matrix* a, u32 col, double num);
bool multiply_with_scalar_r(Matrix* a, double num);
Matrix* multiply_with_scalar(const Matrix* a, double num);

bool add_two_rows_r(Matrix* mat, u32 where, u32 row, double multiplier);
Matrix* add_two_rows(Matrix* mat, u32 where, u32 row, double multiplier);

Matrix* remove_col(Matrix* mat, u32 col);
Matrix* remove_row(Matrix* mat, u32 row);

bool swap_col_r(Matrix* mat, u32 col1, u32 col2);
Matrix* swap_col(Matrix* mat, u32 col1, u32 col2);
bool swap_row_r(Matrix* mat, u32 row1, u32 row2);
Matrix* swap_row(Matrix* mat, u32 row1, u32 row2);

Matrix* vcat_matrices(u32 mnum, Matrix** marr);
int vcat_matrices_inplace(u32 mnum, Matrix** marr, Matrix* r);
Matrix* hcat_matrices(u32 mnum, Matrix** marr);

int find_pivotidx(Matrix* mat, u32 col, u32 row);
Matrix* row_echelon_form(Matrix* mat);
int find_pivotmaxidx(Matrix* mat, u32 col, u32 row);
Matrix* reduced_row_echelon_form(Matrix* mat);

LUP* create_LUP(Matrix* L, Matrix* U, Matrix* P, u32 num_permutations);
void free_LUP(LUP* lu);

int find_absmaxidx(Matrix* mat, u32 k);
LUP* solve_lup(Matrix* mat);

Matrix* solve_ls_forward(Matrix* L, Matrix* b);
Matrix* solve_ls_backward(Matrix* U, Matrix* b);
Matrix* solve_ls_from_lup(LUP* lu, Matrix* b);
Matrix* solve_ls(Matrix* A, Matrix* b);

int solve_ls_backward_inplace(Matrix* U, Matrix* b, Matrix* x);
int solve_ls_from_lup_inplace(LUP* lu, Matrix* b, Matrix* x);
int solve_ls_inplace(Matrix* A, Matrix* b, Matrix* x);

Matrix* inverse(Matrix* A);
double compute_det(Matrix* A);

double inner_product(Matrix* a, Matrix* b);
double l2_vec_norm(const Matrix* mat); // l2 vector norm

QR* create_QR(Matrix* Q, Matrix* R);
void free_QR(QR* qr);
bool normalize_each_col(Matrix* mat);
QR* solve_qr(Matrix* mat);

Matrix* elementwise_product(Matrix* a, Matrix* b);
Matrix* elementwise_quotient(Matrix* a, Matrix* b);
int elementwise_quotient_r(Matrix* a, Matrix* b);
int elementwise_quotient_inplace(Matrix* a, Matrix* b, Matrix* r);


#ifdef __cplusplus
}
#endif

#endif