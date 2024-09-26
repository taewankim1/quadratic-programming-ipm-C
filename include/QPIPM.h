#pragma once
#include "matrix.h"

typedef struct {
    double constraint;
    double gap;
    double cost;
    double iter_ref;
} TOLERANCE;

typedef struct {
    Matrix* x;
    Matrix* s;
    Matrix* z;
    Matrix* y;
} SOLUTION;

SOLUTION* create_solution(u32 nx, u32 ns, u32 nz, u32 ny);
void free_solution(SOLUTION* sol);
void print_solution(SOLUTION* sol);

typedef struct {
    SOLUTION* sol_a;
    SOLUTION* sol_c;
    SOLUTION* sol;
} DELTA;

DELTA* create_delta(u32 nx, u32 ns, u32 nz, u32 ny);
void free_delta(DELTA* del);

typedef struct {
    Matrix* Q;
    Matrix* q;
    Matrix* A;
    Matrix* b;
    Matrix* G;
    Matrix* h;

    u32 nx;
    u32 ns;
    u32 nz;
    u32 ny;
    u32 N;

    f64 eps;

    Matrix* reg;
    Matrix* KKT;
    Matrix* KKT_reg;
    Matrix* rhs_a;
    Matrix* rhs_c;
    Matrix* p_a;
    Matrix* p_c;

    bool converged;
    bool max_iter_reached;

    TOLERANCE tol;
    DELTA* delta;
    SOLUTION* solution;

} IPM;

IPM* create_ipm(Matrix* Q, Matrix* q, Matrix* A, Matrix* b, Matrix* G, Matrix* h);
void free_ipm(IPM* ipm);
bool initialize(IPM* ipm);
void update_kkt(IPM* ipm);
void regularize_kkt(IPM* ipm);
void iterative_refinement(IPM* ipm, Matrix* sol, Matrix* rhs, LUP* lu);
void index_sol_a(IPM* ipm);
f64 linesearch(Matrix* x, Matrix* dx);
void centering_params(IPM* ipm, f64* p_sig, f64* p_mu);


SOLUTION* solve_ipm(IPM* ipm, bool verbose); // NO MALLOC