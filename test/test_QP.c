#include <stdio.h>
#include "matrix.h"
#include "QPIPM.h"

// int main(int argc, char **argv){
int main(){
    int n = 3;
    double gamma = 0.5;
    Matrix* A = create_matrix(1,n);
    Matrix* G = create_identity_matrix(n);
    Matrix* q = create_matrix(n,1);
    Matrix* b = create_matrix(1,1);
    Matrix* h = create_matrix(n,1);
    Matrix* Q_ = create_identity_matrix(n);
    Matrix* Imat = create_identity_matrix(n);

    multiply_with_scalar_r(Q_,2.0);
    set_value(Q_,0,0,3.0);
    set_value(Q_,1,1,5.0);
    Matrix* Q_trans = transpose(Q_);
    Matrix* Q = multiply_matrices(Q_,Q_trans);
    multiply_with_scalar_r(Imat,n);
    add_matrices_r(Q,Imat);
    multiply_with_scalar_r(Q,gamma);

    set_all(q,1.0);
    multiply_with_scalar_r(G,-1.0);
    set_all(h,0.0);
    set_all(A,1.0);
    set_all(b,1.0);

    IPM* ipm = create_ipm(Q,q,A,b,G,h);
    SOLUTION* sol = solve_ipm(ipm,true);

    free_ipm(ipm);
    free_matrix(A);
    free_matrix(Q);
    free_matrix(Q_);
    free_matrix(Q_trans);
    free_matrix(Imat);
    free_matrix(G);
    free_matrix(q);
    free_matrix(b);
    free_matrix(h);

    return 0;
}