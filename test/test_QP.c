#include <stdio.h>
#include "matrix.h"
#include "QPIPM.h"

// int main(int argc, char **argv){
int main(){
    int n = 200;
    double gamma = 0.5;

    // cost
    Matrix* Q_ = create_random_matrix(n,n,-1,1);
    Matrix* Q_trans = transpose(Q_);
    Matrix* Q = multiply_matrices(Q_,Q_trans);
    // Matrix* Q = create_matrix(n,n);
    Matrix* Imat = create_identity_matrix(n);
    add_matrices_with_multiplier_r(Q,Imat,n);
    multiply_with_scalar_r(Q,gamma);
    Matrix* q = create_random_matrix(n,1,-1,1);
    // Matrix* q = create_matrix(n,1);

    // eq
    Matrix* A = create_matrix(1,n);
    set_all(A,1.0);
    Matrix* b = create_matrix(1,1);
    set_all(b,1.0);

    // ineq
    Matrix* G = create_identity_matrix(n);
    multiply_with_scalar_r(G,-1.0);
    Matrix* h = create_matrix(n,1);
    set_all(h,0.0);

    // Solve
    IPM* ipm = create_ipm(Q,q,A,b,G,h);
    solve_ipm(ipm,true);

    // Check
    printf("Converged?: %s\n",ipm->converged == true ? "true" : "false");
    // if (ipm->converged){
    //     printf("Optimal solution: \n");
    //     print_solution(ipm->solution);
    // }

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