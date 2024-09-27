#include "matrix.h"
#include "QPIPM.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

SOLUTION* create_solution(u32 nx, u32 ns, u32 nz, u32 ny){
    SOLUTION* sol = (SOLUTION*) malloc(sizeof(SOLUTION));
    if (sol == NULL){
        fprintf(stderr,"SOLUTION_MEMORY_ALLOCATION_FAILED\n");
        exit(EXIT_FAILURE);
    }
    sol->x = create_matrix(nx,1);
    sol->s = create_matrix(ns,1);
    sol->z = create_matrix(nz,1);
    sol->y = create_matrix(ny,1);
    return sol;
}

void free_solution(SOLUTION* sol){
    free_matrix(sol->x);
    free_matrix(sol->s);
    free_matrix(sol->z);
    free_matrix(sol->y);
    free(sol);
}

void print_solution(SOLUTION* sol){
    printf("--x--\n");
    print_matrix(sol->x);
    printf("--s--\n");
    print_matrix(sol->s);
    printf("--z--\n");
    print_matrix(sol->z);
    printf("--y--\n");
    print_matrix(sol->y);
}

DELTA* create_delta(u32 nx, u32 ns, u32 nz, u32 ny){
    DELTA* del = (DELTA*) malloc(sizeof(DELTA));
    if (del == NULL){
        fprintf(stderr,"DELTA_MEMORY_ALLOCATION_FAILED\n");
        exit(EXIT_FAILURE);
    }
    del->sol = create_solution(nx,ns,nz,ny);
    del->sol_a = create_solution(nx,ns,nz,ny);
    del->sol_c = create_solution(nx,ns,nz,ny);
    return del;
}

void free_delta(DELTA* del){
    free_solution(del->sol);
    free_solution(del->sol_a);
    free_solution(del->sol_c);
    free(del);
}

IPM* create_ipm(Matrix* Q, Matrix* q, Matrix* A, Matrix* b, Matrix* G, Matrix* h){
    IPM* ipm = (IPM*) malloc(sizeof(IPM));
    if (ipm == NULL){
        fprintf(stderr,"IPM_MEMORY_ALLOCATION_FAILED\n");
        exit(EXIT_FAILURE);
    }
    // Variable lengths
    u32 nx = q->rows; // n
    u32 ns = h->rows; // p
    u32 nz = h->rows; // p
    u32 ny = b->rows; // m
    
    ipm->nx = nx;
    ipm->ns = ns;
    ipm->nz = nz;
    ipm->ny = ny;

    u32 N = nx+ns+nz+ny;
    ipm->N = N;

    ipm->Q = Q;
    ipm->q = q;
    ipm->A = A;
    ipm->b = b;
    ipm->G = G;
    ipm->h = h;

    f64 eps = 1e-7;
    ipm->eps = eps;

    Matrix* reg = create_matrix(N,N);
    for(size_t i=0;i<nx+ns;++i){
        set_value(reg,i,i,eps);
    }
    for(size_t i=nx+ns;i<N;++i){
        set_value(reg,i,i,-eps);
    }
    ipm->reg = reg;
    ipm->KKT = create_matrix(N,N);
    ipm->KKT_reg = create_matrix(N,N);
    ipm->rhs_a = create_matrix(N,1);
    ipm->rhs_c = create_matrix(N,1);
    ipm->p_a = create_matrix(N,1);
    ipm->p_c = create_matrix(N,1);

    ipm->converged = false;
    ipm->max_iter_reached = false;
    
    ipm->tol.constraint = 1e-6;
    ipm->tol.gap = 1e-6;
    ipm->tol.cost= 1e-4;
    // ipm->tol.iter_ref = 1e-6;
    ipm->tol.iter_ref = 1e-4;

    ipm->delta = create_delta(nx,ns,nz,ny);
    ipm->solution = create_solution(nx,ns,nz,ny);
    return ipm;
}

void free_ipm(IPM* ipm){
    free_matrix(ipm->reg);
    free_matrix(ipm->KKT);
    free_matrix(ipm->KKT_reg);
    free_matrix(ipm->rhs_a);
    free_matrix(ipm->rhs_c);
    free_matrix(ipm->p_a);
    free_matrix(ipm->p_c);
    free_delta(ipm->delta);
    free_solution(ipm->solution);
    free(ipm);
}

bool initialize(IPM* ipm){
    u32 nx = ipm->nx;
    u32 ns = ipm->ns;
    u32 nz = ipm->nz;
    u32 ny = ipm->ny;
    u32 Ni = nx + nz + ny;

    Matrix* A_init = create_matrix(Ni,Ni);
    for(size_t i = 0;i < nx;++i){
        memcpy(&A_init->data[Ni*i],&ipm->Q->data[nx*i],sizeof(f64)*nx);
        // TODO - unroll?
    }
    if (ipm->G != NULL){
        for(size_t i = 0;i<nx;++i){
            for(size_t j = nx;j<nx+ns;++j){
                set_value(A_init,i,j,get_value(ipm->G,j-nx,i));
            }
        }
        for(size_t i = nx;i<nx+ns;++i){
            memcpy(&A_init->data[Ni*i],&ipm->G->data[nx*(i-nx)],sizeof(f64)*nx);
        }
        for(size_t i = nx;i<nx+ns;++i){
            for(size_t j = nx;j<nx+ns;++j){
                if (i == j) set_value(A_init,i,j,-1.0);
                else set_value(A_init,i,j,0.0);
            }
        }
    }
    if (ipm->A != NULL){
        for(size_t i = 0;i<nx;++i){
            for(size_t j = nx+ns;j<nx+ns+ny;++j){
                set_value(A_init,i,j,get_value(ipm->A,j-nx-ns,i));
            }
        }
        for(size_t i = nx+ns;i<nx+ns+ny;++i){
            memcpy(&A_init->data[Ni*i],&ipm->A->data[nx*(i-nx-ns)],sizeof(f64)*nx);
        }
    }

    Matrix* negative_q = multiply_with_scalar(ipm->q,-1.0);
    Matrix* marr[3] = {negative_q,ipm->h,ipm->b};
    Matrix* rhs = vcat_matrices(3, marr);

    Matrix* sol = solve_ls(A_init,rhs);

    Matrix* px = ipm->solution->x;
    Matrix* ps = ipm->solution->s;
    Matrix* pz = ipm->solution->z;
    Matrix* py = ipm->solution->y;

    memcpy(px->data,sol->data,sizeof(f64)*nx);
    memcpy(pz->data,&sol->data[nx],sizeof(f64)*nz);
    memcpy(py->data,&sol->data[nx+ns],sizeof(f64)*ny);

    if (ipm->G != NULL){
        Matrix* neg_z = copy_matrix(pz);
        multiply_with_scalar_r(neg_z,-1.0);
        double a_p = neg_z->data[0];
        for (size_t i = 1;i<nz;++i){
            if (neg_z->data[i] < a_p){
                a_p = neg_z->data[i];
            }
        }
        a_p *= -1;
        memcpy(ps->data,neg_z->data,sizeof(f64)*ns);
        if (a_p >= 0){
            for(size_t i = 0;i<ns;++i){
                ps->data[i] += 1 + a_p;
            }
        }
        double a_d = pz->data[0];
        for (size_t i = 1;i<nz;++i){
            if (pz->data[i] < a_d){
                a_d = pz->data[i];
            }
        }
        a_d *= -1;
        if (a_d >= 0){
            for(size_t i = 0;i<nz;++i){
                pz->data[i] += 1 + a_d;
            }
        }
        // free
        free_matrix(neg_z);
    }
    free_matrix(A_init);
    free_matrix(negative_q);
    free_matrix(rhs);
    free_matrix(sol);

    // initialize_kkt
    u32 N = nx+ns+nz+ny;
    Matrix* KKT = ipm->KKT;
    for(size_t i = 0;i < nx;++i){
        memcpy(&KKT->data[N*i],&ipm->Q->data[nx*i],sizeof(f64)*nx);
    }
    if (ipm->G != NULL){
        for(size_t i = 0;i<nx;++i){
            for(size_t j = nx+ns;j<nx+ns+nz;++j){
                set_value(KKT,i,j,get_value(ipm->G,j-nx-ns,i));
            }
        }
        for(size_t i = nx+ns;i<nx+ns+nz;++i){
            memcpy(&KKT->data[N*i],&ipm->G->data[nx*(i-nx-ns)],sizeof(f64)*nx);
        }
        for(size_t i = nx+ns;i<nx+ns+nz;++i){
            for(size_t j = nx;j<nx+ns;++j){
                if (i - nx - ns == j - nx) set_value(KKT,i,j,1.0);
                else set_value(KKT,i,j,0.0);
            }
        }
        for(size_t i = nx;i<nx+ns;++i){
            for(size_t j = nx+ns;j<nx+ns+nz;++j){
                if (i - nx == j - nx - ns) set_value(KKT,i,j,1.0);
                else set_value(KKT,i,j,0.0);
            }
        }
    }

    if (ipm->A != NULL){
        for(size_t i = 0;i<nx;++i){
            for(size_t j = nx+ns+nz;j<nx+ns+nz+ny;++j){
                set_value(KKT,i,j,get_value(ipm->A,j-nx-ns-nz,i));
            }
        }
        for(size_t i = nx+ns+nz;i<nx+ns+nz+ny;++i){
            memcpy(&KKT->data[N*i],&ipm->A->data[nx*(i-nx-ns-nz)],sizeof(f64)*nx);
        }
    }
    // print_matrix(KKT);
    return true;
}

void update_kkt(IPM* ipm){
    u32 nx = ipm->nx;
    u32 ns = ipm->ns;

    Matrix* KKT = ipm->KKT;
    Matrix* pz = ipm->solution->z;
    Matrix* ps = ipm->solution->s;
    if (ipm->G != NULL){
        for(size_t i = nx;i<nx+ns;++i){
            for(size_t j = nx;j<nx+ns;++j){
                if (i==j){
                    double tmp = pz->data[i-nx] / ps->data[i-nx];
                    set_value(KKT,i,j,tmp);
                }
                else{
                    set_value(KKT,i,j,0.0);
                }
            }
        }
    }
}

void regularize_kkt(IPM* ipm){
    u32 nx = ipm->nx;
    u32 ns = ipm->ns;
    u32 nz = ipm->nz;
    u32 ny = ipm->ny;
    u32 N = nx+ns+nz+ny;
    Matrix* KKT = ipm->KKT;
    Matrix* KKT_reg = ipm->KKT_reg;
    Matrix* reg = ipm->reg;
    memcpy(KKT_reg->data,KKT->data,sizeof(f64)*N*N);
    add_matrices_r(KKT_reg,reg);
    // print_matrix(KKT_reg);
}

void rhs_kkt_a(IPM* ipm){
    Matrix* rhs_a1 = multiply_matrices(ipm->Q,ipm->solution->x);
    add_matrices_r(rhs_a1,ipm->q);
    if (ipm->A == NULL && ipm->G == NULL){
        multiply_with_scalar_r(rhs_a1,-1.0);
        copy_matrix_inplace(rhs_a1,ipm->rhs_a);
    }
    else if (ipm->A == NULL && ipm->G != NULL){
        Matrix* G_transpose = transpose(ipm->G);
        Matrix* Gz = multiply_matrices(G_transpose,ipm->solution->z);
        add_matrices_r(rhs_a1,Gz);
        multiply_with_scalar_r(rhs_a1,-1.0);

        Matrix* rhs_a2 = copy_matrix(ipm->solution->z);
        multiply_with_scalar_r(rhs_a2,-1.0);

        Matrix* rhs_a3 = multiply_matrices(ipm->G,ipm->solution->x);
        add_matrices_r(rhs_a3,ipm->solution->s);
        subtract_matrices_r(rhs_a3,ipm->h);
        multiply_with_scalar_r(rhs_a3,-1.0);

        Matrix* rhs_a[3] = {rhs_a1,rhs_a2,rhs_a3};
        vcat_matrices_inplace(3,rhs_a,ipm->rhs_a);

        free_matrix(rhs_a1);
        free_matrix(rhs_a2);
        free_matrix(rhs_a3);
        // free_matrix(rhs_a);
        free_matrix(G_transpose);
        free_matrix(Gz);
    }
    else if (ipm->A != NULL && ipm->G == NULL){
        Matrix* A_transpose = transpose(ipm->A);
        Matrix* Ay = multiply_matrices(A_transpose,ipm->solution->y);
        add_matrices_r(rhs_a1,Ay);
        multiply_with_scalar_r(rhs_a1,-1.0);

        Matrix* rhs_a2 = multiply_matrices(ipm->A,ipm->solution->x);
        subtract_matrices_r(rhs_a2,ipm->b);
        multiply_with_scalar_r(rhs_a2,-1.0);

        Matrix* rhs_a[2] = {rhs_a1,rhs_a2};
        vcat_matrices_inplace(2,rhs_a,ipm->rhs_a);

        free_matrix(rhs_a1);
        free_matrix(rhs_a2);
        free_matrix(A_transpose);
        free_matrix(Ay);
    }
    else{
        Matrix* A_transpose = transpose(ipm->A);
        Matrix* Ay = multiply_matrices(A_transpose,ipm->solution->y);
        Matrix* G_transpose = transpose(ipm->G);
        Matrix* Gz = multiply_matrices(G_transpose,ipm->solution->z);
        add_matrices_r(rhs_a1,Ay);
        add_matrices_r(rhs_a1,Gz);
        multiply_with_scalar_r(rhs_a1,-1.0);

        Matrix* rhs_a2 = copy_matrix(ipm->solution->z);
        multiply_with_scalar_r(rhs_a2,-1.0);

        Matrix* rhs_a3 = multiply_matrices(ipm->G,ipm->solution->x);
        add_matrices_r(rhs_a3,ipm->solution->s);
        subtract_matrices_r(rhs_a3,ipm->h);
        multiply_with_scalar_r(rhs_a3,-1.0);

        Matrix* rhs_a4 = multiply_matrices(ipm->A,ipm->solution->x);
        subtract_matrices_r(rhs_a4,ipm->b);
        multiply_with_scalar_r(rhs_a4,-1.0);

        Matrix* rhs_a[4] = {rhs_a1,rhs_a2,rhs_a3,rhs_a4};
        vcat_matrices_inplace(4,rhs_a,ipm->rhs_a);

        free_matrix(rhs_a1);
        free_matrix(rhs_a2);
        free_matrix(rhs_a3);
        free_matrix(rhs_a4);
        free_matrix(G_transpose);
        free_matrix(Gz);
        free_matrix(A_transpose);
        free_matrix(Ay);
    }
    // print_matrix(ipm->solution->x);
    // printf("--\n");
    // print_matrix(ipm->rhs_a);
}

void iterative_refinement(IPM* ipm, Matrix* sol, Matrix* rhs, LUP* lu){
    u32 N = ipm->rhs_a->rows;
    Matrix* dl = create_matrix(N,1) ;

    Matrix* check = multiply_matrices(ipm->KKT,sol);
    subtract_matrices_r(check,rhs);
    multiply_with_scalar_r(check,-1.0);
    while(l2_vec_norm(check) > ipm->tol.iter_ref){
        solve_ls_from_lup_inplace(lu,check,dl);
        add_matrices_r(sol,dl);
        multiply_matrices_inplace(ipm->KKT,sol,check);
        subtract_matrices_r(check,rhs);
        multiply_with_scalar_r(check,-1.0);
    }
    free_matrix(check);
    free_matrix(dl);
}

void index_sol_a(IPM* ipm){
    u32 nx = ipm->nx;
    u32 ns = ipm->ns;
    u32 nz = ipm->nz;
    u32 ny = ipm->ny;

    memcpy(ipm->delta->sol_a->x->data,ipm->p_a->data,sizeof(f64)*nx);
    memcpy(ipm->delta->sol_a->s->data,&ipm->p_a->data[nx],sizeof(f64)*ns);
    memcpy(ipm->delta->sol_a->z->data,&ipm->p_a->data[nx+ns],sizeof(f64)*nz);
    memcpy(ipm->delta->sol_a->y->data,&ipm->p_a->data[nx+ns+nz],sizeof(f64)*ny);
}

f64 linesearch(Matrix* x, Matrix* dx){
    if (x->rows != dx->rows){
        fprintf(stderr,"DIMENSION_MISMATCH\n");
        return 1.0;
    }
    f64 min_val = 1.0; // initialized  as 1
    u32 rows = x->rows;
    for (size_t i = 0;i<rows;++i){
        f64 temp;
        f64 dx_val = dx->data[i];
        if (dx_val < 0){
            temp = - x->data[i] / dx_val;
        }
        else{
            temp = 2.0; // some constant larger than 1
        }
        if (temp < min_val){
            min_val = temp;
        }
    }
    return min_val;
}

void centering_params(IPM* ipm, f64* p_sig, f64* p_mu){
    *p_sig = 0;
    *p_mu = 0;
    if (ipm->G != NULL){
        Matrix* s = ipm->solution->s;
        Matrix* z = ipm->solution->z;
        Matrix* s_a = ipm->delta->sol_a->s;
        Matrix* z_a = ipm->delta->sol_a->z;

        *p_mu = inner_product(s,z) / ipm->ns;
        f64 a1 = linesearch(s,s_a);
        f64 a2 = linesearch(z,z_a);
        f64 a = a1 < a2 ? a1 : a2;
        Matrix* temp1 = multiply_with_scalar(s_a,a);
        add_matrices_r(temp1,s);
        Matrix* temp2 = multiply_with_scalar(z_a,a);
        add_matrices_r(temp2,z);
        f64 sig = inner_product(temp1,temp2) / inner_product(s,z);
        sig = pow(sig,3);
        *p_sig = sig;
        free_matrix(temp1);
        free_matrix(temp2);
    }
}

void rhs_kkt_c(IPM* ipm, f64 sig, f64 mu){
    if (ipm->G != NULL){
        Matrix* s = ipm->solution->s;
        Matrix* s_a = ipm->delta->sol_a->s;
        Matrix* z_a = ipm->delta->sol_a->z;

        Matrix* temp = create_matrix(ipm->ns,1);
        set_all(temp,1.0);
        multiply_with_scalar_r(temp,sig*mu);
        Matrix* sz = elementwise_product(s_a,z_a);
        subtract_matrices_r(temp,sz);
        elementwise_quotient_r(temp,s);

        u32 nx = ipm->nx;
        u32 ns = ipm->ns;
        memcpy(&ipm->rhs_c->data[nx],temp->data,sizeof(f64)*ns);
        free_matrix(temp);
        free_matrix(sz);
    }
}

void index_sol_c(IPM* ipm){
    u32 nx = ipm->nx;
    u32 ns = ipm->ns;
    u32 nz = ipm->nz;
    u32 ny = ipm->ny;

    memcpy(ipm->delta->sol_c->x->data,ipm->p_c->data,sizeof(f64)*nx);
    memcpy(ipm->delta->sol_c->s->data,&ipm->p_c->data[nx],sizeof(f64)*ns);
    memcpy(ipm->delta->sol_c->z->data,&ipm->p_c->data[nx+ns],sizeof(f64)*nz);
    memcpy(ipm->delta->sol_c->y->data,&ipm->p_c->data[nx+ns+nz],sizeof(f64)*ny);
}

void combine_deltas(IPM* ipm){
    add_matrices_inplace(ipm->delta->sol_a->x,ipm->delta->sol_c->x,ipm->delta->sol->x);
    add_matrices_inplace(ipm->delta->sol_a->s,ipm->delta->sol_c->s,ipm->delta->sol->s);
    add_matrices_inplace(ipm->delta->sol_a->z,ipm->delta->sol_c->z,ipm->delta->sol->z);
    add_matrices_inplace(ipm->delta->sol_a->y,ipm->delta->sol_c->y,ipm->delta->sol->y);
}

f64 backtracking_linesearch(IPM* ipm){
    f64 a = 1;
    f64 b = 0.5;
    f64 t = 0.9;

    // aliasing
    Matrix* x = ipm->solution->x;
    Matrix* del_x = ipm->delta->sol->x;
    
    Matrix* xkp1 = copy_matrix(del_x);
    multiply_with_scalar_r(xkp1,a);
    add_matrices_r(xkp1,x);

    Matrix* Qxkp1 = multiply_matrices(ipm->Q,xkp1);
    f64 temp1 = 0.5 * inner_product(xkp1,Qxkp1);
    f64 temp2 = inner_product(ipm->q,xkp1);
    Matrix* Qx = multiply_matrices(ipm->Q,x);
    f64 temp3 = 0.5 * inner_product(Qx,x);
    f64 temp4 = inner_product(ipm->q,x);
    Matrix* temp5 = multiply_matrices(ipm->Q,x);
    add_matrices_r(temp5,ipm->q);
    f64 temp6 = a * t * inner_product(temp5,del_x);
    
    while(temp1 + temp2 > temp3 + temp4 + temp6){
        a *= b;
        copy_matrix_inplace(del_x,xkp1);
        multiply_with_scalar_r(xkp1,a);
        add_matrices_r(xkp1,x);

        multiply_matrices_inplace(ipm->Q,xkp1,Qxkp1);
        temp1 = 0.5 * inner_product(xkp1,Qxkp1);
        temp2 = inner_product(ipm->q,xkp1);
        temp6 = a * t * inner_product(temp5,del_x);
    }
    free_matrix(xkp1);
    free_matrix(Qxkp1);
    free_matrix(Qx);
    free_matrix(temp5);
    return a;
}
void update_vars(IPM* ipm,f64 a){
    add_matrices_with_multiplier_r(ipm->solution->x,ipm->delta->sol->x,a);
    add_matrices_with_multiplier_r(ipm->solution->s,ipm->delta->sol->s,a);
    add_matrices_with_multiplier_r(ipm->solution->z,ipm->delta->sol->z,a);
    add_matrices_with_multiplier_r(ipm->solution->y,ipm->delta->sol->y,a);
}
        // logging(ipm,iter,a,&Jcurr,&eq_res,&ineq_res,&gap);
void logging(IPM* ipm, size_t iter, f64 a,
    f64* pJcurr, f64* peq_res, f64* pineq_res, f64* pgap, bool verbose){
    // aliasing
    Matrix* x = ipm->solution->x;
    Matrix* s = ipm->solution->s;
    Matrix* z = ipm->solution->z;

    Matrix* Qx = multiply_matrices(ipm->Q,x);
    f64 temp1 = 0.5 * inner_product(Qx,x);
    f64 temp2 = inner_product(ipm->q,x);
    free_matrix(Qx);

    *pJcurr = temp1 + temp2;
    *peq_res = 0;
    *pineq_res = 0;
    *pgap = 0;
    if (ipm->A != NULL){
        Matrix* Ax_b = multiply_matrices(ipm->A,x);
        subtract_matrices_r(Ax_b,ipm->b);
        *peq_res = l2_vec_norm(Ax_b);
        free_matrix(Ax_b);
    }
    if (ipm->G != NULL){
        *pgap = inner_product(s,z);
        Matrix* Gxs_h = multiply_matrices(ipm->G,x);
        add_matrices_r(Gxs_h,s);
        subtract_matrices_r(Gxs_h,ipm->h);
        *pineq_res = l2_vec_norm(Gxs_h);
        free_matrix(Gxs_h);
    }
    if (verbose == true){
        printf("%3zu    %10.3e    %9.2e    %9.2e    %9.2e    %6.4f\n", iter, *pJcurr, *pgap, *peq_res, *pineq_res, a);
    }
}

int solve_ipm(IPM* ipm, bool verbose){
    if (verbose){
        printf("iter      objv          gap         |Ax-b|      |Gx+s-h|     step\n");
        printf("------------------------------------------------------------------\n");
    }
    f64 Jprev = 1e3;
    f64 Jcurr = 0;
    f64 ineq_res = 0;
    f64 eq_res = 0;
    f64 gap = 0;

    bool flag_initialize = initialize(ipm);
    if (!flag_initialize){
        fprintf(stderr,"INITIALIZE_FAILED\n");
        return -1;
        // exit(EXIT_FAILURE);
    }

    f64 a = 0;
    u32 max_iter = 10;
    for(size_t iter = 0; iter<max_iter;++iter){
        update_kkt(ipm);
        regularize_kkt(ipm);

        // LU decomposition for KKT_reg
        LUP* lu = solve_lup(ipm->KKT_reg);
        if (lu == NULL){
            fprintf(stderr,"LU_DECOMPOSE_KKT_reg_FAILED\n");
            return -1;
        }

        // Affine scaling step
        rhs_kkt_a(ipm);
        solve_ls_from_lup_inplace(lu,ipm->rhs_a,ipm->p_a);
        // print_matrix(ipm->p_a);
        iterative_refinement(ipm,ipm->p_a,ipm->rhs_a,lu);
        // printf("--\n");
        // print_matrix(ipm->p_a);
        index_sol_a(ipm);

        // Centering and correction step
        f64 sig,mu;
        centering_params(ipm,&sig,&mu);
        rhs_kkt_c(ipm,sig,mu);
        solve_ls_from_lup_inplace(lu,ipm->rhs_c,ipm->p_c);
        iterative_refinement(ipm,ipm->p_c,ipm->rhs_c,lu);
        index_sol_c(ipm);

        // free lu
        free_LUP(lu);

        // combine deltas
        combine_deltas(ipm);

        // Update solution iterate after linesearch
        if (ipm->G != NULL){
            f64 as = linesearch(ipm->solution->s,ipm->delta->sol->s);
            f64 az = linesearch(ipm->solution->z,ipm->delta->sol->z);
            a = as < az ? as : az;
            a *= 0.99;
            a = a < 1.0 ? a : 1.0;
        }
        else{
            a = backtracking_linesearch(ipm);
        }
        update_vars(ipm,a);
        logging(ipm,iter,a,&Jcurr,&eq_res,&ineq_res,&gap,verbose);

        ipm->converged = (fabs(Jcurr - Jprev) < ipm->tol.cost) && (eq_res < ipm->tol.constraint) && (ineq_res < ipm->tol.constraint) && (gap < ipm->tol.gap);
        ipm->max_iter_reached = (iter >= max_iter);
        Jprev = Jcurr;
        if (ipm->converged){
            break;
        }
    }
    return 1;
}