//
// Created by taratt on 2/7/23.
//

/**
 * This library contains functions to perform Matrix Multiplication and Triangular Parallel Solve
*/

#pragma once
#include <iostream>
#include <cmath>

#include "omp.h"
#include "inputFormatter.h"
#include "DependencyGraphCalulations.h"
#include "chrono"

using namespace std::chrono;

/**
 * Implements a naive approach for Sparse Triangular Solve
 *
 * @param n The dimension of the matrix
 * @param Lp Column pointers in the CSC format of the matrix
 * @param Li Row indices of the nonzero elements in the CSC format of the matrix
 * @param Lx Nonezero elements of the matrix in the CSC format of the matrix
 * @param x The right hand side vector of the system
 * @return Whether the operation was successful
 */

int lsolve (int n, size_t * Lp, int* Li, long double* Lx,long double *x){

    int p, j;
    if (!Lp || !Li || !x) {return (0) ;} // check inputs
    auto start = high_resolution_clock::now();

    //Performing the calculations on all columns
    for (j = 0; j<n; j++){
        x[j] /= Lx[Lp[j]];
        for (p = Lp[j]+1; p<Lp[j+1]; p++) {
            if (!(abs(Lx[p])<pow(10,-5) || abs(x[j])<pow(10,-5)))
                x[Li[p]] -= Lx[p] * x[j];
        }
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    //Printing the time taken to solve the system
    cout<<"Solve time:  ";
    cout << duration.count() <<"ms"<< endl;
    return (1);
}

/**
 * Implements an optimized single-threaded approach for Sparse Triangular Solve that optimizes
 * the naive algorithm by only using the columns that take part in the computations
 *
 * @param L CSC format of the matrix
 * @param x The right hand side vector of the system
 * @return Whether the operation was successful
 */

int optimized_single_lsolve (CSC L, Vector_b x){
    int p, i;
    if (!L.col_ptr || !L.row_ind || !x.val) {return (0) ;} // check inputs
    auto start1 = high_resolution_clock::now();

    //Getting the columns that will take part in the computations
    set <int> mattered_columns = get_reach(L, x);

    auto stop1 = high_resolution_clock::now();
    auto start2 = high_resolution_clock::now();

    set<int >::iterator it ;
    //Performing the computations only on the columns that take part in the computations of the naive approach
    for (it = mattered_columns.begin() ; it != mattered_columns.end() ; it++ ){
        int j = *it;
        x.val[j] /= L.val[L.col_ptr[j]];
        for (p = L.col_ptr[j]+1; p<L.col_ptr[j+1]; p++) {
            if (!(abs(L.val[p])<pow(10,-5) || abs(x.val[j])<pow(10,-5)))
                x.val[L.row_ind[p]] -= L.val[p] * x.val[j];
        }
    }
    auto stop2 = high_resolution_clock::now();
    auto duration1 = duration_cast<milliseconds>(stop1 - start1);
    auto duration2 = duration_cast<milliseconds>(stop2 - start2);

    //Printing the time taken to find the columns that take part
    cout<<"Preprocessing time:  ";
    cout << duration1.count() <<"ms"<< endl;

    //Printing the time taken to solve the system
    cout<<"Solve time:  ";
    cout << duration2.count() <<"ms"<< endl;
    return (1);
}

/**
 * Implements a parallel approach for Sparse Triangular Solve namely the level-scheduling approach.
 * The Columns that take part in the computations are divided into levels. Columns in each level are computed in parallel.
 * @param L CSC format of the matrix
 * @param x The right hand side vector of the system
 * @return Whether the operation was successful
 */

int parallel_lsolve_level_scheduling(CSC L, Vector_b x){
    if (!L.col_ptr || !L.row_ind || !x.val) {return (0) ;} // check inputs
    auto start3 = high_resolution_clock::now();

    //Getting the columns that take part in the computations grouped by their levels

    vector<vector<int>> levels= create_levels(L, x);
    auto stop3 = high_resolution_clock::now();
    auto start4 = high_resolution_clock::now();

    //Looping on the levels to solve the system.
    // Lower levels are computed first because the computations of the columns in higher levels rely on them.
    for (int i = 0; i < levels.size(); i++) {
        omp_set_num_threads(6);

    //Computations of columns in each level is parallelized because they are independent of one another.
#pragma omp parallel for
        for(int it = 0 ; it <levels[i].size() ; it++){
            int k = levels[i][it];
            x.val[k] /= L.val[L.col_ptr[k]];
            for (int p = L.col_ptr[k]+1; p<L.col_ptr[k+1]; p++) {
                if (!(abs(L.val[p]) < pow(10, -5) || abs(x.val[k]) < pow(10, -5))) {
                    long double calculated_val = L.val[p] * x.val[k];
#pragma omp atomic //Much lower overhead than the critical operation
                    x.val[L.row_ind[p]] -= calculated_val;
                }

            }
        }
    }
    auto stop4 = high_resolution_clock::now();
    auto duration3 = duration_cast<milliseconds>(stop3 - start3);
    auto duration4 = duration_cast<milliseconds>(stop4 - start4);
    //Printing the time taken to calculate the levels
    cout<<"Preprocessing time:  ";
    cout << duration3.count() <<"ms"<< endl;
    //Printing the time taken to solve the system
    cout<<"Solve time:  ";
    cout << duration4.count() <<"ms"<< endl;
    return (1);
}

/**
 * Implements a parallel approach for Sparse Triangular Solve namely the self-scheduling approach
 * Columns that take part in the computations are calculated in parallel but the thread calculating the calculations
 * must wait for the columns it depends on to start the computations.
 * @param L CSC format of the matrix
 * @param x The right hand side vector of the system
 * @return Whether the operation was successful
 */
int parallel_lsolve_self_scheduling(CSC L, Vector_b x){
    if (!L.col_ptr || !L.row_ind || !x.val) {return (0) ;} // check inputs

    set <int> reached_columns;
    int * dependencies = new int[L.dim]();
    auto start5 = high_resolution_clock::now();

    //Getting the columns each column is depended on
    get_dependencies(L, x, reached_columns,dependencies);
    vector<int> reached_vec (reached_columns.begin(), reached_columns.end());
    auto stop5 = high_resolution_clock::now();
    auto start6 = high_resolution_clock::now();
    omp_set_num_threads(6);

    //Looping through the columns that take part in the computations in parallel
#pragma omp parallel for
    for (int i = 0; i < reached_vec.size(); i++) {
        int col = reached_vec[i];

        //Waiting for all columns that this column's calculations are depended on to finish
        while(dependencies[col]!=0){}
        x.val[col] /= L.val[L.col_ptr[col]];
        for (int p = L.col_ptr[col]+1; p<L.col_ptr[col+1]; p++) {
            long double calculated = 0;
            if (!(abs(L.val[p]) < pow(10, -6) || abs(x.val[col]) < pow(10, -6))) {
                calculated = L.val[p] * x.val[col];
            }
#pragma omp atomic
            x.val[L.row_ind[p]] -= calculated;
//Decreasing the number of unfinished column calculations
#pragma omp atomic
            dependencies[L.row_ind[p]]--;
        }
    }

    auto stop6 = high_resolution_clock::now();
    auto duration5 = duration_cast<milliseconds>(stop5 - start5);
    auto duration6 = duration_cast<milliseconds>(stop6 - start6);

    //Printing the time taken to calculate the dependencies
    cout<<"Preprocessing time:  ";
    cout << duration5.count() <<"ms"<< endl;

    //Printing the time taken to solve the system
    cout<<"Solve time:  ";
    cout << duration6.count() <<"ms"<< endl;
    return (1);
}


/**
 * Implements a single-threaded naive approach to matrix-vector multiplication
 *
 * @param n The dimension of the matrix
 * @param Lp Column pointers in the CSC format of the matrix
 * @param Li Row indices of the nonzero elements in the CSC format of the matrix
 * @param Lx Nonezero elements of the matrix in the CSC format of the matrix
 * @param x The vector
 * @param y The ouput vector
 * @return Whether the operation was successful
 */

int spmv_csc(int n, size_t * Ap,int * Ai, long double* Ax, long double* x, long double* y){
    int p, j;
    if (!Ap || !x || !y) return (0); //Checking inputs
    for (j = 0; j<n; j++){
        for (p = Ap[j]; p<Ap[j+1]; p++){
            if (!(abs(Ax[p])<pow(10,-5) || abs(x[j])<pow(10,-5)))
                y[Ai[p]] += Ax[p]*x[j];
        }

    }
    return (1);
}

