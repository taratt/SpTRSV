//
// Created by taratt on 2/7/23.
//

/**
 * This library contains functions to verify the result of the triangular solve algorithm
*/
#pragma once
#include <iostream>
#include <cmath>

#include "matrixCalculations.h"

/**
 * Checks if the matrix multiplied by the result of the triangular solve algorithm is equal to the right hand side vector
 *
 * @param n The dimension of the matrix
 * @param Lp Column pointers in the CSC format of the matrix
 * @param Li Row indices of the nonzero elements in the CSC format of the matrix
 * @param Lx Nonezero elements of the matrix in the CSC format of the matrix
 * @param x The result of the triangular solve algorithm
 * @param b The right hand side vector of the system
 * @param y A vector that saves the multiplication results
 * @return Whether the algorithm is completely verifies
 */

int verify_original_results( int n, size_t *Ap, int * Ai, long double* Ax,  long double* x, long double * b, long double* y){
    int res;
        y = new long double[n]();
        //Calculating the multiplication of the matrix and the result of the triangular solve algorithm
        res = spmv_csc(n, Ap, Ai, Ax, x, y);
    if (!res) return (0);
    int counter = 0;
    bool nonz = false;
    //Looping over elements of the vectors and see if they are equal
    for (int i = 0; i < n; ++i) {
        if( abs(y[i]- b[i])>0.01f){
            if (b[i]!=0)
                nonz = true;
//            cout<<i<<" "<<y[i]<<" "<<b[i]<<" "<<endl;
            counter++;
        }

    }
    if (counter>0) {
        //Checking if the nonzero elements of the right-hand side were the exact same
        if(!nonz)
            cout<<"All nonzeros were correct."<<endl;
        float percent = round((float(n-counter)/n)*1000)/10;
        //Outputting the percentage of elements that were correct giving a sense of how close the answer is to the exact answer
        cout<<"Approximately "<<percent<<"% of the indices were equal to the result of multiplication."<<endl;

        return (0);
    }
    cout<<"The multiplication completely verifies the results"<<endl;
    return (1);
}

/**
 * Checks if the result of the triangular solve algorithm is equal to the result calculated by the naive approach
 *
 * @param n The dimension of the result
 * @param b The result of the triangular solve algorithm
 * @param y The result of the naive triangular solve algorithm
 * @return Whether the algorithm is confidently verified
 */
int verify_equal_to_naive(int n, long double *b, long double *y){
    bool nonz = false;
    int counter = 0;
    //comparing the elements of the result of the naive approach to the given results
    for (int i = 0; i < n; ++i) {
        if( abs(y[i]- b[i])>0.01){
            if (b[i]!=0)
                nonz = true;
            counter++;
        }
    }
    if (counter>0) {
        //Checking if the nonzero elements of each result were the exact same
        if(!nonz)
            cout<<"All nonzeros were correct."<<endl;
        float percent = round((float(n-counter)/n)*1000)/10;
        //Outputting the percentage of elements that were same giving a sense of how close the answer is to the exact answer
        cout<<"Approximately "<<percent<<"% of the indices were equal to the naive results."<<endl;
        return (0);
    }
    cout<<"The optimized approach gives the same result as the naive approach"<<endl;
    return (1);
}
