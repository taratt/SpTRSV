//
// Created by taratt on 2/7/23.
//
#pragma once
#include <iostream>
#include <cmath>

#include "matrixCalculations.h"


int verify_original_results( int n, size_t *Ap, int * Ai, long double* Ax,  long double* x, long double * b, long double* y){
    int res;
        y = new long double[n]();
        res = spmv_csc(n, Ap, Ai, Ax, x, y);
    if (!res) return (0);
    int counter = 0;
    bool nonz = false;
    for (int i = 0; i < n; ++i) {
        if( abs(y[i]- b[i])>0.01f){
            if (b[i]!=0)
                nonz = true;
//            cout<<i<<" "<<y[i]<<" "<<b[i]<<" "<<endl;
            counter++;
        }

    }
    if (counter>0) {
        if(!nonz)
            cout<<"All nonzeros were correct."<<endl;
        float percent = round((float(n-counter)/n)*1000)/10;
        cout<<"Approximately "<<percent<<"% of the indices were equal to the result of multiplication."<<endl;
//        cout<<counter<<endl;
        return (0);
    }
    cout<<"The multiplication completely verifies the results"<<endl;
    return (1);
}

int verify_equal_to_naive(int n, long double *b, long double *y){
    bool nonz = false;
    int counter = 0;
    for (int i = 0; i < n; ++i) {
        if( abs(y[i]- b[i])>0.01){
            if (b[i]!=0)
                nonz = true;
            //           cout<<i<<" "<<y[i]<<" "<<b[i]<<endl;
            counter++;
        }
    }
    if (counter>0) {
        if(!nonz)
            cout<<"All nonzeros were correct."<<endl;
        float percent = round((float(n-counter)/n)*1000)/10;
        cout<<"Approximately "<<percent<<"% of the indices were equal to the naive results."<<endl;
        return (0);
    }
    cout<<"The optimized approach gives the same result as the naive approach"<<endl;
    return (1);
}
