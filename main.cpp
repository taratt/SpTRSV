#include <iostream>
#include <cmath>


#include "inputFormatter.h"
#include "outputVerifier.h"
#include "matrixCalculations.h"


int main(int argc, char *argv[]) {
//checking the number of inputs
    if (argc != 3){
        cout<< "File paths to the matrix and the vector have to be provided"<<endl;
        return (1);
    }

    string matr_name = argv[1];
    string right_vec_name = argv[2];

//Turning the matrix and the vector into CSC format
    CSC matr = format_csc(matr_name);
    Vector_b a = read_b(right_vec_name);
    Vector_b b = read_b(right_vec_name);
    Vector_b c = read_b(right_vec_name);
    Vector_b d = read_b(right_vec_name);
    Vector_b ground = read_b(right_vec_name);

    int n = matr.dim;

//solving using the single-threaded naive approach
cout << "Reporting on the single-threaded naive implementation:"<<endl;
lsolve(n, matr.col_ptr, matr.row_ind, matr.val,a.val);
long double* y;
verify_original_results(n, matr.col_ptr, matr.row_ind, matr.val, a.val, ground.val,y);
cout<<endl;

//solving using the single-threaded optimized approach
cout << "Reporting on the single-threaded optimized implementation:"<<endl;
optimized_single_lsolve(matr, b);
long double* p;
verify_original_results(n, matr.col_ptr, matr.row_ind, matr.val,  b.val, ground.val,p);
verify_equal_to_naive(n, b.val, a.val);
    cout<<endl;

//solving using the parallel level-scheduling approach
cout << "Reporting on the parallel level scheduling implementation:"<<endl;
parallel_lsolve_level_scheduling(matr, c);
long double* m;
verify_original_results(n, matr.col_ptr, matr.row_ind, matr.val,  c.val, ground.val,m);
verify_equal_to_naive(n, c.val, a.val);
    cout<<endl;

//solving using the parallel self-scheduling approach
//This section is commented because it does not preform well with large matrices without great processing resources (a large number of cores)
//cout << "Reporting on the parallel self scheduling implementation:"<<endl;
//parallel_lsolve_self_scheduling(matr, d);
//long double* q;
//verify_original_results(n, matr.col_ptr, matr.row_ind, matr.val, d.val, ground.val,q);
//verify_equal_to_naive(n, d.val, a.val);

}
