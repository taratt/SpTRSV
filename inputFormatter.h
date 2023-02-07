//
// Created by taratt on 2/3/23.
//
/**
 * This library holds structs and functions to convert the context of .mtx files containing sparse matrices to CSC format.
*/

#pragma once
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <math.h>

using namespace std;

/**
 * This struct represents the CSC format of a matrix.
 * fields:
 *      dim: The dimension of the square matrix
 *      val: An array of the nonzero values
 *      row_ind: An array of the nonzero elements' rows
 *      col_ptr: An array of pointer to the first element of each column
 *      nonzeros: The number of nonzero elements of the matrix
*/
struct CSC {
    int dim;
    long double * val;
    int * row_ind;
    size_t * col_ptr;
    int nonzeros;
};

/**
 * This struct represents a vector.
 * fields:
 *      dim: The number of rows in the vector
 *      val: An array of the nonzero values
 *      nonzeros: The number of nonzero elements of the matrix
 *      nonzero_rows: An array of the nonzero elements' rows
*/

struct Vector_b{
    int dim;
    long double * val;
    int nonzeros;
    int * nonzero_rows;
};

/**
 * Reads a .mtx file and turns it into a matrix with CSC format
 *
 * @param filePath The path to the matrix file
 * @return The CSC of the matrix
 */

CSC format_csc(std::string filePath){
    std::ifstream fin(filePath);
    // Skipping the headers and comments
    while (fin.peek() == '%') {
        fin.ignore(2048, '\n');
    }
    // Reading number of rows, columns and nonzero elements
    int m, n, nonzeros;
    fin >> m >> n >> nonzeros;
    CSC input_matrix;
    input_matrix.dim = m;
    input_matrix.val = new long double [nonzeros];
    input_matrix.col_ptr = new size_t [n + 1]();
    input_matrix.row_ind = new int [nonzeros];
    int col_pointer = 0;
    int counter = 0;

    // Reading nonzero values and filling the val, col_ptr and row_ind arrays
    for (int l = 0; l < nonzeros-counter; l++){
        int row, col;
        long double data;
        fin >> row >> col >> data;

    // Only reading the lower triangular elements of the matrix
        if (row >= col) {
            input_matrix.row_ind[l] = row - 1;
            input_matrix.val[l] = data;
            while (col_pointer < col){
                input_matrix.col_ptr[col_pointer] = l;
                col_pointer++;
            }
        } else{
            l--;
            counter++;
        }
    }
    input_matrix.nonzeros = nonzeros - counter;
    while (col_pointer <= n) {
        input_matrix.col_ptr[col_pointer] = input_matrix.nonzeros;
        col_pointer++;
    }
    fin.close();
    return input_matrix;
}

/**
 * Reads a .mtx file and turns it into a Vector_b
 *
 * @param filePath The path to the vector file
 * @return The Vector_b of the vector
 */
Vector_b read_b(std::string filePath){
    std::ifstream fin(filePath);
    // Skipping the headers and comments
    while (fin.peek() == '%') {
        fin.ignore(2048, '\n');
    }
    // Reading number of rows, columns and nonzero elements
    int m, n, nonzeros;
    fin >> m >> n >> nonzeros;
    Vector_b b;
    b.val = new long double [m]();
    b.dim = m;
    b.nonzeros = nonzeros;
    b.nonzero_rows = new int[nonzeros];

    // Reading nonzero values and filling the val and nonzerorows arrays
    for (int l = 0; l < nonzeros; l++) {
        int row, col;
        long double data;
        fin >> row >> col >> data;
        b.nonzero_rows[l] = row - 1;
        b.val[row - 1] = data;
    }

    fin.close();
    return (b);
}

