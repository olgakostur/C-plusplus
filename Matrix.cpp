#include<iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <map>
#include <cmath>
#include <cassert>
#include "Matrix.h"

using namespace std;

// constructor
template <class T>
Matrix<T>::Matrix(int rows, int cols, bool preallocate): rows(rows), cols(cols), size_of_values(rows * cols), preallocated(preallocate)
{
   // we own this memory
   if (this->preallocated)
   {
      this->values = new T[this->size_of_values];
   }
}

template <class T>
Matrix<T>::Matrix(int rows, int cols, T *values_ptr): rows(rows), cols(cols), size_of_values(rows * cols), values(values_ptr)
{}

//destructor
template <class T>
Matrix<T>::~Matrix()
{
   // only delete this if we own the memory
   if (this->preallocated)
   {
      delete[] this->values;
   }
}

//print values of matrix one after the other in one line
template <class T>
void Matrix<T>::printValues()
{
    std::cout << "Printing values: " << std::endl;
   for (int i = 0; i < this->size_of_values; i++)
   {
      std::cout << this->values[i] << " ";
   }
   std::cout << std::endl;
}

//print the x vector results into a text file
template <class T>
void Matrix<T>::printTextFile(){
    ofstream vectorResult("x_vector_result.txt");
    for (int i = 0; i < this->size_of_values; i++){
        vectorResult << this->values[i] << endl;
    }
    cout << endl;
}

//print the matrix values row by row
template <class T>
void Matrix<T>::printMatrix()
{
    std::cout << "Printing:";
   for (int j = 0; j < this->rows; j++)
   {
      std::cout << std::endl;
      for (int i = 0; i < this->cols; i++)
      {
         // this is row-major ordering
         std::cout << this->values[i + j * this->cols] << " ";
      }
   }
   std::cout << std::endl;
}

//computes matrix-vector multiplication based on 
//Ax = b
// assumption that output already has the correct sizes
// and that it has memory allocated
template<class T>
void Matrix<T>::matVecMult(T *input, T *output)
{
    //Check that our vectors have been allocated memory space
    if (input == nullptr || output == nullptr)
    {
        std::cerr << "Input or output haven't been created" << std::endl;
        return;
    }

    // Set the output to zero
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 0.0;
    }

    // Loop over each row
    for (int i = 0; i < this->rows; i++)
    {
        // Loop over all the entries in this row for matrix and over each row in vector
        // Multiply each entry in row i of matrix by each row in vector 
        // Then sum them up to get value in row i of output vector
        for (int k = 0; k < this->rows; k++)
        {
            output[i] += this->values[i * this->rows + k] * input[k];
        }
    }
}

 // computes matrix-matrix multiplication
 // output = mat_left * this
 // n * m = n * k * k * m
 // assumption that output already has the correct sizes
 // and that it has memory allocated
template <class T>
void Matrix<T>::matMatMult(Matrix<T>& mat_left, Matrix<T>& output)
{

   // Check our dimensions match
   if (this->cols != output.cols)
   {
      std::cerr << "Input dimensions for matrices don't match" << std::endl;
      return;
   }

   // Check if our output matrix has had space allocated to it
   if (output.values != nullptr) 
   {
      // Check our dimensions match
      if (this->rows != mat_left.cols || mat_left.rows != output.rows)
      {
         std::cerr << "Input dimensions for matrices don't match" << std::endl;
         return;
      }      
   }
   // The output hasn't been preallocated, so we are going to do that
   else
   {
      output.values = new T[this->cols * mat_left.rows];
      output.preallocated = true;
   }

   // Set values to zero before hand
   for (int i = 0; i < output.size_of_values; i++)
   {
      output.values[i] = 0;
   }

    //this uses row-major ordering
   for(int i = 0; i < mat_left.rows; i++)
   {
      for(int j = 0; j < this->cols; j++)
      {
         for(int k = 0; k < mat_left.cols; k++)
         {            
               output.values[i * this->cols + j] += mat_left.values[i * mat_left.cols + k] * this->values[k * this->cols + j];
         }
      }
   }
}

//function that adds to matrices
template <class T>
void Matrix<T>::matMatAdd(Matrix<T>& mat_left, Matrix<T>& output)
{

    // Check our dimensions match
    if (this->cols != output.cols)
    {
        cerr << "Input dimensions for matrices don't match" << endl;
        return;
    }

    // Check if our output matrix has had space allocated to it
    if (output.values != nullptr)
    {
        // Check our dimensions match
        if (this->rows != mat_left.cols || mat_left.rows != output.rows)
        {
            cerr << "Input dimensions for matrices don't match" << endl;
            return;
        }
    }
    // The output hasn't been preallocated, so we are going to do that
    else
    {
        output.values = new T[this->rows * mat_left.cols];
        output.preallocated = true;
    }

    // Set values to zero before hand
    for (int i = 0; i < output.size_of_values; i++)
    {
        output.values[i] = 0;
    }

    //this uses row-major ordering
    for (int i = 0; i < mat_left.rows; i++)
    {
        for (int j = 0; j < this->cols; j++)
        {
            output.values[i * this->cols + j] =  mat_left.values[i * mat_left.cols + j] + this->values[i * this->cols + j];
        }
    }
}

//function that subtracts one matrix from another matrix
template <class T>
void Matrix<T>::matMatSubtract(Matrix<T>& mat_left, Matrix<T>& output)
{

    // Check our dimensions match
    if (this->cols != output.cols)
    {
        cerr << "Input dimensions for matrices don't match" << endl;
        return;
    }

    // Check if our output matrix has had space allocated to it
    if (output.values != nullptr)
    {
        // Check our dimensions match
        if (this->rows != mat_left.cols || mat_left.rows != output.rows)
        {
            cerr << "Input dimensions for matrices don't match" << endl;
            return;
        }
    }
    // The output hasn't been preallocated, so we are going to do that
    else
    {
        output.values = new T[this->rows * mat_left.cols];
        output.preallocated = true;
    }

    // Set values to zero before hand
    for (int i = 0; i < output.size_of_values; i++)
    {
        output.values[i] = 0;
    }

    //this uses row-major ordering
    for (int i = 0; i < mat_left.rows; i++)
    {
        for (int j = 0; j < this->cols; j++)
        {
            output.values[i * this->cols + j] = mat_left.values[i * mat_left.cols + j] - this->values[i * this->cols + j];
        }
    }
}

//this function performs a Gaussian elimination with partial pivoting
//partial pivoting takes care of division by zero and round off errors
//returns vector x by reference using output
template <class T>
void Matrix<T>::gauss_elimination(Matrix<T>& mat_left, Matrix<T>& output)
{
   //check our dimensions match
   //checking that rows of matrix A = cols of matrix A, rows of matrix A = rows of vector b, rows of matrix A = rows of vector x
   
   if (mat_left.rows != mat_left.cols || mat_left.rows != this->rows || mat_left.rows != output.rows)
   {
      std::cerr << "Input dimensions for matrices and vectors don't match" << std::endl;
      return;
   }

   //perform the calculations for obtaining upper-triangular forms for the matrix A and vector b
   //this includes row swapping to avoid division by zero and round off errors
   //algorithm performed in L5_linear_solvers
   for (int k = 0; k < (mat_left.rows - 1); k++)
   {
       //initialize kmax = k
       int kmax = k;
       //iterate over all the elements in the column and see whether there is a larger value
       for (int m = (k+1); m < mat_left.rows; m++){
           if (abs(mat_left.values[mat_left.rows * kmax + k]) < abs(mat_left.values[m * mat_left.rows + k])){
               //set kmax to be the position of the largest element
               kmax = m;
           }
       }
        //if we have different rows, we will perform row swapping
        if (kmax != k){
        //initialize two arrays to temporarily store the original row elements when performing the swap
        double *cop_A = new double[mat_left.rows];
        double cop_b;

        //we copy into the copy array for matrix A the new row we want to replace
        for (int cop1 = 0; cop1 < mat_left.rows; cop1++){
            cop_A[cop1] = mat_left.values[kmax * mat_left.rows + cop1];
        }
        //we copy into the copy array for vector b the new row we want to replace
        cop_b = this->values[kmax];

        //we copy the old row into the new row for matrix A
        for (int cop2 = 0; cop2 < mat_left.rows; cop2++){
            mat_left.values[kmax * mat_left.rows + cop2] = mat_left.values[k * mat_left.rows + cop2];
        }
        //we copy the old row into the new row for vector b
        this->values[kmax] = this->values[k];

        //we copy the copy into the old row for matrix A
        for (int cop3 = 0; cop3 < mat_left.rows; cop3++){
            mat_left.values[k * mat_left.rows + cop3] = cop_A[cop3];
        }

        //we copy the copy into the old row for vector b
        this->values[k] = cop_b;
        }

        //this converts both the matrix A and the vector b into upper-triangular form
        //i.e. matrix A will have 0 in its lower triangle based on row operations
        //those same operations will be applied to vector b
      for (int i = (k+1); i < (mat_left.rows); i++)
      {
          double s = mat_left.values[i * mat_left.rows + k] / mat_left.values[(mat_left.rows + 1) * k];
         for (int j = k; j < mat_left.rows; j++)
         {
            mat_left.values[i * mat_left.rows + j] -= s * mat_left.values[k * mat_left.rows + j];
         }
         this->values[i] -= s * this->values[k];
      }
   }

    //perform the calculation for back substitution and obtaining vector x
    for (int k = (mat_left.rows-1); k > -1; k--){
        double s = 0.0;
        for (int j = (k+1); j < mat_left.rows; j++){
            s += mat_left.values[mat_left.rows * k + j] * output.values[j];
        }
        output.values[k] = (this->values[k] - s) / mat_left.values[(mat_left.rows + 1) * k];
    }
}

//function that tests the algorithm is correct (i.e. computes Ax - b and checks that it is close to zero)
template <class T>
void Matrix<T>::testResultGauss(Matrix<T> &mat_left, Matrix<T> &vector_right, Matrix<T> &output){
    //create a new array to store Ax (using matMatMult function)
    auto *result = new Matrix<double>(mat_left.rows, 1, true);
    //perform matrix multiplication between matrix A and vector x and store in result
    this->matMatMult(mat_left, *result);
    //calculate (Ax - b) and store in the output array
    for (int i = 0; i < mat_left.rows; i++){
        output.values[i] = result->values[i] - vector_right.values[i];
    }
    //check that the output array is close to zero array
    //if it is not, then we will have an assertion error on the screen
    //if no assertion error, then the test passes
    for (int i = 0; i < mat_left.rows; i++){
        assert(abs(round(output.values[i])) == 0);
    }

    cout << "Method is correct and algorithm works" << endl;
    //delete result so that we don't leak memory
    delete result;
}

//function that tests the algorithm is correct
//computes Ax - b and checks that it less than or equal to tolerance level
template <class T>
void Matrix<T>::testResultJacobi(Matrix<T> &mat_left, Matrix<T> &vector_right, Matrix<T> &output, double tolerance){
    //create a new array to store Ax (using matMatMult function)
    auto *result = new Matrix<double>(mat_left.rows, 1, true);
    //perform multiplication of matrix A and vector x and store in result
    this->matMatMult(mat_left, *result);
    //calculate (Ax - b) and store in the output array
    for (int i = 0; i < mat_left.rows; i++){
        output.values[i] = result->values[i] - vector_right.values[i];
    }

    //initialize residual squared, which will store the square root of sum squares of each element of the residual
    //this residual squared needs to be less or equal to tolerance level
    double resid_squared = 0;
    for (int i = 0; i < mat_left.rows; i++)
    {
        //compute the square of each element of residual and sum
        resid_squared += pow(output.values[i], 2);
    }
    //take the square root of the sum squares
    resid_squared = sqrt(resid_squared);

    //assert whether the residual squared is less than or equal to tolerance
    //if the test passes, there is no assertion error on the screen
    //if the test doesn't pass, there will be assertion error on the screen
    assert(resid_squared <= tolerance);

    cout << "Method is correct and algorithm works" << endl;
    delete result;
}

//function that computes the residual for the Jacobi method (Ax - b)
//this function needs to be called on vector b
//returns the residual
template <class T>
double Matrix<T>::computeResidual(Matrix<T>& mat_left, Matrix<T>& vector_x, Matrix<T> &output) {


    //create a new array to store Ax (using matMatMult function)
    auto* result = new Matrix<double>(mat_left.rows, 1, true);
    double resid_squared = 0;
    double resid;
    //compute the multiplication between matrix A and vector x and store in result
    vector_x.matMatMult(mat_left, *result);
    //calculate (Ax - b) and store in the output array
    for (int i = 0; i < mat_left.rows; i++) {
        //the output array is residual
        output.values[i] = this->values[i] - result->values[i];
    }
    //iterate over the number of rows
    for (int i = 0; i < mat_left.rows; i++)
    {
        //compute the square of each element of residual and sum
        resid_squared += pow(output.values[i], 2);
    }
    //take the square root of the square residual and return it for the Jacobi method
    resid = sqrt(resid_squared);

    //return the square root of the square residual
    return resid;

    //delete the result so that we don't leak memory
    delete result;
}


//function that calculates the Jacobi method
//based on L5_linear_solvers notes
//it return the vector x by reference via output
template <class T>
void Matrix<T>::jacobi(Matrix<T> &mat_left, Matrix<T> &output, double tolerance){

// reminder of syntax from user input syntax file:

// right_side_A: sub array of external input matrix A on right side of the diagonal
// left_side_A: sub array of external input matrix A on the left side of the diagonal
// mat_left: Matrix A
// output: vector x
// vector_right : vector b

    int rows = mat_left.rows;
    int cols = mat_left.cols;

    // x_new gets filled in at each iteration
    vector<double> x_new(rows);

    //splitting the matrix and vector into different parts of the diagonal: left and right
    vector<double> right_side_A;
    vector<double> left_side_A;
    vector<double> right_side_x;
    vector<double> left_side_x;
    double diagonal;

    //initialize maximum number of iterations to be 1000
    int iter_max = 1000;
    int iter = 0;

    //initialize vector x to be 0
    for (int i = 0; i < rows; i++){
        x_new.push_back(0.0);
    }

    //initialize test_tol to be false, it will be true when the residual will be less than or equal to tolerance level
    bool test_tol = false;

    //the while loop ends when we reach maximum number of iterations or when the residual is less than tolerance level
    while (iter < iter_max || test_tol == false)
    {
        //iterate over number of matrix rows
        for (int i = 0; i < rows; i++)
        {
            //find the diagonal element of the row
            diagonal = mat_left.values[i * (rows + 1)];
            for (int j = 1; j < (rows - i); j++)
            {
                //fill the right side vector (elements to the right of diagonal)
                right_side_A.push_back(mat_left.values[i * (rows + 1) + j]);
                //fill the right side vector
                right_side_x.push_back(output.values[j + i]);
            }
            for (int k = rows; k < (rows + i); k++)
            {
                //fill the left side vector (elements to the left of diagonal)
                left_side_A.push_back(mat_left.values[rows * i + k - rows]);
            }
            for (int s = 0; s < i; s++)
            {
                //fill the left side vector
                left_side_x.push_back(output.values[s]);
            }
            
            //initialize dot product results to be zero
            double dot_product_right = 0;
            double dot_product_left = 0;

             for (int i = 0; i < left_side_A.size(); i++)
            {
                //multiply the left side of matrix A with left side of vector x
                dot_product_left += left_side_A[i] * left_side_x[i];
            }
            for (int i = 0; i < right_side_A.size(); i++)
            {
                //multiply the right side of matrix A with left side of vector x
                dot_product_right += right_side_A[i] * right_side_x[i];
            }
            //assign to the x vector the temporary value
            x_new[i] = (1 / diagonal) * (this->values[i] - dot_product_left - dot_product_right);
            //reset the diagonal to be 0 and clear the vectors
            diagonal = 0;
            left_side_x.clear();
            right_side_x.clear();
            left_side_A.clear();
            right_side_A.clear();
        }
        //increment the iteration counter
        iter++;

        //calculate the residual (Ax - b) and check against the tolerance level
        //calling computeResidual function from above
        auto* difference = new Matrix<double>(mat_left.rows, 1, true);
        if (tolerance >= this->computeResidual(mat_left, output, *difference))
        {
            cout << "Tolerance Level reached at iter: " << iter << endl;
            //if residual is less than or equal to tolerance level, then we set test_tol = true
            //this will exit the while loop and the Jacobi method
            test_tol = true;
            break;
        }

        //copy the values of vector x into output -- will be passed by reference to the main file
        for (int i = 0; i < rows; i++){
            output.values[i] = x_new[i];
        }

        //delete difference so that we don't leak memory
        delete difference;
    }
}