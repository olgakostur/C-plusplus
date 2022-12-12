#pragma once

template <class T>
class Matrix
{
public:

   //constructor where we want to preallocate ourselves
   Matrix(int rows, int cols, bool preallocate);
   //constructor where we already have allocated memory outside
   Matrix(int rows, int cols, T *values_ptr);
   //destructor
   virtual ~Matrix();

   //print out the values of matrix
   void printValues();
   virtual void printMatrix();

   //print out the x vector in a text file
   void printTextFile();

   //function that verifies that Ax - b < tolerance level (algorithm works)
   void testResultJacobi(Matrix<T> &mat_left, Matrix<T> &vector_right, Matrix<T> &output, double tolerance);

   //function that verifies that Ax - b = 0 (Gauss algorithm has residual 0)
   void testResultGauss(Matrix<T> &mat_left, Matrix<T> &vector_right, Matrix<T> &output);

   //functions that performs matrix vector multiplication
   void matVecMult(T *input, T *output);

   //function that performs gauss elimination
   //converts the matrix A and vector b into upper-triangular forms 
   //computes back substitution to compute vector x
   //based on L5_linear_solvers algorithm
   void gauss_elimination(Matrix<T> &mat_left, Matrix<T> &output);

   //function that performs matrix matrix multiplication
   void matMatMult(Matrix<T> &mat_left, Matrix<T> &output);

   //function that performs addition of two matrices
   virtual void matMatAdd(Matrix<T> &mat_left, Matrix<T> &output);

   //function that performs subtraction of one matrix from another matrix
   virtual void matMatSubtract(Matrix<T> &mat_left, Matrix<T> &output);

   //function that performs Jacobi method
   //based on L5_linear_solvers algorithm
   void jacobi(Matrix<T> &mat_left, Matrix<T> &output, double tolerance);

   //function that computes the residual i.e. Ax - b, used for Jacobi method
   double computeResidual(Matrix<T>& mat_left, Matrix<T>& vector_x, Matrix<T> &output);

   //explicitly using the C++11 nullptr here
   T *values = nullptr;   
   int rows = -1;
   int cols = -1;

//accessible by the subclass
protected:
	bool preallocated = false;

//private variables - there is no need for other classes 
//to know about these variables
private:

   int size_of_values = -1;
   
};