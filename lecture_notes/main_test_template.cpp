#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdlib.h> 
#include <typeinfo>
#include "Matrix.h"
#include "Matrix.cpp"

using namespace std;

int readDoubles(double **input_doubles)
{
   //define the strings that will be inputted by the user for the type of input and text file name
   string val, text_file_name;
   cout << "Please type in F for text/csv file input or K for keyboard input: " << endl;
   //input either F or K based on input type
   cin >> val;

   //if the user selects F, we will use a text/csv file input
   if (val == "F"){
      //needs to be in the format: number of rows, number of cols, then the matrix A, then vector b
      //assume that vector b number of rows equals matrix A number of rows
      //assume that matrix A rows = matrix A cols - square matrix
      //assume that the matrix A inputted actually matches the rows and cols
      cout << "File needs to be in the same location as the .cpp and .h files. " << endl;
      cout << "File needs to have comma separated values. The format of file is: " << endl;
      cout << "First row has the number of rows and cols of the matrix: rows, cols " << endl;
      cout << "Next rows have the elements of the matrix " << endl;
      cout << "The final row of the file has the elements of the vector " << endl;
      cout << "For more information, please see documentation file. " << endl;
      cout << "Please type the file name to input (example: test.txt or test.csv): " << endl;
      //ask user to input file name (this file needs to be in the same location as the .cpp and .h files)
      cin >> text_file_name;
      //open the text/csv file
      fstream myfile;
      myfile.open(text_file_name, ios::in);
      vector<string> result;
      //check we opened the file successfully, otherwise show error
      if (myfile.fail())
	   {
		   cout<<"Error opening file"<<endl;
		   return(-1);
	   }

      //read each line from the file and place the values in the string vector
      while (!myfile.eof()){
            string input;
            getline(myfile, input);
            stringstream stream(input);
            while (stream.good()){
               string substr;
               getline(stream, substr, ',' );
               result.push_back(substr);
            }
      }
      //build a double array to store the converted strings
      (*input_doubles) = new double[result.size()];
      //get a pointer to the string data in vector
      string *result_data = result.data();
      //go through and convert the strings to doubles
      for (int i = 0; i < result.size(); i++)
      {
         (*input_doubles)[i] = atof(result_data[i].c_str());
      }       
      //store how many input values we have
      int no_input_vals = result.size();
      //clear the string vector
      result.clear();
      //close the file
	   myfile.close();  
      //return how many values we have in the array
      return no_input_vals;
   }

   //if the user selects K, we will use keyboard input
   else if (val == "K"){
      string rows, cols;
      //user inputs number of rows
      cout << "Please input number of rows for matrix A (as an int): " << endl;
      cin >> rows;

      //user inputs number of columns
      cout << "Please input number of columns for matrix A (as an int): " << endl;
      cin >> cols;

      //check that columns equals numbers of rows - otherwise ERROR
      if (cols != rows){
         cerr << "Number of columns needs to be equal to number of rows. Please start again." << endl;
      }

      //initialize a string vector
      vector<string> result;

      //push back the rows and cols as first elements of the vector
      result.push_back(rows);
      result.push_back(cols);

      //input the matrix elements by row
      //assume that number of elements matches rows x cols
      cout << "Please input the matrix A elements by row: " << endl;
      for (int i = 1; i <= (stoi(rows)+1); i++){
         //reading the matrix elements and pushing them into the vector
         string input;
         getline(cin, input);
         stringstream stream(input);
         while (stream.good()){
               string substr;
               getline(stream, substr, ' ');
               result.push_back(substr);
         }
      }

      //input the vector b in one single row
      //assume that number of elements matches rows of matrix A
      cout << "Please input the vector b on one single row: " << endl;
      //reading the vector elements and pushing them into the vector
      string input;
      getline(cin, input);
      stringstream stream(input);
      while (stream.good()){
         string substr;
         getline(stream, substr, ' ');
         result.push_back(substr);
      }

      //erasing third element of the vector because it is blank
      result.erase(next(result.begin(), 2));

      //build a double array to store the converted strings
      (*input_doubles) = new double[result.size()];
      //get a pointer to the string data in vector
      string *result_data = result.data();

      //go through and convert the strings to doubles
      for (int i = 0; i < result.size(); i++)
      {
         (*input_doubles)[i] = atof(result_data[i].c_str());
      }       
      //store how many input values we have
      int no_input_vals = result.size();
      //clear the string vector
      result.clear();
      //return the size of the array
      return no_input_vals;
   }
   //else return 0
   return 0;
}

int main(){
   //build double arrays to store the converted strings
   double *input_doubles;
   double *input_matrix;
   double *input_vector;
   int no_input_vals;
   int rows, cols;
   //retrieve the number of elements of the large array and the modified array
   no_input_vals = readDoubles(&input_doubles);

   //set rows to be first element of the array
   rows = input_doubles[0];
   //set cols to be second element of the array
   cols = input_doubles[1];

   //build the matrix array - from the large array
   input_matrix = new double[rows * rows];
   //start from the third element from the large array, has size of rows x rows
   //assume that rows = cols (square matrix)
   for (int i = 0; i < rows * rows; i++){
      input_matrix[i] = input_doubles[2 + i];
   }

   //build the vector array - from the large array
   input_vector = new double[rows];
   //insert the last *rows* elements from the large array
   for (int i = 0; i < rows; i++){
      input_vector[i] = input_doubles[no_input_vals - rows + i];
   }

   //create a dense matrix, using class Matrix
   //initialize it with the matrix inputted from user
   auto *dense_mat = new Matrix<double>(rows, cols, input_matrix);

   //now let's test printing our matrix
   dense_mat->printValues();
   //now let's test printing our matrix with our other function
   dense_mat->printMatrix();   

   //calling destructors so that we don't leak memory
   delete[] input_matrix;
   delete[] input_vector;
   delete[] input_doubles;
   delete dense_mat;

   return 0;
}