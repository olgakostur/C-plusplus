
#  :brain:  Meeting Notes – Wednesday 19th January 2022 🧠
*3-4pm*


### General Review of the Group Assignment Instructions



**General Project Summary:** 

-	4 files: header file (.h), a report file (.pdf), two c++ files (dense.cpp, sparse.cpp)
-	Solve the linear system: Ax=b
-	A: matrix, x: vector, b: vector
-	Input: A (matrix) and b (vector)
-	Has to be greater than 10 x 10
-	Symmetric matrix (square)
-	Output= x
-	Main method takes input of A, x and b
-	1 method per solver, no class/per solver nor new main/solver
-	Multiple solvers can be direct solver or iterative

---


### Group Discussion

**General Project Etiquette:**

-	No use of liveshare -> gets too messy with compiling and building
-	Using branches instead with one branch per person
-	Reminder to pull the github repository often to avoid conflicts
-	Use of Microsoft Teams group for direct messages and questions + use of file sharing in Teams
-	Commented pull requests
-	Comment code and explain design 

**Project Design Ideas**
-	Class: for matrix x vector
-	Constructor to create vector

**Assessing Performance**
-	? tbd
-	Using ctime in some areas that have potential for performance improvements
-	NB: don’t overdo it as this also has a computational cost

**User Experience**
-	Make it smooth by checking for errors (eg. Dimensions, user inputs, format) and outputting cerr messages 


**Actions for 20/01/22:**

-	Work on dense matrix first 
-	Pick linear system eg. Simple large square matrix with random b vector 
-	Work on user inputs for matrix in C++
-	Look at inverse multiplication, need for upper triangular method (and lower triangular method)
-	Questions for TAs

**Questions for TA**

-	Ways to assess performance and compliance? Is there anything similar to PEP8 in Python for C++
-	Review project design ideas and especially structure with TAs: eg. Class, constructors, methods 
-	Meeting at 2:40pm ideally??


