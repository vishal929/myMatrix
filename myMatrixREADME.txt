README!


myMatrix is my first java project class for matrix operations!!

how to declare a matrix:
myMatrix test = new myMatrix(3,3);
^^ the above code constructs a new matrix called test that has 3 rows and 3 columns^^
myMatrix test2 = new myMatrix(3);
^^Alternatively, the above code constructs a new square matrix with 3 rows and 3 columns^^


how to add entries to a matrix:
test.addEntry(1,2,3,4,5,6,7,8,9);
^^values are entered left to right to fill the matrix^^
double[] entries = new double[9];
for (int i=0;i<9;i++) {
  entries[i]=i+1;
}
test2.addEntry(entries);
^^If you already have an array with the entries you wish to add to the matrix, you can do so like the above code^^
test.fillColumn(double[]...ds);
^^If you have arrays which represent the columns of the matrix, you can use this method to fill the columns ^^
test.fillRows(double[]..ds);
^^If you have arrays with represent the rows of the matrix, you can use this method to fill the rows^^






The class currently supports the following operations:

1)Adding entries to an empty matrix (addEntry method)(or user can use fillColumns() or fillRows())
2)Changing entries on a matrix (changeEntry method)
3)Getting a row or column of a matrix as an array (getColumn/getRow methods)
4)Getting a single entry of a matrix (getEntry method)
3)Adding/Subtracting matrices together (addMatrix/subtractMatrix methods)
4)Scaling a matrix by some factor (scaleMatrix method)
5)Calculating the matrix exponential (matrixExponential method)
6)Matrix multiplication (multiplyMatrix method)
7)Scaling a row by some factor in a matrix (multiplyRow method)
8)adding rows of a matrix together(addRow method)
9)Swapping rows of a matrix (swapRow method)
10) Getting the row echelon form of a matrix(rowEchelon method)
11) Getting the reduced Row echelon form (reducedRowEchelon method)
12) Solving a system given an augmented column (solveSystem)
13) Getting determinant and inverse of a matrix (getDeterminant/getInverse methods)
14)Getting projection of one vector onto another (getProjection method)
15)Transposing a matrix (getTranspose() method)
16)Getting orthogonal basis of a matrix via gramSchmidt (gramSchmidt method)
17)Getting the QR factorization of a matrix (QR method)(note: my code only has functionality for square matrices)
18)Getting the real and complex eigenvalues of any REAL matrix(eigenValues method)(this is via QR algorithm)
19)Normalizing columns of a matrix (normalize() method)
20)Getting the inverse of an invertible matrix (getInverse() method)
