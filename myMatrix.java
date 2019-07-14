import java.util.Arrays;
import java.util.ArrayList;
public class myMatrix {

	private double[][] matrix;

					//constructor initializes the matrix
	public myMatrix(int numberofrows, int numberofcolumns) {
		this.matrix=new double[numberofrows][numberofcolumns];

		//column have the same number of entries as the number of rows


	}

					//overloading constructor for square matrix
	public myMatrix(int square) {
		this.matrix= new double[square][square];
	}
		//checks if the matrix is a valid matrix
	public boolean isMatrix() {
		if (this==null || (this.getRows()==0 && this.getColumns()==0)) {
			return false;
		} else {
			return true;
		}
	}

	public int getRows() {
		//error case----empty matrix
			if (this.matrix.length!=0 && this.matrix[0].length!=0) {
		return this.matrix.length;
			} else {
				System.out.println("The matrix is empty!");
				return -1;
			}
	}

	public int getColumns() {
		if (this.matrix.length!=0 && this.matrix[0].length!=0) {
			return this.matrix[0].length;
		} else {
			System.out.println("The matrix is empty!");
			return -1;
		}
	}



		//adds entries from left to right to our matrix
	public void addEntry(double ...ds ) {
		// error case--- dimensions do not match matrix dimensions
		if (ds.length!=this.matrix.length*this.matrix[0].length) {
			System.out.println("The Matrix has "+this.matrix.length+" rows and "+ this.matrix[0].length+ " columns.");
			System.out.println("Please call the function again with a valid number of entries.");
		}
		else {
			//iterates through each index in the matrix from left to right and adds the input
			// we need a count for the parameter
			int count=0;
			for (int i=0;i<this.matrix.length;i++) {
				for (int j=0;j<this.matrix[0].length;j++) {

				this.matrix[i][j]=ds[count];
				count++;
				}
			}
		}
	}


	//if a user has row arrays defined already, they can add them to an empty matrix with this method
	public void fillRows(double[]... ds){
		if (ds[0].length!=this.getColumns() || ds.length!=this.getRows()) {
			System.out.println("Dimension mismatch! Please check your rows and the matrix and try again!");
		} else {
			double[] entries = new double[ds.length*ds[0].length];
			//adding entries into one big array to use addEntry
			int count=0;
			for (int i=0;i<ds.length;i++) {
				for (int j=0;j<ds[0].length;j++) {
					entries[count]=ds[i][j];
					count++;
				}
			}
			this.addEntry(entries);
		}
	}

	//if a user has column arrays defined already,they can add them to an empty matrix with this method
	public void fillColumns(double []... ds) {
		if (ds.length!=this.getColumns() || ds[0].length!=this.getRows()) {
			System.out.println("Dimension mismatch! Please check your rows and the matrix and try again!");
		} else {
			double [] entries = new double[ds.length*ds[0].length];
			//adding entries into one big array to use addEntry
			int count=0;
			for (int i=0;i<ds[0].length;i++) {
				for (int j=0;j<ds.length;j++) {
					entries[count]=ds[j][i];
					count++;
				}
			}
			this.addEntry(entries);
		}
	}


		//gets the desired entry in the matrix
	public double getEntry(int row, int column) {
		//error case ---dimension mismatch
		if (row<=this.matrix.length && column<=this.matrix[0].length) {
		return this.matrix[row-1][column-1];
	} else {
			System.out.println("Dimension mismatch. Please call the function again with valid parameters.");
			return -123;
		}


	}

	public void setEntry(int row,int column, double entry) {
		this.matrix[row-1][column-1]=entry;
	}
		// gets the desired row in the matrix
	public double[] getRow(int row) {
		return this.matrix[row-1];
	}

		//gets the column as an array
	public double[] getColumn(int column) {
		//error case---- desired column doesnt exist
		double[] thiscolumn=new double[this.matrix.length];
		if (column <= this.matrix[0].length && column >0) {

		for (int i=0;i<this.matrix.length;i++) {
			//first column is 0, last column is length-1, etc
			thiscolumn[i]=this.matrix[i][column-1];
		}
		return thiscolumn;
		} else {
			System.out.println("Please call the function again with a valid column number.");
			return thiscolumn;
		}
	}
		//returns a copy of the matrix the method is called on
	public myMatrix getCopy() {
		myMatrix copy = new myMatrix(this.getRows(),this.getColumns());
		double [] entries = new double[this.getRows()*this.getColumns()];
		int count=0;
		for (int i=1;i<=this.getRows();i++) {
			for (int j=1;j<=this.getColumns();j++) {
				entries[count]=this.getEntry(i,j);
				count++;
			}
		}
		copy.addEntry(entries);
		return copy;
	}

		//adds the first array parameter to the second and replaces the second
		//multiplier is the first parameter for the first array
	public void addRow(double multiplier, int row1, int row2) {
		//we need to multiply each entry in a duplicate of row1 by multiplier and add to row 2

		//rowtoreplace references the row2 row in the matrix object
		double [] rowtoreplace= this.getRow(row2);

		//replacing row desired by sum of first row multiplied by some multiplier

			for (int i=0;i<this.matrix[0].length;i++) {
				// we use i-1 for column in getEntry because that is how the method is written above
				rowtoreplace[i]=(this.getEntry(row1, i+1)*multiplier)+rowtoreplace[i];
			}



	}

	public void multiplyRow(double multiplier, int row1) {
		double [] rowtochange = this.getRow(row1);
		for (int i=0;i<this.getColumns();i++) {
			rowtochange[i]=rowtochange[i]*multiplier;
		}
	}

		//since we now have row addition, we should also add a method that swaps rows

	public void swapRow(int row1,int row2) {
		double [] placeholder = this.getRow(row1);
		this.matrix[row1-1]=this.getRow(row2);
		this.matrix[row2-1]=placeholder;
	}

	public String toString() {
		if (this.matrix.length==0 || this.matrix[0].length==0) {
			System.out.println("The matrix is empty!!. Please add values!");
			return null;
		} else {
			String stringvalue ="";
			for (int i=0;i<this.matrix.length;i++) {
				for (int j=0;j<this.matrix[0].length;j++) {
						if (j!=this.matrix[0].length-1) {
							stringvalue+=this.matrix[i][j]+", ";
						} else {
							stringvalue+=this.matrix[i][j];
						}
				}
				//line break between matrix rows
				stringvalue+="\n";
			}
			return stringvalue;
		}
	}


	public myMatrix addMatrix( myMatrix othermatrix ) {
		if (othermatrix.getRows()==this.getRows() && othermatrix.getColumns()==this.getColumns()) {
			double[] result = new double[othermatrix.getRows()*othermatrix.getColumns()];
				//adding sum to the result
			int count=0;
			for (int i=1;i<=othermatrix.getRows();i++) {
				for (int j=1;j<=othermatrix.getColumns();j++) {
					result[count]=othermatrix.getEntry(i, j)+this.getEntry(i, j);
					count++;
				}
			}
			myMatrix results = new myMatrix(othermatrix.getRows(), othermatrix.getColumns());
			results.addEntry(result);
			return results;
		} else {

			return null;
		}
	}

		//returns the result of traditional matrix multiplication;
	public myMatrix multiplyMatrix(myMatrix othermatrix) {
		if (this.getColumns()==othermatrix.getRows()) {
			//resulting matrix has as many rows as the first and as many columns as the second
			myMatrix result = new myMatrix(this.getRows(), othermatrix.getColumns());
			double[] entries = new double[this.getRows()*othermatrix.getColumns()];
			double colSum=0;
			int count=0;

			for (int i=1;i<=this.getRows();i++) {
				for (int z=1;z<=othermatrix.getColumns();z++) {
					for (int j=1;j<=this.getColumns();j++) {
						//using fact that column of first is the row of the second
						colSum+=(this.getEntry(i,j)*othermatrix.getEntry(j,z));

					}
					entries[count]=colSum;
					colSum=0;
					count++;
				}
			}

			result.addEntry(entries);
			return result;
		} else {
			System.out.println("This matrix does not have the same number of columns as the rows of the other matrix! (Dimension mismatch)");
			return null;
		}
	}

		//method scales the whole matrix by some scaling factor;
	public void scaleMatrix(double scale) {
		if (this.matrix.length==0 || this.matrix[0].length==0 || this==null) {
			System.out.println("The matrix is empty!");
		} else {
			for (int i=0;i<this.getRows();i++) {
				for (int j=0;j<this.getColumns();j++) {
					this.matrix[i][j]=scale*this.matrix[i][j];
				}
			}
		}

	}
		//method that returns a matrix of 1st entry - 2nd entry;
	public myMatrix subtractMatrix(myMatrix otherMatrix) {
		if (otherMatrix.getRows()==this.getRows() && otherMatrix.getColumns()==this.getColumns()) {
			otherMatrix.scaleMatrix(-1);
			myMatrix result =this.addMatrix(otherMatrix);

			return result;


		} else {
			System.out.println("Dimensions do not match!");
			return null;
		}
	}
		//method to return the transpose of a matrix
	public myMatrix getTranspose() {
		myMatrix transpose = new myMatrix(this.getColumns(),this.getRows());
		int count=0;
		double[] transposed = new double[this.getColumns()*this.getRows()];
		for (int j=1;j<=this.getColumns();j++) {
			for (int i=1;i<=this.getRows();i++) {
				transposed[count]=this.getEntry(i,j);
				count++;
			}
		}
		transpose.addEntry(transposed);
		return transpose;
	}


	public myMatrix matrixExponential(double exponential) {
		if (this.matrix.length==0 || this.matrix[0].length==0 || this==null) {
			System.out.println("The matrix is empty!");
			return null;
		} else {
			double [] result = new double[this.getColumns()*this.getRows()];
			int count=0;
			for (int i=1;i<=this.getRows();i++) {
				for (int j=1;j<=this.getColumns();j++) {

					result[count]=Math.pow(this.getEntry(i,j), exponential);
					count++;
				}
			}
				myMatrix finalresult=new myMatrix(this.getRows(),this.getColumns());
				finalresult.addEntry(result);
				return finalresult;
			}
		}

		//helper function
	public boolean restOfColumnZero(int rowstart, int column) {
		boolean result=true;
		for (int j=rowstart;j<=this.getRows();j++) {
			if (this.getEntry(j, column)!=0) {
				result=false;
				break;
			}
		}
		return result;
	}




	//now we have all the methods we need , so we can make complex matrix operations

	//helper function for rowEchelon
	//moves zero rows to bottom in a specific column
	public void moveZerosToBottom(int column) {
		//once we hit a zero that we put towards the bottom, we can break
		int breakrow=-1;
		for (int i=1;i<=this.getRows();i++) {
			if (this.getEntry(i,column)==0) {
				if (i==breakrow) {
					break;
				}
				for (int j=this.getRows();j>i;j--) {
					if (this.getEntry(j,column)!=0) {
						//swap zero entry row with nonzero entry row closest to the bottom
						this.swapRow(i,j);
						breakrow=j;
						break;
					}
				}
			}
		}
	}

	//tries to account for rounding error
	public void checkZeros() {
		for (int i=1;i<=this.getRows();i++) {
			for (int j=1;j<=this.getColumns();j++) {
				if (Math.abs(this.getEntry(i,j))<=0.000001) {
					this.setEntry(i,j,0.0);
				}
			}
		}
	}

	//function to return rowEchelon form
	//KEEP IN MIND THAT THERE ARE SEVERAL ROW ECHELON FORMS FOR A SINGLE MATRIX (UNLIKE REDUCED ROW ECHELON FORM)
	public myMatrix rowEchelon(){
		//special cases
		if (this==null || this.getRows()==0 || this.getColumns()==0) {
			System.out.println("Please check your matrix and try again!");
			return null;
		}

		if (this.getRows()==1) {
			return this;
		}

		if (this.getColumns()==1) {
			myMatrix copy = this.getCopy();
			copy.moveZerosToBottom(1);
			return copy;
		}
		//creating copy of given matrix to perform operations on
		myMatrix copy=this.getCopy();
		//keeping track of possible pivots
		int pivots =1;
		for (int j=1;j<=copy.getColumns();j++) {
			//first move zeros to bottom
			copy.moveZerosToBottom(j);

			int counter=0;

			for (int i=pivots;i<=copy.getRows();i++) {
				//need to iterate through row and find number nonzero
				if (copy.getEntry(i,j)!=0 ) {
					counter++;
					//making leading entry 1

					copy.multiplyRow(1/copy.getEntry(i,j),i);

				}
			}
				//operation based on number of nonzero entries in row
				if (counter==1) {
					pivots++;
				} else if (counter>1) {
					//need to not count the pivot row in this case
					counter--;
					for (int i=pivots+1;i<=pivots+counter;i++) {
						copy.addRow(-1, pivots,i);
						//accounting for rounding error (assuming the user input valid first small entries to start)
						copy.checkZeros();



					}
					pivots++;

				}


		}
		return copy;
	}

		//there is only ONE reduced row echelon form for a matrix
		public myMatrix reducedRowEchelon() {


			//need to use the fact that the first nonzero entry in the row is a pivot
			myMatrix copy = this.rowEchelon();
			for (int i=1;i<=copy.getRows();i++) {
				int pivot=0;
				int column=0;
				for (int j=1;j<=copy.getColumns();j++) {
					if (copy.getEntry(i,j)!=0 && pivot==0) {
						//recording row of pivot
						pivot = i;
						column=j;
						break;
					}
				}
					//subtracting rows above the pivot if they have nonzero entries in that column
					if (pivot!=1 && pivot!=0) {
						for (int z =1;z<pivot;z++) {
							if (copy.getEntry(z,column)!=0) {
								//pivot entry should be equal to 1
								copy.addRow(-copy.getEntry(z,column),pivot,z);
								copy.checkZeros();

							}
						}

					}
				}


			return copy;
		}


		//method that solves the associated augemented matrix
		//parameter is the augmented part of the matrix
		//this matrix must return some string represenation
		//this is because things like free variables must be represented
	public ArrayList<String> solveSystem(myMatrix augmented) {

		if (this==null || augmented==null || this.getRows()!=augmented.getRows()) {
			System.out.println("Please check that the matrix is correct and try again.");
			return null;
		} else {

			//copying columns of variable coefficients
			myMatrix result = new myMatrix(augmented.getRows(),1);
			double [] entries = new double[augmented.getRows()];


				for (int j=1;j<=augmented.getRows();j++) {
					entries[j-1]=augmented.getEntry(j,1);
				}
			result.addEntry(entries);


			int currentpivots=0;
			ArrayList<Integer> zerorow = new ArrayList<Integer>();




			for (int i=1;i<=this.getColumns();i++) {

				for (int j=currentpivots+1;j<=this.getRows();j++) {

					if (this.getEntry(j, i)==0 && !this.restOfColumnZero(j, i)) {

						zerorow.add(j);
						for (int c=j;j<=this.getRows();c++) {
							if (this.getEntry(c, i)!=0) {
								result.multiplyRow(1/this.getEntry(c,i),c);
								result.swapRow(j,c);
								this.multiplyRow(1/this.getEntry(c, i), c);
								this.swapRow(j, c);

								break;
							}
						}

					} else if (this.getEntry(j, i)==0) {
						break;
					} else {
						result.multiplyRow(1/this.getEntry(j,i),j);
						this.multiplyRow(1/this.getEntry(j, i), j);
						if (currentpivots!=i) {
							currentpivots++;
						}
					}


				}
					//this part subtracts the pivot row from the nonzero rows
					//that are below itself
					//we do not count the end because then we will transform
					//a zero row into a nonzero one
				for (int j=currentpivots;j<=this.getRows();j++) {
					if (this.getEntry(j, i)!=0 && j!=currentpivots) {
					result.addRow(-1,currentpivots,j);
					this.addRow(-1, currentpivots, j);
					}
				}
				zerorow.clear();
			}



			//we are checking if each entry in the row is not zero
			//if so, then we divide that row by itself and move on to the next row
			for (int i=1;i<=this.getRows();i++) {
				for (int j=1;j<=this.getColumns();j++) {
					if (this.getEntry(i, j)!=0) {
						result.multiplyRow(1/this.getEntry(i,j),i);
						this.multiplyRow(1/this.getEntry(i, j), i);
						break;
					}
				}
			}



			//finally converting to row echelon

			//we do not have to check the bottom row because if we get to that row
			//without detecting any pivots, then the other rows must be zero rows
				for (int i=1;i<=this.getRows();i++) {
					//first for loop to iterate across row to find pivot
					for (int j=1;j<=this.getColumns();j++) {
						if (this.getEntry(i, j)!=0) {


							for (int z = i-1;z>=1;z--) {
								if (this.getEntry(z, j)!=0) {
									result.addRow(-this.getEntry(z,j),i,z);
									this.addRow(-this.getEntry(z, j), i, z);

								}
							}

							break;
						}
					}



				}
					//checking if row is zero and augemented entry is nonzero
					//in that case, the system is inconsistent
				for (int i=1;i<=this.getRows();i++) {
					double rowSum=0;
					for (int j=1;j<=this.getColumns();j++) {
						rowSum+=this.getEntry(i,j);
					}
					if (rowSum==0 && result.getEntry(i,1)!=0) {
						System.out.println("The system is inconsistent!");
						return null;
					}
				}

					//putting the answer vectors in an array list

					//creating possible possible variables
				ArrayList<String> variables = new ArrayList<String>();
				for (int i=1;i<=this.getColumns();i++) {
					variables.add("x"+i);
				}
					//creating the answers to add to the ArrayList
				ArrayList<String> answers = new ArrayList<String>();
				for (int i=0;i<this.getRows();i++) {
					String answer="";
					int oneCount=0;
					for (int j=0;j<this.getColumns();j++) {
						if (this.getEntry(i+1,j+1)!=0 && oneCount==0) {
							answer=""+variables.get(j)+":"+result.getEntry(i+1,1);
							oneCount++;
						} else if (this.getEntry(i+1,j+1)!=0) {
							if (this.getEntry(i+1,j+1)>0) {
								if (this.getEntry(i+1,j+1)==1) {
									answer+="-"+variables.get(j);
								} else {
								answer+="-"+this.getEntry(i+1,j+1)+variables.get(j);
								}
							} else {
								if (this.getEntry(i+1,j+1)==-1) {
									answer+="+"+variables.get(j);
								} else {
								answer+=""+(-this.getEntry(i+1,j+1))+variables.get(j);
								}
							}
						}
					}
					if (!answers.equals("")) {
					answers.add(answer);
					}
				}
				return answers;
			}
		}


		//method that returns the determinant of the given matrix
	public double getDeterminant() {

		//check to see if matrix isnt square && not null
		if (this==null || this.getRows()!=this.getColumns()) {
			System.out.println("Please enter a valid matrix");
			return -1;
		} else if (this.getRows()==1) {
				//returns the single entry in the matrix
				//special case
				return this.getEntry(1,1);
		} else if (this.getRows()==2 && this.getColumns()==2) {
			//base definition of the determinant
			return (this.getEntry(1,1)*this.getEntry(2,2))-(this.getEntry(1,2)*this.getEntry(2,1));
		} else {
				double det=0;

				myMatrix toDo;
			//need a recursive definition now for bigger matrices
			for (int i=1;i<=this.getColumns();i++) {
				//if the entry is zero, the term added to the
				//determinant sum is just zero, so it doesnt matter
				if (this.getEntry(1,i)!=0) {

					toDo = new myMatrix(this.getRows()-1,this.getColumns()-1);

					double [] entries = new double[(this.getRows()-1)*(this.getColumns()-1)];

						//adding values to new matrix
					int count=0;
						for (int j=2;j<=this.getRows();j++) {
							for (int z=1;z<=this.getColumns();z++) {
								if(z!=i) {
								entries[count]=this.getEntry(j,z);
								count++;
								}
							}
						}
						toDo.addEntry(entries);

						if (i%2!=0) {
							det+=this.getEntry(1,i)*toDo.getDeterminant();
						} else {
							det-=this.getEntry(1,i)*toDo.getDeterminant();
					}

				}
			}
				return det;

		}
	}
		//method that returns the matrix of minors with checkerboard negative pattern
	public myMatrix getCoFactorMatrix() {

		myMatrix cofactors=new myMatrix(this.getRows(),this.getColumns());
		double[] entries = new double[this.getRows()*this.getColumns()];
		int count=0;
		//i need to find the minors of each entry
		//then i need to add it to this matrix
		//and then I need to make according entries positive and negative
		for (int i=1;i<=this.getRows();i++) {
			for (int j=1;j<=this.getColumns();j++) {
				//i need to form the minor matrix here and call my
				//getDeterminant() function on it
				myMatrix minor = new myMatrix(this.getRows()-1,this.getColumns()-1);
				double[] placer = new double[(this.getRows()-1)*(this.getColumns()-1)];
				int index=0;
				for (int c=1;c<=this.getRows();c++) {
					for (int d=1;d<=this.getColumns();d++) {
						if (c!=i && d!=j) {
							placer[index]=this.getEntry(c,d);
							index++;
						}
					}
				}
				minor.addEntry(placer);

				entries[count]=minor.getDeterminant();
				count++;

			}
		}
		cofactors.addEntry(entries);
		//now i have to add negatives in a "checkerboard" fashion
		for (int i=1;i<=cofactors.getRows();i++) {
			for (int j=1;j<=cofactors.getColumns();j++) {
				if ((i+j)%2!=0) {
					cofactors.matrix[i-1][j-1]=-cofactors.matrix[i-1][j-1];
				}
			}
		}
		return cofactors;

	}
		//method that returns the adjoint matrix-cofactor matrix with elements "reflected" over the diagonal
	public myMatrix getAdjointMatrix() {
		myMatrix adjoint = this.getCoFactorMatrix();
		//need to swap entries over the diagonal to form adjoint
		for (int i=1;i<=this.getRows();i++) {
			for (int j=1;j<=this.getColumns();j++) {
				if (j>i) {
					double temp =adjoint.matrix[i-1][j-1];
					adjoint.matrix[i-1][j-1]=adjoint.matrix[j-1][i-1];
					adjoint.matrix[j-1][i-1]=temp;
				}
			}
		}
		return adjoint;
	}


	public myMatrix getInverse() {
		//if determinant is 0, the inverse is not defined
		//otherwise, it is the adjoint matrix scaled by 1/determinant
		if (this.getDeterminant()==0 || this.getColumns()!=this.getRows()) {
			System.out.println("The given matrix does not have an inverse!");
			return null;
		} else {
				myMatrix inverse = this.getAdjointMatrix();
				//scaling adjoint by 1/determinant
				inverse.scaleMatrix(1.0/this.getDeterminant());

				return inverse;
		}

	}


		//returns the projection of this vector onto the other vector
	public myMatrix getProjection(myMatrix otherVector) {
		//assuming the user actually inputs vectors
		//need to transform given arguments for matrix multiplication
		//I am using the fact that a dot product between vectors is just matrix multiplication of a row vector and column vector
		myMatrix firstVector;
		myMatrix secondVector;
		//TRANSFORMING FIRST VECTOR INTO A ROW VECTOR
		if (this.getRows()>1) {
		
			firstVector = this.getTranspose();
			} else {
				//making a copy of "this" so that we dont alter "this"
			double[] entries = new double[this.getColumns()];
			for (int i=1;i<=this.getColumns();i++) {
				entries[i-1]=this.getEntry(1,i);
			}
			firstVector=new myMatrix(1,this.getColumns());
			firstVector.addEntry(entries);
		}

		//TRANSFORMING SECOND VECTOR INTO A COLUMN VECTOR

		if (otherVector.getColumns()>1) {
			secondVector=otherVector.getTranspose();
		} else {
			//making a copy of otherVector
			double[] entries = new double[otherVector.getRows()];
			for (int i=1;i<=otherVector.getRows();i++) {
				entries[i-1]=otherVector.getEntry(i,1);
			}
			secondVector=new myMatrix(otherVector.getRows(),1);
			secondVector.addEntry(entries);
		}


		if (firstVector.getColumns()!=secondVector.getRows()) {
			System.out.println("Please enter valid vectors!");
			return null;
		} else {
			//general formula is ((v dot s)/(s dot s))*s (scalar times a vector)
			//s dot s is just every entry multiplied by itself and added together
			//in our case s is already a column vector
			double doubleDot=0;
				for (int i=1;i<=secondVector.getRows();i++) {
					double toAdd=secondVector.getEntry(i,1);
					doubleDot+=(toAdd*toAdd);
				}
			//my multiplication function will result in 1x1 matrix, so i need to "extract" that entry
			myMatrix topPart = firstVector.multiplyMatrix(secondVector);
			double top = topPart.getEntry(1,1);
			//need to scale the column vector according to the equation
			secondVector.scaleMatrix(top/doubleDot);
			return secondVector;
		}
	}

	//method that normalizes the columns of a matrix (or vectors)
	public void normalize() {
		//need to iterate through column, then divide each entry by the square root of the sum of each entry squared (definition of norm)
		for (int j=1;j<=this.getColumns();j++) {
			double colSum=0;
			for (int i=1;i<=this.getRows();i++) {
				colSum+=Math.pow(this.getEntry(i,j),2);
			}
			//now colSum is the norm of the vector squared so we need to take the square root to get norm
			colSum=Math.pow(colSum,0.5);
			//now we are dividing each entry
			for (int i=1;i<=this.getRows();i++) {
				this.setEntry(i,j,this.getEntry(i,j)/colSum);
			}
		}
	}
	//method returns an orthogonal basis based on gram schmidt algorithm of the columns of the matrix
	public myMatrix gramSchmidt() {
		//user enters matrix where the columns are vectors
		//<u1,u2,u3,...> is the entered matrix
		//we first say v1=u1
		//then v2=u2-u2projv1
		//v3= u3-u3projv2-u3projv1 and so on
		myMatrix result=null;
		ArrayList<myMatrix> orthogonal = new ArrayList<myMatrix>();
		for (int i=1;i<=this.getColumns();i++) {
			myMatrix vector = new myMatrix(this.getRows(),1);
			vector.addEntry(this.getColumn(i));
				if (i==1) {
					//making a copy of vector
					result = new myMatrix(this.getRows(),1);
					result.addEntry(this.getColumn(i));
						//normalizing the vector

					orthogonal.add(result);


				} else {
						//need to access every vector in orthogonal so far
						result=vector;
					for (int j=0;j<orthogonal.size();j++) {
						result=result.subtractMatrix(vector.getProjection(orthogonal.get(j)));
					}

					orthogonal.add(result);


				}
		}
		myMatrix realResult = new myMatrix(this.getRows(),this.getColumns());
		//adding vectors as columns in the new matrix
		double [] realEntries = new double[this.getRows()*this.getColumns()];
		int count=0;
			for (int i=1;i<=this.getRows();i++) {
				for (int j=0;j<orthogonal.size();j++) {
					realEntries [count]=orthogonal.get(j).getEntry(i,1);
					count++;
				}
			}
		realResult.addEntry(realEntries);
		realResult.checkZeros();
		return realResult;
	}
		//this method returns an array list with QR decomposition in that order as elements of the list
		//my method only works on square matrices
	public ArrayList<myMatrix> QR() {
		if (this.getColumns()!=this.getRows() || this==null) {
			System.out.println("Please check that the matrix is square and try again!");
			return null;
		}
		myMatrix Q = this.gramSchmidt();
		Q.normalize();
		//making the upper triangular matrix
		myMatrix R = new myMatrix(this.getRows(),this.getColumns());
		double[] entries = new double[this.getRows()*this.getColumns()];
		int count=0;
		for (int i=1;i<=this.getRows();i++) {
			double []	normalColumn = Q.getColumn(i);
			for (int j=1;j<=this.getColumns();j++) {
				if (i<=j) {
					double[] thisColumn = this.getColumn(j);
					//dot product
					double dot =0;
					for (int z=0;z<normalColumn.length;z++) {
						dot+=(normalColumn[z])*(thisColumn[z]);
					}
					entries[count]=dot;
					count++;
				} else {
					entries[count]=0;
					count++;
				}
			}
		}
		R.addEntry(entries);
		ArrayList<myMatrix> QR = new ArrayList<myMatrix>();
		QR.add(Q);
		QR.add(R);
		return QR;
	}
		//helper function for checking if matrix is upper triangular (zeros below diagonal)
	public boolean isUpperTriangular() {
		myMatrix copy=this.getCopy();
		copy.checkZeros();
		for (int i=1;i<=this.getRows();i++) {
			for (int j=1;j<=this.getColumns();j++) {
				if (i>j) {
					if (copy.getEntry(i,j)!=0) {
						return false;
					}
				}
			}
		}
		return true;
	}
		//checks if matrix is in blockdiagonal form (blocks or 1x1 entries along the diagonal)
		//this will help us extract complex eigenvalues from QR factorization
		//only valid for square matrices
	public boolean isBlockDiagonal() {
		myMatrix copy = this.getCopy();
		copy.checkZeros();
		for (int i=1;i<=copy.getRows();i++) {
			for (int j=1;j<=copy.getColumns();j++) {
				if (i==j) {
					if (i!=copy.getRows()) {
						//if its not the last entry in the diagonal
						if (copy.getEntry(i+1,j)!=0) {
							//checking for 2x2 block matrix
							if (i+1!=copy.getRows()) {
								if (! restOfColumnZero(i+2,j)) {
									return false;
								}
							}
						} else {
							//checking for rest of column to be zero
							if (!restOfColumnZero(i+1,j)) {
								return false;
							}
						}
					}
				}
			}
		}
		return true;
	}

	//method that returns the non-trivial eigenvalues of the given matrix as a string array
	public ArrayList<String> eigenValues() {

	//uses QR algorithm in order to numerically retrieve eigenValues
	//keeps iterating A = RQ and then finding the QR factorization of That
	if (this.getRows()!=this.getColumns() || this==null) {
		System.out.println("Either the matrix is empty, or it is not square!");
		return null;
	}
	//copying the matrix that the function is called on
	myMatrix copy = this.getCopy();
	//now copy is a duplicate object of THIS

	//Supports Complex and Real eigenvalues for inputted matrices that are REAL
	ArrayList<String> eigenvalues = new ArrayList<String>();
	//logic for 1x1 matrix: the eigenvalue is just the single entry
		if (copy.getRows()==1 && copy.getColumns()==1) {
			eigenvalues.add(Double.toString(copy.getEntry(1,1)));
			return eigenvalues;
		}
	//logic for 2x2 matrix: need to use quadratic formula
	if (copy.getRows()==2 && copy.getColumns()==2) {
		double trace = -(copy.getEntry(1,1)+copy.getEntry(2,2));
		double det = copy.getDeterminant();
		//need to see if discriminant is 0, <0, or >0
		if (Math.pow(trace,2)-4*det>0) {
			//we have two distinct real roots
			double eigen_one= (-trace+Math.pow(Math.pow(trace,2)-4*det,0.5))/(2);
			double eigen_two=(-trace-Math.pow(Math.pow(trace,2)-4*det,0.5))/(2);
			eigenvalues.add(Double.toString(eigen_one));
			eigenvalues.add(Double.toString(eigen_two));
			return eigenvalues;
		} else if (Math.pow(trace,2)-4*det<0) {
			//we have two complex roots
			//i need a string for i
			ArrayList<String> complex = new ArrayList<String>();
			String eigen_one=""+(-trace/2)+"+"+(Math.pow(4*det-Math.pow(trace,2),0.5)/2)+"i";
			String eigen_two=""+(-trace/2)+"-"+(Math.pow(4*det-Math.pow(trace,2),0.5)/2)+"i";
			complex.add(eigen_one);
			complex.add(eigen_two);
			return complex;

		} else {
			//we have a double root
			double eigen = -trace/2;
			eigenvalues.add(Double.toString(eigen));
			eigenvalues.add(Double.toString(eigen));
			return eigenvalues;
		}
	}


	int itcount=0;

		while (!copy.isUpperTriangular() || !copy.isBlockDiagonal()) {
			//getting the QR factorization

			ArrayList<myMatrix> factorization= copy.QR();
			copy = factorization.get(1).multiplyMatrix(factorization.get(0));
			itcount++;

				if (itcount==500) {
					System.out.println("Too many iterations!");
					break;
				}

		}
		ArrayList<String> eigs = new ArrayList<String>();
		System.out.println(copy);
//		myMatrix finalR = copy.QR().get(1);
//		System.out.println(finalR);
		copy.checkZeros();
	//	System.out.println(finalR);
		for (int i=1;i<=copy.getRows();i++) {
			for (int j=1;j<=copy.getColumns();j++) {
				if (i==j) {
					if (i!=copy.getRows()) {
						if (copy.getEntry(i+1,j)!=0) {
							//we have a block
							myMatrix block = new myMatrix(2);
							double[] complex = new double[4];
							complex[0]=copy.getEntry(i,j);
							complex[1]=copy.getEntry(i,j+1);
							complex[2]=copy.getEntry(i+1,j);
							complex[3]=copy.getEntry(i+1,j+1);
							block.addEntry(complex);
							//recursive definition for 2x2 blocks for which we know the eigenvalues!
							ArrayList<String> temp=block.eigenValues();
							for (String value: temp) {
								eigs.add(value);
							}
							i++;
							j++;
							//after the for loop increments, we will be at the second furthest diagonal position
						} else {
							//we have a 1x1 matrix
							eigs.add(Double.toString(copy.getEntry(i,j)));
						}
					} else {
						//if the last entry in the diagonal is not part of a block, then we add the value
						if (copy.getEntry(i,j-1)==0) {
							eigs.add(Double.toString(copy.getEntry(i,j)));
						}
					}
				}
			}
		}
		return eigs;



	}


}
