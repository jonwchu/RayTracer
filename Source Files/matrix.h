#ifndef _MATRIX_H_
#define _MATRIX_H_

class matrix
{
	public:
		//Constructor
		matrix();
		matrix(int rows, int cols);

		~matrix();

		//Matrix Multiply returns this * m
		//Assures checked that they are valid to multiply
		matrix* multiply(matrix* m);

		//Matrix add returns this + m
		//Assures checked that they are valid to add
		matrix* add(matrix* m);

		bool validToAdd(matrix* m);
		bool validToMultiply(matrix* m);

		void setCell(int row, int col, double value);

		int getRows();
		int getCols();
		double getCell(int row, int col);

	private:
		int ROWS;
		int COLS;

		double** array_ptr;
};

#endif // _MATRIX_H_