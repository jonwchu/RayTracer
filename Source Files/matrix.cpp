#include "matrix.h"

matrix::matrix()
{
	ROWS = 4;
	COLS = 4;

	array_ptr = new double* [ROWS];
	for(int i = 0; i < ROWS; i++)
	{
		array_ptr[i] = new double[COLS];
	}

	//initialize all values to 0

	for(int j = 0; j < ROWS; j++)
	{
		for (int k = 0; k < COLS; k++)
		{
			array_ptr[j][k] = 0.0;
		}
	}
}

matrix::matrix(int rows, int cols)
{
	ROWS = rows;
	COLS = cols;

	array_ptr = new double* [ROWS];
	for(int i = 0; i < ROWS; i++)
	{
		array_ptr[i] = new double[COLS];
	}

	for(int j = 0; j < ROWS; j++)
	{
		for (int k = 0; k < COLS; k++)
		{
			array_ptr[j][k] = 0.0;
		}
	}
}

matrix::~matrix()
{
	for(int i = 0; i < ROWS; i++)
	{
		delete[] array_ptr[i]; 
	}
	delete[] array_ptr;
	array_ptr = 0;
}

matrix* matrix::multiply(matrix* m)
{
	matrix* toReturn = new matrix(ROWS,m->getCols());

	for(int i = 0; i < toReturn->getRows(); i++)
	{
		for(int j = 0; j < toReturn->getCols(); j++)
		{
			toReturn->setCell(i,j,0.0f);

			for(int k = 0; k < getCols(); k++)
			{
				double toAdd = toReturn->getCell(i,j) + (getCell(i,k) * m->getCell(k,j));
				toReturn->setCell(i,j,toAdd);
			}
		}
	}
	return toReturn;
}

matrix* matrix::add(matrix* m)
{
	matrix* toReturn = new matrix(ROWS,COLS);

	for(int j = 0; j < toReturn->getRows(); j++)
	{
		for (int k = 0; k < toReturn->getCols(); k++)
		{
			array_ptr[j][k] = m->getCell(j,k) + getCell(j,k);
		}
	}

	return toReturn;
}

bool matrix::validToAdd(matrix* m)
{
	if((m->getCols() == COLS) && (m->getRows() == ROWS))
	{
		return true;
	}
	return false;
}

bool matrix::validToMultiply(matrix* m)
{
	if(COLS == m->getRows())
	{
		return true;
	}
	return false;
}

void matrix::setCell(int row, int col, double value)
{
	array_ptr[row][col] = value;
}

int matrix::getRows()
{
	return ROWS;
}

int matrix::getCols()
{
	return COLS;
}

double matrix::getCell(int row, int col)
{
	return array_ptr[row][col];
}