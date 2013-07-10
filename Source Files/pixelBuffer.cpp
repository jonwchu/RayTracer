#include "pixelBuffer.h"

pixelBuffer::pixelBuffer()
{
	this->WIDTH = 640;
	this->HEIGHT = 480;

	array_ptr = new vector** [HEIGHT];
	for(int i = 0; i < HEIGHT; i++)
	{
		array_ptr[i] = new vector*[WIDTH];
	}

	//initialize all values to 0

	for(int j = 0; j < HEIGHT; j++)
	{
		for (int k = 0; k < WIDTH; k++)
		{
			array_ptr[j][k] = new vector(255.0,255.0,255.0);
		}
	}
}
pixelBuffer::pixelBuffer(int WIDTH, int HEIGHT)
{
	this->WIDTH = WIDTH;
	this->HEIGHT = HEIGHT;

	array_ptr = new vector** [HEIGHT];
	for(int i = 0; i < HEIGHT; i++)
	{
		array_ptr[i] = new vector*[WIDTH];
	}

	//initialize all values to 0

	for(int j = 0; j < HEIGHT; j++)
	{
		for (int k = 0; k < WIDTH; k++)
		{
			array_ptr[j][k] = new vector(255.0,255.0,255.0);
		}
	}
}
void pixelBuffer::setPixel(int x, int y, double r, double g, double b)
{
	this->array_ptr[y][x]->setVector(r,g,b);
}
vector pixelBuffer::getPixel(int x, int y)
{
	return *this->array_ptr[y][x];
}