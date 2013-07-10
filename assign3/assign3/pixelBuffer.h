#ifndef _PIXELBUFFER_H_
#define _PIXELBUFFER_H_

#include "vector.h"

class pixelBuffer
{
	public:
		pixelBuffer();
		pixelBuffer(int WIDTH, int HEIGHT);

		void setPixel(int x, int y, double r, double g, double b);
		vector getPixel(int x, int y);

	private:
		int WIDTH;
		int HEIGHT;
		
		vector*** array_ptr;
};

#endif // _PIXELBUFFER_H_