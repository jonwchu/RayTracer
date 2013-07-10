#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <cmath>

class vector
{
	public:
		//Constructor
		vector();
		vector(double x, double y, double z);

		~vector();

		//Dot Product returns this dot v
		double dot(vector* v);

		//Vector add returns this + v
		vector* add(vector* v);

		//Vector subtract v-this (this to v)
		vector* subtract(vector* v);

		//Cross Product of this x v
		vector* cross(vector* v);

		//Vector of reflection about v
		vector* reflect(vector* v);

		//Get Magnitud
		double magnitude(void);

		//Normalize vector
		void normalize(void);
		
		void addTo(vector* v);
		void scalarMultiply(double value);
		void negate();
		void saturate();
		void setX(double value);
		void setY(double value);
		void setZ(double value);
		void setVector(double xValue, double yValue, double zValue);

		double getX();
		double getY();
		double getZ();

	private:
		double x;
		double y;
		double z;
};

#endif // _VECTOR_H_