#include "vector.h"

 vector::vector()
 {
	 this->x = 0.0;
	 this->y = 0.0;
	 this->z = 0.0;
 }

vector::vector(double x, double y, double z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

vector::~vector()
{
}

double vector::dot(vector* v)
{
	double dotProduct = 0.0;

	dotProduct += (getX() * v->getX());
	dotProduct += (getY() * v->getY());
	dotProduct += (getZ() * v->getZ());

	return dotProduct;
}

vector* vector::add(vector* v)
{
	vector* sumVector = new vector(getX() + v->getX(),getY() + v->getY(), getZ() + v->getZ());
	return sumVector;
}

vector* vector::subtract(vector* v)
{
	vector* diffVector = new vector(v->getX() - getX(),v->getY() - getY(),v->getZ() - getZ());
	return diffVector;
}

vector* vector::cross(vector* v)
{
	vector* crossProduct = new vector(
		(getY()*v->getZ()) - (getZ()*v->getY()),
		(getZ()*v->getX()) - (getX()*v->getZ()),
		(getX()*v->getY()) - (getY()*v->getX()));
	return crossProduct;
}

vector* vector::reflect(vector* v)
{
	double twoIdotV = 0.0;
	twoIdotV = 2.0 * (this->dot(v));
	vector* toReturn = new vector(
		(twoIdotV * v->getX()) - this->getX(),
		(twoIdotV * v->getY()) - this->getY(),
		(twoIdotV * v->getZ()) - this->getZ());
	return toReturn;
}

double vector::magnitude(void)
{
	return sqrt(pow(getX(),2) + pow(getY(),2) + pow(getZ(),2));
}

void vector::normalize(void)
{
	double mag = magnitude();
	setX((getX()/mag));
	setY((getY()/mag));
	double tempZ = getZ();
	setZ((getZ()/mag));
	double tempZ2 = getZ();
}

void vector::addTo(vector* v)
{
	setX(getX() + v->getX());
	setY(getY() + v->getY());
	setZ(getZ() + v->getZ());
}

void vector::scalarMultiply(double value)
{
	setX(getX() * value);
	setY(getY() * value);
	setZ(getZ() * value);
}

void vector::negate()
{
	setX(getX() * (-1.0));
	setY(getY() * (-1.0));
	setZ(getZ() * (-1.0));
}

void vector::saturate()
{
	if(getX() < 0.0)
	{
		setX(0.0);
	}
	if(getX() > 1.0)
	{
		setX(1.0);
	}

	if(getY() < 0.0)
	{
		setY(0.0);
	}
	if(getY() > 1.0)
	{
		setY(1.0);
	}

	if(getZ() < 0.0)
	{
		setZ(0.0);
	}
	if(getZ() > 1.0)
	{
		setZ(1.0);
	}
}

void vector::setX(double value)
{
	this->x = value;
}

void vector::setY(double value)
{
	this->y = value;
}

void vector::setZ(double value)
{
	this->z = value;
}

void vector::setVector(double xValue, double yValue, double zValue)
{
	setX(xValue);
	setY(yValue);
	setZ(zValue);
}

double vector::getX()
{
	return this->x;
}

double vector::getY()
{
	return this->y;
}

double vector::getZ()
{
	return this->z;
}