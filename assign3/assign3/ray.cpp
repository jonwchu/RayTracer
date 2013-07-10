#include "ray.h"

ray::ray()
{
	this->startPos = new vector(0,0,0);
	this->direction = new vector(0,0,0);
}

ray::ray(double xPos, double yPos, double zPos,
	double xDir, double yDir, double zDir)
{
	this->startPos = new vector(xPos,yPos,zPos);
	this->direction = new vector(xDir,yDir,zDir);
}

ray::ray(double xDir, double yDir, double zDir)
{
	this->startPos = new vector(0,0,0);
	this->direction = new vector(xDir,yDir,zDir);
}

vector ray::getStartPosition()
{
	return *this->startPos;
}

vector ray::getDirection()
{
	return *this->direction;
}

vector* ray::getPositionAtTime(double time)
{
	vector* toReturn = new vector(this->getStartPosition().getX() + (this->getDirection().getX() * time),
		this->getStartPosition().getY() + (this->getDirection().getY() * time),
		this->getStartPosition().getZ() + (this->getDirection().getZ() * time));
	return toReturn;
}

void ray::setStartPosition(double x, double y, double z)
{
	this->startPos->setVector(x,y,z);
}

void ray::setStartPosition(vector vecPos)
{
	this->startPos->setVector(vecPos.getX(),vecPos.getY(),vecPos.getZ());
}

void ray::setDirection(double x, double y, double z)
{
	this->direction->setVector(x,y,z);
}

void ray::setDirection(vector vecDir)
{
	this->direction->setVector(vecDir.getX(),vecDir.getY(),vecDir.getZ());
}