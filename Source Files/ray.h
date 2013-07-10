#ifndef _RAY_H_
#define _RAY_H_

#include "vector.h"

class ray
{
	public:
		ray();
		ray(double xPos, double yPos, double zPos,
			double xDir, double yDir, double zDir);
		ray(double xDir, double yDir, double zDir);

		vector getStartPosition();
		vector getDirection();

		vector* getPositionAtTime(double time);

		void setStartPosition(double x, double y, double z);
		void setStartPosition(vector vecPos);
		void setDirection(double x, double y, double z);
		void setDirection(vector vecDir);

	private:
		vector* startPos;
		vector* direction;
};

#endif // _RAY_H_