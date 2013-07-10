/*
CSCI 480
Assignment 3 Raytracer

Name: Jonathan Chu
*/

#include <pic.h>
#include <windows.h>
#include <stdlib.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <stdio.h>
#include <string>

#include "pixelBuffer.h"
#include "ray.h"

#define MAX_TRIANGLES 10
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

#define MINIPIXELS 1

//the field of view of the camera
#define fov 60.0
#define M_PI 3.14159265358979323846

unsigned char buffer[HEIGHT][WIDTH][3];
pixelBuffer* myImage;

vector* eyePos;

vector* topLeftCornerPos;
vector* topRightCornerPos;
vector* bottomLeftCornerPos;
vector* bottomRightCornerPos;

double screenWidth;
double screenHeight; 

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void trace(ray* currentRay, int x, int y);

bool intersectsLight(ray* currentRay, int lightIndex, double & time);
bool intersectsSphere(ray* currentRay, vector* normal, int sphereIndex, double & time, bool & insideSphere, bool forShadows);
bool intersectsTriangle(ray* currentRay, vector* triV0, vector* triV1, vector* triV2, vector* normal, double & time, double & alpha, double & beta, double & gamma);

void calculateCorners();

void plot_pixel_display(int x,int y,double r,double g,double b);
void plot_pixel_jpeg(int x,int y,double r,double g,double b);
void plot_pixel(int x,int y,double r,double g,double b);
void save_jpg();

//MODIFY THIS FUNCTION
void draw_scene()
{
	double curY = topLeftCornerPos->getY();

	for(int i = 0; i < HEIGHT * MINIPIXELS; i++)
	{
		double curX = topLeftCornerPos->getX();
		for(int j = 0; j < WIDTH * MINIPIXELS; j++)
		{
			//DO RAY THINGS
			vector* currentPos = new vector(curX,curY,-1.0);
			vector* direction = eyePos->subtract(currentPos);
			direction->normalize();

			ray* curRay = new ray(direction->getX(),direction->getY(),direction->getZ());

			trace(curRay,j,i);

			delete curRay;
			delete currentPos;
			delete direction;

			curX += (screenWidth/(WIDTH * MINIPIXELS));
		}
		curY -= (screenHeight/(HEIGHT * MINIPIXELS));
	}

	int x = 0;
	int y = 0;

	for(int k = 0; k < HEIGHT * MINIPIXELS; k+= MINIPIXELS)
	{
		glPointSize(2.0);
		glBegin(GL_POINTS);
		x = 0;
		for(int l = 0; l < WIDTH * MINIPIXELS; l+= MINIPIXELS)
		{
			double rValue = 0.0;
			double gValue = 0.0;
			double bValue = 0.0;

			for(int m = 0; m < MINIPIXELS; m++)
			{
				for(int n = 0; n < MINIPIXELS; n++)
				{
					rValue += myImage->getPixel(l+m,k+n).getX();
					gValue += myImage->getPixel(l+m,k+n).getY();
					bValue += myImage->getPixel(l+m,k+n).getZ();
				}
			}

			rValue /= (MINIPIXELS*MINIPIXELS);
			gValue /= (MINIPIXELS*MINIPIXELS);
			bValue /= (MINIPIXELS*MINIPIXELS);

			//plot_pixel(l,((HEIGHT - 1) - k),myImage->getPixel(l,k).getX(),myImage->getPixel(l,k).getY(),myImage->getPixel(l,k).getZ());
			plot_pixel(x,((HEIGHT - 1) - y),rValue,gValue,bValue);
			buffer[y][x][0] = rValue;
			buffer[y][x][1] = gValue;
			buffer[y][x][2] = bValue;
			x++;
		}
		y++;
		glEnd();
		glFlush();
	}
	if(mode == MODE_JPEG)
	{
		save_jpg();
	}
}

void trace(ray* currentRay, int x, int y)
{
	//First Case See if it intersects light source
	double time = 10000.0;

	for(int i = 0; i < num_lights; i++)
	{
		double timeToLight = 0.0;
 		if(intersectsLight(currentRay,i, timeToLight))
		{
			if(timeToLight < time)
			{
				time = timeToLight;
				//myImage->setPixel(x,y,
				//	myImage->getPixel(x,y).getX() + (lights[i].color[0] * 255.0),
				//	myImage->getPixel(x,y).getY() + (lights[i].color[1] * 255.0),
				//	myImage->getPixel(x,y).getZ() + (lights[i].color[2] * 255.0));

				myImage->setPixel(x,y,
					(lights[i].color[0] * 255.0),
					(lights[i].color[1] * 255.0),
					(lights[i].color[2] * 255.0));
			}
		}
	}
	for(int j = 0; j < num_spheres; j++)
	{
		double timeToSphere = 0.0;
		bool insideSphere = false;

		vector* normal = new vector;
		//State shows what kind of intersection
		if(intersectsSphere(currentRay, normal, j, timeToSphere, insideSphere, false))
		{
			if(timeToSphere < time)
			{
				bool clear = true;

				time = timeToSphere;
				
				vector* intersectionPoint = currentRay->getPositionAtTime(timeToSphere);

				vector* viewVector = &(currentRay->getDirection());
				viewVector->negate();
				viewVector->normalize();
				normal->normalize();

				vector* ambientComponent = new vector(ambient_light[0],ambient_light[1],ambient_light[2]); 
				vector* diffuseComponent = new vector();
				vector* specularComponent = new vector();

				double dR, dG, dB, sR, sG, sB, shininess;

				dR = (spheres[j].color_diffuse[0]);

				dG = (spheres[j].color_diffuse[1]);

				dB = (spheres[j].color_diffuse[2]);

				sR = (spheres[j].color_specular[0]);

				sG = (spheres[j].color_specular[1]);

				sB = (spheres[j].color_specular[2]);

				shininess = (spheres[j].shininess);

				vector* diffuseColor = new vector(dR,dG,dB);
				vector* specularColor = new vector(sR,sG,sB);

				for(int i = 0; i < num_lights; i++)
				{
					vector* lightPos = new vector(lights[i].position[0],
						lights[i].position[1], lights[i].position[2]);
					
					vector* lightVector = intersectionPoint->subtract(lightPos);
					lightVector->normalize();
					vector* rVector = lightVector->reflect(normal);
					rVector->normalize();

					ray* toLight = new ray(intersectionPoint->getX(),intersectionPoint->getY(),intersectionPoint->getZ(),
						lightVector->getX(),lightVector->getY(),lightVector->getZ()); 
					vector* tempNormal = new vector();

					double lightTime = 0.0;

					intersectsLight(toLight, i, lightTime);

					bool tempState = false;

					for(int i = 0; i < num_spheres; i++)
					{
						double tempTime = 0.0;
						if(intersectsSphere(toLight, tempNormal, i, tempTime, tempState, true))
						{
							if((tempTime < lightTime) && tempTime > 0.001)
							{
								clear = false;
							}
						}
					}
					for(int j = 0; j < num_triangles; j++)
					{
						vector* triV0 = new vector(triangles[j].v[0].position[0], triangles[j].v[0].position[1], triangles[j].v[0].position[2]);
						vector* triV1 = new vector(triangles[j].v[1].position[0], triangles[j].v[1].position[1], triangles[j].v[1].position[2]);
						vector* triV2 = new vector(triangles[j].v[2].position[0], triangles[j].v[2].position[1], triangles[j].v[2].position[2]);

						double tempAlpha, tempBeta, tempGamma;

						double tempTime = 0.0;

						if(intersectsTriangle(toLight, triV0, triV1, triV2, tempNormal, tempTime, tempAlpha, tempBeta, tempGamma))
						{
							if(tempTime < lightTime)
							{
								clear = false;
							}
						}

						delete triV0;
						delete triV1;
						delete triV2;
					}

					//Diffuse Component
					double diffuseDotProduct = lightVector->dot(normal);

					if(diffuseDotProduct < 0.0)
					{
						diffuseDotProduct = 0.0;
					}
					//if(diffuseDotProduct > 1.0)
					//{
					//	diffuseDotProduct = 1.0;
					//}

					diffuseColor->scalarMultiply(diffuseDotProduct);
					
					diffuseComponent->addTo(diffuseColor);

					//Specular Component
					double specularDotProduct = rVector->dot(viewVector);

					if(specularDotProduct < 0.0)
					{
						specularDotProduct = 0.0;
					}
					//if(specularDotProduct > 1.0)
					//{
					//	specularDotProduct = 1.0;
					//}

					specularColor->scalarMultiply(pow(specularDotProduct,shininess));

					specularComponent->addTo(specularColor);

					delete lightPos;
					delete lightVector;
					delete rVector;
					delete toLight;
				}

				if(clear)
				{
					ambientComponent->saturate();
					diffuseComponent->saturate();
					specularComponent->saturate();

					vector* phongModel = ambientComponent->add(diffuseComponent);
					phongModel->addTo(specularComponent);

					phongModel->saturate();

					myImage->setPixel(x,y,phongModel->getX()*255.0,phongModel->getY()*255.0,phongModel->getZ()*255.0);
					
					delete phongModel;
				}
				else
				{
					myImage->setPixel(x,y,0.0,0.0,0.0);
				}
					
				//delete viewVector;
				delete diffuseColor;
				delete specularColor;
				delete specularComponent;
				delete diffuseComponent;
				delete ambientComponent;
			} //END IF HERE
		}
		delete normal;
	}
	for(int k = 0; k < num_triangles; k++)
	{
		vector* triV0 = new vector(triangles[k].v[0].position[0], triangles[k].v[0].position[1], triangles[k].v[0].position[2]);
		vector* triV1 = new vector(triangles[k].v[1].position[0], triangles[k].v[1].position[1], triangles[k].v[1].position[2]);
		vector* triV2 = new vector(triangles[k].v[2].position[0], triangles[k].v[2].position[1], triangles[k].v[2].position[2]);
		vector* normal = new vector();

		double alpha = 0.0;
		double beta = 0.0;
		double gamma = 0.0;
		double timeToTriangle = 0.0;

		if(intersectsTriangle(currentRay,triV0,triV1,triV2,normal,timeToTriangle,alpha,beta,gamma))
		{
			if(timeToTriangle < time)
			{
				bool clear = true;
				time = timeToTriangle;
				//Do Intersection things
				normal->normalize();
				vector* intersectionPoint = currentRay->getPositionAtTime(timeToTriangle);
				vector* viewVector = &(currentRay->getDirection());
				viewVector->negate();
				viewVector->normalize();

				vector* ambientComponent = new vector(ambient_light[0],ambient_light[1],ambient_light[2]); 
				vector* diffuseComponent = new vector();
				vector* specularComponent = new vector();

				double dR, dG, dB, sR, sG, sB, shininess;

				dR = (triangles[k].v[0].color_diffuse[0] * alpha) +
					(triangles[k].v[1].color_diffuse[0] * beta) +
					(triangles[k].v[2].color_diffuse[0] * gamma);

				dG = (triangles[k].v[0].color_diffuse[1] * alpha) +
					(triangles[k].v[1].color_diffuse[1] * beta) +
					(triangles[k].v[2].color_diffuse[1] * gamma);

				dB = (triangles[k].v[0].color_diffuse[2] * alpha) +
					(triangles[k].v[1].color_diffuse[2] * beta) +
					(triangles[k].v[2].color_diffuse[2] * gamma);

				sR = (triangles[k].v[0].color_specular[0] * alpha) +
					(triangles[k].v[1].color_specular[0] * beta) +
					(triangles[k].v[2].color_specular[0] * gamma);

				sG = (triangles[k].v[0].color_specular[1] * alpha) +
					(triangles[k].v[1].color_specular[1] * beta) +
					(triangles[k].v[2].color_specular[1] * gamma);

				sB = (triangles[k].v[0].color_specular[2] * alpha) +
					(triangles[k].v[1].color_specular[2] * beta) +
					(triangles[k].v[2].color_specular[2] * gamma);

				shininess = (triangles[k].v[0].shininess * alpha) +
					(triangles[k].v[1].shininess * beta) +
					(triangles[k].v[2].shininess * gamma);

				vector* diffuseColor = new vector(dR,dG,dB);
				vector* specularColor = new vector(sR,sG,sB);

				for(int i = 0; i < num_lights; i++)
				{
					vector* lightPos = new vector(lights[i].position[0],
						lights[i].position[1], lights[i].position[2]);
					
					vector* lightVector = intersectionPoint->subtract(lightPos);
					lightVector->normalize();
					vector* rVector = lightVector->reflect(normal);
					rVector->normalize();

					ray* toLight = new ray(intersectionPoint->getX(),intersectionPoint->getY(),intersectionPoint->getZ(),
						lightVector->getX(),lightVector->getY(),lightVector->getZ()); 
					vector* tempNormal = new vector();

					double lightTime = 0.0;

					intersectsLight(toLight, i, lightTime);

					double tempTime = 0.0;
					bool tempState;

					for(int i = 0; i < num_spheres; i++)
					{
						tempTime = 0.0;
						if(intersectsSphere(toLight, tempNormal, i, tempTime, tempState, false))
						{
							if(tempTime < lightTime)
							{
								clear = false;
							}
						}
					}

					tempNormal = new vector();

					for(int j = 0; j < num_triangles; j++)
					{
						tempTime = 0.0;
						vector* triVZ = new vector(triangles[j].v[0].position[0], triangles[j].v[0].position[1], triangles[j].v[0].position[2]);
						vector* triVO = new vector(triangles[j].v[1].position[0], triangles[j].v[1].position[1], triangles[j].v[1].position[2]);
						vector* triVT = new vector(triangles[j].v[2].position[0], triangles[j].v[2].position[1], triangles[j].v[2].position[2]);

						double tempAlpha, tempBeta, tempGamma;

						if(k != j && intersectsTriangle(toLight, triVZ, triVO, triVT, tempNormal, tempTime, tempAlpha, tempBeta, tempGamma))
						{
							if(tempTime < lightTime)
							{
								clear = false;
							}
						}

						delete triVZ;
						delete triVO;
						delete triVT;
					}

					//Diffuse Component
					double diffuseDotProduct = lightVector->dot(normal);

					if(diffuseDotProduct < 0.0)
					{
						diffuseDotProduct = 0.0;
					}
					//if(diffuseDotProduct > 1.0)
					//{
					//	diffuseDotProduct = 1.0;
					//}

					diffuseColor->scalarMultiply(diffuseDotProduct);
					
					diffuseComponent->addTo(diffuseColor);

					//Specular Component
					double specularDotProduct = rVector->dot(viewVector);

					if(specularDotProduct < 0.0)
					{
						specularDotProduct = 0.0;
					}
					//if(specularDotProduct > 1.0)
					//{
					//	specularDotProduct = 1.0;
					//}

					specularColor->scalarMultiply(pow(specularDotProduct,shininess));

					specularComponent->addTo(specularColor);

					delete lightPos;
					delete lightVector;
					delete rVector;
				}

				if(clear)
				{
					ambientComponent->saturate();
					diffuseComponent->saturate();
					specularComponent->saturate();

					vector* phongModel = ambientComponent->add(diffuseComponent);
					phongModel->addTo(specularComponent);

					phongModel->saturate();

					myImage->setPixel(x,y,phongModel->getX()*255.0,phongModel->getY()*255.0,phongModel->getZ()*255.0);
					
					delete phongModel;
				}
				else
				{
					myImage->setPixel(x,y,0.0,0.0,0.0);
				}

				delete intersectionPoint;
				//delete viewVector;
				delete diffuseColor;
				delete specularColor;
				delete specularComponent;
				delete diffuseComponent;
				delete ambientComponent;
			}
		}

		delete triV0;
		delete triV1;
		delete triV2;
		delete normal;
	}
}

bool intersectsLight(ray* currentRay, int lightIndex, double & time)
{
	double tX, tY, tZ;

	tX = (lights[lightIndex].position[0] - currentRay->getStartPosition().getX())/currentRay->getDirection().getX();
	tY = (lights[lightIndex].position[1] - currentRay->getStartPosition().getY())/currentRay->getDirection().getY();
	tZ = (lights[lightIndex].position[2] - currentRay->getStartPosition().getZ())/currentRay->getDirection().getZ();

	if(tX == tY && tY == tZ && tX == tZ)
	{
		time = tX;
		return true;
	}
	else
	{
		time = tX;
		return false;
	}
}

bool intersectsSphere(ray* currentRay, vector* normal, int sphereIndex, double & time, bool & insideSphere, bool forShadows)
{
	double a, b, c;

	a = 1.0;
	b = 2.0 * ((currentRay->getDirection().getX() * (currentRay->getStartPosition().getX() - spheres[sphereIndex].position[0])) +
		(currentRay->getDirection().getY() * (currentRay->getStartPosition().getY() - spheres[sphereIndex].position[1])) +
		(currentRay->getDirection().getZ() * (currentRay->getStartPosition().getZ() - spheres[sphereIndex].position[2]))); 
	//b = 2.0 * (((currentRay->getDirection().getX() - currentRay->getStartPosition().getX()) * (currentRay->getStartPosition().getX() - spheres[sphereIndex].position[0])) +
	//	((currentRay->getDirection().getY() - currentRay->getStartPosition().getY()) * (currentRay->getStartPosition().getY() - spheres[sphereIndex].position[1])) +
	//	((currentRay->getDirection().getZ() - currentRay->getStartPosition().getZ()) * (currentRay->getStartPosition().getZ() - spheres[sphereIndex].position[2])));
	c = ((pow(currentRay->getStartPosition().getX() - spheres[sphereIndex].position[0], 2) + 
		pow(currentRay->getStartPosition().getY() - spheres[sphereIndex].position[1], 2) +
		pow(currentRay->getStartPosition().getZ() - spheres[sphereIndex].position[2], 2)) - pow(spheres[sphereIndex].radius, 2));

	double tZero, tOne;

	tZero = ((-1.0 * b) + sqrt(pow(b, 2) - (4.0 * c)))/(2.0);
	tOne = ((-1.0 * b) - sqrt(pow(b, 2) - (4.0 * c)))/(2.0);

	if(tZero > 0.0 && tOne > 0.0)
	{
		if(!forShadows)
		{
			time = min(tZero,tOne);
		}
		else
		{
			time = max(tZero,tOne);
		}

		vector* intersectionPoint = currentRay->getPositionAtTime(time);
		double inverseRadius = (1.0/spheres[sphereIndex].radius);

		normal->setX(intersectionPoint->getX() - spheres[sphereIndex].position[0]);
		normal->setY(intersectionPoint->getY() - spheres[sphereIndex].position[1]);
		normal->setZ(intersectionPoint->getZ() - spheres[sphereIndex].position[2]);
		normal->scalarMultiply(inverseRadius);
		normal->normalize();

		insideSphere = false;
		delete intersectionPoint;
		return true;
		//return 1.0;
	}
	else
	{
		if(tZero > 0.0 ^ tOne > 0.0)
		{
			time = max(tZero,tOne);

			vector* intersectionPoint = currentRay->getPositionAtTime(time);
			double inverseRadius = (1.0/spheres[sphereIndex].radius);

			normal->setX(intersectionPoint->getX() - spheres[sphereIndex].position[0]);
			normal->setY(intersectionPoint->getY() - spheres[sphereIndex].position[1]);
			normal->setZ(intersectionPoint->getZ() - spheres[sphereIndex].position[2]);
			normal->scalarMultiply(inverseRadius);
			normal->negate();
			normal->normalize();

			insideSphere = true;
			delete intersectionPoint;
			return true;
		}
		else
		{
			time = -1.0;
			insideSphere = false;
			return false;
		}
	}
}

bool intersectsTriangle(ray* currentRay, vector* triVZero, vector* triVOne, vector* triVTwo, vector* triNormal, double & time, double & alpha, double & beta, double & gamma)
{
	//Algorithm for ray-triangle intersection detailed here: http://www.cs.washington.edu/education/courses/cse457/07sp/lectures/triangle_intersection.pdf

	vector* edge1 = triVZero->subtract(triVOne);
	vector* edge2 = triVZero->subtract(triVTwo);

	vector* normal = edge1->cross(edge2);
	normal->normalize();

	triNormal->setX(normal->getX());
	triNormal->setY(normal->getY());
	triNormal->setZ(normal->getZ());

	double determinant = normal->dot(&(currentRay->getDirection()));

	if(determinant == 0.0)
	{
		//Parallel, does not intersect
		delete edge1;
		delete edge2;
		delete normal;

		return false;
	}
	else
	{
		double dVal = triVZero->dot(normal);
		double nDotP = normal->dot(&(currentRay->getStartPosition()));
		time = (dVal - nDotP)/determinant;

		vector* intersectionPoint = currentRay->getPositionAtTime(time);

		vector* AtoB = triVZero->subtract(triVOne);
		vector* BtoC = triVOne->subtract(triVTwo);
		vector* CtoA = triVTwo->subtract(triVZero);
		vector* AtoC = triVZero->subtract(triVTwo); //Used for Barycentric

		vector* AreaVec = AtoB->cross(AtoC);

		vector* AtoIntersect = triVZero->subtract(intersectionPoint);
		vector* BtoIntersect = triVOne->subtract(intersectionPoint);
		vector* CtoIntersect = triVTwo->subtract(intersectionPoint);

		vector* ABCross = AtoB->cross(AtoIntersect);
		vector* BCCross = BtoC->cross(BtoIntersect);
		vector* CACross = CtoA->cross(CtoIntersect);

		double ABval, BCval, CAval;
		ABval = ABCross->dot(normal);
		BCval = BCCross->dot(normal);
		CAval = CACross->dot(normal);

		double alphaNum, betaNum, gammaNum, area;
		alphaNum = BCCross->dot(normal);
		betaNum = CACross->dot(normal);
		gammaNum = ABCross->dot(normal);
		area = AreaVec->dot(normal);

		alpha = alphaNum/area;
		beta = betaNum/area;
		//gamma = gammaNum/area;
		gamma = 1.0 - (alpha + beta);

		delete edge1;
		delete edge2;
		delete normal;
		delete intersectionPoint;
		delete AtoB;
		delete BtoC;
		delete CtoA;
		delete AtoC;
		delete AreaVec;
		delete AtoIntersect;
		delete BtoIntersect;
		delete CtoIntersect;
		delete ABCross;
		delete BCCross;
		delete CACross;

		if(ABval >= 0.0 && BCval >= 0.0 && CAval >= 0.0 && time >= 0.0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}

void calculateCorners()
{
	double aspectRatio = (double)(WIDTH)/(double)(HEIGHT);
	double x, y, z;
	z = -1;

	x = aspectRatio * (tan((fov/2) * M_PI/180));
	y = tan((fov/2) * M_PI/180);

	topLeftCornerPos->setVector((-1.0)*x,y,z);
	topRightCornerPos->setVector(x,y,z);
	bottomLeftCornerPos->setVector((-1.0)*x,(-1.0)*y,z);
	bottomRightCornerPos->setVector(x,(-1.0)*y,z);

	screenWidth = 2*x;
	screenHeight = 2*y;
}

void plot_pixel_display(int x,int y,double r,double g,double b)
{
	glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
	glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,double r,double g,double b)
{
  buffer[HEIGHT-y-1][x][0]=((double)r)/256.f;
  buffer[HEIGHT-y-1][x][1]=((double)g)/256.f;
  buffer[HEIGHT-y-1][x][2]=((double)b)/256.f;
}

void plot_pixel(int x,int y,double r,double g,double b)
{
	while(r >= 256.0)
	{
		r -= 256.0;
	}
	while(g >= 256.0)
	{
		g -= 256.0;
	}
	while(b >= 256.0)
	{
		b -= 256.0;
	}
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

void parse_check(char *expected,char *found)
{
  if(stricmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(stricmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(stricmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(stricmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void display()
{

}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);

  myImage = new pixelBuffer(WIDTH * MINIPIXELS,HEIGHT * MINIPIXELS);
  eyePos = new vector(0.0,0.0,0.0);
  topLeftCornerPos = new vector();
  topRightCornerPos = new vector();
  bottomLeftCornerPos = new vector();
  bottomRightCornerPos = new vector();
  calculateCorners();
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
      draw_scene();
      if(mode == MODE_JPEG)
	save_jpg();
    }
  once=1;
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      filename = argv[2];
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
