#ifndef PLANE__H
#define PLANE__H

#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ray.h"
#include "Angel.h"

#ifndef MATERIALSTRUCTURE_H
#define MATERIALSTRUCTURE_H

#include "materialStructure.h"

#endif


class Plane{

public:

	// struct Material{
	// 	vec3 kdColor;
	// 	vec3 ksColor;	
	// 	float specular;
	// 	float reflection;
	// 	float refraction;
	// };

	float A;
	float B;
	float C;
	float D;
	Material matPlane;

	Plane();
	Plane(float a, float b, float c, float d, Material newMat);
	vec3 plane_normal();
	float intersectPlane(Ray &ray);
	
	
	

};




#endif