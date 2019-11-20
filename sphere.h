#ifndef SPHERE_H
#define SPHERE_H


#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Angel.h"

class Sphere
{
public:


	struct Material{
		vec3 kdColor;
		vec3 ksColor;	
		float specular;
		float reflection;
		float refraction;
	};
	
	vec3 center;
	float radius;
	Material mat;
	

	Sphere();
	Sphere(vec3 point, float rad, Material newMat);
	void draw();
};


#endif