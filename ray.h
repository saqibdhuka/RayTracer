#ifndef RAY_H
#define RAY_H
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Angel.h"



class Ray
{

	public:

		
		vec3 origin;
		vec3 direction;

		Ray();
		Ray(vec3 o, vec3 d);
		
	
};

#endif