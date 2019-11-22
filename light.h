#ifndef LIGHT_H
#define LIGHT_H

#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Angel.h"
#include "ray.h"

class Light{

public:

	Ray ray;
	vec3 colorIntensity;

	Light();
	Light(Ray rayObj, vec3 c);
	
};




#endif