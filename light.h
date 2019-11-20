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

class Light{

public:

	vec3 point;
	vec3 colorIntensity;

	Light();
	Light(vec3 p, vec3 c);
	
	
	

};




#endif