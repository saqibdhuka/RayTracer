#ifndef LIGHT_H
#define LIGHT_H

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

	vec3 origin;
	vec3 colorIntensity;

	Light();
	Light(vec3 position, vec3 c);

};




#endif
