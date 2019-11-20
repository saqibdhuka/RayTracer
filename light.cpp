#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Angel.h"


Light :: Light(){
	point.x = 0.0;
	point.y = 0.0;
	point.z = 0.0;

	colorIntensity.x = 0.0;
	colorIntensity.y = 0.0;
	colorIntensity.z = 0.0;

}


Light :: Light(vec3 p, vec3 c){
	point.x = p.x;
	point.y = p.y;
	point.z = p.z;

	colorIntensity.x = c.x;
	colorIntensity.y = c.y;
	colorIntensity.z = c.z;

}