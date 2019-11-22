#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Angel.h"
#include "light.h"


Light :: Light(){
	ray.origin.x = 0.0;
	ray.origin.y = 0.0;
	ray.origin.z = 0.0;

	ray.direction.x = 0.0;
	ray.direction.y = 0.0;
	ray.direction.z = -1.0;

	colorIntensity.x = 0.0;
	colorIntensity.y = 0.0;
	colorIntensity.z = 0.0;

}


Light :: Light(Ray rayObj, vec3 c){
	ray = rayObj;
	colorIntensity.x = c.x;
	colorIntensity.y = c.y;
	colorIntensity.z = c.z;

}

