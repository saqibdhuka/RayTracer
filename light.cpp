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
	origin.x = 0.0;
	origin.y = 0.0;
	origin.z = 0.0;

	colorIntensity.x = 2.5;
	colorIntensity.y = 2.5;
	colorIntensity.z = 2.5;

}


Light :: Light(vec3 position, vec3 c){
	origin.x = position.x;
	origin.y = position.y;
	origin.z = position.z;

	colorIntensity.x = c.x;
	colorIntensity.y = c.y;
	colorIntensity.z = c.z;

}
