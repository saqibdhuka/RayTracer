#include "ray.h"
#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Angel.h"

Ray :: Ray(){
	origin.x = 0.0;
	origin.y = 0.0;
	origin.z = 0.0;

	direction.x = 0.0;
	direction.y = 0.0;
	direction.z = -1.0;
}
Ray :: Ray(vec3 o, vec3 d){
	origin.x = o.x;
	origin.y = o.y;
	origin.z = o.z;

	direction.x = d.x;
	direction.y = d.y;
	direction.z = d.z;
}

void Ray::shootRay(const float &t){


	glBegin(GL_LINES);
		glVertex3f(origin.x, origin.y, origin.z);
		glVertex3f(origin.x + (t * direction.x),
				   origin.y + (t * direction.y),
				   origin.z + (t * direction.z));
	glEnd();



}
