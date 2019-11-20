#include "ray.h"
#include <GL/glew.h>
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
