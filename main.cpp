#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <time.h>
#include "Angel.h"
#include "ray.h"
#include "plane.h"
#include "sphere.h"
#include "light.h"

#ifndef MATERIALSTRUCTURE_H
#define MATERIALSTRUCTURE_H

#include "materialStructure.h"

#endif



#define WINDOW 400
#define FOVY 45
#define zNear 0.01f
#define zFar 50.0f




float t = 100.0; //t for the use of Ray
int rayDepth=1;
vec3 backColor(0.0, 0.0, 0.0);

const float LIMIT = tan(FOVY/2);



Sphere sph_objects[100];
Plane pla_objects[100];

Light light_object[100];

// vec3 shadingSphere(Ray ray, Sphere sphere){

// 	vec3 color(1.0, 0.0, 0.0);

// 	if(sphere.insersect(ray, t)){
// 		vec3 V = ray.direction;
// 		vec3 P = ray.origin + V * t;
// 		vec3 N = sphere.get_normal(P);

// 		float ratio = dot(N,V);

// 		color = sphere.getkdColor() * (ratio *0.5);
// 	}

// 	return color;


// }

vec3 vecPow(vec3 v, float power){

	vec3 ans(v.x, v.y, v.z);

	for(int i =2; i < power; i++){
		ans *= v;
	}
	// ans.x = pow((float)v.x, power);
	// ans.y = pow((float)v.y, power);
	// ans.z = pow((float)v.z, power);
	// printf("answer.x:%f, answer.y:%f, answer.z:%f \n", ans.x, ans.y, ans.z);
	return ans;
	// return vec3(pow(v.x, power), pow(v.y, power), pow(v.z, power));
}

vec3 trace(Ray ray, Sphere sphere, Light light){ //Work on this

	vec3 I; //Light Intesity to calculate
	if(rayDepth = 0){
		return backColor;
	}else{
		// printf("Inside rayDepth > 0\n\n");
		// printf("Ia.x:%f, Ia.y:%f, Ia.z:%f \n", Ia.x, Ia.y, Ia.z);
		vec3 kd = sphere.mat.kdColor;
		// printf("kd.x:%f, kd.y:%f, kd.z:%f \n", kd.x, kd.y, kd.z);
		vec3 Ii = light.colorIntensity;
		vec3 Ia =  kd * Ii;
		vec3 ks = sphere.mat.ksColor;
		float  kr = sphere.mat.reflection;
		float kt = sphere.mat.refraction;
		int q = sphere.mat.specular;
		float S =1;
		vec3 intersectionPoint = ray.origin + (t * ray.direction);
		vec3 V = normalize(intersectionPoint);
		vec3 N = normalize(sphere.get_normal(intersectionPoint));
		vec3 Li = normalize(intersectionPoint - light_object[0].origin);
		vec3 R = normalize(2*(N * Li)*N - Li);

		// printf("R.x:%f, R.y:%f, R.z:%f \n", R.x, R.y, R.z);

		vec3 rv = R*V;
		// printf("rv.x:%f, rv.y:%f, rv.z:%f \n", rv.x, rv.y, rv.z);

		vec3 tmp = vecPow(rv, q);
		// printf("tmp.x:%f, tmp.y:%f, tmp.z:%f \n", tmp.x, tmp.y, tmp.z);
		
		vec3 ik = vec3(Ia.x * kd.x, Ia.y*kd.y, Ia.z*kd.z);
		// printf("ik.x:%f, ik.y:%f, ik.z:%f \n", ik.x, ik.y, ik.z);
		vec3 ktemp = ks * tmp;
		// printf("ktemp.x:%f, ktemp.y:%f, ktemp.z:%f \n", ktemp.x, ktemp.y, ktemp.z);
		vec3  knl = kd*N*Li;
		// printf("knl.x:%f, knl.y:%f, knl.z:%f \n", knl.x, knl.y, knl.z);
		vec3 sum = S*Ii*(knl + ktemp);
		I = ik + sum;
		// printf("I.x:%f, I.y:%f, I.z:%f \n", I.x, I.y, I.z);

		return (I);
	}

}


// vec3 shadingPlane(Ray ray, Plane plane){

// 	vec3 color(0.5, 0.5, 0.5);

// 	if(plane.insersect(ray, t)){
// 		vec3 V = ray.direction;
// 		vec3 P = ray.origin + V * t;
// 		vec3 N = plane.get_normal(P);

// 		float ratio = dot(N,V);

// 		color = sphere.getkdColor() * (ratio *0.5);
// 	}

// 	return color;


// }

// int count =0;
void render(){

	for(int x=-WINDOW/2; x < WINDOW/2; x++){
		float xN=x * LIMIT* LIMIT*20;
		// printf("xN: %f\n", xN); 
		for (int y = -WINDOW/2; y < WINDOW/2; y++)
		{
			float yN=y * LIMIT* LIMIT*20;			// glBegin(GL_POINTS);
			// 	glVertex3f(xN, yN, -1.0);
			// glEnd();

			vec3 orig(0.0, 0.0, 0.0);
			vec3 dir(xN, yN, -1.0);
			Ray ray(orig, dir);

			// glColor3f(0.0,0.0,0.0);
			// 
			// glColor3f(1.0,0.0,0.0);
			// glBegin(GL_POINTS);
			// 	glVertex3f(xN, yN, -50);
			// glEnd();
			// ray.shootRay(t);
			t = sph_objects[0].intersect(ray); 
			// count++;
			// printf("t %f\n", t);
			if(t!=-1){
				// glColor3f(sph_objects[0].mat.ksColor.x, sph_objects[0].mat.ksColor.y, sph_objects[0].mat.ksColor.z);
				// glColor3f(1.0, 0.5, 1.0);
				vec3 colorLight = trace(ray, sph_objects[0], light_object[0]);
				printf("t is R: %f, G: %f, B: %f\n", colorLight.x, colorLight.y, colorLight.z);
				glColor3f(colorLight.x, colorLight.y, colorLight.z);
				// glEnable(GL_CULL_FACE);
				// glCullFace(GL_FRONT);
				
				glBegin(GL_POINTS);
					glVertex3f(ray.origin.x + t*ray.direction.x , 
						ray.origin.y + t*ray.direction.y , 
						ray.origin.z + t*ray.direction.z );
				glEnd();
			}
			// else{
			// 	glColor3f(0.0,0.0,0.0);
			// 	glBegin(GL_POINTS);
			// 		glVertex3f(xN, yN, -1.0);
			// 	glEnd();
			// }
			
			// if(t != -1){
				
			// }
			// vec3 kd = shadingSphere(ray, sph_objects[0]);
			// sph_objects[0].mat.kdColor = kd;
			// sph_objects[0].draw();
		}
	}
	// printf("count: %d\n", sph_objects[0].count);
}

void init(){
	
    glViewport(0, 0, WINDOW, WINDOW);

    // glClearColor(backColor.x, backColor.y, backColor.z, 0.0);
	glEnable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective((float) FOVY, 1.0f, zNear, zFar);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0);
	// glOrtho(-100.0, 100.0, -100.0,100.0,-5.0,20.0);
	// glLoadIdentity();
}


void callbackDisplay() {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH);
	glMatrixMode(GL_MODELVIEW);

	// glPushMatrix();
		
	// glPopMatrix();

	vec3 lightPos(0.1, 0.1, 0.0);
	vec3 lightColor(50.0, 50.0, 50.0);
	Light light(lightPos, lightColor);
	light_object[0] = light;

	// glTranslatef(0.0, 0.0, -3.0);
	vec3 center_ball1(-1.0, 1.0, -25.0);
	float radius_ball1 = 3;
	vec3 kd_ball1(0.1, 0.8, 0.8);//KD
	vec3 ks_ball1(0.9, 0.9, 0.9);//KS
	float spec = 32;
	float refl = 0.5;
	float refr = 0.0;
	//Setting up the material
	Material newMat;
	newMat.kdColor = kd_ball1;
	newMat.ksColor = ks_ball1;
	newMat.specular = spec;
	newMat.reflection = refl;
	newMat.refraction = refr;
	
	Sphere ball1(center_ball1, radius_ball1, newMat);
	sph_objects[0] = ball1;
	
	render();

	glutSwapBuffers();

}



int main(int argc, char **argv) {

	// init GLUT and create Window
	glutInit(&argc, argv);
	glutInitDisplayMode( GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowPosition(1080, 450);
	glutInitWindowSize(WINDOW, WINDOW);
	glutCreateWindow("Ray Tracer");
	
	init();
	// register callbacks
	glutDisplayFunc(callbackDisplay);

	// glutSpecialFunc(specialKeyboard);
 	//glutKeyboardFunc(keyboard);
    // glutMouseFunc( mouse );

	// enter GLUT event processing cycle
	// glutSwapBuffers();
	glutMainLoop();
	return 0;
}
