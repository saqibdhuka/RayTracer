#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <time.h>
#include <iostream>
#include <array>
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
#define zNear 0.1f
#define zFar 50.0f




float tS; //t for the use of Ray
float t;
float tP; //t for the use of Ray

float xN, yN;
int rayDepth=1;
vec3 backColor(0.0, 0.0, 0.0);

const float LIMIT = tan(FOVY/2);



Sphere sph_objects[2];
Plane pla_objects[1];

Light light_object[4];

// vec3 shadingSphere(Ray ray, Sphere sphere){

// 	vec3 color(1.0, 0.0, 0.0);

// 	if(sphere.intersect(ray, t)){
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


vec3 tracePlane(Ray ray, Plane plane, Light light){ //Work on this

	vec3 I; //Light Intesity to calculate
	if(rayDepth = 0){
		return backColor;
	}else{
		if(yN > yN/2){
			return backColor;
		}
		// printf("Inside rayDepth > 0\n\n");
		// printf("Ia.x:%f, Ia.y:%f, Ia.z:%f \n", Ia.x, Ia.y, Ia.z);
		vec3 kd = plane.matPlane.kdColor;
		// printf("kd.x:%f, kd.y:%f, kd.z:%f \n", kd.x, kd.y, kd.z);
		vec3 Ii = light.colorIntensity;
		vec3 Ia =  kd * Ii;
		vec3 ks = plane.matPlane.ksColor;
		float  kr = plane.matPlane.reflection;
		float kt = plane.matPlane.refraction;
		int q = plane.matPlane.specular;
		float S =1;
		vec3 intersectionPoint = (ray.origin + (t * ray.direction));
		vec3 V = normalize(intersectionPoint);
		vec3 N = normalize(plane.plane_normal());
		vec3 Li = normalize(intersectionPoint - light_object[0].origin);
		vec3 R = normalize(2*((N * Li)*N) - Li);

		// printf("R.x:%f, R.y:%f, R.z:%f \n", R.x, R.y, R.z);

		vec3 rv = R*V;
		// printf("rv.x:%f, rv.y:%f, rv.z:%f \n", rv.x, rv.y, rv.z);

		vec3 tmp = vecPow(rv, q);
		// printf("tmp.x:%f, tmp.y:%f, tmp.z:%f \n", tmp.x, tmp.y, tmp.z);
		
		vec3 ik = Ia * kd;
		// printf("ik.x:%f, ik.y:%f, ik.z:%f \n", ik.x, ik.y, ik.z);
		vec3 ktemp = ks * tmp;
		// printf("ktemp.x:%f, ktemp.y:%f, ktemp.z:%f \n", ktemp.x, ktemp.y, ktemp.z);
		vec3  knl = kd*N*Li;
		// printf("knl.x:%f, knl.y:%f, knl.z:%f \n", knl.x, knl.y, knl.z);
		vec3 sum = S*Ii*(knl + ktemp);
		// printf("sum.x:%f, sum.y:%f, sum.z:%f \n", sum.x, sum.y, sum.z);

		I = ik + sum;
		// printf("I.x:%f, I.y:%f, I.z:%f \n", I.x, I.y, I.z);
		// if(I.x > 1.0)
		// 	I.x = 1.0;
		// if(I.y > 1.0)
		// 	I.y = 1.0;
		// if(I.z > 1.0)
		// 	I.z = 1.0;
		return normalize(I);
	}

}

vec3 trace(Ray ray, Sphere sphere, Light light){ //Sphere Trace 

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
		vec3 R = normalize(2*((N * Li)*N) - Li);
		

		// printf("R.x:%f, R.y:%f, R.z:%f \n", R.x, R.y, R.z);

		vec3 rv = R*V;
		// printf("rv.x:%f, rv.y:%f, rv.z:%f \n", rv.x, rv.y, rv.z);

		vec3 tmp = vecPow(rv, q);
		// printf("tmp.x:%f, tmp.y:%f, tmp.z:%f \n", tmp.x, tmp.y, tmp.z);
		
		vec3 ik = Ia * kd;
		// printf("ik.x:%f, ik.y:%f, ik.z:%f \n", ik.x, ik.y, ik.z);
		vec3 ktemp = ks * tmp;
		// printf("ktemp.x:%f, ktemp.y:%f, ktemp.z:%f \n", ktemp.x, ktemp.y, ktemp.z);
		vec3  knl = kd*N*Li;
		// printf("knl.x:%f, knl.y:%f, knl.z:%f \n", knl.x, knl.y, knl.z);
		vec3 sum = S*Ii*(knl + ktemp);
		I = ik + sum;
		// printf("I.x:%f, I.y:%f, I.z:%f \n", I.x, I.y, I.z);
		// if(I.x > 1.0)
		// 	I.x = 1.0;
		// if(I.y > 1.0)
		// 	I.y = 1.0;
		// if(I.z > 1.0)
		// 	I.z = 1.0;
		return normalize(I);
	}

}


// vec3 shadingPlane(Ray ray, Plane plane){

// 	vec3 color(0.5, 0.5, 0.5);

// 	if(plane.intersect(ray, t)){
// 		vec3 V = ray.direction;
// 		vec3 P = ray.origin + V * t;
// 		vec3 N = plane.get_normal(P);

// 		float ratio = dot(N,V);

// 		color = sphere.getkdColor() * (ratio *0.5);
// 	}

// 	return color;


// }

int sizeSphere(Sphere sphere[]){
	int count=0;
	for(int i=0; i< 4; i++){
		if(sphere[i].radius > 0)
			count++;
	}

	return count;
}

int countS =0;
int countP =0;
float tArray[3] = {500, 500, 500};//Array of t

vec2 minT(Ray ray){

	vec2 min(10000,-1);//min[MinTValue, Index]
	for(int i=0; i < 3;i++){
		
		tArray[i] = sph_objects[i].intersect(ray);
		if(tArray[i] != -1 && tArray < min){
			min.x = tArray[i];
			min.y = i;
		}
		// if(tArray[i] != -1)
		// 	printf("tArray[%d]: %f\n", i,tArray[i] );
		// if(tArray[i] < 0)
		// 	tArray[i] = 500;
	}
	// if(tArray[4] > 0)
	// 	tArray[4] = 500;

	// float min = tArray[0];
	// for(int k=0; k < 3; k++){
	// 	if(tArray[k] == -1){
	// 		min = tArray[k+1];
	// 		continue;
	// 	}
	// 	if(tArray[k] < min){
	// 		min = tArray[k];
	// 	}
	// }

	return min;

}

// float minTindex(float tCheck){
// 	for(int i =0; i<3; i++){
// 		if(tArray[i] == tCheck)
// 			return i;
// 	}
// }

void render(){

	float scale = 26.3;

	for(int x=-WINDOW/2; x < WINDOW/2; x++){
		 xN=x * LIMIT* LIMIT*scale;
		// printf("xN: %f\n", xN); 
		for (int y = -WINDOW/2; y < WINDOW/2; y++)
		{
			 yN=y * LIMIT* LIMIT*scale;			// glBegin(GL_POINTS);
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
			// // ray.shootRay(t);
			// if(y < y/2){
			// 	t = sph_objects[0].intersect(ray);
			// 	sphereHit = 1;
			// }else{
			
			// tP = pla_objects[0].intersectPlane(ray);

			// if(sph_objects[0].intersect(ray) < pla_objects[0].intersectPlane(ray)){
			// 	t = sph_objects[0].intersect(ray); 
			// 	sphereHit = 1;
			// }else{
			// 	t = pla_objects[0].intersectPlane(ray);
			// }
			// }
			// t = tS < tP ? tS : tP;
			
			// count++;
			// if(tS >0)
			// 	printf("ts: %f\n",tS);
			// if(tP >0)
			vec2 T = minT(ray);
			t = T.x;
			if(t!=-1)
				printf("t: %f\n",t);
			int index = T.y;
			printf("Index: %d\n",index);
			// t = sph_objects[0].intersect(ray);
			if(index != -1){
					countS++;
					vec3 colorLight = trace(ray, sph_objects[index], light_object[0]);
					glColor3f(colorLight.x, colorLight.y, colorLight.z);
						
					glBegin(GL_POINTS);
						glVertex3f(ray.origin.x + t*ray.direction.x , 
							ray.origin.y + t*ray.direction.y , 
							ray.origin.z + t*ray.direction.z );
					glEnd();
					// glutSwapBuffers();

			}else{
				t = pla_objects[0].intersectPlane(ray);

				// printf("Plane t: %f\n", t);
				countP++;
				vec3 colorLight = tracePlane(ray, pla_objects[0], light_object[0]);
				// printf(" PLANE R: %f, PLANE G: %f, PLANE B: %f\n", colorLight.x, colorLight.y, colorLight.z);
				glColor3f(colorLight.x, colorLight.y, colorLight.z);
				// glEnable(GL_CULL_FACE);
				// glCullFace(GL_FRONT);
				glPushMatrix();
					// glTranslatef(0.0, 0.0, 0.046);
					
					glBegin(GL_POINTS);
						glVertex3f(ray.origin.x + t*ray.direction.x , 
							ray.origin.y + t*ray.direction.y , 
							ray.origin.z + t*ray.direction.z );
					glEnd();
				glPopMatrix();
				// glutSwapBuffers();
				
			}
			// for(int i =0 ; i< 2; i++){
			// 	tS = sph_objects[i].intersect(ray);
			// 	if(tS != -1){
			// 		t = tS;
			// 		countS++;
			// 		vec3 colorLight = trace(ray, sph_objects[i], light_object[0]);
			// 		glColor3f(colorLight.x, colorLight.y, colorLight.z);
						
			// 		glBegin(GL_POINTS);
			// 			glVertex3f(ray.origin.x + t*ray.direction.x , 
			// 				ray.origin.y + t*ray.direction.y , 
			// 				ray.origin.z + t*ray.direction.z );
			// 		glEnd();
			// 		glutSwapBuffers();

			// 	}
			// 	else{
			// 		t = tP;
			// 		// printf("Plane t: %f\n", t);
			// 		countP++;
			// 		vec3 colorLight = tracePlane(ray, pla_objects[0], light_object[0]);
			// 		// printf(" PLANE R: %f, PLANE G: %f, PLANE B: %f\n", colorLight.x, colorLight.y, colorLight.z);
			// 		glColor3f(colorLight.x, colorLight.y, colorLight.z);
			// 		// glEnable(GL_CULL_FACE);
			// 		// glCullFace(GL_FRONT);
			// 		glPushMatrix();
			// 			glTranslatef(0.0, 0.0, -2.45);
						
			// 			glBegin(GL_POINTS);
			// 				glVertex3f(ray.origin.x + t*ray.direction.x , 
			// 					ray.origin.y + t*ray.direction.y , 
			// 					ray.origin.z + t*ray.direction.z );
			// 			glEnd();
			// 		glPopMatrix();
			// 		// glutSwapBuffers();
					
			// 	}
			// }
				

			

			// if(tS != -1 && tP != -1){ // Hit the sphere
				// glColor3f(sph_objects[0].mat.ksColor.x, sph_objects[0].mat.ksColor.y, sph_objects[0].mat.ksColor.z);
				// glColor3f(1.0, 0.5, 1.0);
				// t = tS;
				// if(tS < tP){
				// 	t = tS;
				// 	countS++;
				// 	vec3 colorLight = trace(ray, sph_objects[0], light_object[0]);
				// 	// printf(" R: %f, G: %f, B: %f\n", colorLight.x, colorLight.y, colorLight.z);
				// 	// if(t != -1){
				// 	// 	glColor3f(colorLight.x, colorLight.y, colorLight.z);
				// 	// }else{
				// 	// 	glColor3f(0.0, 0.0, 0.0);
				// 	// }
				// 	// glEnable(GL_CULL_FACE);
				// 	// glCullFace(GL_FRONT);
				// 	glColor3f(colorLight.x, colorLight.y, colorLight.z);
					
				// 	glBegin(GL_POINTS);
				// 		glVertex3f(ray.origin.x + t*ray.direction.x , 
				// 			ray.origin.y + t*ray.direction.y , 
				// 			ray.origin.z + t*ray.direction.z );
				// 	glEnd();
				// }
				// else if(yN < (400*LIMIT*LIMIT*25)/2){
				// 	t = tP;
				// 	countP++;
				// 	vec3 colorLight = tracePlane(ray, pla_objects[0], light_object[0]);
				// 	printf(" PLANE R: %f, PLANE G: %f, PLANE B: %f\n", colorLight.x, colorLight.y, colorLight.z);
				// 	glColor3f(colorLight.x, colorLight.y, colorLight.z);
				// 	// glEnable(GL_CULL_FACE);
				// 	// glCullFace(GL_FRONT);
					
				// 	glBegin(GL_POINTS);
				// 		glVertex3f(ray.origin.x + t*ray.direction.x , 
				// 			ray.origin.y + t*ray.direction.y , 
				// 			ray.origin.z + t*ray.direction.z );
				// 	glEnd();
				// }
			// }
			// else {
			// 	t = pla_objects[0].intersectPlane(ray);
			// 	countP++;
			// 	vec3 colorLight = tracePlane(ray, pla_objects[0], light_object[0]);
			// 	printf(" PLANE R: %f, PLANE G: %f, PLANE B: %f\n", colorLight.x, colorLight.y, colorLight.z);
			// 	glColor3f(colorLight.x, colorLight.y, colorLight.z);
			// 	// glEnable(GL_CULL_FACE);
			// 	// glCullFace(GL_FRONT);
				
			// 	glBegin(GL_POINTS);
			// 		glVertex3f(ray.origin.x + t*ray.direction.x , 
			// 			ray.origin.y + t*ray.direction.y , 
			// 			ray.origin.z + t*ray.direction.z );
			// 	glEnd();
			// }

			// else{
			// 	glColor3f(0.0,0.0,0.0);
			// 	glBegin(GL_POINTS);
			// 		glVertex3f(ray.origin.x + t*ray.direction.x , 
			// 			ray.origin.y + t*ray.direction.y , 
			// 			ray.origin.z + t*ray.direction.z );
			// 	glEnd();
			// }
			
			// if(t != -1){
				
			// }
			// vec3 kd = shadingSphere(ray, sph_objects[0]);
			// sph_objects[0].mat.kdColor = kd;
			// sph_objects[0].draw();
		}
	}
	printf("countS: %d\n", countS);
	printf("countP: %d\n", countP);
}

void init(){
	
    glViewport(0, 0, WINDOW, WINDOW);

    glClearColor(backColor.x, backColor.y, backColor.z, 0.0);
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
	// glPushMatrix();
		// glLoadIdentity();

	vec3 lightPos(0.0, 0.0, 0.0);
	vec3 lightColor(5.0, 5.0, 5.0);
	Light light(lightPos, lightColor);
	light_object[0] = light;

	

	
// glPopMatrix();
	vec3 center_ball1(-1.0, 0.5, -20.0);
	vec3 center_ball2(-2.0,-1.0,-15.0);
	vec3 center_ball3(1.0, 1.0, -15.0);

	float radius_ball1 = 1.5;
	float radius_ball2 = 1.0;
	float radius_ball3 = 2.0;

	//Setting up the material
	Material Blue = Material(vec3(0.0, 0.651, 0.851), vec3(1.0, 1.0, 1.0),24, 0.5f, 0.0f, 1.5f);

	Material DarkRed = Material(vec3(0.502, 0.0353, 0.0353), vec3(1.0, 1.0, 1.0),47, 0.5f, 0.0f, 1.5f);

	Material Orange = Material(vec3(0.9804, 0.4588, 0.0), vec3(1.0, 1.0, 1.0),90, 0.5f, 0.0f, 1.5f);

	Material PlaneMaterial = Material(vec3(0.0, 1.0, 1.0), vec3(0.0, 1.0, 1.0), 5.0, 1.0f, 0.0f, 1.5f);

	Sphere ball1(center_ball1, radius_ball1, Blue);
	Sphere ball2(center_ball2, radius_ball2, DarkRed);
	Sphere ball3(center_ball3, radius_ball3, Orange);

	sph_objects[0] = ball1;
	sph_objects[1] = ball2;
	sph_objects[2] = ball3;




	Plane plane(0, 0, 1, -1, PlaneMaterial);
	pla_objects[0] = plane;
	
	render();
	// glFlush();
	glutSwapBuffers();

}



int main(int argc, char **argv) {

	// init GLUT and create Window
	glutInit(&argc, argv);
	// glutInitDisplayMode( GLUT_SINGLE | GLUT_RGB);
	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
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
