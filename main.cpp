// #include <GL/glew.h>
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




// float tS; //t for the use of Ray
// float t;
// float tP; //t for the use of Ray

float xN, yN;

int rayDepth = 4;

vec3 backColor(0.0, 0.0, 0.0);
char scene; //Scene called

const float LIMIT = tan(DegreesToRadians *FOVY/2)/WINDOW;

int countS =0;
int countP =0;
float tArray[3] = {500, 500, 500};//Array of t

Sphere sph_objects[3];
Plane pla_objects[2];

Light light_object[4];


void setRayDepth(int depth){
	rayDepth = depth;
}
vec3 I_plane = backColor; //Light Intesity to calculate


vec2 minT(Ray ray){
	vec2 min(10000,-1);//min[MinTValue, Index]
		for(int i=0; i < 3 ;i++){
			tArray[i] = sph_objects[i].intersect(ray);
			if(tArray[i] != -1 && tArray[i] < min.x){
				min.x = tArray[i];
				min.y = i;
		}
	}
	return min;

}



vec3 refraction(Ray ray, vec3 normal, vec3 incident, float ratio){ //WORKING

		double dotProd = -dot(normal, incident);
		double tmp = ratio * ratio * (1.0 - dotProd * dotProd);
		if(tmp > 1.0)
			return vec3(-1, -1, -1);
		double root = sqrt(1.0 - tmp);
		return ratio * incident + (ratio * dotProd - root) * normal;

}



vec3 reflection(vec3 normal, vec3 incident){
		double tmp = -dot(normal, incident);
		return incident + 2 * (tmp) * normal;

}




	// float planeIntersect = pla_objects[0].intersectPlane(ray);
	// if(planeIntersect < min.x){
	// min.x = planeIntersect;
	// min.y = -1;
	// }


vec3 tracePlane(Ray ray, Plane plane, Light light, int depth){ //Work on this

	if(depth == 0){
		return backColor;
	}else{
		vec3 kd = (plane.matPlane.kdColor);
		vec3 Ii = (light.colorIntensity);
		// float ambientStrength = 0.75f;
		vec3 ks = (plane.matPlane.ksColor);
		float kr = plane.matPlane.reflection;
		float kt = plane.matPlane.refraction;
		// vec3  kr = vec3(sphere.mat.reflection, sphere.mat.reflection , sphere.mat.reflection);
		// vec3 kt = vec3(sphere.mat.refraction, sphere.mat.refraction, sphere.mat.refraction);
		// vec3 Ia =  Ii;
		vec3 Ia(0.5, 0.5, 0.5);
		// vec3 Ia =  Ii;
		// vec3 Ia = Ii * kr;
		if(scene == '1' || scene =='2' || scene =='3'){
			if(yN > yN/2){
				float diff = (yN - yN/2) * 10;
				Ia = vec3(0.5 - diff,0.5 - diff,0.5 - diff);
			}else{
				Ia = vec3(0.5, 0.5, 0.5);
			}
		}
		int q = plane.matPlane.specular;
		// float S =1;

		// vec2 T = minT(ray);
		// double t = T.x;
		// int index = T.y;
		// if(t < 0 ){
		// 	return backColor;
		// }

		double t = plane.intersectPlane(ray);
		if(t != -1){
			plane = pla_objects[0];

		vec3 intersectionPoint = (ray.origin + (t * ray.direction));


		// vec3 incidentRayVector = normalize(ray.origin - intersectionPoint);

		vec3 V = normalize(-ray.direction);
		vec3 N = normalize(plane.plane_normal());
		vec3 Li = normalize(light.origin - intersectionPoint);
		// vec3 R = normalize(((2*(dot(N,incidentRayVector))*N) - incidentRayVector));
		vec3 R = normalize((2*(dot(N,Li))*N) - Li);

		// vec3 Ir = 0;
		// vec3 It = 0;

		int S = 1;

		if(scene == '2'){
			// Ia = (0.5,0.5,0.5);
			Ray rayShadow(intersectionPoint, Li + intersectionPoint);//ray.direction );
			for(int s=0; s < 3; s++){
				// float t_shadow = minT(rayShadow).x;
				double t_shadow = sph_objects[s].intersect(rayShadow);
				// printf("t_shadow: %f\n", t_shadow);
				// vec3 point_intersect = rayShadow.origin + (t_shadow * rayShadow.direction);
				if(/*length(light.origin - point_intersect) > length(light.origin - sphere.center)*/ t_shadow > 0.01 /*&& length(light.origin - point_intersect) <= 0.95*/ /*&& point_intersect.z < sph_objects[s].center.z*/){
					S = 0;
					break;
				}
			}
		}else if(scene =='3'){
			// Ia = (0.5,0.5,0.5);

			Ray rayShadow(intersectionPoint, Li - intersectionPoint);//ray.direction );
			for(int s=0; s < 3; s++){
				// float t_shadow = minT(rayShadow).x;
				double t_shadow = sph_objects[s].intersect(rayShadow);
				// printf("t_shadow: %f\n", t_shadow);
				// vec3 point_intersect = rayShadow.origin + (t_shadow * rayShadow.direction);
				if(/*length(light.origin - point_intersect) > length(light.origin - sphere.center)*/ t_shadow > 0.01 /*&& length(light.origin - point_intersect) <= 0.95*/ /*&& point_intersect.z < sph_objects[s].center.z*/){
					S = 0;
					break;
				}
			}
		}else{
			// Ia = (0.2,0.2,0.2);
			Ray rayShadow(intersectionPoint, Li - intersectionPoint);//ray.direction );
			for(int s=0; s < 3; s++){
				// float t_shadow = minT(rayShadow).x;
				double t_shadow = sph_objects[s].intersect(rayShadow);
				// printf("t_shadow: %f\n", t_shadow);
				// vec3 point_intersect = rayShadow.origin + (t_shadow * rayShadow.direction);
				if(/*length(light.origin - point_intersect) > length(light.origin - sphere.center)*/ t_shadow > 0.01 /*&& length(light.origin - point_intersect) <= 0.95*/ /*&& point_intersect.z < sph_objects[s].center.z*/){
					S = 0;
					break;
				}
			}
		}
		// if(pla_objects[0].intersectPlane(rayShadow) != -1){
		// S=1;
		// }

		// I_planes = ik + (sum);

		// I = ik + (Ir * kr) + (It * kt) + (sum);
		// vec3 Ir = normalize(I);//Reflected Ray Direction
		// vec3 It = normalize(I); //Refracted Ray Direction
		// vec3 reflectDir = ((2*N)*(dot(N,V))) - V;

		vec3 incidenceRayDir = normalize(ray.direction);
		// float tmpDiv = length(N) / length(incidenceRayDir);
		// N = normalize(N);
		// vec3 ik = kd*N;
		// float knl = dot(ik,Li);
		// if(knl < 0){
		// knl = 0;
		// }
		// float theta1 = acos(tmpDiv);

		// float theta1 = dot(N, Li);

		// printf("cos: %f\n", theta1);
		double refraction_index_air = 1.00029;
		double ratio = refraction_index_air / plane.matPlane.index_refraction;
		// ratio = refraction_index_air / sphere.mat.index_refraction;
		vec3 Nref = (plane.plane_normal());
		if(dot(incidenceRayDir, Nref) < 0){
			ratio = refraction_index_air / plane.matPlane.index_refraction; //n1/n2
		}else {
			ratio = plane.matPlane.index_refraction / refraction_index_air; //n2/n1
			Nref = -Nref;
		}
		Nref = normalize(Nref);

		// vec3 ReflectionVec = normalize((2 * (dot((sphere.get_normal(intersectionPoint)),(-ray.direction)) * (sphere.get_normal(intersectionPoint)))) - (-ray.direction));
		// vec3 ReflectionVec = normalize((2 * (dot(N,V)) * N)- (V));
		vec3 ReflectionVec = (reflection(Nref, incidenceRayDir));
		// vec3 ReflectionVec = normalize((2*(dot(N,V))*N) - V);

		// float c1 = dot(Nref, incidenceRayDir);
		// if(c1 < 0)
		// c1 = 0;
		// // float c1 = cos(theta1);
		// float c2temp = 1-(pow(ratio,2)*((pow((sin(theta1)), 2))));
		// float c2;

		// if(c2temp > 0)
		// c2 = sqrt(c2temp);
		// else
		// c2 = 0;
		// float temp = ratio * sin(DegreesToRadians * theta1);
		// float theta2 = asin(temp);

		vec3 refractDir = (refraction(ray, Nref, incidenceRayDir, ratio));
		if(refractDir.x == -1 && refractDir.y == -1 && refractDir.z == -1)
		{
			refractDir = ReflectionVec;
		}

		Ray ray2(intersectionPoint, ReflectionVec);
		Ray ray3(intersectionPoint, refractDir);
		vec3 Ir_into_kr(0,0,0), It_into_kt(0,0,0);
		// Ir = trace(ray2, sphere, light, depth-1);
		// It = trace(ray3, sphere, light, depth-1);

		// if(kr > 0){
			Ir_into_kr += kr * tracePlane(ray2, plane, light, depth-1);
		// }
		// else{
		// Ir_into_kr = 0;
		// }


		// if(kt > 0){
			It_into_kt += kt * tracePlane(ray3, plane, light, depth-1);
		// }
		// else{
		// It_into_kt = 0;
		// }

		if(scene == '1'){
			// Ray rayShadow(intersectionPoint, Li);
			// for(int s=0; s< 3; s++){
			// float t_shadow = sph_objects[s].intersect(rayShadow);
			// vec3 point_intersect = rayShadow.origin + (t_shadow * rayShadow.direction);
			// if(t_shadow != -1 /*&& sph_objects[s].center.z > point_intersect.z*/){
			// S = 1;
			// break;
			// }
			// }
			S=1;
			Ir_into_kr =0;
			It_into_kt =0;
		}

		if(scene == '2'){
			// S=1;
			It_into_kt = 0;//NP REGRACTION
		}

		if(scene == '3'){
			Ir_into_kr = 0;
		}


		// if(scene == '4'){
		// // S =1;
		// }

		double rv = dot(R,V);

		if(rv < 0){
			rv = 0;
		}

		double tmp = pow(rv, q);
		// vec3 ik = kd*N;
		vec3 ktemp = (ks * tmp);
		vec3 nli = dot(N,Li);

		if(nli < 0){
			nli = 0;
		}

		vec3 knl = kd * nli;

		vec3 sum = S*Ii*(knl + ktemp);

		I_plane = (Ia * kd) + (Ir_into_kr) + (It_into_kt) + (sum);
		// if(I.x > 1)
		// I.x = 1;
		// if(I.y > 1)
		// I.y = 1;
		// if(I.z > 1)
		// I.z = 1;
		return (I_plane);
		}
	}

}


vec3 I = backColor; //Light Intesity to calculate


vec3 trace(Ray ray, Sphere sphere, Light light, int depth){ //Sphere Trace


	if(depth == 0){
		return backColor;
	}else{
		vec3 kd = sphere.mat.kdColor;
		vec3 Ii = light.colorIntensity;
		// float ambientStrength = 0.75f;
		vec3 ks = sphere.mat.ksColor;
		float kr = sphere.mat.reflection;
		float kt = sphere.mat.refraction;
		// vec3  kr = vec3(sphere.mat.reflection, sphere.mat.reflection , sphere.mat.reflection);
		// vec3 kt = vec3(sphere.mat.refraction, sphere.mat.refraction, sphere.mat.refraction);
		// vec3 Ia =  Ii;
		vec3 Ia;
		// vec3 Ia =  Ii;
		// vec3 Ia = Ii * kr;
		int q = sphere.mat.specular;
		// float S =1;

		vec2 T = minT(ray);
		double t = T.x;
		int index = T.y;
		// if(t < 0 ){
		// 	return backColor;
		// }

		// double t = sphere.intersect(ray);
		if(t != -1){
			sphere = sph_objects[index];


		vec3 intersectionPoint = (ray.origin + (t * ray.direction));


		// vec3 incidentRayVector = normalize(ray.origin - intersectionPoint);

		vec3 V = normalize(-ray.direction);
		vec3 N = normalize(sphere.get_normal(intersectionPoint));
		vec3 Li = normalize(light.origin - intersectionPoint);
		// vec3 R = normalize(((2*(dot(N,incidentRayVector))*N) - incidentRayVector));
		vec3 R = normalize((2*(dot(N,Li))*N) - Li);

		// vec3 Ir = 0;
		// vec3 It = 0;

		int S = 1;


		if(scene == '2'){
			S = 1;
			Ia = 0.5,0.5,0.5;

			Ray rayShadow(intersectionPoint, light.origin);//ray.direction );
			for(int s=0; s < 3; s++){
				// float t_shadow = minT(rayShadow).x;
				double t_shadow = sph_objects[s].intersect(rayShadow);
				// printf("t_shadow: %f\n", t_shadow);
				// vec3 point_intersect = rayShadow.origin + (t_shadow * rayShadow.direction);
				if(/*length(light.origin - point_intersect) > length(light.origin - sphere.center)*/ t_shadow > 0.01 /*&& length(light.origin - point_intersect) <= 0.95*/ /*&& point_intersect.z < sph_objects[s].center.z*/){
					S = 0;
					break;
				}
			}
		}else if(scene =='3'){
			S=1;
			Ia = (0.4,0.4,0.4);
			// Ia = Ii * ks;

			Ray rayShadow(intersectionPoint, Li + intersectionPoint);//ray.direction );
			for(int s=0; s < 3; s++){
				// float t_shadow = minT(rayShadow).x;
				double t_shadow = sph_objects[s].intersect(rayShadow);
				// printf("t_shadow: %f\n", t_shadow);
				// vec3 point_intersect = rayShadow.origin + (t_shadow * rayShadow.direction);
				if(/*length(light.origin - point_intersect) > length(light.origin - sphere.center)*/ t_shadow > 0.01 /*&& length(light.origin - point_intersect) <= 0.95*/ /*&& point_intersect.z < sph_objects[s].center.z*/){
					S = 0;
					break;
				}
			}
		}else{
			S=1;
			Ia = (0.4,0.4,0.4);
			// Ia = Ii * ks;

			Ray rayShadow(intersectionPoint, Li + intersectionPoint);//ray.direction );
			for(int s=0; s < 3; s++){
				// float t_shadow = minT(rayShadow).x;
				double t_shadow = sph_objects[s].intersect(rayShadow);
				// printf("t_shadow: %f\n", t_shadow);
				// vec3 point_intersect = rayShadow.origin + (t_shadow * rayShadow.direction);
				if(/*length(light.origin - point_intersect) > length(light.origin - sphere.center)*/ t_shadow > 0.01 /*&& length(light.origin - point_intersect) <= 0.95*/ /*&& point_intersect.z < sph_objects[s].center.z*/){
					S = 0;
					break;
				}
			}
		}



		// if(pla_objects[0].intersectPlane(rayShadow) != -1){
		// S=1;
		// }

		// I_planes = ik + (sum);

		// I = ik + (Ir * kr) + (It * kt) + (sum);
		// vec3 Ir = normalize(I);//Reflected Ray Direction
		// vec3 It = normalize(I); //Refracted Ray Direction
		// vec3 reflectDir = ((2*N)*(dot(N,V))) - V;

		vec3 incidenceRayDir = normalize(ray.direction);
		// float tmpDiv = length(N) / length(incidenceRayDir);
		// N = normalize(N);
		// vec3 ik = kd*N;
		// float knl = dot(ik,Li);
		// if(knl < 0){
		// knl = 0;
		// }
		// float theta1 = acos(tmpDiv);

		// float theta1 = dot(N, Li);

		// printf("cos: %f\n", theta1);
		double refraction_index_air = 1.00029;
		double ratio = refraction_index_air / sphere.mat.index_refraction;
		// ratio = refraction_index_air / sphere.mat.index_refraction;
		vec3 Nref = (sphere.get_normal(intersectionPoint));
		if(dot(incidenceRayDir, Nref) < 0){
			ratio = refraction_index_air / sphere.mat.index_refraction; //n1/n2
		}else {
			ratio = sphere.mat.index_refraction / refraction_index_air; //n2/n1
			Nref = -Nref;
		}
		Nref = normalize(Nref);

		// vec3 ReflectionVec = normalize((2 * (dot((sphere.get_normal(intersectionPoint)),(-ray.direction)) * (sphere.get_normal(intersectionPoint)))) - (-ray.direction));
		// vec3 ReflectionVec = normalize((2 * (dot(N,V)) * N)- (V));
		vec3 ReflectionVec = (reflection(Nref, -incidenceRayDir));
		// vec3 ReflectionVec = normalize((2*(dot(N,V))*N) - V);

		// float c1 = dot(Nref, incidenceRayDir);
		// if(c1 < 0)
		// c1 = 0;
		// // float c1 = cos(theta1);
		// float c2temp = 1-(pow(ratio,2)*((pow((sin(theta1)), 2))));
		// float c2;

		// if(c2temp > 0)
		// c2 = sqrt(c2temp);
		// else
		// c2 = 0;
		// float temp = ratio * sin(DegreesToRadians * theta1);
		// float theta2 = asin(temp);

		vec3 refractDir = (refraction(ray, Nref, incidenceRayDir, ratio));
		if(refractDir.x == -1 && refractDir.y == -1 && refractDir.z == -1)
		{
			refractDir = ReflectionVec;
		}

		Ray ray2(intersectionPoint, ReflectionVec);
		Ray ray3(intersectionPoint, refractDir);
		vec3 Ir_into_kr(0,0,0), It_into_kt(0,0,0);
		// Ir = trace(ray2, sphere, light, depth-1);
		// It = trace(ray3, sphere, light, depth-1);

		// if(kr > 0){
			Ir_into_kr += kr * trace(ray2, sphere, light, depth-1);
		// }
		// else{
		// Ir_into_kr = 0;
		// }


		// if(kt > 0){
			It_into_kt += kt * trace(ray3, sphere, light, depth-1);
		// }
		// else{
		// It_into_kt = 0;
		// }

		if(scene == '1'){
			// Ray rayShadow(intersectionPoint, Li);
			// for(int s=0; s< 3; s++){
			// float t_shadow = sph_objects[s].intersect(rayShadow);
			// vec3 point_intersect = rayShadow.origin + (t_shadow * rayShadow.direction);
			// if(t_shadow != -1 /*&& sph_objects[s].center.z > point_intersect.z*/){
			// S = 1;
			// break;
			// }
			// }

			S=1;
			Ir_into_kr =0;
			It_into_kt =0;
		}

		if(scene == '2'){
			// S=1;
			It_into_kt = 0;//NP REGRACTION
		}

		// if(scene == '4'){
		// // S =1;
		// }

		double rv = dot(R,V);

		if(rv < 0){
			rv = 0;
		}

		double tmp = pow(rv, q);
		// vec3 ik = kd*N;
		vec3 ktemp = (ks * tmp);
		vec3 nli = dot(N,Li);

		if(nli < 0){
			nli = 0;
		}

		vec3 knl = kd * nli;

		vec3 sum = S*Ii*(knl + ktemp);

		I = (Ia * kd) + (Ir_into_kr) + (It_into_kt) + (sum);
		// if(I.x > 1)
		// I.x = 1;
		// if(I.y > 1)
		// I.y = 1;
		// if(I.z > 1)
		// I.z = 1;
		return (I);
	}
		// return trace(ray2, sphere, light, depth-1);
		// return trace(ray3, sphere, light, depth-1);
		// for(int i = depth; i >= 1; i--){
		// I = ik + (Ir * kr) + (It * kt) + (sum);
		// vec3 Ir = I;//Reflected Ray Direction
		// vec3 It = I; //Refracted Ray Direction
		// }


		// printf("Ix: %f,  Iy: %f,  Iz: %f\n", I.x, I.y, I.z);


		// return normalize(I);
		// return trace(ray, sphere, light, depth-1);
		}
}


void render(){

	float scale = 2;

	for(int x=-WINDOW/2; x < WINDOW/2; x++){
		xN=x * LIMIT*scale;
		// printf("xN: %f\n", xN);
		for (int y = -WINDOW/2; y < WINDOW/2; y++)
		{
			yN=y * LIMIT*scale; // glBegin(GL_POINTS);
			// glVertex3f(xN, yN, -1.0);
			// glEnd();

			vec3 orig(0.0, 0.0, 0.0);
			vec3 dir(xN, yN, -1.0);
			Ray ray(orig, dir);

			vec2 T = minT(ray);
			double t = T.x;
			// if(t!=-1)
			// printf("t: %f\n",t);
			int index = T.y;
			// printf("Index: %d\n",index);
			// t = sph_objects[0].intersect(ray);
			if(index != -1){
				countS++;
				// Ir = vec3(0,0,0);
				// It = vec3(0,0,0);
				// I = vec3(0,0,0);
				vec3 colorLight = trace(ray, sph_objects[index], light_object[0], rayDepth);
				// I = vec3(0,0,0);
				glColor3f(colorLight.x, colorLight.y, colorLight.z);
				// printf("  R: %f,  G: %f,  B: %f\n", colorLight.x, colorLight.y, colorLight.z);
				glBegin(GL_POINTS);
				glVertex3f(ray.origin.x + t*ray.direction.x,
								ray.origin.y + t*ray.direction.y,
								ray.origin.z + t*ray.direction.z);
				glEnd();
				// glutSwapBuffers();


			}
			else {
				t = pla_objects[0].intersectPlane(ray);

				countP++;
				// Ir = vec3(0,0,0);
				// It = vec3(0,0,0);
				// I = vec3(0,0,0);
				vec3 colorLight = tracePlane(ray, pla_objects[0], light_object[0], rayDepth);
				// I = vec3(0,0,0);
				glColor3f(colorLight.x, colorLight.y, colorLight.z);

				glBegin(GL_POINTS);
				glVertex3f(ray.origin.x + t*ray.direction.x,
					ray.origin.y + t*ray.direction.y ,
					ray.origin.z + t*ray.direction.z );
				glEnd();
				// glutSwapBuffers();

			}
		}
	}
	// printf("countS: %d\n", countS);
	// printf("countP: %d\n", countP);
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
}

void scene1(){
	vec3 center_ball1(-2.0, 0.0, -8.0);
	vec3 center_ball2(-0.15,0.0,-10.0);
	vec3 center_ball3(3.0, 0.0, -15.0);

	float radius_ball1 = 1.0;
	float radius_ball2 = 1.0;
	float radius_ball3 = 1.5;

	//Setting up the material
	Material Yellow = Material(vec3(0.93, 1.0, 0.0), vec3(1.0, 1.0, 1.0),20, 0.5, 0.0, 1.5);

	Material DarkRed = Material(vec3(0.502, 0.0353, 0.0353), vec3(1.0, 1.0, 1.0),50, 0.5, 0.0, 1.5);

	Material Purple = Material(vec3(0.4314, 0.0589, 0.6196), vec3(1.0, 1.0, 1.0),100, 0.5, 0.0, 1.5);

	Sphere ball1(center_ball1, radius_ball1, Yellow);
	Sphere ball2(center_ball2, radius_ball2, Purple);
	Sphere ball3(center_ball3, radius_ball3, DarkRed);

	sph_objects[0] = ball1;
	sph_objects[1] = ball2;
	sph_objects[2] = ball3;

	//Setting up the Plane
	Material PlaneMaterial = Material(vec3(0.0, 1.0, 1.0), vec3(1.0, 1.0, 1.0), 50, 0.0, 0.0, 1.5);
	Plane plane(0, 0, 1, -5, PlaneMaterial);
	pla_objects[0] = plane;


	//Setting up the Light
	vec3 lightPos(-20.0, -20.0, 0.0);
	vec3 lightColor(0.5, 0.5, 0.5);//Ambient Light Color
	Light light(lightPos, lightColor);
	light_object[0] = light;

	render();
}


void scene2(){
	vec3 center_ball1(0.0, 0.0, -20.0);
	vec3 center_ball2(-4.5 ,1.0,-17.0);
	vec3 center_ball3(4.0, 3.0, -20.0);

	float radius_ball1 = 2.0;
	float radius_ball2 = 2.0;
	float radius_ball3 = 2.0;

	//Setting up the material
	Material Orange = Material(vec3(0.8, 0.4, 0.0), vec3(1.0, 1.0, 1.0),70, 1.0, 0.0, 1.52);

	Material DarkRed = Material(vec3(0.502, 0.0353, 0.0353), vec3(1.0, 1.0, 1.0),50, 1.0, 0.0, 1.52);

	Material Purple = Material(vec3(0.4314, 0.0589, 0.6196), vec3(1.0, 1.0, 1.0),100, 1.0, 0.0, 1.52);

	Sphere ball1(center_ball1, radius_ball1, Orange);
	Sphere ball2(center_ball2, radius_ball2, Purple);
	Sphere ball3(center_ball3, radius_ball3, DarkRed);

	sph_objects[0] = ball1;
	sph_objects[1] = ball2;
	sph_objects[2] = ball3;

	//Setting up the Plane
	Material PlaneMaterial = Material(vec3(0.0, 0.2, 0.2), vec3(1.0, 1.0, 1.0), 50, 0.0, 0.0, 1.33);
	Plane plane(0, 0, 1, -5, PlaneMaterial);
	pla_objects[0] = plane;

	//Setting up the Light
	vec3 lightPos(0.0, 50.0, 50.0);
	vec3 lightColor(0.5, 0.5, 0.5);//Ambient Light Color
	Light light(lightPos, lightColor);
	light_object[0] = light;

	render();
}

void scene3(){
	vec3 center_ball1(0.0, -3.0, -15.0);
	vec3 center_ball2(-4.0,-3.0,-20.0);
	vec3 center_ball3(2.5, -2.0, -25.0);

	float radius_ball1 = 3.0;
	float radius_ball2 = 2.0;
	float radius_ball3 = 2.0;

	//Setting up the material
	Material BlueTransparent = Material(vec3(0.3, 0.3, 0.3), vec3(1.0, 1.0, 1.0),50, 0.0, 1.0, 1.52);

	Material DarkRed = Material(vec3(0.502, 0.0353, 0.0353), vec3(0.7, 0.7, 0.7),80, 0.0, 1.0, 1.52);

	Material Purple = Material(vec3(0.4314, 0.0589, 0.6196), vec3(0.7, 0.7, 0.7),40, 0.5, 0.5, 1.52);

	Sphere ball1(center_ball1, radius_ball1, BlueTransparent);
	Sphere ball2(center_ball2, radius_ball2, Purple);
	Sphere ball3(center_ball3, radius_ball3, DarkRed);

	sph_objects[0] = ball1;
	sph_objects[1] = ball2;
	sph_objects[2] = ball3;

	//Setting up the Plane
	Material PlaneMaterial = Material(vec3(0.0, 0.5, 0.5), vec3(1.0, 1.0, 1.0), 80.0, 0.5, 0.0, 1.5);
	Plane plane(0, 0, 1, -5, PlaneMaterial);
	pla_objects[0] = plane;

	//Setting up the Light
	vec3 lightPos(0.0, 10.0, 100.0);
	vec3 lightColor(0.5, 0.5, 0.5);//Ambient Light Color
	Light light(lightPos, lightColor);
	light_object[0] = light;

	render();
}


void scene4(){
	vec3 center_ball1(0.0, -1.5, -10.0);
	vec3 center_ball2(-2.5,-1.0,-15.0);
	vec3 center_ball3(4.0, -0.5, -25.0);

	float radius_ball1 = 1.75;
	float radius_ball2 = 1.5;
	float radius_ball3 = 2.0;

	//Setting up the material
	Material BlueTransparent = Material(vec3(0.3, 0.3, 0.3), vec3(1.0, 1.0, 1.0),25, 0.0, 1.0, 1.52);

	Material DarkRed = Material(vec3(0.502, 0.0353, 0.0353), vec3(1.0, 1.0, 1.0),80, 1.0f, 0.0, 1.52);

	Material Purple = Material(vec3(0.4314, 0.0589, 0.6196), vec3(1.0, 1.0, 1.0),128, 1.0f, 0.0, 1.52);

	Sphere ball1(center_ball1, radius_ball1, BlueTransparent);
	Sphere ball2(center_ball2, radius_ball2, Purple);
	Sphere ball3(center_ball3, radius_ball3, DarkRed);

	sph_objects[0] = ball1;
	sph_objects[1] = ball2;
	sph_objects[2] = ball3;

	//Setting up the Plane
	Material PlaneMaterial = Material(vec3(0.1, 0.1, 0.1), vec3(1.0, 1.0, 1.0), 80.0, 0.0, 0.0, 1.52);
	Plane plane(0, 0, 1, -5, PlaneMaterial);
	pla_objects[0] = plane;


	//Setting up the Light
	vec3 lightPos(0.0, -1.5, 100.0);
	vec3 lightColor(0.5, 0.5, 0.5);//Ambient Light Color
	Light light(lightPos, lightColor);
	light_object[0] = light;

	render();
}

void callbackDisplay() {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH);
	glMatrixMode(GL_MODELVIEW);

	switch (scene) {
	case '1':
		scene1();
		break;
	case '2':
		scene2();
		break;
	case '3':
		scene3();
		break;
	case '4':
		scene4();
		break;
	default:
		printf("Please input a number between 1 and 4 (inclusive):\n./RayTracer <(1-4)>\n");
		exit(EXIT_FAILURE);
	break;

	}


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
	setRayDepth(4);
	if(argc == 1){
		printf("Please input a number between 1 and 4 (inclusive):\n./RayTracer <(1-4)>\n");
		exit(EXIT_FAILURE);
	}
	scene = argv[1][0];
	glutDisplayFunc(callbackDisplay);

	// enter GLUT event processing cycle
	// glutSwapBuffers();
	glutMainLoop();
	return 0;
}
