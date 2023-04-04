#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <iostream>
#include <xlocale>
#include <cstdio>
#include <cmath>
#include <corecrt_math_defines.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

static inline double sqr(double x) { return x * x; }

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}
	double& operator[](int i) { return coord[i]; }
	double operator[](int i) const { return coord[i]; }

	Vector& operator+=(const Vector& v) {
		coord[0] += v[0];
		coord[1] += v[1];
		coord[2] += v[2];
		return *this;
	}

	double norm2() const {
		return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
	}
	void normalize() {
		double vecNorm = sqrt(sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]));
		coord[0] /= vecNorm;
		coord[1] /= vecNorm;
		coord[2] /= vecNorm;
	}

	double coord[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const Vector& a, double b) {
	return Vector(a[0] * b, a[1] * b, a[2] * b);
}
Vector operator*(double a, const Vector& b) {
	return Vector(a * b[0], a * b[1], a * b[2]);
}

double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}


class Ray {
public:
	Ray(const Vector& O, const Vector& u) : O(O), u(u) {};
	Vector O;
	Vector u;
};


class Sphere {
public:
	Sphere(const Vector& C, const Vector& albedo, double R) : C(C), albedo(albedo), R(R) {};
	Vector C;
	Vector albedo;
	double R{};

};


// Bool�en permettant de d�terminer si il y a intersection entre une sph�re et un rayon
// et stocke le param�tre t dans param
bool intersect(const Sphere& s, const Ray& r, double& param) {
	param = std::numeric_limits<double>::infinity();
	// equation at� + bt + c = 0
	double a = 1;
	Vector Inter = r.O - s.C;
	double b = dot(r.u, Inter);
	double c = (Inter).norm2() - s.R * s.R;
	double delta = (sqr(b) - a * c);
	double param1;
	double param2;
	if (delta < 0) {
		return false;
	}
	else {
		param1 = dot(r.u, s.C - r.O) - sqrt(delta);
		param2 = dot(r.u, s.C - r.O) + sqrt(delta);
		if (param2 < 0) {
			return false;
		}
		else if (param1 >= 0) {
			param = param1;
			return true;
		}
		else {
			param = param2;
			return true;
		}
	};
};

class Scene {
public:
	vector<Sphere> spheres;

	Scene() {}

	bool first_intersection(const Ray& r, double& first_param, Vector& first_center, Vector& first_albedo, double& first_r) {
		first_param = std::numeric_limits<double>::infinity();
		bool intersect_exists = false;
		for (const auto& sphere : spheres) {
			double t;
			bool intersection = intersect(sphere, r, t);
			if (intersection) {
				if (t < first_param) {
					first_param = t;
					first_center = sphere.C;
					first_albedo = sphere.albedo;
					first_r = sphere.R;
				}
				intersect_exists = true;
			}
		}
		return intersect_exists;
	};
};



//Vector Scene::getColor(const Ray& ray, int ray_depth) {
//	if (ray_depth < 0) return Vector(0., 0., 0.);  // terminates recursion at some point
//	if (intersect(ray, P, N, sphere_id)) {
//		if (spheres[sphere_id].mirror) {
//			Ray reflected_ray = ....;
//			return getColor(reflected_ray, ray_depth - 1);
//		}
//		else {
//			// handle diffuse surfaces
//		}
//	}
//}



int main() {
	int W = 512;
	int H = 512;
	double fovAlpha = 60 * M_PI / 180;
	double fov = 2 * tan(fovAlpha / 2);
	double epsilon = 0.0001;

	Vector C(0, 0, -55);
	Vector O(0, 0, 0);
	double R = 10;
	Vector alb(0.5, 0.5, 0.5);
	Sphere s(C, alb, R);

	// Source de lumi�re
	Vector source_lumiere(20, 20, -20);
	double intensite = 10000000;

	// Sc�ne
	Scene scene;

	// create a red sphere centered at (0, 1000, 0) with radius 940
	Sphere red_sphere(Vector(0, 1000, 0), Vector(1, 0, 0), 940);
	scene.spheres.push_back(red_sphere);

	// create a green sphere centered at (0, 0, 1000) with radius 940
	Sphere green_sphere(Vector(0, 0, 1000), Vector(0, 1, 0), 940);
	scene.spheres.push_back(green_sphere);

	// create a blue sphere centered at (0, -1000, 0) with radius 990
	Sphere blue_sphere(Vector(0, -1000, 0), Vector(0, 0, 1), 990);
	scene.spheres.push_back(blue_sphere);

	// create a yellow sphere centered at (0, 0, -1000) with radius 940
	Sphere yellow_sphere(Vector(0, 0, -1000), Vector(1, 1, 0), 940);
	scene.spheres.push_back(yellow_sphere);

	scene.spheres.push_back(s);

	std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			Vector u(j - W / 2 + 0.5, -i + H / 2 - 0.5, -W / fov);
			u.normalize();
			Ray r(O, u);

			//double t;
			//bool intersection = intersect(s, r, t);

			double first_param = 100000000000000000000.;
			Vector first_center;
			Vector first_albedo;
			double first_r;
			bool first_intersect = scene.first_intersection(r, first_param, first_center, first_albedo, first_r);

			Vector P = O + u * first_param;
			Vector N = P - first_center;
			N.normalize();
			Vector omega = source_lumiere - P;
			omega.normalize();
			double d = sqrt((source_lumiere - P).norm2());

			// Calcul de la visibilit�
			//Ray v(P, omega);
			// On lance le rayon v pas de P mais d'un point un peu �lev� de P suivant N
			Vector P_eps = P + N * epsilon;
			Ray v(P_eps, omega);

			//Sphere source_lumiere(S, Vector (1, 1, 1), d);

			double paramv = 0.;
			bool intersect_visibilite = false;
			for (const auto& sphere : scene.spheres) {
				double t1;
				if (intersect(sphere, v, t1)) {
					paramv = t1;
					intersect_visibilite = true;
				}
			}
			double visibilite;
			if (!intersect_visibilite || paramv > d) {
				visibilite = 1;
			}
			else {
				visibilite = 0;
			};

			double luminosite = intensite * visibilite * max(dot(N, omega), 0.) / (4 * sqr(M_PI * d));
			Vector luminosite_couleur = first_albedo * luminosite;
			//printf(intersection);

			// Correction Gamma
			//double gamma = 2.2;

			image[(i * W + j) * 3 + 0] = first_intersect ? min(luminosite_couleur[0], 255.) : 0;   // RED
			image[(i * W + j) * 3 + 1] = first_intersect ? min(luminosite_couleur[1], 255.) : 0;  // GREEN
			image[(i * W + j) * 3 + 2] = first_intersect ? min(luminosite_couleur[2], 255.) : 0;  // BLUE
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}