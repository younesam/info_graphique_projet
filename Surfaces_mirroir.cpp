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
	Sphere(const Vector& C, const Vector& albedo, double R, bool mirror = false, bool transparent = false) : C(C), albedo(albedo), R(R), mirror(mirror), transparent(transparent) {};
	Vector C;
	Vector albedo;
	double R{};
	bool mirror;
	bool transparent;

};


// Booléen permettant de déterminer si il y a intersection entre une sphère et un rayon
// et stocke le paramètre t dans param
bool intersect(const Sphere& s, const Ray& r, double& param) {
	//param = std::numeric_limits<double>::infinity();
	// equation at² + bt + c = 0
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

	bool first_intersection(const Ray& r, double& first_param, Vector& first_center, Vector& first_albedo, double& first_r, bool& miroir, bool& transparence) const {
		//first_param = std::numeric_limits<double>::infinity();
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
					miroir = sphere.mirror;
					transparence = sphere.transparent;
				}
				intersect_exists = true;
			}
		}
		return intersect_exists;
	};
};

Vector getColor(const Ray& ray, const Scene& scene, int ray_depth, const Vector& S) {
	if (ray_depth < 0) return Vector(0., 0., 0.);  // terminates recursion at some point
	double t = 100000000000;
	Vector luminosite_couleur(0., 0., 0.);
	Vector center, albedo;
	double rayon;
	bool miroir;
	bool transparent;
	double epsilon = 0.0001;
	bool intersection = scene.first_intersection(ray, t, center, albedo, rayon, miroir, transparent);
	if (intersection) {
		Vector P = ray.O + ray.u * t;
		Vector N = P - center;
		Vector P_eps = P + epsilon * N;
		N.normalize();
		if (miroir == true) {
			// handle mirror surfaces
			Vector omega_reflected = ray.u - 2 * dot(ray.u, N) * N;
			omega_reflected.normalize();
			Ray reflected_ray(P_eps, omega_reflected);
			luminosite_couleur = getColor(reflected_ray, scene, ray_depth - 1, S);
		}
		else if (transparent == true) {
			// handle transparent surfaces
			double ni = 1.0, nt = 1.5;
			double cosi = -dot(N, ray.u);
			if (cosi < 0) {
				std::swap(ni, nt);
				cosi = -cosi;
				N = 2 * N - N;
			}
			double eta = ni / nt;
			double k = 1 - sqr(eta) * (1 - sqr(cosi));
			if (k >= 0) {
				Vector omega_refracted = eta * ray.u + (eta * cosi - sqrt(k)) * N;
				omega_refracted.normalize();
				Ray refracted_ray(P - epsilon * N, omega_refracted);
				luminosite_couleur = getColor(refracted_ray, scene, ray_depth - 1, S);
			}
		}
		else {
			// handle diffuse surfaces
			Vector omega = S - P;
			omega.normalize();
			double d = sqrt((S - P).norm2());
			double intensite = 10000000000;
			Ray v(P_eps, omega);
			Sphere source_lumiere(S, Vector(1, 1, 1), d, false);
			double paramv;
			bool intersect_visibilite = intersect(source_lumiere, v, paramv);
			double visibilite;
			if (!intersect_visibilite || paramv > d) {
				visibilite = 1;
			}
			else {
				visibilite = 0;
			}
			double lumiere = intensite * max(dot(N, omega), 0.) / (4 * sqr(M_PI * d));
			Vector lumiere_couleur = albedo * lumiere;
			luminosite_couleur = visibilite * lumiere_couleur;
		}
	}
	return luminosite_couleur;
};




int main() {
	int W = 1024;
	int H = 1024;
	double fovAlpha = 60 * M_PI / 180;
	double fov = 2 * tan(fovAlpha / 2);

	Vector C(15, 0, -55);
	Vector O(0, 0, 0);
	double R = 10;
	Vector alb(0.5, 1, 0);
	Sphere s(C, alb, R, true);

	Sphere s2(Vector(-15, 0, -55), alb, R, false, true);

	// Source de lumière S
	Vector Source(-10, 20, -15);

	// Scène
	Scene scene;

	// create a red sphere centered at (0, 1000, 0) with radius 940
	Sphere red_sphere(Vector(0, 1000 + 50, 0), Vector(1, 0, 0), 1000, false);
	scene.spheres.push_back(red_sphere);

	// create a green sphere centered at (0, 0, 1000) with radius 940
	Sphere green_sphere(Vector(0, 0, -1000 - 50), Vector(0, 1, 0), 1000, false);
	scene.spheres.push_back(green_sphere);

	// create a blue sphere centered at (0, -1000, 0) with radius 990
	Sphere blue_sphere(Vector(0, -1000 - 10, 0), Vector(0, 0, 1), 1000, false);
	scene.spheres.push_back(blue_sphere);

	// create a yellow sphere centered at (0, 0, -1000) with radius 940
	Sphere yellow_sphere(Vector(0, 0, 1000 + 50), Vector(1, 1, 0), 1000, false);
	scene.spheres.push_back(yellow_sphere);

	// create a pink  sphere centered at (1000, 0, 0) with radius 940
	Sphere pink_sphere(Vector(1000 + 25, 0, 0), Vector(1, 0, 1), 1000, false);
	scene.spheres.push_back(pink_sphere);

	// create a purple  sphere centered at (-1000, 0, 0) with radius 940
	Sphere purple_sphere(Vector(-1000 - 25, 0, 0), Vector(0.5, 0, 1), 1000, false);
	scene.spheres.push_back(purple_sphere);

	scene.spheres.push_back(s);
	scene.spheres.push_back(s2);

	std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			Vector u(j - W / 2 + 0.5, -i + H / 2 - 0.5, -W / fov);
			u.normalize();
			Ray r(O, u);

			//double t;
			//bool intersection = intersect(s, r, t);


			//printf(intersection);

			// Correction Gamma
			double gamma = 2.2;

			Vector lumiere_couleur = getColor(r, scene, 100, Source);

			image[(i * W + j) * 3 + 0] = min(std::pow(lumiere_couleur[0], 1 / gamma), 255.);   // RED
			image[(i * W + j) * 3 + 1] = min(std::pow(lumiere_couleur[1], 1 / gamma), 255.);  // GREEN
			image[(i * W + j) * 3 + 2] = min(std::pow(lumiere_couleur[2], 1 / gamma), 255.);  // BLUE
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}