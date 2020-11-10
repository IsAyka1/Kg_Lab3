#pragma once
#include <cmath>
struct Vec {
	double x = 0;
	double y = 0;
	double z = 0;
	Vec() = default;
	Vec(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {};
	Vec operator+(const Vec& a) const {
		return { this->x + a.x, this->y + a.y, this->z + a.z };
	}
	Vec operator*(const Vec& a) const {
		return { this->x * a.x, this->y * a.y, this->z * a.z };
	}
	Vec operator*(double a) const {
		return { this->x * a, this->y * a, this->z * a };
	}
	Vec operator-(const Vec& a) const {
		return { this->x - a.x, this->y - a.y, this->z - a.z };
	}
	Vec operator/(const Vec& a) const {
		return { this->x / a.x, this->y / a.y, this->z / a.z };
	}
	Vec operator/(double a) const {
		return { this->x / a, this->y / a, this->z / a };
	}
	Vec operator*=(double a) {
		return { this->x *= a, this->y *= a, this->z *= a };
	}
	Vec operator+=(const Vec& a) {
		return { this->x += a.x, this->y += a.y, this->z += a.z };
	}
	double* ToArray() const{
		return (double*)this;
	}
	Vec Normaliz() const {
		return { this->x / Len(), this->y / Len(), this->z / Len() };
	}
	double Len() const {
		return sqrt(pow(this->x, 2) + pow(this->y, 2) + pow(this->z, 2));
	}
	double Scal(const Vec& a) const {
		return this->x * a.x + this->y * a.y + this->z * a.z;
	}

};