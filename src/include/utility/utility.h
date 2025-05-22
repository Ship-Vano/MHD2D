//
// Created by ivan on 19.05.2025.
//

#ifndef MHD2D_UTILITY_H
#define MHD2D_UTILITY_H

#include <cmath>

struct Vec2 {double x, y;};

struct Vec3 {double x, y, z;};


Vec3 operator+(const Vec3& a, const Vec3& b);

Vec3 operator-(const Vec3& a, const Vec3& b);

double operator*(const Vec3& a, const Vec3& b);

Vec3 operator*(const Vec3& a, const double b);

Vec3 operator*(double b, const Vec3& a);

bool operator==(const Vec3&a, const Vec3&b);

Vec2 operator+(const Vec2& a, const Vec2& b);

Vec2 operator-(const Vec2& a, const Vec2& b);

double operator*(const Vec2& a, const Vec2& b);

Vec2 operator*(const Vec2& a, const double b);

Vec2 operator*(const double b, const Vec2& a);

bool operator==(const Vec2&a, const Vec2&b);

#endif //MHD2D_UTILITY_H