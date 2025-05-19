//
// Created by ivan on 19.05.2025.
//
#include "utility.h"



Vec3 operator+(const Vec3& a, const Vec3& b) {
    return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vec3 operator-(const Vec3& a, const Vec3& b) {
    return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vec3 operator*(const Vec3& a, const double b) {
    return Vec3(a.x * b, a.y * b, a.z * b);
}

Vec3 operator*(const double b, const Vec3& a) {
    return a * b;
}

double operator*(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

bool operator==(const Vec3&a, const Vec3&b){
    return (std::fabs(a.x - b.x) < 1e-16) && (std::fabs(a.y - b.y) < 1e-16) && (std::fabs(a.z - b.z) < 1e-16);
}

Vec2 operator+(const Vec2& a, const Vec2& b) {
    return Vec2(a.x + b.x, a.y + b.y);
}

Vec2 operator-(const Vec2& a, const Vec2& b) {
    return Vec2(a.x - b.x, a.y - b.y);
}

double operator*(const Vec2& a, const Vec2& b) {
    return a.x * b.x + a.y * b.y;
}

Vec2 operator*(const Vec2& a, const double b) {
    return Vec2(a.x * b, a.y * b);
}

Vec2 operator*(const double b, const Vec2& a) {
    return a * b;
}

bool operator==(const Vec2&a, const Vec2&b){
    return (std::fabs(a.x - b.x) < 1e-16) && (std::fabs(a.y - b.y) < 1e-16);
}