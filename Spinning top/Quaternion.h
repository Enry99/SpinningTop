#ifndef QUATERNION_H
#define QUATERNION_H

#include <iostream>
#include <cmath>
#include <stdexcept>

class Quaternion {

    double a, b, c, d;

public:

    Quaternion() : a(0.), b(0.), c(0.), d(0.) {}

    Quaternion(double a, double b = 0., double c = 0., double d = 0.) : a(a), b(b), c(c), d(d) {}

    //constructor for vector-representing quaternion
    Quaternion(const double * const & values) : a(0.), b(values[0]), c(values[1]), d(values[2]) {}

    Quaternion operator +(const Quaternion& quat) const {
        return Quaternion(a + quat.a, b + quat.b, c + quat.c, d + quat.d);
    }

    Quaternion& operator +=(const Quaternion& quat) {
        a += quat.a, b += quat.b, c += quat.c, d += quat.d;
        return *this;
    }

    Quaternion operator -(const Quaternion& quat) const {
        return Quaternion(a - quat.a, b - quat.b, c - quat.c, d - quat.d);
    }

    Quaternion& operator -=(const Quaternion& quat) {
        a -= quat.a, b -= quat.b, c -= quat.c, d -= quat.d;
        return *this;
    }

    Quaternion operator +(double number) const {
        return Quaternion(a + number, b, c, d);
    }

    Quaternion& operator +=(double number) {
        a += number;
        return *this;
    }

    Quaternion operator -(double number) const {
        return Quaternion(a - number, b, c, d);
    }

    Quaternion& operator -=(double number) {
        a -= number;
        return *this;
    }

    Quaternion operator -() const {//unary, equals *(-1)
        return Quaternion(-a, -b, -c, -d);
    }

    Quaternion operator *(const Quaternion& quat) const {
        return Quaternion(
            a * quat.a - b * quat.b - c * quat.c - d * quat.d,
            a * quat.b + quat.a * b + c * quat.d - d * quat.c,
            a * quat.c - b * quat.d + c * quat.a + d * quat.b,
            a * quat.d + b * quat.c - c * quat.b + d * quat.a);
    }

    Quaternion& operator *=(const Quaternion& quat) {
        Quaternion result(
            a * quat.a - b * quat.b - c * quat.c - d * quat.d,
            a * quat.b + quat.a * b + c * quat.d - d * quat.c,
            a * quat.c - b * quat.d + c * quat.a + d * quat.b,
            a * quat.d + b * quat.c - c * quat.b + d * quat.a);
        *this = result;
        return *this;
    }

    Quaternion operator *(double number) const {
        return Quaternion(a * number, b * number, c * number, d * number);
    }

    Quaternion& operator *=(double number) {
        a *= number, b *= number, c *= number, d *= number; return *this;
    }

    Quaternion operator / (double number) const {
        return Quaternion(a / number, b / number, c / number, d / number);
    }

    Quaternion& operator /=(double number) {
        a /= number, b /= number, c /= number, d /= number; return *this;
    }

    double operator [] (int n) {
        switch (n) {
        case 0: return a;
        case 1: return b;
        case 2: return c;
        case 3: return d;
        default: throw std::runtime_error("Out of range in Quaternion operator []");
        }
    }

    explicit operator bool() const { return !(a == 0 && b == 0 && c == 0 && d == 0); }

    //friend functions
    friend double modulus(const Quaternion&);
    friend double modulus_squared(const Quaternion&);
    friend Quaternion conjugate(const Quaternion&);
    friend Quaternion inverse(const Quaternion&);
    friend Quaternion rotation(const Quaternion&, const Quaternion&, const double&);// ruota attorno a quaternione f il vettore e
    friend Quaternion normalize(Quaternion);
    friend Quaternion quaternionFromVector(const double* const&);
    friend std::ostream& operator << (std::ostream& o, const Quaternion&);
};

// Friend definitions
std::ostream& operator << (std::ostream& o, const Quaternion& quat) {
    return o << quat.a << " + (" << quat.b << ")i + (" << quat.c << ")j + (" << quat.d << ")k";
}

double modulus(const Quaternion& quat) {

    return sqrt(quat.a * quat.a + quat.b * quat.b + quat.c * quat.c + quat.d * quat.d);
}

double modulus_squared(const Quaternion& quat) {

    return quat.a * quat.a + quat.b * quat.b + quat.c * quat.c + quat.d * quat.d;
}

Quaternion conjugate(const Quaternion& quat) {
    return Quaternion(quat.a, -quat.b, -quat.c, -quat.d);
}

Quaternion inverse(const Quaternion& quat) {//inverse/modulus_squared
    return Quaternion(quat.a, -quat.b, -quat.c, -quat.d) / (quat.a * quat.a + quat.b * quat.b + quat.c * quat.c + quat.d * quat.d);
}

Quaternion normalize(Quaternion quat) {
    return quat / modulus(quat);
}
void normalize_this(Quaternion& quat) {
    quat = quat / modulus(quat);
}

//rotates quaternion "rotant" around the axis defined by the 3D vector-like quaternion "axis"
Quaternion rotation(const Quaternion& axis, const Quaternion& rotant) {
    Quaternion normalized_axis = normalize(axis);
    return normalized_axis * rotant * conjugate(normalized_axis);
}

//rotation angle in radians
Quaternion rotation(const Quaternion& axis, const Quaternion& rotant, const double& angle) {

    Quaternion normalized_axis = normalize(axis);
    double alpha = angle / 2;

    Quaternion q(
        cos(alpha),
        normalized_axis.b * sin(alpha),
        normalized_axis.c * sin(alpha),
        normalized_axis.d * sin(alpha)
    );

    return q * rotant * inverse(q);
}

Quaternion quaternionFromVector(const double* const & v) { //builds a quaternion from vector v in the form (0, x, y, z)
    Quaternion q;
    double modulus = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (modulus == 0)  return q;
    q.a = cos(modulus / 2);
    q.b = sin(modulus / 2) / modulus * v[0];
    q.c = sin(modulus / 2) / modulus * v[1];
    q.d = sin(modulus / 2) / modulus * v[2];
    return q;
}

#endif