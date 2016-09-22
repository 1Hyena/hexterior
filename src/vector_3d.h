#ifndef VECTOR_3D_H_15_09_2016
#define VECTOR_3D_H_15_09_2016

#include <math.h>

class VECTOR_3D {
    public:
    VECTOR_3D() : x(0.0), y(0.0), z(0.0) {}
    VECTOR_3D(float x, float y, float z) : x(x), y(y), z(z) {}
    VECTOR_3D(const VECTOR_3D &v) : x(v.x), y(v.y), z(v.z) {}
   ~VECTOR_3D() {}

    float x, y, z;

    float &operator[] (size_t i) {
        Again:
        switch (i) {
            case  0: return x;
            case  1: return y;
            case  2: return z;
            default: break;
        }
        i %= 3;
        goto Again;
    }

    VECTOR_3D& operator  = (const VECTOR_3D& v) {
        x = v.x; y = v.y; z = v.z; return *this;
    }

    VECTOR_3D& operator += (const VECTOR_3D& v) {
        x += v.x; y += v.y; z += v.z; return *this;
    }

    VECTOR_3D& operator -= (const VECTOR_3D& v) {
        x -= v.x; y -= v.y; z -= v.z; return *this;
    }

    VECTOR_3D& operator *= (float m) {
        mul(m); return *this;
    }

    VECTOR_3D& operator /= (float m) {
        mul(1.0/m); return *this; 
    }

    void mul(float);
    void normalize();
    float norm();
    static float dot_product(const VECTOR_3D&, const VECTOR_3D&);
    static VECTOR_3D cross_product(const VECTOR_3D&, const VECTOR_3D&);

    private:
};

inline VECTOR_3D VECTOR_3D::cross_product(const VECTOR_3D& a, const VECTOR_3D& b) {
    // Calculate the cross product of two vectors. This produces a normal to the
    // plane containing the operands.
    return VECTOR_3D(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

inline float VECTOR_3D::dot_product(const VECTOR_3D& a, const VECTOR_3D& b) {
    // Calculate the dot product between two vectors. This corresponds to the
    // angle between them times their lengths.
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline float VECTOR_3D::norm() {
    // Return the vector norm (length).
    return sqrt(dot_product(*this, *this));
}

inline void VECTOR_3D::normalize() {
    // Return a normalized version of the given vector.
    float s = norm();
    if (s == 0.0) return;
    mul(1.0/s);
}

inline void VECTOR_3D::mul(float s) {
    // Return a vector multiplied by a scalar.
    x*=s; y*=s; z*=s;
}

inline VECTOR_3D operator - (const VECTOR_3D& a) {
    return VECTOR_3D(-a.x, -a.y, -a.z);
}

inline VECTOR_3D operator + (const VECTOR_3D& a, const VECTOR_3D& b) {
    return VECTOR_3D(a.x+b.x, a.y+b.y, a.z+b.z);
}

inline VECTOR_3D operator - (const VECTOR_3D& a, const VECTOR_3D& b) {
    return VECTOR_3D(a.x-b.x, a.y-b.y, a.z-b.z);
}

inline VECTOR_3D operator * (const VECTOR_3D& a, const float sc) {
    return VECTOR_3D(a.x*sc, a.y*sc, a.z*sc);
}

inline float operator * (const VECTOR_3D& a, const VECTOR_3D& b) {
    return VECTOR_3D::dot_product(a, b);
}

inline VECTOR_3D operator ^ (const VECTOR_3D& a, const VECTOR_3D& b) {
    return VECTOR_3D::cross_product(a, b);
}

inline int operator == (const VECTOR_3D& a, const VECTOR_3D& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

inline int operator != (const VECTOR_3D& a, const VECTOR_3D& b) {
    return !(a == b);
}

#endif

