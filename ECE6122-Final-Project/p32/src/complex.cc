//
// Created by brian on 11/20/18.
//

#include <cmath>
#include <math.h>
#include <iostream>
#include <string>

#include "complex.h"

const float PI = 3.14159265358979f;
using namespace std;

Complex::Complex() : real(0.0f), imag(0.0f) {}

Complex::Complex(float r) : real(r), imag(0.0f) {}

Complex::Complex(float r, float i) : real(r), imag(i) {}

//add complex numbers
Complex Complex::operator+(const Complex &b) const {
    return Complex(real + b.real, imag + b.imag);
}
// subtract complex numbers
Complex Complex::operator-(const Complex &b) const {
    return Complex(real - b.real, imag - b.imag);
}

Complex Complex::operator*(const Complex &b) const {
    return Complex(real*b.real-imag*b.imag, real*b.imag+imag*b.real);
}

Complex Complex::mag() const {
    return Complex(sqrt(real*real+imag*imag));
}
// tan^-1(b/a) in degrees
Complex Complex::angle() const {
    return Complex(atan2(imag,real) * 360/(2*PI));

}

Complex Complex::conj() const {
    return Complex(real, -imag);

}

void Complex::Print() const {
    float r = real;
    float i = imag;
    if (fabs(i) < 1e-10) i = 0;
    if (fabs(r) < 1e-10) r = 0;

    if (i == 0) {
        cout << real << " ";
    }
    else {
        cout << "(" << r << "," << i << ")";
    }
}

std::ostream& operator<< (std::ostream& os, const Complex& rhs) {
    Complex c(rhs);
    if(fabsf(rhs.imag) < 1e-10) c.imag = 0.0f;
    if(fabsf(rhs.real) < 1e-10) c.real = 0.0f;

    if(c.imag == 0) {
        os << c.real;
    }
    else {
        os << "(" << c.real << "," << c.imag << ")";
    }
    return os;
}