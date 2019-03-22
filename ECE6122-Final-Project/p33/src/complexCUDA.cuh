#include <cmath>
#include <cuda.h>
#include <iostream>

class Complex {

public:
    float real;
    float imag;

    __device__ __host__ Complex();
    __device__ __host__ Complex(float r, float i);
    __device__ __host__ Complex(float r);
    __device__ __host__ Complex operator+(const Complex& b) const;
    __device__ __host__ Complex operator-(const Complex& b) const;
    __device__ __host__ Complex operator*(const Complex& b) const;

    Complex mag() const;
    Complex angle() const;
    Complex conj() const;


};

std::ostream& operator<< (std::ostream& os, const Complex& rhs);