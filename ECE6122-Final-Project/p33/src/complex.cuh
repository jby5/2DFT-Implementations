#include <cmath>

//const float PI = 3.14159265358979f;

class Complex {
    const float PI = 3.14159265358979f;
public:
    __device__ __host__ Complex() : real(0.0f), imag(0.0f) {}

    __device__ __host__ Complex(float r) : real(r), imag(0.0f) {}

    __device__ __host__ Complex(float r, float i) : real(r), imag(i) {}

//add complex numbers
    __device__ __host__ Complex operator+(const Complex &b) const {
        return Complex(real + b.real, imag + b.imag);
    }
// subtract complex numbers
    __device__ __host__ Complex operator-(const Complex &b) const {
        return Complex(real - b.real, imag - b.imag);
    }

    __device__ __host__ Complex operator*(const Complex &b) const {
        return Complex(real*b.real-imag*b.imag, real*b.imag+imag*b.real);
    }
    __device__ __host__ Complex operator=(const Complex &b) const {
        return Complex(real, imag);
    }

    __device__ __host__ Complex mag() const {
        return Complex(sqrt(real*real+imag*imag));
    }
// tan^-1(b/a) in degrees
    __device__ __host__ Complex angle() const {
        return Complex(atan2(imag,real) * 360/(2*PI));

    }

    __device__ __host__ Complex conj() const {
        return Complex(real, -imag);

    }

    float real;
    float imag;
};

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