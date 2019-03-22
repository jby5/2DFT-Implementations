//
// Created by brian on 11/20/18.
//
#include "complexCUDA.cuh"
#include <cuda.h>
#include <fstream>
#include <sstream>
#include <iostream>
//class Complex;

class InputImage {
public:

    __host__ InputImage(const char* filename) {
        std::ifstream ifs(filename);
        if(!ifs) {
            std::cout << "Can't open image file " << filename << std::endl;
            exit(1);
        }

        ifs >> w >> h;
        data = new Complex[w * h];
        for(int r = 0; r < h; ++r) {
            for(int c = 0; c < w; ++c) {
                float real;
                ifs >> real;
                data[r * w + c] = Complex(real);
            }
        }
    }

    __host__ int get_width() const {
        return w;
    }

    __host__ int get_height() const {
        return h;
    }

    __host__ Complex* get_image_data() const {
        return data;
    }

    __host__ void save_image_data(const char *filename, Complex *d, int w, int h) {
        std::ofstream ofs(filename);
        if(!ofs) {
            std::cout << "Can't create output image " << filename << std::endl;
            return;
        }

        ofs << w << " " << h << std::endl;

        for(int r = 0; r < h; ++r) {
            for(int c = 0; c < w; ++c) {
                ofs << d[r * w + c] << " ";
                //std::cout<<d[r * w + c] << " ";
            }
            ofs << std::endl;
            //std::cout<<std::endl;
        }
    }

    __host__ void save_image_data_real(const char* filename, Complex* d, int w, int h) {
        std::ofstream ofs(filename);
        if(!ofs) {
            std::cout << "Can't create output image " << filename << std::endl;
            return;
        }

        ofs << w << " " << h << std::endl;

        for (int r = 0; r < h; ++r) {
            for (int c = 0; c < w; ++c) {
                ofs << d[r * w + c].real << " ";
            }
            ofs << std::endl;
        }
    }
private:
    int w;
    int h;
    Complex* data;
};

