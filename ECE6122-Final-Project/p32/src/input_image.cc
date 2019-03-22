//
// Created by brian on 11/20/18.
//
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>

#include "input_image.h"
#include "complex.h"

using namespace std;

InputImage::InputImage(const char* filename) {
    ifstream ifs(filename);
    if(!ifs) {
        cout << "Can't open image file " << filename << endl;
        exit(1);
    }

    ifs >> w >> h;
    data = new Complex[w * h];
    for(int r = 0; r < h; ++r) {
        for(int c = 0; c < w; ++c) {
            float real;
            ifs >> real;
            data[r * w + c] = Complex((float)real);
        }
    }
}

int InputImage::get_width() const {
    return w;
}

int InputImage::get_height() const {
    return h;
}

Complex* InputImage::get_image_data() const {
    return data;
}

void InputImage::save_image_data(const char *filename, Complex *d, int w, int h) {
    ofstream ofs(filename);
    if(!ofs) {
        cout << "Can't create output image " << filename << endl;
        return;
    }

    ofs << w << " " << h << endl;
    for(int r = 0; r < h; ++r) {
        for(int c = 0; c < w; ++c) {
            ofs << d[r * w + c] << " ";
        }
        ofs << endl;
    }
}

void InputImage::save_image_data_real(const char* filename, Complex* d, int w, int h) {
    std::ofstream ofs(filename);
    if(!ofs) {
        out << "Can't create output image " << filename << endl;
        return;
    }

    ofs << w << " " << h << std::endl;

    for (int r = 0; r < h; ++r) {
        for (int c = 0; c < w; ++c) {
            ofs << d[r * w + c].real << " ";
        }
        ofs << endl;
    }
}
