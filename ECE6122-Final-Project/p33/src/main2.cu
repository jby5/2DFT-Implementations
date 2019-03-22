#include <stdio.h>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <bits/stdc++.h>
#include <vector>
#include "input_image.cuh"
//#include "complex.cuh"
#include "math.h"
using namespace std;
#define pi 3.14159265358979323846

__global__ void makeCol(Complex *d_col,Complex *d_img, int which, int N){
    int k=threadIdx.x+blockIdx.x*blockDim.x;
//    d_col[k]->real=d_img[k*N+which]->real;
//    d_col[k]->imag=d_img[k*N+which]->imag;
    d_col[k]=d_img[k*N+which];
}
__global__ void makeRow(Complex *d_row,Complex *d_img, int which, int N){
    int k=threadIdx.x+blockIdx.x*blockDim.x;
        d_row[k]=d_img[which*N+k];
}



__global__ void fft(Complex *d_img, int N){ //fft for each row
    int k=threadIdx.x+blockIdx.x*blockDim.x; //which row
//    Complex** d_vec=new Complex*[N];

            //Get your column/row
//            for (int i=0;i<N;i++){
//                if (dim==1){ //columns
//                    d_vec[i]=d_img[i*N+k];
//                }else {  //rows
//                    d_vec[i]=d_img[k*N+i];
//                }
//            }

            //split even and odd halves
//            Complex temp[N/2];
//            for(int i=0; i<N/2; i++)   //copy all odd elements to upper half
//                temp[i] = column[i*2+1];
//            for(int i=0; i<N/2; i++)    // copy all even elements to lower half
//                column[i] = column[i*2];
//            for(int i=0; i<N/2; i++)    // copy all odd back to upper-half
//                column[i+N/2] = temp[i];

    int half, first, last;
    int current=N;
    int num=1;
    while(current>1) { //starting from N down to 2
        half = current/2;
        for (int i = 0; i < num; i++) {
            Complex W(cos(2 * pi / current * i), -sin(2 * pi / current * i ));
            //Jfirst=i*current;
            //Jlast=i*current+half-1;
            for(int J=i*current;J<i*current+half-1;J++) {
                Complex odd = d_img[k * N + J + half + i];
                Complex even = d_img[k * N + J + i];

                d_img[k * N + J + half + i] = even - W * odd;
                d_img[k * N + J + i] = even + W * odd;
            }
        }
        num*=2;
        Pcurrent=half;
    }

}

//__global__ void fft1d(Complex *d_vec, int N){
//    int k=threadIdx.x+blockIdx.x*blockDim.x;
//    if (N < 2) { //base case - done
//    } else {
//
//        //Complex W = exp(Complex(0, -2. * pi * k / N));
//
//        //split even and odd halves
////            Complex temp[N/2];
////            for(int i=0; i<N/2; i++)   //copy all odd elements to upper half
////                temp[i] = column[i*2+1];
////            for(int i=0; i<N/2; i++)    // copy all even elements to lower half
////                column[i] = column[i*2];
////            for(int i=0; i<N/2; i++)    // copy all odd back to upper-half
////                column[i+N/2] = temp[i];
//
//        //Recurse for each half of column/row
//        fft1d<<<1, N/2>>>(d_vec,N/2);
//        fft1d<<<1, N/2>>>(d_vec+N/2,N/2);
//
//        for (int i = 0; i < N / 2; i++) {
//            Complex W(cos(2*pi/N*i*k), -sin(2*pi/N*i*k));
//            Complex even = d_vec[i * 2];
//            Complex odd = d_vec[i * 2 + 1];
//
//            d_vec[i*2] = even + W * odd;
//            d_vec[i*2+1] = even - W * odd;
//
//        }
//
//    }
//
//
//}

__global__ void putImage(Complex *d_vec, Complex *d_img, int which,int N,int mode){
    int k=threadIdx.x+blockIdx.x*blockDim.x;
    if (mode==1) {
        d_img[k * N + which] = d_vec[k];
    } else{
        d_img[which*N+k] = d_vec[k];
    }
}

__global__ void ifft(Complex *imgvec, int N){

}


int main(int argc, char *argv[]) {
    //read inputs
    string direction = argv[1];
    InputImage input = InputImage(argv[2]);
    int N=input.get_width(); //N = width = height
    string outputStr = argv[3];
    const char *outputfile = outputStr.c_str();

    Complex *d_vec, *d_row, *d_img;
    int blockNum=ceil((float)N/1024);
    int colsize=N*sizeof(Complex);
    int imgsize=(colsize^2);
    int threadNum;
    if (N==2048){
        threadNum=1024;
    } else{
        threadNum=N;
    }
    cudaMalloc((void **)&d_vec,colsize);
    cudaMalloc((void **)&d_img,imgsize);
    cudaMemcpy(d_img,input.get_image_data(),imgsize,cudaMemcpyHostToDevice);
    //do the fft
    if (direction=="forward"){
        //columns first
        for (int i=0;i<N;i++) {
            //makeCol<<<blockNum,threadNum>>>(d_vec,d_img, i, N);
            //cudaDeviceSynchronize();
            fft<<<blockNum,threadNum>>>(d_img, N);
            //putImage<<<blockNum,threadNum>>>(d_vec,d_img,i,N,1);
            //fft<<<blockNum,threadNum>>>(d_img,N,2);
            cudaDeviceSynchronize();

        }
        //then rows
//        for (int i=0;i<N;i++) {
//            makeRow<<<blockNum,threadNum>>>(d_vec,d_img, i, N);
//            cudaDeviceSynchronize();
//            fft1d<<<blockNum,threadNum>>>(d_vec, N);
//            putImage<<<blockNum,threadNum>>>(d_vec,d_img,i,N,2);
//            //fft<<<blockNum,threadNum>>>(d_img,N,2);
//            cudaDeviceSynchronize();
//
//        }
        //save image to file
        cudaMemcpy(input.get_image_data(),d_img,imgsize,cudaMemcpyDeviceToHost);
        input.save_image_data(outputfile,input.get_image_data(),N,N);
        cout<<"FFTs saved."<<endl;

    }
    if (direction=="reverse"){
        //ifft(input.data,width);
    }

    return 0;
}