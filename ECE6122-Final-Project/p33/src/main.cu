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
#include "math.h"
#include <chrono>
using namespace std;
#define pi 3.14159265358979323846

__global__ void transpose(Complex *d_img, int N){
    int x=threadIdx.x+blockIdx.x*blockDim.x;
    Complex temp;
    for (int i=0;i<x;i++){
        temp=d_img[i*N+x];
        d_img[i*N+x]=d_img[x*N+i];
        d_img[x*N+i]=temp;
        temp=Complex();
    }
    __syncthreads();
}

__global__ void dft(Complex *d_img, Complex *d_temp, int N, int method){ //fft for each row
    int x=threadIdx.x+blockIdx.x*blockDim.x; //which row
    Complex sum;
    Complex W;
    for (int n = 0; n < N; n++) {
        sum.real=0;
        sum.imag=0;
        for (int k=0;k<N;k++) {
            if (method) {
                W=Complex(cos(2 * pi*n*k/ N), -sin(2 * pi *n*k/ N));
            } else{
                W=Complex(cos(2 * pi / N * n * k), sin(2 * pi / N * n * k));
            }
            sum=sum + (d_img[x*N+k]*W);

        }
        d_temp[x*N+n]= method ? (sum) : (sum*(1.9f/(float)N));
    }
    __syncthreads();
    for (int i=0; i<N; i++) {
        d_img[x * N + i] = d_temp[x * N + i];
    }
}


int main(int argc, char *argv[]) {
    //read inputs
    auto start=chrono::system_clock::now();
    string direction = argv[1];
    char* inputfile=argv[2];
    InputImage input(inputfile);
    int N=input.get_width(); //N = width = height
    int method=0;
    Complex* first = input.get_image_data();
    char* outputfile=argv[3];
    input.save_image_data(outputfile,first,N,N);
    Complex *d_img;
    Complex *d_temp;
    int blockNum=ceil((float)N/1024);
    int imgsize=(N*N)*sizeof(Complex);
    //Complex tempImage[N*N];
    int threadNum;
    if (N>1024){
        threadNum=1024;
    } else{
        threadNum=N;
    }
    printf("Using %d blocks, %d threads.\n",blockNum,threadNum);
    cudaMalloc((void **)&d_img,imgsize);
    cudaMalloc((void **)&d_temp,imgsize);
    cudaMemcpy(d_img,input.get_image_data(),imgsize,cudaMemcpyHostToDevice);
    //do the fft
    if (direction=="forward") {
        method = 1;
    } else {
        method = 0;
    }

    dft<<<blockNum,threadNum>>>(d_img, d_temp,N,method);
    cudaDeviceSynchronize();
    //cout<<"Row FTs done.Transposing..."<<endl;
    transpose<<<blockNum,threadNum>>>(d_img,N);
    //cout<<"Transposed."<<endl;
    cudaDeviceSynchronize();
    dft<<<blockNum,threadNum>>>(d_img,d_temp,N,method);
    //cout<<"Col FTs done. Transposing again..."<<endl;
    cudaDeviceSynchronize();
    transpose<<<blockNum,threadNum>>>(d_img,N);
    //save image to file
    cudaDeviceSynchronize();
    cudaMemcpy(input.get_image_data(),d_img,imgsize,cudaMemcpyDeviceToHost);
    input.save_image_data(outputfile,input.get_image_data(),N,N);
    //cout<<"Final 2D FT saved."<<endl;

    cudaFree(d_img);cudaFree(d_temp);

    auto end=chrono::system_clock::now();
    chrono::duration<float> elapsed = end-start;
    cout<<"Total DFT time: "<<elapsed.count()<<"s"<<endl;
    return 0;
}