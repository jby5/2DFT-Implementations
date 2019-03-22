// 2D FFT using MPI- Danielson-Lanczos Algorithm

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <mpi.h>
#include <math.h>

#include "complex.h"
#include "complex.cc"
#include "input_image.h"
#include "input_image.cc"

using namespace std;

void FourierTransform1D(Complex* tDomainInput, int width, Complex* opFFT);
void InverseFourierTransform1D(Complex* freqDomainInput, int width, Complex* opInverseFFT);
void CollectData(Complex* ipData, int numProc, int procRank, int totalSize);
void TransposeMatrix(int rank, Complex* inputData, Complex* outputData, int widthMatrix, int heightMatrix);
void DistributeData(Complex* ipData, int numProc, int procRank, int totalSize);
void FourierTransform2D(const char* ipFileName, string typeFT, char* outputFileName);

void FourierTransform2D(const char* ipFileName, string typeFT, char* outputFileName) {
    
    InputImage ipimage(ipFileName); //Pass filename into InputImage constructor
    
    int numTasks, myRank;

    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    printf("Process Number: %d / %d \n", myRank, numTasks);

    int image_w = ipimage.get_width();
    int image_h = ipimage.get_height();
    int matrixSize = image_w*image_h;

    Complex* image_data = ipimage.get_image_data();

    int rowPerRank = image_h/numTasks;
    int rowStart = rowPerRank*myRank;
    if (typeFT == "forward") 
    {
        Complex* fft1D_result = new Complex[matrixSize]; //output1D
        Complex* transpose1DMatrix = new Complex[matrixSize];
        Complex* fft2D_result = new Complex[matrixSize];
        Complex* transpose2DMatrix = new Complex[matrixSize];

        //1D FFT
        for (int rowItr = 0; rowItr < rowPerRank; rowItr++) {
            FourierTransform1D(image_data+((rowItr+rowStart)*image_w), image_w, fft1D_result + image_w*rowItr);
        }
        CollectData(fft1D_result, numTasks, myRank, matrixSize);
        if(myRank == 0)
        {
            ipimage.save_image_data("After1D.txt", fft1D_result, image_w, image_h);
            cout << "Created After1D.txt: First Fourier 1D Transform" << endl;
        }
        //2D FFT Process
        TransposeMatrix(myRank, fft1D_result, transpose1DMatrix, image_w, image_h);
        DistributeData(transpose1DMatrix, numTasks, myRank, matrixSize);

        for (int rowitr = 0; rowitr < rowPerRank; rowitr++)
        {
            FourierTransform1D(transpose1DMatrix+((rowitr+rowPerRank)*image_w), image_w, fft2D_result + image_w*rowitr);
        }

        CollectData(fft2D_result, numTasks, myRank, matrixSize);
        TransposeMatrix(myRank, fft2D_result, transpose2DMatrix, image_w, image_h);
        if(myRank == 0)
        {
            ipimage.save_image_data(outputFileName, transpose2DMatrix, image_w, image_h);
            cout << "Created Output file: Second Fourier 1D Transform" << endl;
        }
    }
    
    if(typeFT == "reverse") 
    {
        Complex* firstInverse = new Complex[matrixSize]; //1D FT
        Complex* firstInverserTranspose = new Complex[matrixSize]; // 1D FT transpose
        Complex* secondInverse = new Complex[matrixSize]; //2D FT after transpose, final result
        Complex* secondInverseTranspose = new Complex[matrixSize]; // 2D FT 

        for (int ii = 0; ii < rowPerRank; ii++) {
            InverseFourierTransform1D(image_data + ((rowStart+ii)*image_w), image_w, firstInverse+ii*image_w);
        }
        CollectData(firstInverse, numTasks, myRank, matrixSize);
        if(myRank == 0)
        {
            ipimage.save_image_data("After1D_inverse.txt", firstInverse, image_w, image_h);
            cout << "Created After1D_inverse.txt: First Fourier 1D Inverse Transform" << endl;
        }
        TransposeMatrix(myRank, firstInverse, firstInverserTranspose, image_w, image_h);
        DistributeData(firstInverserTranspose, numTasks, myRank, matrixSize);

        for( int jj = 0; jj < rowPerRank; jj++) {
            InverseFourierTransform1D(firstInverserTranspose+((rowStart+jj)*image_w), image_w, secondInverseTranspose+jj*image_w);
        }
        CollectData(secondInverseTranspose, numTasks, myRank, matrixSize);
        TransposeMatrix(myRank, secondInverseTranspose, secondInverse, image_w, image_h);

        if (myRank == 0) {
            ipimage.save_image_data_real(outputFileName, secondInverse, image_w, image_h);
            cout << "Created Output file: Second Fourier 1D Inverse Transform" << endl;
        }
    }
}

void CollectData(Complex* ipData, int numProc, int procRank, int totalSize) 
{
    //Send data if not rank 0
    int error;
    if(procRank !=0)
    {
        Complex* sendData = ipData;
        MPI_Request request;
        error = MPI_Isend(sendData, totalSize*2/numProc, MPI_COMPLEX, 0, 0, MPI_COMM_WORLD, &request);
    }
    if(procRank == 0)
    {
        for( int itr = 1; itr < numProc; itr++) {
            MPI_Status status;
            Complex* getData = ipData + totalSize*itr/numProc;
            error = MPI_Recv(getData, totalSize*2/numProc, MPI_COMPLEX, itr, 0, MPI_COMM_WORLD, &status);
            if(error != MPI_SUCCESS) {
                MPI_Abort(MPI_COMM_WORLD,1);
            }
        }
    }

}

void DistributeData(Complex* ipData, int numProc, int procRank, int totalSize)
{
    int error;
    if(procRank == 0)
    {
        for (int itr2 = 1; itr2 < numProc; itr2++) {
            Complex* buffData = ipData;
            MPI_Request request;
            error = MPI_Isend(buffData, totalSize*2, MPI_COMPLEX, itr2, 0, MPI_COMM_WORLD, &request);
            }
    }
    if(procRank != 0)
    {
        MPI_Status status;
        Complex* recvData = ipData;
        error = MPI_Recv(recvData, totalSize*2, MPI_COMPLEX, 0, 0, MPI_COMM_WORLD, &status);
    }
}

void TransposeMatrix(int rank, Complex* inputData, Complex* outputData, int w, int h) 
{
    if(rank == 0)
    {
        int currindex = 0;
        for (int rows=0; rows < h; rows++) {
            for (int cols=0; cols< w; cols++) {
                outputData[currindex] = inputData[rows + cols*w];
                currindex++;
            }
        }   
    }
}

void FourierTransform1D(Complex* hInput, int width, Complex* opH) 
{
    Complex totalH;
    Complex W; 
    double angle;

    int k, n;
    for (k=0; k<width; k++) {
        for (n=0; n<width; n++) {
            angle= 2*PI*n*k/width;
            W = Complex(cos(angle), -sin(angle));
            totalH = totalH + hInput[n]*W;
        }
        opH[k] = totalH;
        totalH = Complex(0,0); 
    }
}

void InverseFourierTransform1D(Complex* h, int width, Complex* H)
{
    Complex totalh;
    Complex W;
    double angle;

    int k, n;
    for (n = 0; n < width; n++) {
        for (k = 0; k < width; k++) {
            angle= 2*PI*n*k/(float)width;
            W = Complex(cos(angle), sin(angle));
            totalh = totalh + h[k]*W;
        }
        H[n] = totalh;
        H[n].real = H[n].real/(float)width;
        H[n].imag = H[n].imag/(float)width;
        totalh = Complex(0.0f,0.0f); 
    }
}

int main(int argc, char* argv[]) 
{
    int error;
    if (argc != 4) {
        cout << "Please enter all parameters\n";
    }
    string typeFT = argv[1]; //first commandline input
    string filename(argv[2]); //input file name
    char* outputFileName(argv[3]); //output file name
    error = MPI_Init(&argc, &argv);
    if(error != MPI_SUCCESS) {
        printf("Error starting MPI. End.\n");
        MPI_Abort(MPI_COMM_WORLD, error);
    }
    //cout << "Direction " << typeFT << endl;
    FourierTransform2D(filename.c_str(), typeFT, outputFileName);

    MPI_Finalize();
}
