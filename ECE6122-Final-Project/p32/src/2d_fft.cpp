// 2D FFT using MPI- Danielson-Lanczos Algorithm

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <mpi.h>
#include <math.h>
#include <chrono>
#include <ctime>

#include "complex.h"
#include "input_image.h"
#include "complex.cc"
#include "input_image.cc"

//const float PI = 3.14159265358979f;

using namespace std;

void FourierTransform1D(Complex* tDomainInput, int width, Complex* opFFT);
void InverseFourierTransform1D(Complex* freqDomainInput, int width, Complex* opInverseFFT);
void CollectData(Complex* ipData, int numProc, int procRank, int totalSize);
void TransposeMatrix(int rank, Complex* inputData, Complex* outputData, int widthMatrix, int heightMatrix);
void DistributeData(Complex* ipData, int numProc, int procRank, int totalSize);
void FourierTransform2D(const char* ipFileName, string typeFT, char* outputFileName);
void separateInput(Complex* rowArray, int N);
void fftFunc(Complex* imageData, int N);

void FourierTransform2D(const char* ipFileName, string typeFT, char* outputFileName) {
    auto start = std::chrono::system_clock::now();
    double start_MPI = MPI_Wtime();
    
    InputImage ipimage(ipFileName); //Pass filename into InputImage constructor
    
    int numTasks, myRank;

    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    printf("Process Number: %d / %d \n", myRank, numTasks);

    int image_w = ipimage.get_width();
    int image_h = ipimage.get_height();
    int matrixSize = image_w*image_h;

    Complex* image_data = ipimage.get_image_data();

    Complex* fft1D_result = new Complex[matrixSize]; //output1D
    Complex* transpose1DMatrix = new Complex[matrixSize];
    Complex* fft2D_result = new Complex[matrixSize];
    Complex* transpose2DMatrix = new Complex[matrixSize];

    int rowPerRank = image_h/numTasks;
    int from = rowPerRank*myRank; //starting 
    int to= from + rowPerRank;
    if (typeFT == "forward") 
    {
        //1D FFT
        for (int rowItr = from; rowItr < to; rowItr++) { //iterate through rows
        // copy the row into the buffer

            fftFunc(image_data+(rowItr*image_w), image_w);
        }
        CollectData(image_data, numTasks, myRank, matrixSize);
        if(myRank == 0)
        {
            ipimage.save_image_data("After1D.txt", image_data, image_w, image_h);
            cout << "Created After1D.txt: First Fourier 1D Transform" << endl;
        }
        //2D FFT Process
        TransposeMatrix(myRank, image_data, transpose1DMatrix, image_w, image_h);
        DistributeData(transpose1DMatrix, numTasks, myRank, matrixSize);

        for (int rowitr = from; rowitr < to; rowitr++)
        {
            fftFunc(transpose1DMatrix+(rowitr*image_w), image_w);
        }

        CollectData(transpose1DMatrix, numTasks, myRank, matrixSize);
        TransposeMatrix(myRank, transpose1DMatrix, transpose2DMatrix, image_w, image_h);
        if(myRank == 0)
        {
            ipimage.save_image_data(outputFileName, transpose2DMatrix, image_w, image_h);
            cout << "Created Output file: Second Fourier 1D Transform" << endl;
        }
    }

    if (myRank == 0) {
        double end_MPI = MPI_Wtime();
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        cout << "MPI: Elapsed time: " << end_MPI-start_MPI << "s\n";
        cout << "Chrono: Elapsed time: " << elapsed_seconds.count() << endl;
    }
    // if(typeFT == "reverse") 
    // {
    //     Complex* firstInverse = new Complex[matrixSize]; //1D FT
    //     Complex* firstInverseTranspose = new Complex[matrixSize]; // 1D FT transpose
    //     Complex* secondInverse = new Complex[matrixSize]; //2D FT after transpose, final result
    //     Complex* secondInverseTranspose = new Complex[matrixSize]; // 2D FT 

    //     for (int ii = 0; ii < rowPerRank; ii++) {
    //         InverseFourierTransform1D(transpose2DMatrix + ((rowStart+ii)*image_w), image_w, firstInverse+ii*image_w);
    //     }
    //     CollectData(firstInverse, numTasks, myRank, matrixSize);
    //     if(myRank == 0)
    //     {
    //         ipimage.save_image_data("After1D_inverse.txt", firstInverse, image_w, image_h);
    //         cout << "Created After1D_inverse.txt: First Fourier 1D Inverse Transform" << endl;
    //     }
    //     printf("After 1D,before transpose\n");
    //     TransposeMatrix(myRank, firstInverse, firstInverseTranspose, image_w, image_h);
    //     printf("After transpose\n");
    //     DistributeData(firstInverseTranspose, numTasks, myRank, matrixSize);

    //     for( int jj = 0; jj < rowPerRank; jj++) {
    //         InverseFourierTransform1D(firstInverseTranspose+((rowStart+jj)*image_w), image_w, secondInverseTranspose+jj*image_w);
    //     }
    //     CollectData(secondInverseTranspose, numTasks, myRank, matrixSize);
    //     TransposeMatrix(myRank, secondInverseTranspose, secondInverse, image_w, image_h);
    //     printf("Before saying inverse data\n");
    //     if (myRank == 0) {
    //         ipimage.save_image_data("InverseFT.txt", secondInverse, image_w, image_h);
    //         cout << "Created Output file: Second Fourier 1D Inverse Transform" << endl;
    //     }
    // }
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
}

void fftFunc(Complex* X, int N) {
    if (N<2) {

    }
    else {
        double angle;
        Complex W;
        separateInput(X, N);
        fftFunc(X, N/2); //even numbers
        fftFunc(X+N/2, N/2); //odd numbers
        for (int k=0; k<N/2; k++) {
            Complex e = X[k];  //even
            Complex o = X[k+N/2];
            angle= 2*PI*k/N;
            W = Complex(cos(angle), -sin(angle));
            X[k] = e + W*o;
            X[k+N/2] = e-W*o;
        }
    }
}

void separateInput(Complex* in, int n) {
    Complex* buff = new Complex[n/2]; //temp heap
    int l, c,z;
    for (l=0; l<n/2; l++) //copy odd elements into first half buff[]
    buff[l] = in[l*2+1];
    for (c=0; c<n/2; c++) //copy even to lower half of in[]
    in[c] = in[c*2];
    for (z=0; z<n/2; z++)
    in[z+n/2] = buff[z];
    delete[] buff;
    //in[] has all the even first, then odd elements
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
    FourierTransform2D(filename.c_str(), typeFT, outputFileName);
    
    MPI_Finalize();
    
}
