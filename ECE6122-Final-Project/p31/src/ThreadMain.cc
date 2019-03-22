#include <thread>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include "input_image.h"
#include "complex.h"
#include <cmath>

using namespace std;

const float PI = 3.14159265358979f;

void rowComputation(int threadID, bool forward, int dftPerThread, int width, int height, Complex const *firstArray, Complex *secondArray)
{
	// A temporary storage of complex values 
	Complex tempComplex;

	// determine the starting row for this thread
	const int start = threadID * dftPerThread;

	// iterate through the array
	for(int i = start; i < (start+dftPerThread); i++)
	{
		for(int j = 0; j < width; j++)
		{
			// this complex variable will store the sum of complex variables
			Complex sumComplex;
			// iterate through the current row to sum up the complex variables 
			for(int n = 0; n < width; n++)
			{
				tempComplex.real = cos( (2 * PI * n * j) / width );
				tempComplex.imag = forward ? (- sin( (2 * PI * n * j) / width )) : ( sin( (2 * PI * n * j) / width ) );
				sumComplex = sumComplex + ( firstArray[n + i * width] * tempComplex );
			}
			// store the result in an array
			secondArray[j + i * width] = forward ? (sumComplex) : ( sumComplex * ( 1.0f / (float)width ) );
		}
	}
}	

void transposeMatrixThread(Complex *matrix, Complex *tempMatrix, int width, int height, int threadID, int elementsPerThread)
{
	// determine the start row for the thread
	const int startRow = threadID * elementsPerThread;

        // transpose and store the result in a temp matrix
        for( int i = startRow; i < (startRow+elementsPerThread); i++ )
        {
                for( int j = 0;  j < width; j++ )
                {
                        tempMatrix[i + j * height] = matrix[j + i * width];
                }
        }
}

void transposeMatrix(Complex *matrix, Complex *tempMatrix, int width, int height, int threadCount)
{
	int elementsPerThread = height / threadCount;

	thread *threadArray = new thread[threadCount];

	for(int i = 0; i < (threadCount-1); i++)
		threadArray[i] = thread(transposeMatrixThread, matrix, tempMatrix, width, height, i+1, elementsPerThread);

	transposeMatrixThread(matrix, tempMatrix, width, height, 0, elementsPerThread);

	for(int i = 0; i < (threadCount-1); i++)
		threadArray[i].join();

	delete[] threadArray;
}

int main(int argc, char **argv)
{
	// start the time measurement
	auto start = chrono::system_clock::now();

	// check for the right amount of arguments to the program
	if( argc != 4 )
	{
		cerr << "Wrong amount of arguments.";
		return -1;
	}

	// read first argument to determine if forward or reverse method
	bool forwardMethod;
	string direction = argv[1];
	if( direction == "forward" )
		forwardMethod = true;
	else
		forwardMethod = false;

	// read in the second argument
	InputImage input_image(argv[2]);

	int maxThreads = 8;	
	int width = input_image.get_width();
	int height = input_image.get_height();
	int totalMatrixElements = width * height;
	
	int threadCount;
	if( width == 2 )
                threadCount = 2;
        else if( width == 4 )
                threadCount = 4;
        else
                threadCount = 8;
	
	Complex *originalSampleArray = input_image.get_image_data();
	Complex *halfCompletedDFTArray = new Complex[totalMatrixElements];
	Complex *completedDFTArray = new Complex[totalMatrixElements];
	Complex *transposedArray = new Complex[totalMatrixElements];
	Complex *completeTransposedArray = new Complex[totalMatrixElements];
	thread threadArray[maxThreads];
	int dftPerThread = width / threadCount;
	
	// start the time for the DFT computation
	auto dft_start = chrono::system_clock::now();

	// row computation
	for( int i = 0; i < threadCount-1; i++ )
		threadArray[i] = thread(rowComputation,i+1,forwardMethod,dftPerThread, width, height, originalSampleArray, halfCompletedDFTArray);
	
	// this current thread should do the row computation as well
	rowComputation(0, forwardMethod, dftPerThread, width, height, originalSampleArray, halfCompletedDFTArray);

	// wait for all threads to complete before beginning column computation
	for( int i = 0; i < threadCount-1; i++ )
		threadArray[i].join();

	// transpose
	transposeMatrix(halfCompletedDFTArray, transposedArray, width, height, threadCount);

	// column computation
	for( int i = 0; i < threadCount-1; i++ )
		threadArray[i] = thread(rowComputation, i+1, forwardMethod, dftPerThread, width, height, transposedArray, completeTransposedArray);

	rowComputation(0, forwardMethod, dftPerThread, width, height, transposedArray, completeTransposedArray);
	
	// wait for all threads to finish
	for( int i = 0; i < threadCount-1; i++ )
		threadArray[i].join();

	// transpose
	transposeMatrix(completeTransposedArray, completedDFTArray,width, height, threadCount);

	auto dft_end = chrono::system_clock::now();

	// write the DFT result to file
	input_image.save_image_data(argv[3], completedDFTArray, width, height );
	
	// deallocate dynamic memory
	delete[] halfCompletedDFTArray;
	delete[] completedDFTArray;
	delete[] transposedArray;
	delete[] completeTransposedArray;
	
	auto end = chrono::system_clock::now();
	chrono::duration<float> elapsed_seconds = end - start;
	chrono::duration<float> dft_time_seconds = dft_end - dft_start;
	cout << "Finished C++ Thread computation. Total Elapsed Time: " << elapsed_seconds.count() << "s\tDFT Elapsed Time: " << dft_time_seconds.count() << "s\n";

	return 0;
}
