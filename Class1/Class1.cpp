// Class1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "omp.h"
#include <xmmintrin.h>
#include <emmintrin.h>
#include <cmath>
#include <cstdio>

using namespace std;

// V - 6
// (y - x)sqrt(xz^2 + xy)
void init(float* x, float* y, float* z, int size) {
    float number = -100;
    for (int i = 0; i < size; i++) {
        x[i] = number;
        number++;
    }
}

void ComputeArrayCpp(float* xArray, float* yArray, float* zArray, float* pResult, int nSize) {
    int i;
    float* xSource = xArray;
    float* ySource = yArray;
    float* zSource = zArray;
    float* pDest = pResult;
    for (i = 0; i < nSize; i++) {
        *pDest - (float)( ((*ySource) -  (*xSource)) * sqrt( (*xSource) * (*zSource) * (*zSource) + (*xSource) * (*ySource) ) ); // Your function
        pDest++; // Move to next element of Source++
    }
}

void ComputeArrayCppSSE(float* xArray, float* yArray, float* zArray, float* pResult, int nSize) {
    int nLoop = nSize / 4;
    // Float to m128
    __m128 m1;  // y - x
    __m128 m2; // z * z
    __m128 m3; // x * m2
    __m128 m4; // x * y
    __m128 m5; // m3 + m4
    __m128 m6; // sqrt(m5)
    __m128 m7; // m1 * m6
    __m128* xSrc = (__m128*)xArray;
    __m128* ySrc = (__m128*)yArray;
    __m128* zSrc = (__m128*)zArray;
    __m128* pDest = (__m128*)pResult;

    for (int i = 0; i < nLoop; i++) {
        m1 = _mm_sub_ps(*ySrc, *xSrc);
        m2 = _mm_mul_ps(*zSrc, *zSrc);
        m3 = _mm_mul_ps(*xSrc, m2);
        m4 = _mm_mul_ps(*xSrc, *ySrc);
        m5 = _mm_add_ps(m3, m4);
        m6 = _mm_sqrt_ps(m5);
        m7 = _mm_mul_ps(m1, m6);
        
        *pDest = m7;
        xSrc++;
        ySrc++;
        zSrc++;
        pDest++;
    }
}


int main()
{
    const int MAX_SIZE = 201;
    // Skip 4float in array
    float* x = (float*)_aligned_malloc(sizeof(float) * MAX_SIZE, 16);
    float* y = (float*)_aligned_malloc(sizeof(float) * MAX_SIZE, 16);
    float* z = (float*)_aligned_malloc(sizeof(float) * MAX_SIZE, 16);
    float* result = (float*)_aligned_malloc(sizeof(float) * MAX_SIZE, 16);
    double startTime, endTime;

    init(x, y, z, MAX_SIZE); 

    startTime = omp_get_wtime();
    ComputeArrayCpp(x, y, z, result, MAX_SIZE);    // Compute using C++ functions
    endTime = omp_get_wtime();
    printf_s("C++: %.16g\n", endTime - startTime);

    startTime = omp_get_wtime();
    ComputeArrayCppSSE(x, y, z, result, MAX_SIZE);    // Compute using SSE
    endTime = omp_get_wtime();
    printf_s("SEE: %.16g\n", endTime - startTime);

    _aligned_free(x);
    _aligned_free(y);
    _aligned_free(z);
    _aligned_free(result);
    return 0;
}