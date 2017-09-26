#pragma once
#include <vector>
#include <random>
#include <time.h>
#include <algorithm>
#include <iostream>
#include <omp.h>
#include "Antsets.h"
#include "Ant.h"

using namespace std;
class CAntScan
{
public:
	CAntScan(CAntSets aSet);
	~CAntScan();
public:
	uniform_int_distribution<unsigned> uint;
	uniform_real_distribution<double> ureal;
	normal_distribution<> unormal;
	default_random_engine e;

	double Alpha;
	double Beta;
	double Rho;
	double Q;

	int rows;
	int ants;
	int iterations;
	int initTime;
	double** eta;
	int* clusterId;

	int** rBest;
	double* lBest;
	int** rPath;
public:
	double** initEta(CAntSets aSet);
	double** init2DArrFor1(int rows, int cols);
	int binSearch(int*Array, int start, int send, int key);
	int* randPerm(int sizeRange, int sizeArray);
	double* arrayCumSum(double** arr, int arrLen);
	int* antCluster(CAntSets aSet);
	double varianceOfArr(int* arr, int arrLen);
	double aveValue(int* arr, int arrLen);
	int maxIndex(double* arr, int arrLen);
	int sumOfArr(int* arr, int arrLen);
	bool accuCalculate(int* clusterI, CAntSets aSet);
	int* lidSort(double* llr);
	double** tauUpdate(CAnt* ant, int* lid, int upt_tau, double info);
	int normLen(double varKK, double kAve);
public:
	int* getClusterId();
};

