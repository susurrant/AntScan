#pragma once
#include <vector>
#include <random>
#include <time.h>
#include <algorithm>
#include "Antsets.h"
using namespace std;

const int N_CITY_COUNT = 100;

class CAnt
{
public:
	CAnt(void);
	CAnt(double alpha, double beta, CAntSets aSet);
	~CAnt(void);

public:
	int* visited;
	int* unvisited;

	int rows;
	int length;
	int model;           //llr¼ÆËãÄ£ÐÍ
	double ALPHA;
	double BETA;
	double randR;
	double llr;

	CAntSets aSet;
public:
	void unVisitedp(int indPnt);
	int sumOfArr(int* arr);
	double sumOfArrDouble(double* arr);
	double LLRofPath();
	int Roulette(double** tau);
	void forage(double** tau);
	void update(int randPos, double** tau, int length, double randR);
};