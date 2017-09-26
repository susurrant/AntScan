#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream> 
#include <sstream>
#include "cases.h"
using namespace std;
class CAntSets
{
public:
	CAntSets();
	CAntSets(string pFileName, string nbFileName);
	~CAntSets();
private:
	// data files
	string pFileName;
	string nbFileName;

	// data contents
	int caseNum;
	int popNum;
	int pAntNum;
public:
	int** pSet;
	double** pCoords;
	int** nbRelation;
	double* ratio;
	double** eta;
private:
	bool readCase();
	bool readNeighbor();
public:
	int getCaseNum();
	int getPopNum();
	int getPSize();
	vector<CCases> ca;
};