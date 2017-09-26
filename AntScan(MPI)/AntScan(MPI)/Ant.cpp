#include "Ant.h"
#include <iostream>

using namespace std;

CAnt::CAnt()
{
}

CAnt::CAnt(double alpha, double beta, CAntSets aSet)
{
	this->aSet = aSet;
	this->ALPHA = alpha;
	this->BETA = beta;
	this->rows = aSet.getPSize() + 1;
	this->model = 0;
	visited = new int[rows];
	unvisited = new int[rows];
}

CAnt::~CAnt()
{

}

void CAnt::unVisitedp(int indPnt)
{
	int i;
	int* id = new int[rows];
	for (int w = 0; w < rows; w++)
	{
		id[w] = 0;
	}
	if (indPnt == 0)
	{
		memcpy(unvisited, aSet.nbRelation[visited[0]], sizeof(int)*rows);
	}
	else
	{
		for (i = 1; i < rows; i++)
		{
			if (unvisited[i] == 1 || aSet.nbRelation[indPnt][i] == 1)
			{
				id[i] = 1;
			}
			if (visited[i] == 1)
			{
				id[i] = 0;
			}
			if (id[i] == 1)
			{
				unvisited[i] = 1;
			}
		}
	}

	delete[]id;
}

int CAnt::sumOfArr(int* arr)
{
	int j = 0;
	for (int i = 1; i < rows; i++)
	{
		j = j + arr[i];
	}
	return j;
}

double CAnt::sumOfArrDouble(double* arr)
{
	double j = 0;
	for (int i = 1; i< rows; i++)
	{
		j = j + arr[i];
	}
	return j;
}

double CAnt::LLRofPath()
{

	int i, l;
	int casAnt = 0;				//the total case value in one path
	int popAnt = 0;				//the total population value in one path
	double llr = 0.0;					//the likelihood ratio value which is calculated by the four values above
	l = rows;
	for (i = 1; i < l; i++)
	{
		if (visited[i] == 1)
		{
			casAnt = casAnt + aSet.pSet[i][0];
			popAnt = popAnt + aSet.pSet[i][1];
		}
	}
	int casAll = aSet.getCaseNum();
	int popAll = aSet.getPopNum();
	double ratioAll = (double)casAll / popAll;
	double ratioInside = (double)casAnt / popAnt;
	double ratioOutside = (double)(casAll - casAnt) / (popAll - popAnt);
	if (ratioInside <= ratioOutside)
	{
		llr = 0.0;
	}
	else
	{
		if (model == 0)
		{
			llr = casAnt * log(ratioInside) + (casAll - casAnt) * log(ratioOutside) - casAll * log(ratioAll);
		}
		else if (model == 1)
		{
			double eCaseInside = popAnt*ratioAll;
			double eCaseOutside = casAll - eCaseInside;
			llr = casAnt * log(casAnt / eCaseInside) + (casAll - casAnt) * log((casAll - casAnt) / eCaseOutside);
		}
	}
	return llr;
}

int CAnt::Roulette(double** tau)
{
	int i, j;
	double** prob = new double*[rows];
	for (int w = 0; w < rows; w++)
	{
		prob[w] = new double[rows];
	}
	double* probCand = new double[rows];
	for (int w = 0; w < rows; w++)
	{
		probCand[w] = 0;
		for (int x = 0; x < rows; x++)
		{
			prob[w][x] = 0;
		}
	}
	int  nextPnt = 0;

	for (i = 1; i < rows; i++)
	{
		if (visited[i] == 1)
		{
			for (j = 1; j < rows; j++)
			{
				if (unvisited[j] == 1)
				{
					prob[i][j] = pow(tau[i][j], ALPHA)*pow(aSet.eta[i][j], BETA);
					if (probCand[i]<prob[i][j])
					{
						prob[i][0] = j;
						probCand[0] = j;
						probCand[i] = prob[i][j];
					}
				}
			}
		}
	}

	if (sumOfArr(visited) == 1)
	{
		nextPnt = (int)probCand[0];
	}
	else
	{
		// calculate the cumulative standardization of those probCandidates, which will be used for Roulette
		double sumProb = sumOfArrDouble(probCand);
		for (i = 1; i < rows; i++)
		{
			probCand[i] = probCand[i] / sumProb;
		}
		for (i = 1; i < rows - 1; i++)
		{
			if (probCand[i] != 0)
			{
				for (j = i + 1; j < rows; j++)
				{
					if (probCand[j] != 0)
					{
						probCand[j] = probCand[i] + probCand[j];
						j = rows;
					}
				}
			}
		}

		//Roulette
		double probRand = randR;
		for (i = 1; i < rows; i++)
		{
			if (probCand[i]>probRand)
			{
				nextPnt = (int)prob[i][0];
				i = rows;
			}
		}
	}

	for (int w = 0; w < rows; w++)
	{
		delete[]prob[w];
	}
	delete[]prob;
	delete[]probCand;

	return nextPnt;
}

//ËÑÑ°
void CAnt::forage(double** tau)
{
	int newPnt = 0;
	for (int j = 2; j <= length; j++)
	{
		unVisitedp(newPnt);
		if (sumOfArr(unvisited) != 0)
		{
			newPnt = Roulette(tau);
			unvisited[newPnt] = 0;
			visited[newPnt] = 1;
		}
	}

	this->llr = LLRofPath();
}

void CAnt::update(int randPos, double** tau, int length, double randR)
{
	for (int i = 0; i < rows; i++)
	{
		visited[i] = 0;
		unvisited[i] = 0;
	}

	visited[0] = randPos;
	visited[randPos] = 1;
	this->length = length;
	this->randR = randR;
	forage(tau);
}
