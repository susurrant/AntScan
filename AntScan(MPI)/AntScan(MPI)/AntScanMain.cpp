#include <iostream>
#include <vector>
#include <random>
#include <time.h>
#include <algorithm>
#include <mpi.h>
#include "AntSets.h"
#include "Ant.h"
using namespace std;

double** init2DArrFor1(int rows, int cols)
{
	int i, j;
	double** arr = new double*[rows];
	for (int w = 0; w < rows; w++)
	{
		arr[w] = new double[cols];
	}
	for (i = 0; i < rows; i++)
	{
		for (j = 0; j < cols; j++)
		{
			arr[i][j] = 1;
		}
	}
	return arr;
}

int binSearch(int*Array, int start, int send, int key)
{
	int left, right;
	int mid;
	left = start;
	right = send;


	while (left <= right)
	{
		mid = (left + right) / 2;
		if (key == Array[mid])
		{
			return mid;
		}
		else if (key<Array[mid])
		{
			right = mid - 1;
		}
		else if (key>Array[mid])
		{
			left = mid + 1;
		}
	}
	return -1;
}

int* randPerm(int sizeRange, int sizeArray, uniform_int_distribution<unsigned> uint, default_random_engine e)
{
	uint = uniform_int_distribution<unsigned>(0, sizeRange - 1);
	int i = 1;
	int* randPos = new int[sizeArray];
	int* tmpArr = new int[sizeArray];
	for (int j = 0; j < sizeArray; j++)
	{
		randPos[j] = 0;
		tmpArr[j] = 0;
	}


	while (i < sizeArray)
	{
		int tmp = uint(e);
		memcpy(tmpArr, randPos, sizeof(int)*sizeArray);
		sort(tmpArr, tmpArr + sizeArray);
		if (binSearch(tmpArr, 0, sizeArray, tmp) < 0)
		{
			randPos[i] = tmp;
			i = i + 1;
		}
	}
	return randPos;
}

double* arrayCumSum(double** arr, int arrLen)
{

	int i, j;
	int len = arrLen;

	double* rowSum = new double[len];
	double* colSum = new double[len];

	for (int w = 0; w < len; w++)
	{
		rowSum[w] = 0;
		colSum[w] = 0;
	}

	for (i = 1; i < len; i++)
	{
		for (j = 1; j < len; j++)
		{
			rowSum[i] = rowSum[i] + arr[i][j];
			colSum[i] = colSum[i] + arr[j][i];
		}
	}

	double* arrSum = new double[len];
	for (i = 1; i < len; i++)
	{
		arrSum[i] = rowSum[i] + colSum[i];
	}

	double totalSum = 0;
	for (i = 1; i < len; i++)
	{
		totalSum = totalSum + arrSum[i];
	}

	double* normSum = new double[len];
	for (i = 1; i < len; i++)
	{
		normSum[i] = arrSum[i] / totalSum;
	}

	double* arrCumSum = new double[len];
	arrCumSum[1] = normSum[1];
	for (i = 2; i < len; i++)
	{
		arrCumSum[i] = arrCumSum[i - 1] + normSum[i];
	}

	delete[]rowSum;
	delete[]colSum;
	delete[]arrSum;
	delete[]normSum;

	return arrCumSum;
}

double aveValue(int* arr, int arrLen)
{
	int i;
	int arrSize = 0;
	double aveVal, sumVal = 0.0;
	for (i = 1; i < arrLen; i++)
	{
		if (arr[i] != 0)
		{
			sumVal = sumVal + arr[i];
			arrSize = arrSize + 1;
		}
	}
	aveVal = sumVal / arrSize;
	return aveVal;
}

double varianceOfArr(int* arr, int arrLen)
{
	int i;
	int arrSize = 0;
	double varVal = 0.0;
	double ave = aveValue(arr, arrLen);
	for (i = 1; i < arrLen; i++)
	{
		if (arr[i] != 0)
		{
			arrSize = arrSize + 1;
			varVal = varVal + (arr[i] - ave)*(arr[i] - ave);
		}
	}
	return varVal;
}

int maxIndex(double* arr, int arrLen)
{
	int i, ind;
	double maxTmp = arr[0];
	ind = 0;
	for (i = 1; i < arrLen; i++)
	{
		if (maxTmp <= arr[i])
		{
			maxTmp = arr[i];
			ind = i;
		}
	}
	return ind;
}

int sumOfArr(int* arr, int arrLen)
{
	int i, j = 0;
	for (i = 1; i < arrLen; i++)
	{
		j = j + arr[i];
	}
	return j;
}

bool accuCalculate(int* clusterI, CAntSets aSet, int rows)
{
	int diffNum = 0;
	int* clusterRl = new int[rows];
	for (int i = 0; i < rows; i++)
	{
		clusterRl[i] = 0;
	}
	//cout << "Miss matching pairs: " << endl;
	for (int i = 1; i < rows; i++)
	{
		if (aSet.pSet[i][1] == 20)
		{
			clusterRl[i] = 1;
		}
		if (clusterRl[i] != clusterI[i])
		{
			diffNum = diffNum + 1;
			//cout << "clusterReal[" << i << "] = " << clusterRl[i] << ", clusterFound[" << i << "] = " << clusterI[i] << endl;
		}
	}
	if (diffNum == 0)
	{
		//cout << "All match!" << endl;
		return true;
	}
	return false;
}

int* lidSort(double* llr, int ants)
{
	int i, j, k;
	int l = ants;
	int* ls = new int[l];
	for (i = 0; i < l; i++)
	{
		ls[i] = i;
	}
	double lt = 0.0;
	for (i = 1; i < l - 1; i++)
	{
		for (j = i + 1; j < l; j++)
		{
			if (llr[i] < llr[j])
			{
				lt = llr[i]; llr[i] = llr[j]; llr[j] = lt;
				k = ls[i]; ls[i] = ls[j]; ls[j] = k;
			}
		}
	}
	return ls;
}

double** tauUpdate(int** visited, int* lid, int upt_tau, double info, int rows, double Q)
{
	int i, j, l;
	int startId;
	l = rows;
	double** delta_tau = new double*[l];
	for (int w = 0; w < l; w++)
	{
		delta_tau[w] = new double[l];
		for (int x = 0; x < l; x++)
		{
			delta_tau[w][x] = 0;
		}
	}
	for (i = 1; i <= upt_tau; i++)
	{
		startId = visited[lid[i]][0];
		for (j = 1; j < l; j++)
		{
			if (visited[lid[i]][j] == 1 && lid[i] != j)
			{
				delta_tau[startId][j] = delta_tau[startId][j] + Q*pow(info, upt_tau - i);
				delta_tau[j][startId] = delta_tau[startId][j];
			}
		}
	}
	return delta_tau;
}

double maxValue(double* arr, int arrLen)
{
	int i;
	double maxTmp = arr[0];
	for (i = 1; i < arrLen; i++)
	{
		if (maxTmp <= arr[i])
		{
			maxTmp = arr[i];
		}
	}
	return maxTmp;
}

int main(int argc, char *argv[])
{
	CAntSets aSet("pSimulate.txt", "nbSimulate.txt");
	int rows = aSet.getPSize() + 1;
	int ants = rows / 2;
	double Alpha = 1.0;
	double Beta = 1.0;
	int iterations = 12;

	int tag = 97;
	int tag1 = 98;
	int tag2 = 99;
	int tag3 = 100;
	int tag4 = 101;
	int tag5 = 102;
	int tag6 = 103;
	int size, rank;

	MPI::Init(argc, argv);
	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();
	if (rank == 0)
	{
		clock_t t1, t2;
		//double start, end;
		clock_t totalClock_t = 0;
		//double totalTime = 0;
		uniform_int_distribution<unsigned> uint;
		uniform_real_distribution<double> ureal = uniform_real_distribution<double>(0, 1);
		normal_distribution<> unormal = normal_distribution<>(0, 1);
		default_random_engine e;
		e.seed(time(0));

		double Rho = 0.6;
		double Q = 1.0;
		int initTime = 5;
		int ex = 1;
		int** rBest;
		double* lBest;
		rBest = new int*[iterations];
		lBest = new double[iterations];
		for (int i = 0; i < iterations; i++)
		{
			rBest[i] = new int[rows];
			lBest[i] = 0;
			for (int j = 0; j < rows; j++)
			{
				rBest[i][j] = 0;
			}
		}

		int nochange = 0;
		int nIteration = 1;
		double** tau = init2DArrFor1(rows, rows);
		int* kk = new int[iterations];
		for (int w = 0; w < iterations; w++)
		{
			kk[w] = 0;
		}
		double varKK = 0.0;
		int kAve = 0;
		int* randPos = new int[ants];

		while (nIteration < iterations)
		{
			if (nIteration <= initTime)
			{
				randPos = randPerm(rows, ants, uint, e);
			}
			else
			{
				randPos[0] = 0;
				double* tauCumSum = arrayCumSum(tau, rows);
				for (int i = 1; i < ants; i++)
				{
					double randTau = fabs(ureal(e));

					int j = 0;
					while (tauCumSum[j] < randTau && j < rows)
					{
						j = j + 1;
					}
					randPos[i] = j;
				}
			}

			if (nIteration <= initTime)
			{
				uint = uniform_int_distribution<unsigned>(0, rows / 4);
				kAve = uint(e);
				varKK = ceil(rows*rows / 16);
			}
			else
			{
				if (lBest[nIteration - 2] == lBest[nIteration - 1])
				{
					int lBestId = maxIndex(lBest, iterations);
					kAve = rBest[lBestId][0];
					varKK = varianceOfArr(kk, iterations) / 4;
				}
				else
				{
					kAve = kk[nIteration - 1];
					varKK = varianceOfArr(kk, iterations) / 9;
				}
			}

			//start = MPI_Wtime();
			t1 = clock();
			//准备消息
			int* p = new int[ants];
			double* aR = new double[ants];
			for (int i = 1; i < ants; i++)
			{
				p[i] = randPos[i];
				aR[i] = ureal(e);
			}
			p[0] = kAve;
			aR[0] = varKK;
			int* alength = new int[ants];
			alength[0] = 0;
			for (int i = 1; i < ants; i++)
			{
				int k = -1;
				int j = 1;
				alength[i] = 0;
				while (k < 0)
				{
					j = (int)(ceil(sqrt(varKK)*unormal(e)) + kAve);
					if (j > 0 && j < (rows - 1) / 2)
					{
						k = 1;
						alength[i] = j;
					}
				}
			}
			double* tau1D = new double[rows*rows];
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < rows; j++)
				{
					tau1D[i*rows + j] = tau[i][j];
				}


			//发送消息
			for (int i = 1; i < size; i++)
			{
				MPI::COMM_WORLD.Send(&ex, 1, MPI::INT, i, tag6);    //必须最先发送，不然阻塞
				MPI::COMM_WORLD.Send(alength, ants, MPI::INT, i, tag);
				MPI::COMM_WORLD.Send(p, ants + 2, MPI::INT, i, tag1);
				MPI::COMM_WORLD.Send(aR, ants, MPI::DOUBLE, i, tag2);
				MPI::COMM_WORLD.Send(tau1D, rows*rows, MPI::DOUBLE, i, tag3);
			}

			double* llrs = new double[ants];
			llrs[0] = 0;
			int** visited = new int*[ants];
			for (int i = 0; i < ants; i++)
			{
				visited[i] = new int[rows];
			}

			//接收结果
			for (int i = 1; i < size; i++)
			{
				double* llrRecv = new double[ants];
				int* visitedRecv = new int[rows*ants];
				MPI::COMM_WORLD.Recv(llrRecv, ants, MPI::DOUBLE, MPI::ANY_SOURCE, tag4);
				MPI::COMM_WORLD.Recv(visitedRecv, ants*rows, MPI::INT, MPI::ANY_SOURCE, tag5);
				
				for (int j = 1; j < ants; j++)
				{
					if (llrRecv[j] != -1)
					{
						llrs[j] = llrRecv[j];
					}
					for (int k = 0; k < rows; k++)
					{
						if (visitedRecv[j*rows + k] == -1)
							break;
						visited[j][k] = visitedRecv[j*rows + k];
					}
				}
				delete[]llrRecv;
				delete[]visitedRecv;
			}

			//end = MPI_Wtime();
			t2 = clock();
			totalClock_t += (t2 - t1);
			//totalTime += (end - start);

			if (nIteration > 1)
			{
				if (llrs[1] < lBest[nIteration - 1])
				{
					for (int j = 1; j < rows; j++)
					{
						{
							visited[1][j] = rBest[nIteration - 1][j];
						}
					}
					llrs[1] = lBest[nIteration - 1];
				}
			}

			lBest[nIteration] = maxValue(llrs, ants);
			int ind = maxIndex(llrs, ants);
			for (int i = 1; i < rows; i++)
			{
				rBest[nIteration][i] = visited[ind][i];
			}
			rBest[nIteration][0] = sumOfArr(rBest[nIteration], rows);
			kk[nIteration] = kAve;
			//cout << "iteration[" << nIteration << "]: ant[" << ind << "]'s LLR = " << lBest[nIteration];
			//cout << ", Length of Path = " << rBest[nIteration][0] << endl;

			bool matchs = accuCalculate(rBest[nIteration], aSet, rows);
			if (matchs)
			{
				nIteration = iterations - 1;
				ex = -1;
			}

			double info = 0.0;
			if (nIteration > 1)
			{
				if (lBest[nIteration] == lBest[nIteration - 1])
				{
					nochange = nochange + 1;
				}
				else{
					nochange = 0;
				}
			}
			if (nochange == 5)
			{		
				info = 1.0;
				tau = init2DArrFor1(rows, rows);
				nochange = 0;
			}
			else
			{
				info = 1.8;
			}
			//Update the pheromones for the path of all elite ants which have LLR values above the 1/3
			int upt_tau = (int)ceil(ants / 3);
			if (upt_tau > 10)
			{
				upt_tau = 10;
			}
			int* lid = lidSort(llrs, ants);
			double** delta_tau = tauUpdate(visited, lid, upt_tau, info, rows, Q);
			for (int i = 1; i < rows; i++)
			{
				for (int j = 1; j<rows; j++)
				{
					tau[i][j] = (1 - Rho)*tau[i][j] + delta_tau[i][j];
				}
			}
			nIteration = nIteration + 1;

			delete[]llrs;
			for (int i = 0; i < ants; i++)
			{
				delete[]visited[i];
			}
			delete[]visited;
			delete[]tau1D;
			delete[]p;
			delete[]alength;
			
		}  //结束一次迭代

		ex = -1;
		for (int i = 1; i < size; i++)
		{
			MPI::COMM_WORLD.Send(&ex, 1, MPI::INT, i, tag6);
		}

		int* cid = new int[rows];
		int pos = maxIndex(lBest, iterations);
		for (int i = 1; i < rows; i++)
		{
			cid[i] = rBest[pos][i];
		}
		cid[0] = sumOfArr(cid, rows);

		cout << "view of result:" << endl;
		for (int i = 1; i < rows; i++)
		{
			cout << cid[i] << "  ";
			if (i % 20 == 0)
			{
				cout << endl;
			}
		}

		delete[]lBest;
		for (int i = 0; i < iterations; i++)
		{
			delete[]rBest[i];
		}
		delete[]rBest;
		delete[]kk;
		delete[]randPos;
		delete[]cid;

		//cout << "Total Time(MPI) = " << totalTime << endl;
		cout << "Total Time(Clock) = " << totalClock_t << endl;
	}
	else
	{
		int nIteration = 1;
		clock_t t1, t2;
		clock_t searchTime = 0;
		while (nIteration < iterations)
		{
			int ex;
			MPI::COMM_WORLD.Recv(&ex, 1, MPI::INT, 0, tag6);
			if (ex == -1)
			{
				t2 = clock();
				searchTime += (t2 - t1);
				break;
			}

			int* p = new int[ants + 2];
			double* aR = new double[ants];
			double* tau1D = new double[rows*rows];
			int* length = new int[ants];
			MPI::COMM_WORLD.Recv(length, ants, MPI::INT, 0, tag);
			MPI::COMM_WORLD.Recv(p, ants + 2, MPI::INT, 0, tag1);
			MPI::COMM_WORLD.Recv(aR, ants, MPI::DOUBLE, 0, tag2);
			MPI::COMM_WORLD.Recv(tau1D, rows*rows, MPI::DOUBLE, 0, tag3);

			double** tau = new double*[rows];
			for (int i = 0; i < rows; i++)
			{
				tau[i] = new double[rows];
				for (int j = 0; j < rows; j++)
				{
					tau[i][j] = tau1D[i*rows + j];
				}
			}

			double* llr = new double[ants];
			for (int i = 0; i < ants; i++)
			{
				llr[i] = -1;
			}
			int* visited = new int[ants*rows];
			for (int i = 0; i < ants*rows; i++)
			{
				visited[i] = -1;
			}
			t1 = clock();
			for (int i = rank; i < ants; i += (size - 1))
			{
				CAnt a(Alpha, Beta, aSet);
				a.update(p[i], tau, length[i], aR[i]);
				llr[i] = a.llr;
				for (int j = 0; j < rows; j++)
				{
					visited[i*rows + j] = a.visited[j];
				}
			}
			t2 = clock();
			searchTime += (t2 - t1);
			MPI::COMM_WORLD.Send(llr, ants, MPI::DOUBLE, 0, tag4);
			MPI::COMM_WORLD.Send(visited, rows*ants, MPI::INT, 0, tag5);

			nIteration++;

			delete[]p;
			delete[]aR;
			delete[]tau1D;
			for (int i = 0; i < rows; i++)
			{
				delete[]tau[i];
			}
			delete[]tau;
			delete[]llr;
			delete[]visited;
			delete[]length;

		}
		cout << "proc." << rank << " : searchTime: " << searchTime << endl;
	}

	MPI_Finalize();
	return 0;
}
