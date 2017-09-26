#include "AntScan.h"


CAntScan::CAntScan(CAntSets aSet)
{
	e.seed(time(0));
	ureal = uniform_real_distribution<double>(0, 1);
	unormal = normal_distribution<>(0, 1);
	Alpha = 1.0;
	Beta = 1.0;
	Rho = 0.6;
	Q = 1.0;

	rows = aSet.getPSize() + 1;
	ants = rows / 2;
	iterations = 12;
	initTime = 5;
	eta = initEta(aSet);

	rBest = new int*[iterations];
	rPath = new int*[iterations];
	lBest = new double[iterations];
	for (int i = 0; i < iterations; i++)
	{
		rBest[i] = new int[rows];
		rPath[i] = new int[rows];
		lBest[i] = 0;
		for (int j = 0; j < rows; j++)
		{
			rBest[i][j] = 0;
			rPath[i][j] = 0;
		}
	}

	clusterId = antCluster(aSet);
}


CAntScan::~CAntScan()
{
	for (int i = 0; i < iterations; i++)
	{
		delete[]rBest[i];
		delete[]rPath[i];
	}
	delete[]rBest;
	delete[]lBest;
	delete[]rPath;
}

int* CAntScan::antCluster(CAntSets aSet)
{
	clock_t totalClock_t = 0;
	clock_t begin, end;
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

	CAnt* ant = new CAnt[ants];
	for (int i = 1; i < ants; i++)
	{
		CAnt a(Alpha, Beta, aSet, eta);
		ant[i] = a;
	}

	while (nIteration < iterations)
	{
		if (nIteration <= initTime)
		{
			randPos = randPerm(rows, ants);
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
		//更新蚂蚁,并行部分
//#pragma omp parallel for num_threads(2) firstprivate(varKK, kAve, nIteration)
		begin = clock();
		for (int i = 1; i < ants; i++)
		{
			int length = normLen(varKK, kAve);
			ant[i].update(randPos[i], tau, length, ureal(e));
		}
		end = clock();
		totalClock_t += (end - begin);

		if (nIteration > 1)
		{
			if (ant[1].llr < lBest[nIteration - 1])
			{
				for (int j = 1; j < rows; j++)
				{
					{
						ant[1].visited[j] = rBest[nIteration - 1][j];
						//ant[1].antPath[j] = rPath[nIteration - 1][j];
					}
				}
				ant[1].llr = lBest[nIteration - 1];
			}
		}

		double* llr = new double[ants];
		llr[0] = 0;
		int ind = 1;
		for (int i = 1; i < ants; i++)
		{
			llr[i] = ant[i].llr;
			if (ant[i].llr > lBest[nIteration])
			{
				lBest[nIteration] = ant[i].llr;
				ind = i;
			}
		}
		for (int i = 1; i < rows; i++)
		{
			rBest[nIteration][i] = ant[ind].visited[i];
			//rPath[nIteration][i] = ant[ind].antPath[i];
		}
		rBest[nIteration][0] = sumOfArr(rBest[nIteration], rows);
		kk[nIteration] = kAve;

		//cout << "iteration[" << nIteration << "]: ant[" << ind << "]'s LLR = " << lBest[nIteration];
		//cout << ", Length of Path = " << rBest[nIteration][0] << endl;

		bool matchs = accuCalculate(rBest[nIteration], aSet);
		if (matchs)
		{
			nIteration = iterations - 1;
		}


		/**** Fourth step: Check the precocious of antScan and Update the pheromones ****/
		double info = 0.0;
		if (nIteration > 1)
		{
			if (lBest[nIteration] == lBest[nIteration - 1]){
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
		int* lid = lidSort(llr);
		double** delta_tau = tauUpdate(ant, lid, upt_tau, info);
		for (int i = 1; i < rows; i++)
		{
			for (int j = 1; j < rows; j++)
			{
				tau[i][j] = (1 - Rho)*tau[i][j] + delta_tau[i][j];
			}
		}

		nIteration = nIteration + 1;

		delete[]llr;
	}
	cout << "runtime: " << totalClock_t << "ms" << endl;

	//get the final cluster ids and its llr
	int* cid = new int[rows];
	int pos = maxIndex(lBest, iterations);
	for (int i = 1; i < rows; i++)
	{
		cid[i] = rBest[pos][i];
	}
	cid[0] = sumOfArr(cid, rows);

	delete[]ant;
	return cid;
}

double** CAntScan::initEta(CAntSets aSet)
{

	int i, j;
	double** eta = new double*[rows];
	for (int w = 0; w < rows; w++)
	{
		eta[w] = new double[rows];
		for (int x = 0; x < rows; x++)
		{
			eta[w][x] = 0;
		}
	}

	for (i = 1; i < rows; i++)
	{
		for (j = i; j < rows; j++)
		{
			if (i == j)
			{
				eta[i][j] = 1 / (10e-10);
			}
			else
			{
				eta[i][j] = aSet.ratio[j] / aSet.ratio[i];
				eta[j][i] = aSet.ratio[i] / aSet.ratio[j];
			}
		}
	}
	return eta;
}

double** CAntScan::init2DArrFor1(int rows, int cols)
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

int CAntScan::binSearch(int*Array, int start, int send, int key)
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

int* CAntScan::randPerm(int sizeRange, int sizeArray)
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

double* CAntScan::arrayCumSum(double** arr, int arrLen)
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

double CAntScan::varianceOfArr(int* arr, int arrLen)
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

double CAntScan::aveValue(int* arr, int arrLen)
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

int CAntScan::maxIndex(double* arr, int arrLen)
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

int CAntScan::sumOfArr(int* arr, int arrLen)
{
	int i, j = 0;
	for (i = 1; i < arrLen; i++)
	{
		j = j + arr[i];
	}
	return j;
}

bool CAntScan::accuCalculate(int* clusterI, CAntSets aSet)
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

int* CAntScan::lidSort(double* llr)
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

double** CAntScan::tauUpdate(CAnt* ant, int* lid, int upt_tau, double info)
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
		startId = ant[lid[i]].visited[0];
		for (j = 1; j < l; j++)
		{
			if (ant[lid[i]].visited[j] == 1 && lid[i] != j)
			{
				delta_tau[startId][j] = delta_tau[startId][j] + Q*pow(info, upt_tau - i);
				delta_tau[j][startId] = delta_tau[startId][j];
			}
		}
	}
	return delta_tau;
}

int* CAntScan::getClusterId()
{
	return this->clusterId;
}

int CAntScan::normLen(double varKK, double kAve)
{
	int k = -1;
	int j = 1;
	int length = 0;
	while (k < 0)
	{
		j = (int)(ceil(sqrt(varKK)*unormal(e)) + kAve);
		if (j > 0 && j < (rows - 1) / 2)
		{
			k = 1;
			length = j;
		}
	}
	return length;
}