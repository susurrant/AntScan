#include "AntSets.h"


CAntSets::CAntSets()
{
}


CAntSets::~CAntSets()
{
}

CAntSets::CAntSets(string pFileName, string nbFileName)
{
	this->pFileName = pFileName;
	this->nbFileName = nbFileName;
	caseNum = 0;
	popNum = 0;
	readCase();
	readNeighbor();

	int rows = pAntNum + 1;
	pSet = new int*[rows];
	pCoords = new double*[rows];
	nbRelation = new int*[rows];
	ratio = new double[rows];
	for (int i = 0; i < rows; i++)
	{
		pSet[i] = new int[2];
		pCoords[i] = new double[2];
		nbRelation[i] = new int[rows];
		for (int j = 0; j < rows; j++)
		{
			nbRelation[i][j] = 0;
		}
	}

	int tmpNblen = 0;
	for (int i = 0; i < rows - 1; i++)
	{
		pSet[i + 1][0] = ca[i].cas;
		pSet[i + 1][1] = ca[i].pop;
		ratio[i + 1] = (double)pSet[i + 1][0] / pSet[i + 1][1];
		pCoords[i + 1][0] = ca[i].x;
		pCoords[i + 1][1] = ca[i].y;
		tmpNblen = ca[i].nb.size();
		for (int j = 0; j < tmpNblen; j++)
		{
			nbRelation[i + 1][ca[i].nb[j]] = 1;
		}
	}
}


bool CAntSets::readCase()
{
	ifstream  ifs(pFileName, ios::in);
	if (!ifs)
	{
		cout << "Fail to open case file!" << endl;
		return false;
	}

	int x, y, cas, pop;
	while (ifs >> x >> y >> cas >> pop)
	{
		this->ca.push_back(CCases(x, y, cas, pop));
		this->caseNum += cas;
		this->popNum += pop;
	}

	ifs.close();
	return true;
}

bool CAntSets::readNeighbor()
{
	ifstream ifs(nbFileName, ios::in);
	if (!ifs)
	{
		cout << "Fail to open neighbor file!" << endl;
		return false;
	}

	vector<int> neighbor;
	int i, nbor;
	for (i = 0; !ifs.eof(); i++)
	{
		string line;
		getline(ifs, line);
		istringstream ss(line);
		while (!ss.eof())
		{
			ss >> nbor;
			this->ca[i].nb.push_back(nbor);
		}
	}
	this->pAntNum = i;
	ifs.close();
	return true;
}

int CAntSets::getCaseNum()
{
	return caseNum;
}

int CAntSets::getPopNum()
{
	return popNum;
}

int CAntSets::getPSize()
{
	return pAntNum;
}