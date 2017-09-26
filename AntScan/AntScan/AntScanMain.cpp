#include <iostream>
#include "AntSets.h"
#include "AntScan.h"
using namespace std;

int main()
{
	//clock_t begin, end;
	//begin = clock();
	CAntSets aSet("pSimulate.txt", "nbSimulate.txt");
	CAntScan aScan(aSet);
	int* cluster = aScan.getClusterId();
	cout << "view of result:" << endl;
	for (int i = 1; i < aScan.rows; i++)
	{
		cout << cluster[i] << "   ";
		if (i % 20 == 0)
		{
			cout << endl;
		}
	}
	//end = clock();
	//cout << "runtime: " << (end - begin) << "ms" << endl;
	return 0;
}
