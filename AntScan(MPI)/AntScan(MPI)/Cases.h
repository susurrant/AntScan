#pragma once

#include <vector>
using namespace std;
class CCases
{
public:
	CCases();
	CCases(int x, int y, int cas, int pop);
	~CCases();
public:
	int x;
	int y;
	int cas;
	int pop;
	vector<int> nb;
};



