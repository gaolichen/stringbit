#pragma once
#include <cassert>
#include <vector>
#include <map>
#include "BitUtility.h"

class BitManager
{
private:
	vector<vector<vector<int> > > bitCollection;
	vector<vector<vector<int> > > bitchecks;
public:
	void Init(int maxBit);

	vector<int>& GetNumbers(int maxBit, int bits)
	{
		return bitCollection[maxBit][bits];
	}

	vector<int>& GetChecks(int maxBit, int bits)
	{
		return bitchecks[maxBit][bits];
	}
};

class Partitioner
{
private:
	int M;
	int s;
	int bitUnit;
	int allOnes;
	vector<i64> modes;
	vector<vector<int> > res;
	BitManager bm; 
public:
	Partitioner(int M_, int s_);
	void Doit();
	void distribute(i64 mode, int bit, vector<int>& curr);
};
