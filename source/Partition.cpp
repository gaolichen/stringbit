#pragma warning(disable:4018)
#include <iostream>
#include <vector>
#include "BitUtility.h"
#include "Partition.h"
using namespace std;

void BitManager::Init(int maxBit)
{
	for (int i = 0; i <= maxBit; i++)
	{
		bitCollection.push_back(vector<vector<int> >(i + 1));
		bitchecks.push_back(vector<vector<int> >(i + 1));
	}

	for (int i = 0; i < (1 << maxBit); i++)
	{
		int check = 0;
		for (int j = 0; j < maxBit; j++)
		{
			if (IsBitSet(i, j) && !IsBitSet(i, j + 1))
			{
				check |= (1 << (j + 1));
			}
		}

		int bits = BitCount(i);
		for (int j = maxBit; (1<<j) > i; j--)
		{
			bitCollection[j][bits].push_back(i);
			bitchecks[j][bits].push_back(PickBits(check, j));
		}
	}
}

Partitioner::Partitioner(int M_, int s_): M(M_), s(s_)
{
	bitUnit = 1;
	while ((1 << bitUnit) <= s) bitUnit++;
	allOnes = (1 << bitUnit) - 1;
	bm.Init(s + 1);
}

void Partitioner::Doit()
{
	i64 n = 1;
	for (int i = 0; i < M - 1; i++) n *= s + 1;
	//cout << "bitUnit=" << bitUnit << endl << " n=" << n << endl;

	for (int i = 0; i < n; i++)
	{
		i64 toadd = 0;
		int ii = i;
		int count = 0;
		int k = 0;
		//cout << "i=" << i << endl;
		while (ii > 0)
		{
			int j = ii % (s + 1);
			ii /= s + 1;
			toadd |= (j<<(k * bitUnit));
			k++;
			count += j * k;
		}

		if (((M & 1) == 1) || ((s & 1) == 0))
		{
			if (count % M == 0) 
			{
				//cout << "a " << count << " " << toadd << endl; 
				modes.push_back(toadd);
			}
		}
		else if (count % M == M / 2)
		{
			//cout << "b " << cout << " " << toadd << endl;
			modes.push_back(toadd);
		}
	}

	for (int i = 0; i < modes.size(); i++)
	{
		i64 op = modes[i];
		cout << "modes " << i << " " << modes[i] << ": ";
		for (int k = 0; k < M - 1; k++)
		{
			cout << (op & allOnes) << " ";
			op >>= bitUnit;
		}

		cout << endl;
	}
		
	for (int i = 0; i < modes.size(); i++)
	{
		//distribute(modes[i], M - 1);
	}
}

void Partitioner::distribute(i64 mode, int bit, vector<int>& curr)
{
	if (bit == 0)
	{
		// found a valid one.
		res.push_back(curr);
		return;
	}
		
	int n = (mode >> (bitUnit * (bit - 1))) & allOnes;
	vector<int> numbers = bm.GetNumbers(s, n);
	vector<int> checks = bm.GetChecks(s, n);
		
	for (int i = 0; i < numbers.size(); i++)
	{
		bool ok = true;
		for (int j = 1; j < s; j++)
		{
			if (!IsBitSet(checks[i], j)) continue;
			if (curr[j] < curr[j-1] + bit) { ok = false; break;}
		}
			
		if (!ok) continue;
		for (int j = 0; j < s; j++)
		{
			if (IsBitSet(i, j))
			{
				curr[j] |= (1 << bit);
			}
		}
			
		distribute(mode, bit - 1, curr);
			
		for (int j = 0; j < s; j++)
		{
			if (IsBitSet(i, j))
			{
				curr[j] ^= (1 << bit);
			}
		}
	}
}