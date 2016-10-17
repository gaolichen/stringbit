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

vector<i64>& Partitioner::AllPartitions()
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

	return modes;	
}

vector<int> Partitioner::ToArray(i64 mode)
{
	vector<int> ret(M);
	for (int i = 0; i < M - 1; i++)
	{
		ret[i] = (mode & allOnes);
		mode >>= bitUnit;
	}

	return ret;
}

vector<vector<int> >& Partitioner::AllPartitionsBruteForce()
{
	vector<int> partition(M);
	AllPartitionsRecursive(partition, 0, M - 1);
	return allPartitionsTest;
}

void Partitioner::AllPartitionsRecursive(vector<int>& partition, int total, int bit)
{
	if (bit == 0)
	{
		if (s * (M - 1) % 2 == 0 && total % M == 0)
		{
			allPartitionsTest.push_back(partition);
		}
		else if (s * (M - 1) % 2 == 1 && total % M == M / 2)
		{
			allPartitionsTest.push_back(partition);
		}

		return;
	}
	
	for (int i = 0; i <= s; i++)
	{
		partition[bit - 1] = i;
		AllPartitionsRecursive(partition, total + i * bit, bit - 1);
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
