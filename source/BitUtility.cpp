#pragma warning(disable:4018)
#include "BitUtility.h"
#include <time.h>
using namespace std;

int BitCount(int i)
{
	 i = i - ((i >> 1) & 0x55555555);
     i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
     return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

int CyclicRotation(int n, int bitNumber)
{
	return ((n & 1) << (bitNumber -1)) | (n >> 1);
}

i64 BinomialCoefficient(i64 a,i64 b)
{
	if (b > a) return 0;
	if (b + b > a) b = a - b;
	i64 ret = 1;
	for (int i = 0; i < b; i++)
	{
		ret *= (a - i);
		ret /= (i+1);
	}

	return ret;
}

int PickBits(int n, int bitNumber)
{
	return n & ((1 << bitNumber) - 1);
}

int InverseNumber(const vector<int>& v)
{
	int ret = 0;
	for (int i = 1; i < v.size(); i++)
	{
		for (int j = 0; j < i; j++)
			if (v[i] < v[j]) ret++;
	}

	return ret;
}

int BuildMask(int trace, int bitNumber)
{
	return (bitNumber << MAX_TRACE_BITS) | trace;
}

int Gcd(int a, int b)
{
	if (b == 0) return a;
	return Gcd(b, a % b);
}

Stopwatch::Stopwatch()
{
	start = 0;
}

void Stopwatch::Start()
{
	start = clock();
}

/// return how much time in second since last Start() function call
double Stopwatch::Stop()
{
	double ret = (clock() - start) / (double)CLOCKS_PER_SEC;
	start = 0;

	return ret;
}