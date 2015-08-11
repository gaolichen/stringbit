#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <time.h>
using namespace std;

//#define HAM_PARAMETER

#if WIN32
typedef __int64 i64;
#else
typedef long long i64;
#endif

//typedef double snum;
typedef i64 snum;

#define MAX_TRACE_BITS 25
#define FULL_TRACE_BITS ((1<<MAX_TRACE_BITS)-1)

int BitCount(int);
int CyclicRotation(int, int);
snum BinomialCoefficient(snum a, i64 b);
int PickBits(int n, int bitNumber);
int InverseNumber(const vector<int>& v);
int BuildMask(int trace, int bitNumber);
int Gcd(int a, int b);
string CombinePath(string path1, string path2);
string Bits2String(int bits, int bitNumber);
string ToUpper(string s);

template<class T> string ToString(T a)
{
	ostringstream oss;
	oss << a;
	return oss.str();
};

class Stopwatch
{
private:
	clock_t start;
public:
	Stopwatch();
	void Start();
	double Stop();
};
