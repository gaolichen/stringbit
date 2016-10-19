#pragma warning(disable:4018)
#include<iostream>
#include<complex>
#include <cmath>
#include <sstream>
#include <cassert>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <map>
#include "BitUtility.h"
#include "StringBitMatrices.h"
#include "EnergyCalculator.h"
#include "Partition.h"
#include <Eigen/Dense>

using namespace std;

void ChopInplace(MatrixSB &m)
{
	for (int i = 0; i < m.innerSize(); i++)
	{
		for (int j = 0; j < m.innerSize(); j++)
		{
			m(i, j) = Chop(m(i, j));
		}
	}
}

MatrixSB Chop(MatrixSB m)
{
	MatrixSB ret = m;
	for (int i = 0; i < ret.innerSize(); i++)
	{
		for (int j = 0; j < ret.innerSize(); j++)
		{
			ret(i, j) = Chop(ret(i, j));
		}
	}

	return ret;
}

void TestMatrixCS()
{
	StringBitMatrices sbm;
	for (int M = 3; M <= 10; M++)
	{
		MatrixSB realZero = MatrixSB::Zero(M, M);
		MatrixSB id = MatrixSB::Identity(M, M);
		for (int L = 1; L < M - 1; L++)
		{
			MatrixSB c = sbm.MatrixC(M, L);
			MatrixSB s = sbm.MatrixS(M, L);
			MatrixSB zero = c * s.transpose() + s * c.transpose();
			MatrixSB one = c * c.adjoint() + s * s.adjoint();
			ChopInplace(zero);

			if (zero != realZero)
			{
				cout << "error! M=" << M << ", L=" << L << ". Expect zero but returned " << endl << zero << endl;
				return;
			}
			if (Chop(one - id) != realZero)
			{
				cout << "error! M=" << M << ", L=" << L << ". Expect identity but returned " << endl << one << endl;
				return;
			}
		}
	}

	cout << "TestMatrixCS passed!" << endl;
}

void TestMatrixA()
{
	StringBitMatrices sbm;
	for (int M = 3; M < 10; M++)
	{
		MatrixSB zero = MatrixSB::Zero(M, M);
		for (int L = 1; L < M; L++)
		{
			MatrixSB z = sbm.MatrixA(M, L + 1, M) - sbm.MatrixAV(M, L);
			ChopInplace(z);

			if (z != zero)
			{
				cout << "Test MatrixAV error! M=" << M << ", L=" << L << endl;
				return;
			}

			z = sbm.MatrixA(L, L + 1, M) + sbm.MatrixA(M, 1, M) - sbm.MatrixAW(M, L);
			ChopInplace(z);

			if (z != zero)
			{
				cout << "Test MatrixAW error! M=" << M << ", L=" << L << endl;
				return;
			}
		}
	}

	cout << "TestMatrixA passed!" << endl;
}

void TestAllStates()
{
	vector<StateInfo> states = EnergyCalculator::AllStates(5);
	for (int i = 0; i < states.size(); i++)
	{
		cout << states[i].Ops << " " << states[i].Energy << endl;
	}
}

void TestVevCalculator()
{
	StringBitMatrices sbm;
	MatrixSB matM = sbm.MatrixM(3, 1);
	VevCalculator calc(matM);
	CDT res = calc.CalculateVev(3);
	cout << res << endl;
}

void TestEnergyCorrection(int maxM = 17)
{
	EnergyCalculator calc;
	Stopwatch watch;
	DT expected[8] = {-0.2942, 1.5161, 8.8474, 24.8336, 52.487, 94.7313, 154.4182, 234.55};
	for (int i = 3; i <= maxM; i += 2)
	{
		watch.Start();
		DT res = calc.EnergyCorrection(i);
		cout << "M= " << i << ", totalStates = " << calc.TotalStates() << ", calculateTime=" << calc.CalculateTime();
		cout << ", totalTime=" << watch.Stop() << ", E=" << res << endl;

		// TODO: the error looks a little too big. May require investigate later.
		if (i <= 17 && abs(res - expected[(i-3)/2])/abs(res) > 1e-2)
		{
			cout << "But expected energy correction is " << expected[(i-3)/2] << endl;
			cout << "!!!!!!!!!!!TestEnergyCorrection failed!!!!!!!!" << endl;
			return;
		}
	}

	cout << "TestEnergyCorrection passed!" << endl;
}

void TestEnergyCorrection2(int s, int maxM)
{
	EnergyCalculator calc(s);
	for (int M = 2; M <= maxM; M++)
	{
		if ((s * (M - 1)) % 2 == 1) continue;
		Stopwatch watch;
		watch.Start();
		cout << "M=" << M << " s=" << s << " E=" << calc.EnergyCorrection(M);
		cout << " time=" << watch.Stop() << " seconds" << endl;
	}
}

void TestMatrices(int M, int L)
{
	StringBitMatrices sbm;
	MatrixSB omegaV = sbm.OmegaV(M, L);
	MatrixSB omegaW = sbm.OmegaW(M, L);
	MatrixSB matM = sbm.MatrixM(M, L);
	ChopInplace(omegaV);
	ChopInplace(omegaW);
	ChopInplace(matM);
	CDT gammaV = sbm.GammaPV(M, L);
	CDT gammaW = sbm.GammaPW(M, L);
	gammaV = Chop(gammaV);
	gammaW = Chop(gammaW);

	cout << "M=" << M << ", L=" << L << endl;
	cout << "Matrix M:" << endl << matM << endl;
	cout << "OmegaV:" << endl << omegaV << endl;
	cout << "OmegaW: " << endl << omegaW << endl;
	cout << "gammaV=" << gammaV << endl;
	cout << "gammaW=" << gammaW << endl;
}

void TestOperatorVev()
{
	int M = 5;
	int L = 2;

	StringBitMatrices sbm;

	MatrixSB omegaV = sbm.OmegaV(M, L);
	MatrixSB omegaW = sbm.OmegaW(M, L);
	MatrixSB matM = sbm.MatrixM(M, L);
	ChopInplace(omegaV);
	ChopInplace(omegaW);
	ChopInplace(matM);
	CDT gammaV = sbm.GammaPV(M, L);
	CDT gammaW = sbm.GammaPW(M, L);
	gammaV = Chop(gammaV);
	gammaW = Chop(gammaW);

	VevCalculator calc(matM);

	for (int i = 0; i < (1<<M); i++)
	{
		CDT res = EnergyCalculator::OperatorVev(i, M, calc, gammaW, omegaW);
		res = Chop(res);
		if (res != .0)
			cout << i << ": " << res << endl;
	}
}

void TestBitManager()
{
	BitManager bm;
	bm.Init(10);
	vector<int>& v = bm.GetNumbers(5,2);
	for (int i = 0; i < v.size(); i++)
	{
		cout << v[i] << ' ';
	}
	cout << endl;
}

bool TestPartition(int M, int s, bool output = false)
{
	ostream* os = &cout;
	ostringstream oss;
	if (!output)
	{
		os = &oss;
	}

	*os << "M=" << M << " s=" << s << endl;
	Partitioner part(M, s);
	Stopwatch watch;
	watch.Start();
	vector<i64>& res1 = part.AllPartitions();
	*os << "returned " << res1.size() << " partitions in " <<  watch.Stop() << "s." << endl;
	watch.Start();
	vector<vector<int> > & res2 = part.AllPartitionsBruteForce();
	*os << "expected " <<  res2.size() << " partitions in " << watch.Stop() << "s." << endl;
	if (res1.size() != res2.size())
	{
		cout << "Numbers of partitions are different. Test failed!!!" << endl;
		return false;
	}

	// cout << "Partitios from AllPartitions(): " << endl;

	for (int i = 0; i < res1.size(); i++)
	{
		vector<int> partition = part.ToArray(res1[i]);
		if (partition != res2[i])
		{
			cout << "partition not match at position " << i << endl;
			cout << "returned: " << partition << endl;
			cout << "expected: " << res2[i] << endl;
			return false;
		}
	}

	/*cout << "partitions from AllPartitionsBruteForce: " << endl;
	for (int i = 0; i < res2.size(); i++)
	{
	cout << "partition " << i << ": " << res2[i] << endl;
	}*/

	return true;
	//cout << "TestPartition passed!" << endl;
}

void TestPartition()
{
	for (int M = 3; M <= 7; M++)
	{
		for (int s = 1; s <= 5; s++)
		{
			if (TestPartition(M, s) == false)
			{
				cout << "TestPartition Failed!! " << endl;
				return;
			}
		}
	}

	cout << "TestPartition passed!" << endl;
}

vector<int> ToVector(int* array, int size)
{
	return vector<int>(array, array + size);
}

void OutputReturnAndExpectedDivision(vector<vector<i64> > &actual, vector<vector<i64> > &expect)
{
	cout << "returned: ";
	for (int i = 0; i < actual.size(); i++) cout << actual[i] << ' ';
	cout << endl;

	cout << "expected: ";
	for (int i = 0; i < expect.size(); i++) cout << expect[i] << ' ';
	cout << endl;
}

bool TestDividePartition(int s, int offset, vector<int>& partition, vector<vector<i64> >& expect)
{
	DividePartition dp(s);
	vector<vector<i64> > res = dp.Divide(partition, offset);
	if (res != expect)
	{
		cout << "TestDividePartition failed! s = " << s << endl;
		OutputReturnAndExpectedDivision(res, expect);
		return false;
	}

	return true;
}

#define LENGTH(a) sizeof(a)/sizeof((a)[0])
bool TestDividePartition2(int s, int offset, vector<int>& partition)
{
	DividePartitionForTest dpt(s);
	vector<vector<i64> > &expect = dpt.Divide(partition, offset);
	sort(expect.begin(), expect.end());

	DividePartition dp(s);
	vector<vector<i64> > &actual = dp.Divide(partition, offset);
	sort(actual.begin(), actual.end());

	if (expect != actual)
	{
		cout << "TestDividePartition2 failed!" << endl;
		OutputReturnAndExpectedDivision(actual, expect);
		return false;
	}
	else
	{
		return true;
	}
}

void TestDividePartition2()
{
	int total = 10;
	srand(time(NULL));

	for (int i = 0; i < total; i++)
	{
		int s = rand() % 7 + 2;
		int offset = rand() % 6;
		int M = rand() % 6 + 2;
		if (s * M > 32) M = 32 / s;

		vector<int> partition(M);
		for (int j = 0; j < M; j++)
		{
			partition[j] = rand() % (s + 1);
		}

		if (!TestDividePartition2(s, offset, partition))
		{
			cout << "s=" << s << ", M=" << M << ", offset=" << offset << endl;
			cout << "partition = " << partition << endl;
			return;
		}
	}
	
	cout << "TestDividePartition2 passed" << endl;
}

bool TestDividePartition3(int s, vector<int> &partition, vector<i64> &division)
{
	DividePartitionForTest dpt(s);
	vector<vector<i64> > &expect = dpt.Divide2(partition, division);
	sort(expect.begin(), expect.end());

	DividePartition dp(s);
	vector<vector<i64> > &actual = dp.Divide2(partition, division);
	sort(actual.begin(), actual.end());
	
	if (actual != expect)
	{
		cout << "TestDividePartition3 failed!" << endl;
		OutputReturnAndExpectedDivision(actual, expect);
		return false;
	}
	
	return true;
}

void TestDividePartition3()
{
	int total = 10;
	srand(time(NULL));

	for (int i = 0; i < total; i++)
	{
		int s = rand() % 7 + 2;
		int M = rand() % 6 + 2;
		if (s * M > 32) M = 32 / s;
		vector<int> partition(M);
		for (int j = 0; j < M;  j++)
		{
			partition[j] = rand() % (s + 1);
		}

		vector<i64> division(s);
		for (int j = 0; j < s; j++)
		{
			division[j] = (rand() % M) << (M + 1);
		}

		sort(division.begin(), division.end());
		reverse(division.begin(), division.end());

		if (!TestDividePartition3(s, partition, division))
		{
			cout << "s=" << s << ", M=" << M << endl;
			cout << "partition: " << partition << endl;
			cout << "division: " << division << endl;
			return;
		}
		
	}

	cout << "TestDividePartition3 passed!" << endl;
}

void TestDividePartition()
{
	int s = 2;
	int offset = 0;
	int arr[] = {2, 1, 2};
	vector<int> partition(arr, arr + 3);
	vector<vector<i64> > expect(1, vector<i64>(s));
	int expect00[] = {1, 2, 3};
	int expect01[] = {1, 3};

	expect[0][0] = Digit2Num(expect00, 3);
	expect[0][1] = Digit2Num(expect01, 2);

	if (!TestDividePartition(s, offset, partition, expect))
	{
		return;
	}

	s = 3;
	expect.clear();
	expect.resize(5, vector<i64>(s));
	int expect40[] = {1, 3};
	int expect41[] = {3};
	int expect42[] = {1, 2};

	int expect30[] = {1, 3};
	int expect31[] = {3, 1};
	int expect32[] = {2};

	int expect20[] = {2, 3};
	int expect21[] = {1,3};
	int expect22[] = {1};

	int expect00b[] = {1, 2, 3};
	int expect01b[] = {1, 3};

	int expect10[] = {1, 2, 3};
	int expect11[] = {3};
	int expect12[] = {1};

	expect[0][0] = Digit2Num(expect00b, LENGTH(expect00));
	expect[0][1] = Digit2Num(expect01b, LENGTH(expect01b));
	expect[0][2] = 0;

	expect[1][0] = Digit2Num(expect10, LENGTH(expect10));
	expect[1][1] = Digit2Num(expect11, LENGTH(expect11));
	expect[1][2] = Digit2Num(expect12, LENGTH(expect12));

	expect[2][0] = Digit2Num(expect20, LENGTH(expect20));
	expect[2][1] = Digit2Num(expect21, LENGTH(expect21));
	expect[2][2] = Digit2Num(expect22, LENGTH(expect22));

	expect[3][0] = Digit2Num(expect30, LENGTH(expect30));
	expect[3][1] = Digit2Num(expect31, LENGTH(expect31));
	expect[3][2] = Digit2Num(expect32, LENGTH(expect32));

	expect[4][0] = Digit2Num(expect40, LENGTH(expect40));
	expect[4][1] = Digit2Num(expect41, LENGTH(expect41));
	expect[4][2] = Digit2Num(expect42, LENGTH(expect42));

	if (!TestDividePartition(s, offset, partition, expect))
	{
		return;
	}

	cout << "TestDividePartition passed!" << endl;
}

void TestModesGenerator()
{
	int M = 5;
	int L = 2;
	int s = 2;
	ModesGenerator gen(M, L, s);
	vector<vector<i64> > res = gen.Generate();
	cout << "Total modes: " << res.size() << endl;
	for (int i = 0; i < res.size(); i++)
	{
		cout << "modes " << i << ":\t";
		for (int j = 0; j < res[i].size(); j++)
		{
			cout << Num2Digit(res[i][j], M) << " ";
		}
		cout << endl;
	}
}

int main()
{
	TestMatrixCS();
	TestMatrixA();
	//TestAllStates();
	//TestVevCalculator();
	TestEnergyCorrection();
	TestEnergyCorrection2(2, 12);
	TestEnergyCorrection2(3, 9);
	//TestMatrices(4,1);
	//TestOperatorVev();
	TestPartition();
	//TestBitManager();
	//TestModesGenerator();
	TestDividePartition();
	TestDividePartition2();
	TestDividePartition3();
	return 0;
}
