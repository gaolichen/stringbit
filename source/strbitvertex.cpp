#pragma warning(disable:4018)
#include<iostream>
#include<complex>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <vector>
#include <map>
#include "BitUtility.h"
#include "StringBitMatrices.h"
#include "EnergyCalculator.h"
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

void TestEnergyCorrection()
{
	EnergyCalculator calc;
	Stopwatch watch;
	for (int i = 3; i <= 25; i += 2)
	{
		watch.Start();
		DT res = calc.EnergyCorrection(i);
		cout << "M= " << i << ", totalStates = " << calc.TotalStates() << ", calculateTime=" << calc.CalculateTime();
		cout << ", totalTime=" << watch.Stop() << ", E=" << res << endl;
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

	// test VevCalculator 
	/*for (int i = 0; i < (1<<M); i++)
	{
		CDT res = calc.CalculateVev(i, omegaV);
		res = Chop(res);
		if (res != .0)
			cout << i << ": " << res << endl;
	}*/

	for (int i = 0; i < (1<<M); i++)
	{
		CDT res = EnergyCalculator::OperatorVev(i, calc, gammaW, omegaW);
		res = Chop(res);
		if (res != .0)
			cout << i << ": " << res << endl;
	}
}

class BitManager
{
private:
	vector<vector<vector<int> > > bitCollection;
	vector<vector<vector<int> > > bitchecks;
public:
	void Init(int maxBit)
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

	vector<int>& GetNumbers(int maxBit, int bits)
	{
		return bitCollection[maxBit][bits];
	}

	vector<int>& GetChecks(int maxBit, int bits)
	{
		return bitchecks[maxBit][bits];
	}
};

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

class ModeDistributor
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
	ModeDistributor(int M_, int s_): M(M_), s(s_)
	{
		bitUnit = 1;
		while ((1 << bitUnit) <= s) bitUnit++;
		allOnes = (1 << bitUnit) - 1;
		bm.Init(s + 1);
	}
	
	void Doit()
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

	void distribute(i64 mode, int bit, vector<int>& curr)
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
};

void TestModeDistributor()
{
	cout << "input M and L: ";
	int M, L;
	cin >> M >> L;
	ModeDistributor dis(M, L);
	dis.Doit();
}

int main()
{
	TestMatrixCS();
	TestMatrixA();
	TestAllStates();
	TestVevCalculator();
	//TestEnergyCorrection();
	TestMatrices(4,1);
	TestOperatorVev();
	//TestModeDistributor();
	TestBitManager();
	return 0;
}
