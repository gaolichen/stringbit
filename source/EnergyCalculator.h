#pragma once
#include<iostream>
#include<complex>
#include <cmath>
#include <cassert>
#include <vector>
#include <map>
#include "BitUtility.h"
#include "StringBitMatrices.h"
#include <Eigen/Dense>

using namespace std;

class VevCalculator
{
private:
	MatrixSB matM;
	int M;
	map<int, CDT> cache;
	CDT DoCalculate(int ops);

public:
	VevCalculator(MatrixSB &mat) : matM(mat)
	{
		M = matM.innerSize();
	}

	CDT CalculateVev(int ops)
	{
		if (BitCount(ops) % 2 == 1) return 0.0;
		return DoCalculate(ops);
	}

	CDT CalculateVev(int ops, MatrixSB &omega);
};

class StateInfo
{
public:
	int Ops;
	DT Energy;
	StateInfo(int ops, DT energy):Ops(ops), Energy(energy)
	{
	}

};

class EnergyCalculator
{
private:
	int totalStates;
	double calculateTime;
	Stopwatch watch;
	DT EnergyCorrection(int M, int L);

public:

	static CDT OperatorVev(int ops, VevCalculator &calc, CDT gamma, MatrixSB &omega)
	{
		return 2.0 * (calc.CalculateVev(ops) * gamma + calc.CalculateVev(ops, omega));
	}

	static vector<StateInfo> AllStates(int M)
	{
		vector<StateInfo> ret;
		DT e0 = -4 / tan(PI/(2 * M));
	
		for (int i = 0; i < (1<<M); i+=2)
		{
			DT deltE = 0.0;
			int modes = 0;
			for (int j = 1; (1<<j) <= i; j++)
			{
				if ((i & (1<<j)) == 0) continue;
				deltE += 8 * sin(j * PI/M);
				modes += j;	
			}
		
			if ((M % 2 == 1 && modes % M == 0) || (M % 2 ==0 && modes % M == M / 2))
			{
				ret.push_back(StateInfo(i, e0 + deltE));
			}
		}

		return ret;
	}

	DT EnergyCorrection(int M);

	int TotalStates()
	{
		return totalStates;
	}
	
	double CalculateTime()
	{
		return calculateTime;
	}
};