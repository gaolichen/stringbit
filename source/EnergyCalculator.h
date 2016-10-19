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
	int s;
	int totalStates;
	double calculateTime;
	Stopwatch watch;
	DT EnergyCorrection(int M, int L);
public:
	EnergyCalculator() : s(1) {}

	EnergyCalculator(int s_) : s(s_) {}

	static CDT OperatorVev(int ops, int M, VevCalculator &calc, CDT &gamma, MatrixSB &omega);

	static CDT OperatorVevAllZeros(int ops, int M, VevCalculator &calc, CDT &gammaW, MatrixSB &omegaW, CDT &gammaV, MatrixSB &omegaV);

	static vector<StateInfo> AllStates(int M);

	static DT Operators2Energy(vector<int>& ops, int M, int s);

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
