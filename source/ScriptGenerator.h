#pragma once
#include <string>
#include <vector>
#include "StateType.h"
#include "Hamiltonian.h"
using namespace std;

class ScriptGenerator
{
private:
	string rootFolder;
	StateType type;
	void GroupStates(int bits, vector<vector<int> >& states);
	double NormMatrixScale(int bits, int index);
public:
	ScriptGenerator(string rootFolder, StateType stateType);

	void OutputHamToMatlab(int bits, Hamiltonian& ham);
	void OutputHamToLaTeX(vector<int> bits, vector<double> scales, string filename);
	void OutputNormToLaTeX(int minBits, int maxBits, string filename);
	void OutputNormToMatlab(int bits);
	void OutputStateToLaTeX(int minBit, int maxBit, string filename);
};