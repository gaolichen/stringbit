#pragma once
#include <string>
#include <vector>
#include <fstream>
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
	Hamiltonian h0;
	Hamiltonian deltaH;
	Hamiltonian hPrime;
	void InitHamiltonians();
	void HamToMatlab(int bits, Hamiltonian& ham, ofstream& os);
public:
	ScriptGenerator(string rootFolder, StateType stateType);

	void OutputHamToMatlab(int bits, Hamiltonian& ham);
	void OutputHamToMatlab(int bits, bool inverted);
	void OutputHamToLaTeX(vector<int> bits, vector<double> scales, string filename);
	void OutputNormToLaTeX(int minBits, int maxBits, string filename);
	void OutputNormToMatlab(int bits);
	void OutputStateToLaTeX(int minBit, int maxBit, string filename);
};