#pragma once
#include <string>
#include <vector>
using namespace std;

class ScriptGenerator
{
private:
	void GroupStates(int bits, vector<vector<int> >& states);
	double NormMatrixScale(int bits, int index);
public:
	ScriptGenerator();

	void OutputHamToMatlab(int bits, bool invert, string folder);
	void OutputHamToLaTeX(vector<int> bits, vector<double> scales, string filename);
	void OutputNormToLaTeX(int minBits, int maxBits, string filename);
	void OutputNormToMatlab(int bits, string folder, bool output);
	void OutputStateToLaTeX(int minBit, int maxBit, int type, string filename);
};