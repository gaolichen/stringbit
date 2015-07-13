#pragma warning(disable:4018)
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "SingleTrace.h"
#include "BitUtility.h"
#include "StateCollection.h"
#include "TraceState.h"
#include "StateGenerator.h"
#include "Hamiltonian.h"
#include "NormCalculator.h"
#include "ScriptGenerator.h"

using namespace std;

void TestHamiltonianMultiTrace(TraceState& state)
{
	Hamiltonian ham;
	MixState re, im;
	ham.Apply(state, re, im);
	cout << "H" << state << "|0> = ";
	cout << re.ToString();
	cout <<" + i{" << im.ToString() << "}" << endl;
}

void GenerateStates()
{
	StateGenerator generator;
	for (int i = 1; i <= StateGenerator::MAX_BIT_TO_COUNT; i++)
	{
		cout << i << "\t" << generator.SingleStateNumber(i) << "\t" << generator.BosonNumber(i) << endl;
	}

	generator.GenerateAllStates();
	generator.InitStateCollection(StateCollection::Inst());
}

void OutputHamiltonianMatrix(int bits)
{
	Hamiltonian ham;
	vector<vector<Coefficient> > rem, imm;
	ham.Matrix(bits, rem, imm);
	cout << bits << " bits Hamiltonian: " <<endl;
	//for (int i = 0; i < rem.size(); i++)
	//{
	//	for (int j = 0; j < rem[i].size(); j++)
	//	{
	//		if (!imm[i][j].IsZero())
	//		{
	//			cout << imm[i][j] << "i ";
	//		}
	//		else
	//		{
	//			cout << rem[i][j] << " ";
	//		}
	//		
	//	}
	//	cout << endl;
	//}
}

void OutputHamToLatex()
{
	int arr[6] = {2, 3, 4, 5, 6, 7};
	double scale[6] = {1.0, 1.0, 1.0, 1, 0.57, 0.27};
	vector<int> bits(arr, arr + 6);
	ScriptGenerator sc;
	sc.OutputHamToLaTeX(bits, vector<double>(scale, scale + 6), "E:\\Dropbox\\proj\\stringbit\\hamiltonian.tex");
}

void CalculateNorm(int bits, bool output)
{
	BruteForceCalculator calc;
	StateCollection* sc = StateCollection::Inst();
	Stopwatch watch;
	
	watch.Start();
	cout << "Calculating norm matrix for " << bits << " bits states..." << endl;
	for (int i = 0; i < sc->StateNumber(bits); i++)
	{
		for (int j = 0; j < sc->StateNumber(bits); j++)
		{
			Polynomial res = calc.Calculate(sc->GetBosonState(bits, i), sc->GetBosonState(bits, j));
			if (output)
			{
				cout << res << "\t";
			}
		}

		if (output)
		{
			cout << endl;
		}
	}
	cout << "It takes " << watch.Stop() << " seconds." << endl;
}

void TestNormCalculator()
{
	int bits = 5;
	BruteForceCalculator calc;
	StateCollection* sc = StateCollection::Inst();

	const TraceState& a = sc->GetBosonState(bits, 7);
	const TraceState& b = sc->GetBosonState(bits, 8);

	cout << "(" << a << "," << b << ")=";
	cout << calc.Calculate(a, b) << endl;;
}

void GenerateMatlabScript()
{
	//string scriptFolder = "E:\\Dropbox\\proj\\stringbit\\script";
	string scriptFolder = "D:\\proj\\stringbit\\temp";
	ScriptGenerator sc;

	Stopwatch watch;
	
	watch.Start();
	bool invert = false;
	sc.OutputHamToMatlab(3, invert, scriptFolder);
	sc.OutputHamToMatlab(4, invert, scriptFolder);
	sc.OutputHamToMatlab(5, invert, scriptFolder);
	sc.OutputHamToMatlab(6, invert, scriptFolder);
	sc.OutputHamToMatlab(7, invert, scriptFolder);
	sc.OutputHamToMatlab(8, invert, scriptFolder);
	sc.OutputHamToMatlab(9, invert, scriptFolder);
	sc.OutputHamToMatlab(10, invert, scriptFolder);
	sc.OutputHamToMatlab(11, invert, scriptFolder);

	//sc.OutputNormToMatlab(3, scriptFolder, false);
	//sc.OutputNormToMatlab(4, scriptFolder, false);
	//sc.OutputNormToMatlab(5, scriptFolder, false);
	//sc.OutputNormToMatlab(6, scriptFolder, false);
	//sc.OutputNormToMatlab(7, scriptFolder, false);
	//sc.OutputNormToMatlab(8, scriptFolder, false);
	//sc.OutputNormToMatlab(9, scriptFolder, false);
	//sc.OutputNormToMatlab(10, scriptFolder, true);
	//sc.OutputNormToMatlab(11, scriptFolder, true);

	cout << "Time: " << watch.Stop() << " seconds." << endl;
}

void GenerateLaTeX()
{
	//ScriptGenerator sg;
	//string filename = "E:\\Dropbox\\proj\\stringbit\\boson-states-1-7.tex";
	//sg.OutputStateToLaTeX(1, 7, 1, filename);
	//filename = "E:\\Dropbox\\proj\\stringbit\\fermion-states-1-7.tex";
	//sg.OutputStateToLaTeX(1, 7, 2, filename);
	//OutputHamToLatex();
	//sg.OutputNormToLaTeX(2, 7, "E:\\Dropbox\\proj\\stringbit\\norms-2-7.tex");
}

int main()
{
	GenerateStates();
	//CalculateNorm(3, false);
	//CalculateNorm(4, false);
	//CalculateNorm(5, false);
	//CalculateNorm(6, false);
	//CalculateNorm(7, false);
	//CalculateNorm(8, false);
	//CalculateNorm(9, false);
	//CalculateNorm(10, false);
	//TestNormCalculator();
	
	//TraceState state;
	//state.AddTrace(SingleTrace(1, 2));
	//state.AddTrace(SingleTrace(1, 1));
	//TestHamiltonianMultiTrace(state);
	//OutputHamiltonianMatrix(2);
	//OutputHamiltonianMatrix(3);
	//OutputHamiltonianMatrix(4);
	//OutputHamiltonianMatrix(5);
	//OutputHamiltonianMatrix(6);
	//OutputHamiltonianMatrix(7);
	//OutputHamiltonianMatrix(11);

	//GenerateLaTeX();
	GenerateMatlabScript();
	return 0;
}