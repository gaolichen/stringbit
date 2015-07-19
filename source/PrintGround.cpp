#pragma warning(disable:4018)
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "StateCollection.h"
#include "TraceState.h"
#include "NormCalculator.h"
#include "StateGenerator.h"

using namespace std;

struct StateInfo
{
	int Id;
	double Amplitude;
	StateInfo(int id, double amplitude)
	{
		Id = id;
		Amplitude = amplitude;
	}
};

struct StateList
{
	int Xi;
	int M;
	double E0;
	double E1;
	vector<StateInfo> leadingOrder;
	vector<StateInfo> firstOrder;
	StateList(int xi, int m)
	{
		Xi = xi;
		M = m;
	}

	void AddLeadingOrder(int id, double amplitude)
	{
		leadingOrder.push_back(StateInfo(id, amplitude));
	}

	void AddFirstOrder(int id, double amplitude)
	{
		firstOrder.push_back(StateInfo(id, amplitude));
	}

	void ToLatex(ofstream& os)
	{
		StateCollection* sc = StateCollection::Inst();
		BruteForceCalculator nc;
		if (M >= 8)
		{
			os << "\\newpage{}" << endl;
		}
		os << "\\subsubsection*{$M=" << M << "$}" << endl;
		os << "\\begin{itemize}" << endl;
		os << "\\item Leading order ground energy: $" << E0 << "$. " << endl;
		os << "\\item Leading order ground state" << endl;
		os << "\\end{itemize}" << endl;
		// create table
		os << "\\begin{center}" << endl;
		os << "\\begin{tabular}{|c|c|c|}" << endl;
		os << "\\hline " << endl;
		os << "Trace State & State Norm & Amplitude\\tabularnewline" << endl;
		os << "\\hline " << endl;
		os << "\\hline " << endl;

		for (int i = 0; i < leadingOrder.size(); i++)
		{
			StateInfo& info = leadingOrder[i];
			TraceState ts = sc->GetBosonState(M, info.Id - 1);
			os << "$" << ts.ToLaTeX() << "$ & $" << nc.Calculate(ts, ts).ToLaTeX(M);
			os << "$ & " << info.Amplitude << "\\tabularnewline" << endl;
			os << "\\hline " << endl;
		}

		os << "\\end{tabular}" << endl;
		os << "\\par\\end{center}" << endl << endl;
	}
};

void PrintHeader(ofstream& os)
{
	os << "\\documentclass[english]{article}" << endl;
	os << "\\usepackage[T1]{fontenc}" << endl;
	os << "\\usepackage[latin9]{inputenc}" << endl;

	os << "\\makeatletter" << endl;
	os << "\\providecommand{\\tabularnewline}{\\\\}" << endl;
	os << "\\makeatother" << endl;
	os << "\\usepackage{babel}" << endl;
	os << "\\begin{document}" << endl << endl;

	os << "\\title{Ground States}" << endl;
	os << "\\maketitle" << endl;
}

void PrintEnd(ofstream& os)
{
	os << "\\end{document}" << endl;
}

int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		cout << "Argument required!" << endl;
		return -1;
	}

	StateGenerator generator;
	generator.GenerateAllStates();
	generator.InitStateCollection(StateCollection::Inst());

	ifstream ifs(argv[1]);
	cout << argv[1] << endl;
	cout << argv[2] << endl;
	ofstream ofs(argv[2]);
	
	PrintHeader(ofs);

	char line[256];
	int xi;
	ifs.getline(line, 256);
	while (ifs >> xi)
	{
		int stateId;
		int bit, size1, size2;
		double e0, e1, amplitude;
		cout << xi << endl;
		ofs << "\\subsection*{$\\xi=" << xi << "$}" << endl << endl;
		ofs << "\\[" << endl;
		ofs << "H=aa" << endl;
		ofs << "\\]" << endl;
		for (int i = 3; i <= 11; i++)
		{
			cout << i << endl;
			ifs >> bit >> size1 >> size2;
			ifs >> e0 >> e1;
			StateList list(xi, bit);
			list.E0 = e0;
			list.E1 = e1;

			for (int j = 0; j < size1; j++)
			{
				ifs >> stateId >> amplitude;
				list.AddLeadingOrder(stateId, amplitude);
			}

			for (int j = 0; j < size2; j++)
			{
				ifs >> stateId >> amplitude;
				list.AddFirstOrder(stateId, amplitude);
			}

			if (i <= 11)
			{
				list.ToLatex(ofs);
			}
		}
	}

	PrintEnd(ofs);

	ifs.close();
	ofs.close();

	return 0;
}

