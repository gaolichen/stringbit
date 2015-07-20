#pragma warning(disable:4018)
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "StateCollection.h"
#include "TraceState.h"
#include "NormCalculator.h"
#include "StateGenerator.h"
#include "BitUtility.h"

using namespace std;

map<pair<int, int>,  double> vanishN;
map<pair<int, int>,  double> vanishN2;

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

	void PrintTable(vector<StateInfo>& list, ofstream& os, int order, int limit = 100000)
	{
		StateCollection* sc = StateCollection::Inst();
		BruteForceCalculator nc;
		if (list.size() >= 30)
		{
			os << "\\newpage{}" << endl << endl;
		}

		os << "\\begin{itemize}" << endl;
                os << "\\item Trace states of $\\mathcal{O}\\left(";
		if (order == 0)
		{
			os << 1;
		}
		else
		{
			os << "\\frac{1}{N}";
		}
		os  << "\\right)$ contribution to M=" << M << " ground state";

                if (leadingOrder.size() > 50)
                {
                        os << " (Only top 50 are listed)";
                }
                os << endl;
                os << "\\end{itemize}" << endl;

		// create table
                os << "\\begin{center}" << endl;
                os << "\\begin{tabular}{|c|c|c|c|}" << endl;
                os << "\\hline " << endl;
                os << "Trace State & $\\bar{b}$ Number & State Norm & Amplitude\\tabularnewline" << endl;
                os << "\\hline " << endl;
                os << "\\hline " << endl;
		int count = min(limit, (int)list.size());

                for (int i = 0; i < count; i++)
                {
                        StateInfo& info = list[i];
                        TraceState ts = sc->GetBosonState(M, info.Id - 1);
                        os << "$" << ts.ToLaTeX() << "$ & " << ts.FermionNumber();
			os << " & $" << nc.Calculate(ts, ts).ToLaTeX(M);
                        os << "$ & " << info.Amplitude << "\\tabularnewline" << endl;
                        os << "\\hline " << endl;
                }

                os << "\\end{tabular}" << endl;
                os << "\\par\\end{center}" << endl << endl;
	}

	void ToLatex(ofstream& os)
	{
		os << "\\subsubsection*{$M=" << M << ",\\, E=" << E0 << "$}" << endl;
		os << "Ground state vanishes at $N=";
		double n = vanishN[make_pair(Xi, M)];
		if (abs(n - floor(n)) < 1e-8)
		{
			os <<(int)n;
		}
		else
		{
			os << (int)floor(n) << ".5";
		}
		os << "$." << endl;

		// create table
		PrintTable(leadingOrder, os, 0, 50);
		PrintTable(firstOrder, os, 1, 50);
	}
};

void PrintHeader(int xi, ofstream& os)
{
	os << "\\documentclass[english]{article}" << endl;
	os << "\\usepackage[T1]{fontenc}" << endl;
	os << "\\usepackage[latin9]{inputenc}" << endl;
	os << "\\usepackage{geometry}" << endl;
	os << "\\geometry{verbose,tmargin=0.4cm,bmargin=0.4cm}" << endl;

	os << "\\makeatletter" << endl;
	os << "\\providecommand{\\tabularnewline}{\\\\}" << endl;
	os << "\\makeatother" << endl;
	os << "\\usepackage{babel}" << endl;
	os << "\\begin{document}" << endl << endl;

	os << "\\title{Ground States $\\xi=" << xi << "$}" << endl;
	os << "\\maketitle" << endl;
}

void PrintEnd(ofstream& os)
{
	os << "\\end{document}" << endl;
}

void PrintHam(int xi, ofstream& os)
{
	os << "The Hamiltonian is" << endl << endl;
	os << "\\begin{eqnarray*}" << endl;
	os << "H & = & \\frac{2}{N}\\mathrm{Tr}\\left[\\left(\\bar{a}^{2}-i\\bar{b}^{2}\\right)a^{2}-\\left(\\bar{b}^{2}" << endl;
	os << "i\\bar{a}^{2}\\right)b^{2}+\\left(\\bar{a}\\bar{b}+\\bar{b}\\bar{a}\\right)ba+\\left(\\bar{a}\\bar{b}-" << endl;
	os << "\\bar{b}\\bar{a}\\right)ab\\right]\\\\" << endl;
	if (xi != 0)
	{
		os << " &  &";
		if (xi > 0) os << " + ";
		else os << " - ";
		os << "\\frac{2";
		os << abs(xi);
		os <<  "}{N}\\mathrm{Tr}\\left[\\bar{a}\\bar{b}ba+\\bar{b}\\bar{a}ab+\\bar{a}^{2}a^{2}+\\bar{b}^{2}b^{2}-" << endl;
		os << "M\\right]\\\\" << endl;
	}

	if (xi < -1)
	{
		os << " &  & +\\frac{" << -(2*xi + 2) << "}{N}" << endl;
		os << "\\mathrm{Tr}\\left[\\bar{a}a\\bar{a}a+\\bar{b}b\\bar{a}a-\\bar{a}b\\bar{b}a\\right]" << endl;
	}
	
	os << "\\end{eqnarray*}" << endl;
}


void ReadVanishN(char *file)
{
	ifstream ifs(file);
	char line[256];

	int xi;
	while (ifs >> xi)
	{
		ifs.getline(line, 256);
		ifs.getline(line, 256);
		for (int i = 3; i <= 11; i++)
		{
			int bit;
			double n;
			ifs >> bit >> n;
			vanishN[make_pair(xi, bit)] = n;
		}
		
		ifs.getline(line, 256);
		ifs.getline(line, 256);
		for (int i = 3; i <= 11; i++)
		{
			int bit;
			double n;
			ifs >> bit >> n;
			vanishN2[make_pair(xi, bit)] = n;
		}
	}

	ifs.close();
}

int main(int argc, char *argv[])
{
	if (argc < 4)
	{
		cout << "Argument required!" << endl;
		return -1;
	}

	StateGenerator generator;
	generator.GenerateAllStates();
	generator.InitStateCollection(StateCollection::Inst());

	ifstream ifs(argv[1]);
	ReadVanishN(argv[2]);
	cout << argv[1] << endl;
	cout << argv[2] << endl;
	cout << argv[3] << endl;
	string folder = argv[3];
	
	int xi;
	while (ifs >> xi)
	{
		string file = folder + "/xi=" + ToString(xi) + ".tex";

		ofstream ofs(file.c_str());
		int stateId;
		int bit, size1, size2;
		double e0, e1, amplitude;
		cout << xi << endl;
		PrintHeader(xi, ofs);
		PrintHam(xi, ofs);
		ofs << endl;
		for (int i = 3; i <= 11; i++)
		{
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

			list.ToLatex(ofs);
		}
		
		PrintEnd(ofs);
		ofs.close();
		//break;
	}

	ifs.close();

	return 0;
}

