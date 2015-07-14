#pragma warning(disable:4018)
#include "ScriptGenerator.h"
#include "Hamiltonian.h"
#include "BitUtility.h"
#include "NormCalculator.h"
#include <fstream>
using namespace std;

ScriptGenerator::ScriptGenerator(string rootFolder, StateType stateType)
{
	if (stateType == Boson)
	{
		this->rootFolder = rootFolder + "\\boson";
	}
	else
	{
		this->rootFolder = rootFolder + "\\fermion";
	}

	this->type = stateType;
}

void ScriptGenerator::OutputHamToMatlab(int bits, bool invert)
{
	Hamiltonian ham(1, invert);
	vector<vector<Coefficient> > rem, imm;
	ham.Matrix(bits, type, rem, imm);

	string function;
	if (!invert)
	{
		function = "ham" + ToString(bits);
	}
	else
	{
		function = "pham" + ToString(bits);
	}

	string file = rootFolder + "\\" + function + ".m";
	ofstream ofs(file.c_str());
	ofs << "function f = " << function << "(N)" << endl;
	ofs << "f=[" << endl;
	
	for (int j = 0; j < rem.size(); j++)
	{
		ofs << "   ";
		for (int k = 0; k < rem[j].size(); k++)
		{
			ofs << ' ';
			if (!imm[j][k].IsZero())
			{
				ofs << imm[j][k] << "*1i";
			}
			else
			{
				ofs << rem[j][k];
			}
		}

		if (j != rem.size() - 1)
		{
			ofs << ";" << endl;
		}
	}

	ofs << "];" << endl;
	ofs << "end" << endl;


	ofs.close();
}

void ScriptGenerator::OutputHamToLaTeX(vector<int> bits, vector<double> scales, string filename)
{
	ofstream ofs(filename.c_str());
	ofs << "\\documentclass[english]{article}" << endl;
	ofs << "\\usepackage[T1]{fontenc}" << endl;
	ofs << "\\usepackage{graphicx}" << endl;
	ofs << "\\usepackage{geometry}" << endl;
	ofs << "\\geometry{verbose,tmargin=0.1cm,bmargin=0.1cm,lmargin=0.1cm,rmargin=0.1cm}" << endl;
	ofs << "\\usepackage{pdflscape}" << endl;
	//ofs << "\\usepackage{lscape}" << endl;
	ofs << "\\usepackage[latin9]{inputenc}" << endl;
	ofs << "\\usepackage{babel}" << endl;
	ofs <<"\\begin{document}" << endl << endl;
	ofs << "\\newcommand\\scalemath[2]{\\scalebox{#1}{\\mbox{\\ensuremath{\\displaystyle #2}}}}";
	ofs << endl << endl;
	ofs << "\\begin{landscape}" << endl;

	Hamiltonian ham;
	
	for (int i = 0; i < bits.size(); i++)
	{
		ofs << "\\subsection*{"<< bits[i]<< " bits:}" << endl << endl;
		ofs << "\\[" << endl;

		vector<vector<Coefficient> > rem, imm;
		ham.Matrix(bits[i], type, rem, imm);
		
		string s(rem.size(), 'c');
		ofs << "\\left(\\scalemath{" << scales[i] <<"}" << endl;
		ofs << "{\\begin{array}{"<< string(rem.size(), 'c') << "}" << endl;
		for (int j = 0; j < rem.size(); j++)
		{
			for (int k = 0; k < rem[j].size(); k++)
			{
				if (k > 0)
				{
					ofs << "& ";
				}
				if (!imm[j][k].IsZero())
				{
					ofs << imm[j][k].ToLatex() << "i ";
				}
				else
				{
					ofs << rem[j][k].ToLatex() << " ";
				}
			}

			if (j != rem.size() - 1)
			{
				ofs << "\\\\";
			}
			ofs << endl;
		}

		ofs << "\\end{array}}\\right)" << endl;
		ofs << "\\]" << endl << endl << endl;
	}

	ofs << "\\end{landscape}" << endl;

	ofs << "\\end{document}" << endl;
	ofs.close();
}

void ScriptGenerator::GroupStates(int bits, vector<vector<int> >& states)
{
	StateCollection* sc = StateCollection::Inst();
	int n = sc->StateNumber(bits);
	int b = 0;
	if (type == Fermion)
	{
		b = 1;
	}

	vector<int> v;
	if (bits < 5)
	{
		for (int i = 0; i < n; i++)
		{
			v.push_back(i);
		}
		states.push_back(v);
	}
	else
	{
		for (int j = 0; j < n; j++)
		{
			if (sc->GetState(bits, j, type).FermionNumber() == b)
			{
				v.push_back(j);
			}
			else
			{
				b = sc->GetState(bits, j, type).FermionNumber();
				states.push_back(v);
				v.clear();
				v.push_back(j);
			}
		}

		states.push_back(v);
	}
}

double ScriptGenerator::NormMatrixScale(int bits, int index)
{
	map<int, vector<double> > scales;
	double a[4] = {1.0, 0.6, 0.79, 1.0};
	double b[4] = {0.51, 0.25, 0.26, 1.0};
	scales[6] = vector<double>(a, a + 4);
	scales[7] = vector<double>(b, b + 4);

	if (scales.find(bits) == scales.end())
	{
		return 1.0;
	}
	else
	{
		return scales[bits][index];
	}
}

void ScriptGenerator::OutputNormToLaTeX(int minBits, int maxBits, string filename)
{
	StateCollection* sc = StateCollection::Inst();
	BruteForceCalculator calc;

	ofstream ofs(filename.c_str());

	ofs << "\\documentclass[english]{article}" << endl;
	ofs << "\\usepackage[T1]{fontenc}" << endl;
	ofs << "\\usepackage{graphicx}" << endl;
	ofs << "\\usepackage{geometry}" << endl;
	ofs << "\\geometry{verbose,tmargin=0.1cm,bmargin=0.1cm,lmargin=0.1cm,rmargin=0.1cm}" << endl;
	ofs << "\\usepackage{pdflscape}" << endl;
	//ofs << "\\usepackage{lscape}" << endl;
	ofs << "\\usepackage[latin9]{inputenc}" << endl;
	ofs << "\\usepackage{babel}" << endl;
	ofs <<"\\begin{document}" << endl << endl;
	ofs << "\\newcommand\\scalemath[2]{\\scalebox{#1}{\\mbox{\\ensuremath{\\displaystyle #2}}}}";
	ofs << endl << endl;
	ofs << "\\begin{landscape}" << endl;

	ofs << "\\title{States norm matrices}" << endl << endl;
	ofs << "\\maketitle" << endl << endl;

	for (int i = minBits; i <= maxBits; i++)
	{
		ofs << "\\subsection*{" << i << " bits:}" << endl;
		vector<vector<int> > states;
		GroupStates(i, states);
		if (states.size() > 1)
		{
			ofs << "\\begin{itemize}" << endl;
		}

		for (int j = 0; j < states.size(); j++)
		{
			if (states.size() > 1)
			{
				ofs << "\\item State " << states[j][0] + 1 << " to State " << states[j].back() + 1 << ": " << endl; 
			}
			
			ofs << "\\[" << endl;
			ofs << "\\left(\\scalemath{" << NormMatrixScale(i, j) <<"}" << endl;
			ofs << "{\\begin{array}{" << string(states[j].size(), 'c') << "}" << endl;

			for (int k = 0; k < states[j].size(); k++)
			{
				for (int h = 0; h < states[j].size(); h++)
				{
					if (h > 0)
					{
						ofs << " & "; 
					}
					ofs << calc.Calculate(sc->GetState(i, states[j][k], type), sc->GetState(i, states[j][h], type)).ToLaTeX(i);
				}

				if (k + 1 != states[j].size())
				{
					ofs << "\\\\";
				}

				ofs << endl;
			}

			ofs << "\\end{array}}\\right)" << endl;
			ofs << "\\]" << endl;
		}

		if (states.size() > 1)
		{
			ofs << "\\end{itemize}" << endl;
		}
	}

	ofs << "\\end{landscape}" << endl;
	ofs << "\\end{document}" << endl;
	ofs.close();
}

void ScriptGenerator::OutputStateToLaTeX(int minBit, int maxBit, string filename)
{
	ofstream ofs(filename.c_str());
	ofs << "\\documentclass[english]{article}" << endl;
	ofs << "\\usepackage[T1]{fontenc}" << endl;

	ofs << "\\usepackage[latin9]{inputenc}" << endl;
	ofs << "\\usepackage{babel}" << endl;
	ofs <<"\\begin{document}" << endl;

	for (int i = minBit; i <= maxBit; i++)
	{
		ofs << "\\subsection*{" << i << " bits}" << endl << endl;
		if (type == Boson)
		{
			ofs << StateCollection::Inst()->StateNumber(i) << " bosonic states" << endl;
			ofs << "\\begin{itemize}" << endl;
			for (int j = 0; j < StateCollection::Inst()->StateNumber(i); j++)
			{
				ofs << "\\item $\\left|" << j + 1 << "\\right\\rangle=";
				ofs << StateCollection::Inst()->GetBosonState(i, j).ToLaTeX() <<"$" <<  endl;
			}
			ofs << "\\end{itemize}" << endl;
		}

		if (type == Fermion)
		{
			ofs << StateCollection::Inst()->StateNumber(i) << " fermionic states" << endl;
			ofs << "\\begin{itemize}" << endl;
			for (int j = 0; j < StateCollection::Inst()->StateNumber(i); j++)
			{
				ofs << "\\item $\\left|" << j + 1 << "\\right\\rangle=";
				ofs << StateCollection::Inst()->GetFermionState(i, j).ToLaTeX() <<"$" <<  endl;
			}
			ofs << "\\end{itemize}" << endl;
		}
	}
	ofs << "\\end{document}" << endl;

	ofs.close();
}

void ScriptGenerator::OutputNormToMatlab(int bits)
{
	StateCollection* sc = StateCollection::Inst();
	BruteForceCalculator calc;
	vector<vector<int> > states;
	GroupStates(bits, states);

	bool output = (bits >= 10);
	string function = "norm" + ToString(bits);
	string file = rootFolder + "\\" + function + ".m";
	ofstream ofs(file.c_str());
	ofs << "function f = " << function << "(N)" << endl;
	for (int i = 0; i < states.size(); i++)
	{
		if (output)
		{
			cout << "Calculating state " << i + 1 << "..." << endl;
		}

		ofs << (char)((int)'A' + i) << "=[" << endl;
		for (int k = 0; k < states[i].size(); k++)
		{
			ofs << "    ";
			for (int h = 0; h < states[i].size(); h++)
			{
				if (h > 0)
				{
					ofs << ", ";
				}
				ofs << calc.Calculate(sc->GetState(bits, states[i][k], type), sc->GetState(bits, states[i][h], type)).ToMatlab(bits);
			}

			if (k + 1 != states[i].size())
			{
				ofs << ";" << endl;
			}
		}

		ofs << endl << "];" << endl << endl;
	}

	if (states.size() == 1)
	{
		ofs << "f=A;" << endl;
	}
	else
	{
		ofs << "f=blkdiag(";
		for (int i = 0; i < states.size(); i++)
		{
			if (i > 0)
			{
				ofs << ", ";
			}

			ofs << (char)((int)'A' + i);
		}

		ofs << ");" << endl;
	}
	
	ofs << "end" << endl;
}
