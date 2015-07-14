#pragma warning(disable:4018)
#include "Hamiltonian.h"

Hamiltonian::Hamiltonian()
{
	Init(1, false);
}

Hamiltonian::Hamiltonian(int xi1, bool invert)
{
	Init(xi1, invert);
}

void Hamiltonian::Init(int xi1, bool invert)
{
	int a = 1;
	if (invert)
	{
		a = -1;
	}

	realOps.push_back(new HamOperatorA(HamOperator::AA, HamOperator::AA));
	rePrefactors.push_back((2 + 2* xi1 - 2) * a);
	realOps.push_back(new HamOperatorA(HamOperator::BB, HamOperator::BB));
	rePrefactors.push_back((-2 + 2* xi1 - 2) * a);
	realOps.push_back(new HamOperatorA(HamOperator::AB, HamOperator::BA));
	rePrefactors.push_back(2 * xi1 * a);
	realOps.push_back(new HamOperatorA(HamOperator::BA, HamOperator::BA));
	rePrefactors.push_back(2 * a);
	realOps.push_back(new HamOperatorA(HamOperator::AB, HamOperator::AB));
	rePrefactors.push_back(2 * a);
	realOps.push_back(new HamOperatorA(HamOperator::BA, HamOperator::AB));
	rePrefactors.push_back(2 * (xi1 - 2) * a);

	realOps.push_back(new BitNumberHamOperator());
	rePrefactors.push_back(-(2 * xi1 - 2) * a);

	imaginaryOps.push_back(new HamOperatorA(HamOperator::BB, HamOperator::AA));
	imPrefactor.push_back(-2 * a);
	imaginaryOps.push_back(new HamOperatorA(HamOperator::AA, HamOperator::BB));
	imPrefactor.push_back(2 * a);

	if (invert)
	{
		int i = (2 + 2* xi1 - 2) * a;
		realOps.push_back(new HamOperatorB(HamOperator::AA, HamOperator::AA));
		rePrefactors.push_back(-i);
		realOps.push_back(new HamOperatorB(HamOperator::BA, HamOperator::BA));
		rePrefactors.push_back(-i);
		realOps.push_back(new HamOperatorB(HamOperator::AB, HamOperator::BA));
		rePrefactors.push_back(i);
	}
}

Hamiltonian::~Hamiltonian()
{
	for (int i = 0; i < realOps.size(); i++)
	{
		delete realOps[i];
	}

	for (int i = 0; i < imaginaryOps.size(); i++)
	{
		delete imaginaryOps[i];
	}
}

void Hamiltonian::Apply(const TraceState& state, MixState& real, MixState& imaginary)
{
	for (int i = 0; i < realOps.size(); i++)
	{
		MixState res;
		realOps[i]->ApplyOn(state, res);
		real.Merge(res, rePrefactors[i]);
	}

	for (int i = 0; i < imaginaryOps.size(); i++)
	{
		MixState res;
		imaginaryOps[i]->ApplyOn(state, res);
		imaginary.Merge(res, imPrefactor[i]);
	}
}

void Hamiltonian::Matrix(int bits, StateType type, vector<vector<Coefficient> >& rem, vector<vector<Coefficient> >& imm)
{
	int n = StateCollection::Inst()->StateNumber(bits);
	map<StateId, Coefficient>::const_iterator it;
	rem.resize(n, vector<Coefficient>());
	imm.resize(n, vector<Coefficient>());
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			rem[i].push_back(Coefficient(0, 0));
			imm[i].push_back(Coefficient(0, 0));
		}
	}

	for (int i = 0; i < n; i++)
	{
		const TraceState& state = StateCollection::Inst()->GetState(bits, i, type);
		MixState re, im;
		Apply(state, re, im);
		for (it = re.Begin(); it != re.End(); ++it)
		{
			rem[it->first.Index / 2][i] = it->second;
		}

		for (it = im.Begin(); it != im.End(); ++it)
		{
			imm[it->first.Index / 2][i] = it->second;
		}
	}
}