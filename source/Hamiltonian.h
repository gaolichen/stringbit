#pragma once
#include "HamOperator.h"
#include "StateCollection.h"
#include "SingleTrace.h"
#include "TraceState.h"
#include "StateType.h"

class Hamiltonian
{
private:
	vector<HamOperator*> realOps;
	vector<HamOperator*> imaginaryOps;
	vector<int> rePrefactors;
	vector<int> imPrefactor;
	void Init(int xi1, bool invert);
public:
	Hamiltonian();
	Hamiltonian(int xi1, bool invert);
	~Hamiltonian();

	void Apply(const TraceState& state, MixState& real, MixState& imaginary);
	void Matrix(int bits, StateType type, vector<vector<Coefficient> >& rem, vector<vector<Coefficient> >& imm);
};