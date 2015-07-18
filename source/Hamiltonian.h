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
	bool inverted;
	void Init(int xi);
public:
	Hamiltonian();
	Hamiltonian(int xi, bool invert);
	~Hamiltonian();

	inline bool Inverted() { return inverted; };
	void Apply(const TraceState& state, MixState& real, MixState& imaginary);
	void Matrix(int bits, StateType type, vector<vector<Coefficient> >& rem, vector<vector<Coefficient> >& imm);
	void AddReadOp(HamOperator *op, int prefactor);
	void AddImaginaryOp(HamOperator *op, int prefactor);
};