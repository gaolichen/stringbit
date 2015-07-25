#pragma once
#include "HamOperator.h"
#include "StateCollection.h"
#include "SingleTrace.h"
#include "TraceState.h"
#include "StateType.h"

class Hamiltonian
{
protected:
	vector<HamOperator*> realOps;
	vector<HamOperator*> imaginaryOps;
	vector<int> rePrefactors;
	vector<int> imPrefactor;
private:
	bool inverted;
	void Init(int xi);
public:
	Hamiltonian();
	Hamiltonian(int xi, bool invert);
	~Hamiltonian();

	inline bool Inverted() { return inverted; };
	HamOperator* RealOp(int);
	int RealOpSize();
	void Apply(const TraceState& state, MixState& real, MixState& imaginary);
	void Matrix(int bits, StateType type, vector<vector<Coefficient> >& rem, vector<vector<Coefficient> >& imm);
	void AddReadOp(HamOperator *op, int prefactor);
	void AddImaginaryOp(HamOperator *op, int prefactor);
	friend ostream& operator<<(ostream& os, const Hamiltonian& ham);
};

class ZeroHamiltonian : public Hamiltonian
{
public:
	ZeroHamiltonian();
};
