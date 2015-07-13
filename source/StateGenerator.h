#pragma once
#include "SingleTrace.h"
#include "TraceState.h"
#include "BitUtility.h"
#include "StateCollection.h"
#include <vector>

using namespace std;

class StateGenerator
{
private:
	bool* myFlags;
	bool visited[30][30][2];

	vector<vector<vector<i64> > > stateNumbers;
	vector<i64> singleTraceNumbers;
	vector<vector<SingleTrace> > singleFermions;
	vector<vector<SingleTrace> > singleBosons;
	vector<vector<vector<TraceState> > > fermions;
	vector<vector<vector<TraceState> > > bosons;

	void InitSingleTraceNumber();
	void FindSingleStates(int n);
	i64 StateNumbers(int bit, int remain, int b);
	void GeneratSingleStates();
	static void DoPickFermion(vector<SingleTrace>& allstates, int index, int remain, vector<SingleTrace>& curr, vector<TraceState>& ret);
	static void DoPickBoson(vector<SingleTrace>& allstates, int index, int remain, vector<SingleTrace>& curr, vector<TraceState>& ret);
	static vector<TraceState> PickFermionFromSingleState(vector<SingleTrace> &states, int number);
	static vector<TraceState> PickBosonFromSingleState(vector<SingleTrace> &states, int number);
	static TraceState CombineStates(TraceState& a, TraceState& b, TraceState& c);
public:
	const static int MAX_BIT_TO_COUNT = 62;
	const static int MAX_BIT_TO_GENERATE = 11;

	StateGenerator();
	~StateGenerator();
	
	void GenerateAllStates();
	i64 BosonNumber(int n);
	i64 FermionNumber(int n);
	i64 SingleStateNumber(int n);
	TraceState BosonState(int n, int index);
	TraceState FermionState(int n, int index);
	void InitStateCollection(StateCollection* collection);
};