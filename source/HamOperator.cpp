#include "HamOperator.h"
#include "SingleTrace.h"
#include <iostream>
#include "StateId.h"

using namespace std;

HamOperator::HamOperator()
{
	// ..
}

void HamOperator::ApplyOn(const TraceState& state, MixState& res)
{
	int f = 0, g = 0;
	for (int i = 0; i < state.TraceNumber(); i++)
	{
		TraceState toAppend;
		state.CopyTo(toAppend, i);
		MixState ms = ApplyOn(state.Trace(i));
		ms.Extend(toAppend, f * state.Trace(i).FermiNumber());
		res.Merge(ms);
		g = f;
		for (int j = i + 1; j < state.TraceNumber(); j++)
		{
			TraceState toAdd;
			state.CopyTo(toAdd, i, j);
			MixState ms2 = ApplyOn(state.Trace(i), state.Trace(j));
			ms2.Extend(toAdd, f * state.Trace(i).FermiNumber() + g * state.Trace(j).FermiNumber());
			res.Merge(ms2);
			g += state.Trace(j).FermiNumber();
		}
		f += state.Trace(i).FermiNumber();
	}
}

MixState& HamOperator::ApplyOn(const SingleTrace& single)
{
	map<SingleTrace, MixState>::iterator it = cache1.find(single);
	if (it != cache1.end())
	{
		return it->second;
	}

	MixState res;
	ApplyOnSingle(single, res);

	cache1[single] = res;
	return cache1[single];
}

MixState& HamOperator::ApplyOn(const SingleTrace& single1, const SingleTrace& single2)
{
	pair<SingleTrace, SingleTrace> key = make_pair(single1, single2);

	map<pair<SingleTrace, SingleTrace>, MixState>::iterator it = cache2.find(key);
	if (it != cache2.end())
	{
		return it->second;
	}

	MixState res;
	ApplyOnTwoSingle(single1, single2, res);

	cache2[key] = res;
	return cache2[key];
}

void HamOperator::AddState(TraceState& state, int parity, MixState& res, bool decreaseOrder)
{
	Coefficient coef = state.Normalize();
	StateId id = StateCollection::Inst()->GetId(state);
	if (!id.IsValid())
	{
		return;
	}

	if (parity % 2 == 1)
	{
		coef.Opposite();
	}
	
	if (decreaseOrder)
	{
		coef.DecreaseOrder();
	}

	res.AddState(id, coef);
}

BitNumberHamOperator::BitNumberHamOperator()
{
}

void BitNumberHamOperator::ApplyOn(const TraceState& state, MixState& res)
{
	StateId id = StateCollection::Inst()->GetId(state);
	res.AddState(id, Coefficient(0, state.TotalBits()));
}

void BitNumberHamOperator::ApplyOnSingle(const SingleTrace& single, MixState& res)
{
	// for single trace of length 1, the operator do not make contribution
	if (single.BitNumber() == 1)
	{
		return;
	}

	SingleTrace a, b;
	int n = single.BitNumber();
	for (int i = 0; i < n; i++)
	{
		single.Split(i, a, b);
		TraceState state;
		int parity = 0;
		if (single.Bit(i) == 0)
		{
			// if ith bit is bosonic
			state.AddTrace(SingleTrace(1, 0));
		}
		else
		{
			// if ith bit is fermionic
			state.AddTrace(SingleTrace(1, 1));
			parity += a.FermiNumber();
		}

		state.AddTrace(SingleTrace::Merge(a, b));
		// TODO: need to consider 1/N factor.
		AddState(state, parity, res, true);
	}
}

void BitNumberHamOperator::ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res)
{
	//.. do nothing.
}

HamOperatorA::HamOperatorA(int creator, int annihilator)
{
	this->creator = creator;
	this->annihilator = annihilator;
}

void HamOperatorA::ApplyOnSingle(const SingleTrace& single, MixState& res)
{
	int parity;
	int n = single.BitNumber();
	SingleTrace b, c, d;
	SingleTrace a(creator, 2);
	int a2 = (this->annihilator & 1);
	int a1 = this->annihilator / 2;
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			if (single.Bit(i) == a1 && single.Bit(j) == a2)
			{
				single.Split(j, i, b, c, d);
				TraceState state;
				state.AddTrace(SingleTrace::Merge(a, d, b));
				state.AddTrace(c);
				parity = (d.FermiNumber() + a1) * (b.FermiNumber() + c.FermiNumber());
				parity += a2 * b.FermiNumber();
				AddState(state, parity, res);
			}

			if (single.Bit(i) == a2 && single.Bit(j) == a1)
			{
				single.Split(j, i, b, c, d);
				TraceState state2;
				state2.AddTrace(SingleTrace::Merge(a, c));
				state2.AddTrace(SingleTrace::Merge(d, b));
				parity = b.FermiNumber() * (c.FermiNumber() + d.FermiNumber());
				parity +=  b.FermiNumber() * (a1 + a2) + c.FermiNumber() * a2 + a1 * a2;
				AddState(state2, parity, res);
			}
		}
	}
}

void HamOperatorA::ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res)
{
	int parity;
	int a1 = annihilator / 2;
	int a2 = (annihilator & 1);
	SingleTrace a(creator, 2);
	SingleTrace b, c, d, e, f;
	for (int i = 0; i < single1.BitNumber(); i++)
	{
		if (single1.Bit(i) != a1 && single1.Bit(i) != a2) 
		{
			continue;
		}

		single1.Split(i, b, c);
		for (int j = 0; j < single2.BitNumber(); j++)
		{
			if (single1.Bit(i) == a2 && single2.Bit(j) == a1)
			{
				//..
				single2.Split(j, d, e);
				f = SingleTrace::Merge(a, e, d);
				f = SingleTrace::Merge(f, c, b);
				parity = e.FermiNumber() * (b.FermiNumber() + c.FermiNumber() + d.FermiNumber()) +
					d.FermiNumber() * (b.FermiNumber() + c.FermiNumber()) + b.FermiNumber() * c.FermiNumber();
				parity += a2 * b.FermiNumber() + a1 * (b.FermiNumber() + c.FermiNumber() + d.FermiNumber());
				TraceState state;
				state.AddTrace(f);
				AddState(state, parity, res);
			}

			if (single1.Bit(i) == a1 && single2.Bit(j) == a2)
			{
				//..
				single2.Split(j, d, e);
				f = SingleTrace::Merge(a, c, b);
				f = SingleTrace::Merge(f, e, d);
				parity = b.FermiNumber() * c.FermiNumber() + d.FermiNumber() * e.FermiNumber();
				parity += a1 * b.FermiNumber() + a2 * (b.FermiNumber() + c.FermiNumber() + d.FermiNumber() + a1);
				TraceState state;
				state.AddTrace(f);
				AddState(state, parity, res);
			}
		}
	}
}


HamOperatorB::HamOperatorB(int creator, int annihilator)
{
	this->creator = creator;
	this->annihilator = annihilator;
}

void HamOperatorB::ApplyOnSingle(const SingleTrace& single, MixState& res)
{
	int parity;
	int n = single.BitNumber();
	SingleTrace c, d, e;
	SingleTrace a(creator / 2, 1);
	SingleTrace b(creator & 1, 1);
	int a2 = (this->annihilator & 1);
	int a1 = this->annihilator / 2;
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			if (single.Bit(i) == a1 && single.Bit(j) == a2)
			{
				single.Split(j, i, c, d, e);
				TraceState state;
				state.AddTrace(SingleTrace::Merge(a, e, c));
				state.AddTrace(SingleTrace::Merge(b, d));
				parity = (e.FermiNumber() + a1) * (b.FermiNumber() + c.FermiNumber() + d.FermiNumber());
				parity += (a2 + b.FermiNumber()) * c.FermiNumber();
				AddState(state, parity, res);
			}

			if (single.Bit(i) == a2 && single.Bit(j) == a1)
			{
				single.Split(j, i, c, d, e);
				TraceState state2;
				state2.AddTrace(SingleTrace::Merge(a, d));
				state2.AddTrace(SingleTrace::Merge(b, e, c));
				parity = (a1 + d.FermiNumber()) * (b.FermiNumber() + c.FermiNumber());
				parity += e.FermiNumber() * c.FermiNumber();
				parity += a2 * (c.FermiNumber() + d.FermiNumber() + a1);
				AddState(state2, parity, res);
			}
		}
	}

	if (a1 != b.Bit(0))
	{
		return;
	}

	for (int i = 0; i < n; i++)
	{
		if (single.Bit(i) == a2)
		{
			single.Split(i, c, d);
			TraceState state;
			state.AddTrace(SingleTrace::Merge(a, d, c));
			parity = c.FermiNumber() * (a2 + d.FermiNumber());
			AddState(state, parity, res);
		}
	}
}

void HamOperatorB::ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res)
{
	int parity;
	int a1 = annihilator / 2;
	int a2 = (annihilator & 1);
	SingleTrace a(creator / 2, 1);
	SingleTrace b(creator & 1, 1);
	SingleTrace c, d, e, f, g, h;
	for (int i = 0; i < single1.BitNumber(); i++)
	{
		if (single1.Bit(i) != a1 && single1.Bit(i) != a2) 
		{
			continue;
		}

		single1.Split(i, c, d);
		for (int j = 0; j < single2.BitNumber(); j++)
		{
			if (single1.Bit(i) == a2 && single2.Bit(j) == a1)
			{
				//..
				single2.Split(j, e, f);
				g = SingleTrace::Merge(a, f, e);
				h = SingleTrace::Merge(b, d, c);
				TraceState state;
				state.AddTrace(SingleTrace::Merge(g, h));

				parity = (a1 + f.FermiNumber()) * (b.FermiNumber() + c.FermiNumber() + d.FermiNumber() + e.FermiNumber()) +
					e.FermiNumber() * (b.FermiNumber() + c.FermiNumber() + d.FermiNumber());
				parity += c.FermiNumber() * (a2 + d.FermiNumber());
				AddState(state, parity, res);
			}

			if (single1.Bit(i) == a1 && single2.Bit(j) == a2)
			{
				//..
				single2.Split(j, e, f);
				g = SingleTrace::Merge(a, d, c);
				h = SingleTrace::Merge(b, f, e);
				TraceState state;
				state.AddTrace(SingleTrace::Merge(g, h));

				parity = d.FermiNumber() * (b.FermiNumber() + c.FermiNumber()) + b.FermiNumber() * c.FermiNumber();
				parity += e.FermiNumber() * f.FermiNumber();
				parity += a1 * (b.FermiNumber() + c.FermiNumber());
				parity += a2 * (c.FermiNumber() + d.FermiNumber() + e.FermiNumber() + a1);
				
				AddState(state, parity, res);
			}
		}
	}
}
