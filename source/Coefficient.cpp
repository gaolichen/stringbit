#include "Coefficient.h"
#include <sstream>
#include <cmath>
using namespace std;

Coefficient::Coefficient()
	:N(0), One(0), InvN(0)
{
}

Coefficient::Coefficient(int n, int one, int invN)
	:N(n), One(one), InvN(invN)
{
}

Coefficient Coefficient::operator+ (const Coefficient& a) const
{
	return Coefficient(N + a.N, One + a.One, InvN + a.InvN);
}

Coefficient& Coefficient::operator+= (const Coefficient& a)
{
	this->N += a.N;
	this->One += a.One;
	this->InvN += a.InvN;

	return *this;
}

void Coefficient::Opposite()
{
	N = -N;
	One = -One;
	InvN = -InvN;
}

void Coefficient::DecreaseOrder()
{
	if (InvN != 0)
	{
		cout << "DecreaseOrder(): Unexpected: InvN should be 0.";
	}

	InvN = One;
	One = N;
	N = 0;
}

bool Coefficient::IsZero() const
{
	return N == 0 && One == 0 && InvN == 0;
}

Coefficient& Coefficient::operator*= (int n)
{
	N *= n;
	One *= n;
	InvN *= n;

	return *this;
}

Coefficient operator* (const Coefficient& a, int n)
{
	return Coefficient(a.N * n, a.One * n, a.InvN * n);
}

ostream& operator<<(ostream& os, const Coefficient& coef)
{
	// TODO: handle InvN.
	if (coef.IsZero())
	{
		os << "0";
	}
	else if (coef.N != 0 && coef.One != 0)
	{
		os << coef.N;
		if (coef.One > 0)
		{
			os << "+";
		}

		os << coef.One << "/N";
	}
	else if (coef.N != 0)
	{
		os << coef.N;
	}
	else
	{
		os << coef.One << "/N";
	}

    return os;
}

string Coefficient::ToLatex() const
{
	// TODO: handle InvN.
	if (IsZero())
	{
		return "0";
	}
	ostringstream oss;
	if (N != 0)
	{
		oss << N;
	}

	if (One != 0)
	{
		if (One > 0)
		{
			if (N != 0)
			{
				oss << "+";
			}
		}
		else
		{
			oss << "-";
		}

		oss <<"\\frac{" << abs(One) << "}{N}";
	}

	return oss.str();
}
