#include "Coefficient.h"
#include <sstream>
#include <cmath>
using namespace std;

Coefficient::Coefficient()
{
	Coefficient(0, 0);
}

Coefficient::Coefficient(int n, int one)
	:N(n), One(one)
{
}

Coefficient Coefficient::operator+ (const Coefficient& a) const
{
	return Coefficient(N + a.N, One + a.One);
}

Coefficient& Coefficient::operator+= (const Coefficient& a)
{
	this->N += a.N;
	this->One += a.One;

	return *this;
}

void Coefficient::Opposite()
{
	N = -N;
	One = -One;
}

bool Coefficient::IsZero() const
{
	return N == 0 && One == 0;
}

Coefficient& Coefficient::operator*= (int n)
{
	N *= n;
	One *= n;

	return *this;
}

Coefficient operator* (const Coefficient& a, int n)
{
	return Coefficient(a.N * n, a.One * n);
}

ostream& operator<<(ostream& os, const Coefficient& coef)
{
	if (coef.N == 0 && coef.One == 0)
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
	
	if (N == 0 && One == 0)
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
