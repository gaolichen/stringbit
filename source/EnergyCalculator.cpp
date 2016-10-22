#pragma warning(disable:4018)
#include<fstream>
#include "EnergyCalculator.h"
#include "Partition.h"
#include "BitUtility.h"
using namespace std;

CDT VevCalculator::DoCalculate(int ops)
{
	if (ops == 0) return 1.0;

	map<int,CDT>::iterator it = cache.find(ops);
	if (it != cache.end())
	{
		return it->second;
	}

	CDT res = 0.0;
	int lowestBit = 0;
	int sign = 0;
	while (((1 << lowestBit) & ops) == 0) lowestBit++;
	ops ^= (1 << lowestBit);
	for (int i = lowestBit + 1; (1<<i) <= ops; i++)
	{
		if (((1 << i) & ops) == 0) continue;
		if (sign % 2 == 0)
			res += (*matM)(lowestBit, i) * DoCalculate(ops ^ (1 << i));
		else res -= (*matM)(lowestBit, i) * DoCalculate(ops ^ (1 << i));
		sign++;
	}

	cache[ops | (1 << lowestBit)] = res;
	return res;
}

CDT VevCalculator::CalculateVev(int ops, MatrixSB& omega)
{
	CDT ret = 0.0;
	for (int i = 0; i < M; i++)
	{
		if (((1<<i) & ops) == 0) continue;
		int cnt = 0;
		for (int j = i + 1; j < M; j++)
		{
			if (((1 << j) & ops) == 0) continue;
			if (cnt %2 == 0)
				ret += omega(i, j) * CalculateVev(ops - (1<<i) - (1<<j));
			else ret -= omega(i, j) * CalculateVev(ops - (1<<i) - (1<<j));
			cnt++;
		}
	}

	return ret + ret;
}

CDT VevCalculator::VevV(int ops)
{
	return CalculateVev(ops) * gammaV + CalculateVev(ops, *omegaV);
}

CDT VevCalculator::VevW(int ops)
{
        return CalculateVev(ops) * gammaW + CalculateVev(ops, *omegaW);
}

CDT VevCalculator::VevAll(int ops)
{
	map<int, CDT>::iterator it = cacheAll.find(ops);
	if (it != cacheAll.end()) return it->second;

	CDT delta;
        if (BitCount(ops) % 2 == 0)
        {
                delta = conj(VevW(ops)) * VevV(ops);
                int ops2 = ops + (1 << (M - 1)) + 1;
                delta += conj(VevW(ops2)) * VevV(ops2);
        }
        else
        {
                delta = conj(VevW(ops + 1)) * VevV(ops + 1);
                delta += conj(VevW(ops | (1 << (M - 1)))) * VevV(ops | (1 << (M - 1)));
        }

	// the result of each VevW or VevV need to multiple a 2/M, hence here we need 4/M^2.
	delta *= 4.0/(M * M);
	cacheAll[ops] = delta; 
       return delta;
}

DT EnergyCalculator::Operators2Energy(vector<int> &ops, int M, int s)
{
	DT E0 = -4 * s / tan(PI/2/M);
	DT ret = .0;
	for (int i = 0; i < ops.size(); i++)
	{
		ret += 8 * ops[i] * sin((i+1) * PI / M);
	}

	return ret + E0;
}


vector<StateInfo> EnergyCalculator::AllStates(int M)
{
	vector<StateInfo> ret;
	DT e0 = -4 / tan(PI/(2 * M));

	for (int i = 0; i < (1<<M); i+=2)
	{
		DT deltE = 0.0;
		int modes = 0;
		for (int j = 1; (1<<j) <= i; j++)
		{
			if ((i & (1<<j)) == 0) continue;
			deltE += 8 * sin(j * PI/M);
			modes += j;
		}

		if ((M % 2 == 1 && modes % M == 0) || (M % 2 ==0 && modes % M == M / 2))
		{
			ret.push_back(StateInfo(i, e0 + deltE));
		}
	}

	return ret;
}
/*
CDT EnergyCalculator::OperatorVev(int ops, int M,  VevCalculator &calc, CDT &gamma, MatrixSB &omega)
{
	return (calc.CalculateVev(ops) * gamma + calc.CalculateVev(ops, omega)) * (2.0/M);
}


CDT EnergyCalculator::OperatorVevAllZeros(int ops, int M,  VevCalculator &calc, CDT &gmW, MatrixSB &omegaW, CDT &gmV, MatrixSB &omegaV)
{
	CDT delta;
	if (BitCount(ops) % 2 == 0)
	{
		delta = conj(OperatorVev(ops, M, calc, gmW, omegaW)) * OperatorVev(ops, M, calc, gmV, omegaV);
		int ops2 = ops + (1 << (M - 1)) + 1;
		delta += conj(OperatorVev(ops2, M, calc, gmW, omegaW)) * OperatorVev(ops2, M, calc, gmV, omegaV);
		delta += conj(calc.VevW(ops2)) * calc.VevV(ops2);
	}
	else
	{
		delta = conj(OperatorVev(ops + 1, M, calc, gmW, omegaW)) * OperatorVev(ops + 1, M, calc, gmV, omegaV);
		delta += conj(OperatorVev(ops | (1 << (M - 1)), M, calc, gmW, omegaW))
			* OperatorVev(ops | (1 << (M - 1)), M, calc, gmV, omegaV);
	}

	return delta;
}*/

DT EnergyCalculator::EnergyCorrection(int M, int L)
{
	int K = M - L;
	CDT ret = .0;
	DT E0 = -4 / tan(PI/(2*M)) * s;
	StringBitMatrices sbm;
	MatrixSB matM = sbm.MatrixM(M, L);
	DT detC = abs(sbm.MatrixC(M, L).determinant());
	//cout << "detC=" << detC << endl;
	MatrixSB omegaV = sbm.OmegaV(M, L);
	MatrixSB omegaW = sbm.OmegaW(M, L);
	CDT gmV = sbm.GammaPV(M, L);
	CDT gmW = sbm.GammaPW(M, L);
	VevCalculator calc(matM, omegaV, omegaW, gmV, gmW);

	watch.Start();

	if (s == 1)
	{
		vector<StateInfo> states1 = AllStates(L);
		vector<StateInfo> states2 = AllStates(K);
		totalStates += states1.size() * states2.size();

		for (int i = 0; i < states2.size(); i++)
		{
			for (int j = 0; j < states1.size(); j++)
			{
				int ops = (states2[i].Ops << (L - 1)) | states1[j].Ops;
				//CDT delta = OperatorVevAllZeros(ops, M, calc, gmW, omegaW, gmV, omegaV);
				CDT delta = calc.VevAll(ops);
				delta = Chop(delta);
				// delta is not necessary real!!!!
				// assert(delta.imag() == 0);
				ret += delta /(E0 - states1[j].Energy - states2[i].Energy) ;
			}
		}
	}
	else
	{
		ModesGenerator generator(M, L, s);
		vector<vector<i64> >& modes = generator.Generate();
		vector<DT>& allEnergies = generator.AllEnergies();
		//cout << "allEnergies: " << allEnergies << endl;
		for (int i = 0; i < modes.size(); i++)
		{
			CDT delta = 1.0;
			for (int j = 0; j < modes[i].size(); j++)
			{
				int ops = modes[i][j];
				//delta *= OperatorVevAllZeros(ops, M, calc, gmW, omegaW, gmV, omegaV);
				delta *= calc.VevAll(ops);
			}

			ret += delta * (SymmetryFactor(modes[i], s) /(E0 - allEnergies[i]));
			//cout << "ret = " << ret << endl;
		}
	}

	ret = Chop(ret);
	assert(ret.imag() == 0);
	if (ret.imag() != 0) cout << "not real energy: " << ret << endl;
	normalizedE.push_back(ret.real());
	calculateTime += watch.Stop();

	return ret.real() * K * L * M * pow(detC, s);
}

DT EnergyCalculator::EnergyCorrection(int M, bool outputCorrections)
{
	totalStates = 0;
	calculateTime = 0.0;
	DT ret = 0.0;
	normalizedE.clear();
	vector<DT> delta(M - 2);
	for (int L = 1; L < M - 1; L++)
	{
		delta[L - 1] = EnergyCorrection(M, L);
		ret += delta[L - 1];
	}

	if (outputCorrections)
	{
		string file = "s="+ToString(s) + "-M=" + ToString(M) + ".txt";
		ofstream os(file.c_str());
		os << "M\tL\tdeltaE" << endl;
		for (int i = 0; i < delta.size(); i++)
		{
			os << M << "\t" << i + 1 << "\t" << delta[i] << endl;
		}

		os << "=============================" << endl;
		os << "M\tL\tdeltaE" << endl;
		for (int i = 0; i < normalizedE.size(); i++)
		{
			os << M << "\t" << i + 1 << "\t" << normalizedE[i] << endl;
		}

		os.close();
	}

	return ret;
}
