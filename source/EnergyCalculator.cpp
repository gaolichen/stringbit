#pragma warning(disable:4018)
#include "EnergyCalculator.h"
#include "Partition.h"

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
			res += matM(lowestBit, i) * DoCalculate(ops ^ (1 << i));
		else res -= matM(lowestBit, i) * DoCalculate(ops ^ (1 << i));
		sign++;
	}

	cache[ops | (1 << lowestBit)] = res;
	return res;
}

CDT VevCalculator::CalculateVev(int ops, MatrixSB &omega)
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

/*
EnergyCalculator::EnergyCalculator(int s_) : s(s_)
{
}
*/


CDT EnergyCalculator::OperatorVevAllZeros(int ops, int M,  VevCalculator &calc, CDT &gmW, MatrixSB &omegaW, CDT &gmV, MatrixSB &omegaV)
{
	CDT delta;
        if (BitCount(ops) % 2 == 0)
        {
		delta = conj(OperatorVev(ops, calc, gmW, omegaW)) * OperatorVev(ops, calc, gmV, omegaV);
                int ops2 = ops + (1 << (M - 1)) + 1;
                delta += conj(OperatorVev(ops2, calc, gmW, omegaW)) * OperatorVev(ops2, calc, gmV, omegaV);
        }
        else
        {
		delta = conj(OperatorVev(ops + 1, calc, gmW, omegaW)) * OperatorVev(ops + 1, calc, gmV, omegaV);
                delta += conj(OperatorVev(ops | (1 << (M - 1)), calc, gmW, omegaW))
			* OperatorVev(ops | (1 << (M - 1)), calc, gmV, omegaV);
	}

	return delta;
}

DT EnergyCalculator::EnergyCorrection(int M, int L)
{
	int K = M - L;
	CDT ret = .0;
	DT E0 = -4 / tan(PI/(2*M)) * s;
	StringBitMatrices sbm;
	MatrixSB matM = sbm.MatrixM(M, L);
	DT detC = abs(sbm.MatrixC(M, L).determinant());
	MatrixSB omegaV = sbm.OmegaV(M, L);
	MatrixSB omegaW = sbm.OmegaW(M, L);
	CDT gmV = sbm.GammaPV(M, L);
	CDT gmW = sbm.GammaPW(M, L);
	VevCalculator calc(matM);

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
				CDT delta = OperatorVevAllZeros(ops, M, calc, gmW, omegaW, gmV, omegaV);
/*
				CDT delta;
				if (BitCount(ops) % 2 == 0)
				{
					delta = conj(OperatorVev(ops, calc, gmW, omegaW)) * OperatorVev(ops, calc, gmV, omegaV);
					int ops2 = ops + (1 << (M - 1)) + 1;
					delta += conj(OperatorVev(ops2, calc, gmW, omegaW)) * OperatorVev(ops2, calc, gmV, omegaV);
				}
				else
				{	
					delta = conj(OperatorVev(ops + 1, calc, gmW, omegaW)) * OperatorVev(ops + 1, calc, gmV, omegaV);
					delta += conj(OperatorVev(ops | (1 << (M - 1)), calc, gmW, omegaW)) 
						* OperatorVev(ops | (1 << (M - 1)), calc, gmV, omegaV);
				}*/

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
		for (int i = 0; i < modes.size(); i++)
		{
			CDT delta = 1.0;
			for (int j = 0; j < modes[i].size(); j++)
			{
				int ops = modes[i][j];
				delta *= OperatorVevAllZeros(ops, M, calc, gmW, omegaW, gmV, omegaV);
			}
		
			ret += delta * (generator.SymmetryFactor(modes[i]) /(E0 - allEnergies[i]));
		}
	}

	ret = Chop(ret);
	assert(ret.imag() == 0);
	if (ret.imag() != 0) cout << "not real energy: " << ret << endl;
	calculateTime += watch.Stop();

	return ret.real() * K * L * pow(detC, s) / M;
}

DT EnergyCalculator::EnergyCorrection(int M)
{
	totalStates = 0;
	calculateTime = 0.0;
	DT ret = 0;
	for (int L = 1; L < M - 1; L++)
		ret += EnergyCorrection(M, L);

	return ret;
}
