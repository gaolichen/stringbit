#pragma warning(disable:4018)
#include<iostream>
#include<complex>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <vector>
#include <map>
#include "BitUtility.h"
#include <Eigen/Dense>

using namespace std;

#define esp 1e-10

typedef double DT;
typedef std::complex<DT> CDT;
//typedef Eigen::Matrix<CDT, Dynamic, Dynamic> MatrixSB;
typedef Eigen::MatrixXcd MatrixSB;
DT PI = acos(-1.0);
CDT I(0,1);
CDT TPI = I * CDT(2 * PI);

CDT ByPolar(DT r, DT angle)
{
	return CDT(r*cos(angle) + I*r*sin(angle));
}

CDT ByPolar(DT angle)
{
	return CDT(cos(angle) + I*sin(angle));
}

CDT DotV1(int M, int L, int K, int m, int n)
{
	if((n * M - L * m) % (L * M) == 0)
	{
		return sqrt(L/(DT)M);
	}
	else
	{
		return -1/sqrt((DT)M * L) * ((DT)1 - ByPolar(-2 * PI * m * L / M)) /((DT)1 - ByPolar(-2 * PI * (n / (DT)L - m / (DT)M)));
	}
}

CDT DotV2(int M, int L, int K, int m, int n)
{
	if((n * M - K * m) % (K * M) == 0)
	{
		return sqrt(K/(DT)M) * ByPolar(-2 * PI * n * L / K);
	}
	else
	{
		return 1/sqrt((DT)M * K) * ((DT)1 - ByPolar(-2 * PI * m * L / M) ) /((DT)1 - ByPolar(-2 * PI * (n / (DT)K - m / (DT)M)));
	}
}

CDT DotW(int M, int L, int K, int m)
{
	if (m == 0) return 0;

	return -1/sqrt((DT)L * K) * ((DT)1.0 - ByPolar(- 2 * PI * m * L / M)) / ((DT)1.0 - ByPolar(2 * PI * m / M));
}

CDT Cm0(int M, int L, int K, int m)
{
	if (m == 0) return 1;
	return 0;
}

CDT Cmn1(int M, int L, int K, int m, int n)
{
	return DotV1(M, L, K, m, n) * cos(n * PI / (2 * L) - PI * m / (2 * M));
}

CDT Cmn2(int M, int L, int K, int m, int n)
{
	return DotV2(M, L, K, m, n) * cos(n * PI / (2 * K) - PI * m / (2 * M));
}

CDT CmM(int M, int L, int K, int m)
{
	return DotW(M, L, K, m) * cos(m * PI / (2 * M) - PI / 4) * ByPolar(-PI / 4);
}

CDT Smn1(int M, int L, int K, int m, int n)
{
	return DotV1(M, L, K, m, L - n) * cos(n * PI / (2 * L) + PI * m / (2 * M));
}

CDT Smn2(int M, int L, int K, int m, int n)
{
	return DotV2(M, L, K, m, K - n) * cos(n * PI / (2 * K) + PI * m / (2 * M));
}

CDT SmM(int M, int L, int K, int m)
{
	return DotW(M, L, K, m) * cos(m * PI / (2 * M) + PI / 4) * ByPolar(PI / 4);
}

CDT ElementC(int M, int L, int K, int m, int n)
{
	if (m == 0) return Cm0(M, L, K, n);
	if (n == 0) return 0;
	if (n < L) return Cmn1(M, L, K, m, n);
	if (n < M - 1) return Cmn2(M, L, K, m, n - L + 1);
	assert(n == M - 1);
	return CmM(M, L, K, m);
}

CDT ElementS(int M, int L, int K, int m, int n)
{
	if (m == 0) return 0;
	if (n == 0) return 0;
	if (n < L) return Smn1(M, L, K, m, n);
	if (n < M - 1) return Smn2(M, L, K, m, n - L + 1);
	assert(n == M - 1);
	return SmM(M, L, K, m);
}

MatrixSB MatrixC(int M, int L)
{
	MatrixSB ret(M, M);
	for (int i = 0; i < M; i++)
		for (int j = 0; j < M; j++)
			ret(i, j) = ElementC(M, L, M - L, i, j);

	return ret;
}

MatrixSB MatrixS(int M, int L)
{
	MatrixSB ret(M, M);
	for (int i = 0; i < M; i++)
		for (int j = 0; j < M; j++)
			ret(i, j) = ElementS(M, L, M - L, i, j);

	return ret;
}

MatrixSB MatrixM(int M, int L)
{
	return -MatrixC(M, L).inverse() * MatrixS(M, L);
}

CDT ElementA(int k, int l, int M, int n, int m)
{
	return ByPolar(2 * PI * k * (m + n) / M) * sin((m - n) * PI / (2 * M)) 
		+ ByPolar(2 * PI * (k + l) * (m + n) / (2 * M)) * sin((l - k - .5) * (m - n) * PI / M);
}

MatrixSB MatrixA(int k, int l, int M)
{
	MatrixSB ret(M, M);
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < M; j++)
		{
			ret(i, j) = ElementA(k, l, M, i, j);
		}
	}

	return ret;
}

MatrixSB MatrixAV(int M, int L)
{
	MatrixSB ret(M, M);
	for (int n = 0; n < M; n++)
	{
		for (int m = 0; m < M; m++)
		{
			ret(n, m) = sin((m - n) * PI / (2 * M)) 
				+ ByPolar(2 * PI * (L + 1) * (m + n) / (2 * M)) * sin((L + 0.5) * (m - n) * PI / M);
		}
	}

	return ret;
}

MatrixSB MatrixAW(int M, int L)
{
	MatrixSB ret(M, M);
	for (int n = 0; n < M; n++)
	{
		for (int m = 0; m < M; m++)
		{
			ret(n, m) = ((DT)1 + ByPolar(PI * (m + n) / M)) * ((DT)1 
				+ ByPolar(2 * PI * L * (m + n) / M)) * sin((m - n) * PI / (2 * M));
		}
	}

	return ret;
}

MatrixSB OmegaV(int M, int L)
{
	MatrixSB invc = MatrixC(M, L).inverse();
	return invc * MatrixAV(M, L).adjoint() * invc.transpose();
}

MatrixSB OmegaW(int M, int L)
{
	MatrixSB invc = MatrixC(M, L).inverse();
	return invc * MatrixAW(M, L).adjoint() * invc.transpose();
}

CDT GammaPV(int M, int L)
{
	CDT tr = (MatrixS(M, L).conjugate() * MatrixC(M, L).inverse() * MatrixAV(M, L).adjoint()).trace();
	return -1/tan(PI / (2 * M)) - 1 / tan((2 * L + 1) * PI / (2 * M)) - tr;
}

CDT GammaPW(int M, int L)
{
	CDT tr = (MatrixS(M, L).conjugate() * MatrixC(M, L).inverse() * MatrixAW(M, L).adjoint()).trace();
	return -4/tan(PI / (2 * M)) - tr;
}

DT Chop(DT a)
{
	if (abs(a) < esp) return 0;
	return a;
}

CDT Chop(CDT& a)
{
	return CDT(Chop(a.real()), Chop(a.imag()));
}

void ChopInplace(MatrixSB &m)
{
	for (int i = 0; i < m.innerSize(); i++)
	{
		for (int j = 0; j < m.innerSize(); j++)
		{
			m(i, j) = Chop(m(i, j));
		}
	}
}

MatrixSB Chop(MatrixSB m)
{
	MatrixSB ret = m;
	for (int i = 0; i < ret.innerSize(); i++)
	{
		for (int j = 0; j < ret.innerSize(); j++)
		{
			ret(i, j) = Chop(ret(i, j));
		}
	}

	return ret;
}

void TestMatrixCS()
{
	for (int M = 3; M <= 10; M++)
	{
		MatrixSB realZero = MatrixSB::Zero(M, M);
		MatrixSB id = MatrixSB::Identity(M, M);
		for (int L = 1; L < M - 1; L++)
		{
			MatrixSB c = MatrixC(M, L);
			MatrixSB s = MatrixS(M, L);
			MatrixSB zero = c * s.transpose() + s * c.transpose();
			MatrixSB one = c * c.adjoint() + s * s.adjoint();
			ChopInplace(zero);

			if (zero != realZero)
			{
				cout << "error! M=" << M << ", L=" << L << ". Expect zero but returned " << endl << zero << endl;
				return;
			}
			if (Chop(one - id) != realZero)
			{
				cout << "error! M=" << M << ", L=" << L << ". Expect identity but returned " << endl << one << endl;
				return;
			}
		}
	}

	cout << "TestMatrixCS passed!" << endl;
}

void TestMatrixA()
{
	for (int M = 3; M < 10; M++)
	{
		MatrixSB zero = MatrixSB::Zero(M, M);
		for (int L = 1; L < M; L++)
		{
			MatrixSB z = MatrixA(M, L + 1, M) - MatrixAV(M, L);
			ChopInplace(z);

			if (z != zero)
			{
				cout << "Test MatrixAV error! M=" << M << ", L=" << L << endl;
				return;
			}

			z = MatrixA(L, L + 1, M) + MatrixA(M, 1, M) - MatrixAW(M, L);
			ChopInplace(z);

			if (z != zero)
			{
				cout << "Test MatrixAW error! M=" << M << ", L=" << L << endl;
				return;
			}
		}
	}

	cout << "TestMatrixA passed!" << endl;
}

class StateInfo
{
public:
	int Ops;
	DT Energy;
	StateInfo(int ops, DT energy):Ops(ops), Energy(energy)
	{
	}

};

vector<StateInfo> AllStates(int M)
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

void TestAllStates()
{
	vector<StateInfo> states = AllStates(5);
	for (int i = 0; i < states.size(); i++)
	{
		cout << states[i].Ops << " " << states[i].Energy << endl;
	}
}

class VevCalculator
{
private:
	MatrixSB matM;
	int M;
	map<int, CDT> cache;

	CDT DoCalculate(int ops)
	{
		if (ops == 0) return 1.0;

		map<int,CDT>::iterator it = cache.find(ops);
		if (it != cache.end())
		{
			return it->second;
		}

		CDT res = 0.0;
		int lowestBit = 0;
		int sign = 1;
		while (((1 << lowestBit) & ops) == 0) lowestBit++;
		ops ^= (1 << lowestBit);
		for (int i = lowestBit + 1; (1<<i) <= ops; i++)
		{
			if (((1 << i) & ops) == 0) continue;
			res += (DT)sign * matM(lowestBit, i) * DoCalculate(ops ^ (1<<i));
			sign *= -1;
		}

		cache[ops | (1 << lowestBit)] = res;
		return res;
	}
public:
	VevCalculator(MatrixSB mat) : matM(mat)
	{
		M = matM.innerSize();
	}

	CDT CalculateVev(int ops)
	{
		if (BitCount(ops) % 2 == 1) return 0.0;
		return DoCalculate(ops);
	}

	CDT CalculateVev(int ops, MatrixSB &omega)
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

		return 2.0 * ret;
	}
};

CDT OperatorVev(int ops, VevCalculator &calc, CDT gamma, MatrixSB &omega)
{
	return 2.0 * (calc.CalculateVev(ops) * gamma + calc.CalculateVev(ops, omega));
}

DT EnergyCorrection(int M, int L)
{
	int K = M - L;
	CDT ret = .0;
	DT E0 = -4 / tan(PI/(2*M));
	vector<StateInfo> states1 = AllStates(L);
	vector<StateInfo> states2 = AllStates(K);
	MatrixSB matM = MatrixM(M, L);
	DT detC = abs(MatrixC(M, L).determinant());
	MatrixSB omegaV = OmegaV(M, L);
	MatrixSB omegaW = OmegaW(M, L);
	CDT gmV = GammaPV(M, L);
	CDT gmW = GammaPW(M, L);
	VevCalculator calc(matM);
	
	for (int i = 0; i < states2.size(); i++)
	{
		for (int j = 0; j < states1.size(); j++)
		{
			int ops = (states2[i].Ops << (L - 1)) | states1[j].Ops;
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

			delta = Chop(delta);
			// delta is not necessary real!!!!
			//assert(delta.imag() == 0);
			ret += delta /(E0 - states1[j].Energy - states2[i].Energy) ;
		}
	}

	ret = Chop(ret);
	assert(ret.imag() == 0);
	if (ret.imag() != 0) cout << "not real energy: " << ret << endl;

	return ret.real() * K * L * detC / M;
}

DT EnergyCorrection(int M)
{
	DT ret = 0;
	for (int L = 1; L < M - 1; L++)
		ret += EnergyCorrection(M, L);

	return ret;
}

void TestVevCalculator()
{
	MatrixSB matM = MatrixM(3, 1);
	VevCalculator calc(matM);
	CDT res = calc.CalculateVev(3);
	cout << res << endl;
}

void TestEnergyCorrection()
{
	for (int i = 3; i <= 15; i += 2)
	{
		cout << "M= " << i << ": " << EnergyCorrection(i) << endl;
	}
}

void TestMatrices(int M, int L)
{
	MatrixSB omegaV = OmegaV(M, L);
	MatrixSB omegaW = OmegaW(M, L);
	MatrixSB matM = MatrixM(M, L);
	ChopInplace(omegaV);
	ChopInplace(omegaW);
	ChopInplace(matM);
	CDT gammaV = GammaPV(M, L);
	CDT gammaW = GammaPW(M, L);
	gammaV = Chop(gammaV);
	gammaW = Chop(gammaW);

	cout << "M=" << M << ", L=" << L << endl;
	cout << "Matrix M:" << endl << matM << endl;
	cout << "OmegaV:" << endl << omegaV << endl;
	cout << "OmegaW: " << endl << omegaW << endl;
	cout << "gammaV=" << gammaV << endl;
	cout << "gammaW=" << gammaW << endl;
}

void TestOperatorVev()
{
	int M = 5;
	int L = 2;

	MatrixSB omegaV = OmegaV(M, L);
	MatrixSB omegaW = OmegaW(M, L);
	MatrixSB matM = MatrixM(M, L);
	ChopInplace(omegaV);
	ChopInplace(omegaW);
	ChopInplace(matM);
	CDT gammaV = GammaPV(M, L);
	CDT gammaW = GammaPW(M, L);
	gammaV = Chop(gammaV);
	gammaW = Chop(gammaW);

	VevCalculator calc(matM);

	// test VevCalculator 
	/*for (int i = 0; i < (1<<M); i++)
	{
		CDT res = calc.CalculateVev(i, omegaV);
		res = Chop(res);
		if (res != .0)
			cout << i << ": " << res << endl;
	}*/

	for (int i = 0; i < (1<<M); i++)
	{
		CDT res = OperatorVev(i, calc, gammaW, omegaW);
		res = Chop(res);
		if (res != .0)
			cout << i << ": " << res << endl;
	}
}

int main()
{
	TestMatrixCS();
	TestMatrixA();
	//TestAllStates();
	//TestVevCalculator();
	TestEnergyCorrection();
	//TestMatrices(4,1);
	TestOperatorVev();
	return 0;
}
