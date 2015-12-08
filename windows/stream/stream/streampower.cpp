//#include<malloc.h>
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include "streampower.h"

/* old implementation */
float* StreamPower::vector(int nl, int nh)
/*	allocate a float vector with subscript range v[nl..nh] 
	NOTE: changed from malloc to calloc to allow testing of value initialisation with std::vector
*/
{
	float *v;

	v = (float *)calloc((unsigned int)(nh - nl + 1 + NR_END), sizeof(float));
	return v - nl + NR_END;
}

/* allocate a float vector with subscript range v[nl..nh] */
std::vector<float> StreamPower::Vector(int nl, int nh)
{
	int size = nh - nl + 1 + NR_END;
	return std::vector<float>(size);
}

int* StreamPower::ivector(int nl, int nh)
/*	allocate a float vector with subscript range v[nl..nh]
NOTE: changed from malloc to calloc to allow testing of value initialisation with std::vector
*/
{
	int *v;

	v = (int *)calloc((unsigned int)(nh - nl + 1 + NR_END), sizeof(int));
	return v - nl + NR_END;
}

/* allocate a float vector with subscript range v[nl..nh] */
std::vector<int> StreamPower::IVector(int nl, int nh)
{
	int size = nh - nl + 1 + NR_END;
	return std::vector<int>(size);
}

/*
	We don't need to manually free vectors based on std::vector
	so we don't need to duplicate this method
*/
void StreamPower::free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

/*
	We don't need to manually free vectors based on std::vector
	so we don't need to duplicate this method
*/
void StreamPower::free_vector(float *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

int** StreamPower::imatrix(long nrl, long nrh, long ncl, long nch)
/*	allocate a int matrix with subscript range m[nrl..nrh][ncl..nch]
	NOTE: changed from malloc to calloc to allow testing of value initialisation with std::vector
*/
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	int **m;

	m = (int **)calloc((size_t)(nrow + NR_END), sizeof(int*));
	m += NR_END;
	m -= nrl;

	m[nrl] = (int *)calloc((size_t)(nrow*ncol + NR_END), sizeof(int));
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

	return m;
}

std::vector<std::vector<int>> StreamPower::IMatrix(int nrl, int nrh, int ncl, int nch)
{
	int rsize = nrh - nrl + 1 + NR_END;
	int csize = nch - ncl + 1 + NR_END;
	return std::vector<std::vector<int>>(rsize, std::vector<int>(csize));
}

float** StreamPower::matrix(long nrl, long nrh, long ncl, long nch)
/*	allocate a float matrix with subscript range m[nrl..nrh][ncl..nch]
	NOTE: changed from malloc to calloc to allow testing of value initialisation with std::vector
*/
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	float **m;

	m = (float **)calloc((size_t)(nrow + NR_END), sizeof(float*));
	m += NR_END;
	m -= nrl;

	m[nrl] = (float *)calloc((size_t)(nrow*ncol + NR_END), sizeof(float));
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

	return m;
}

std::vector<std::vector<float>> StreamPower::Matrix(int nrl, int nrh, int ncl, int nch)
{
	int rsize = nrh - nrl + 1 + NR_END;
	int csize = nch - ncl + 1 + NR_END;
	return std::vector<std::vector<float>>(rsize, std::vector<float>(csize));
}

float StreamPower::ran3(int* idum)
{
/*
	Assuming this is supposed to produce a uniform random number between 0 and 1
	as in http://cics.umd.edu/~lee/z/h1/w/src_jpl/ran3.F
*/
	static int inext, inextp;
	static long ma[56];
	static int iff = 0;
	long mj, mk;
	int i, ii, k;

	if (*idum < 0 || iff == 0) {
		iff = 1;
		mj = MSEED - (*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for (i = 1; i <= 54; i++) {
			ii = (21 * i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if (mk < MZ) mk += MBIG;
			mj = ma[ii];
		}
		for (k = 1; k <= 4; k++)
			for (i = 1; i <= 55; i++) {
				ma[i] -= ma[1 + (i + 30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext = 0;
		inextp = 31;
		*idum = 1;
	}
	if (++inext == 56) inext = 1;
	if (++inextp == 56) inextp = 1;
	mj = ma[inext] - ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext] = mj;
	return mj*FAC;
}

float StreamPower::Ran3(std::default_random_engine& generator, std::uniform_real_distribution<float>& distribution)
{
	return distribution(generator);
}

float StreamPower::gasdev(int* idum)
{
	static int iset = 0;
	static float gset;
	float fac, r, v1, v2;
	//float ran3();

	if (iset == 0)
	{
		do
		{
			v1 = 2.0*StreamPower::ran3(idum) - 1.0;
			v2 = 2.0*StreamPower::ran3(idum) - 1.0;
			r = v1*v1 + v2*v2;
		} while (r >= 1.0);
		fac = sqrt(-2.0*log(r) / r);
		gset = v1*fac;
		iset = 1;
		return v2*fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}

float StreamPower::Gasdev(std::default_random_engine& generator, std::normal_distribution<float>& distribution)
{
/*
	Assuming this is the same code from here: http://www.stat.berkeley.edu/~paciorek/diss/code/regression.binomial/gasdev.C
	We need to return a standard, normally distributed gaussian random number
*/

	return distribution(generator);

}

void StreamPower::indexx(int n, float* arr, int* indx)
{
	unsigned long i, indxt, ir = n, itemp, j, k, l = 1;
	int jstack = 0, *istack;
	float a;

	istack = ivector(1, NSTACK);
	for (j = 1; j <= n; j++) indx[j] = j;
	for (;;)
	{
		if (ir - l < M)
		{
			for (j = l + 1; j <= ir; j++)
			{
				indxt = indx[j];
				a = arr[indxt];
				for (i = j - 1; i >= 1; i--)
				{
					if (arr[indx[i]] <= a) break;
					indx[i + 1] = indx[i];
				}
				indx[i + 1] = indxt;
			}
			if (jstack == 0) break;
			ir = istack[jstack--];
			l = istack[jstack--];
		}
		else
		{
			k = (l + ir) >> 1;
			SWAP(indx[k], indx[l + 1]);
			if (arr[indx[l + 1]] > arr[indx[ir]])
			{
				SWAP(indx[l + 1], indx[ir]);
			}
			if (arr[indx[l]] > arr[indx[ir]])
			{
				SWAP(indx[l], indx[ir]);
			}
			if (arr[indx[l + 1]] > arr[indx[l]])
			{
				SWAP(indx[l + 1], indx[l]);
			}
			i = l + 1;
			j = ir;
			indxt = indx[l];
			a = arr[indxt];
			for (;;)
			{
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i], indx[j])
			}
			indx[l] = indx[j];
			indx[j] = indxt;
			jstack += 2;
			if (ir - i + 1 >= j - l)
			{
				istack[jstack] = ir;
				istack[jstack - 1] = i;
				ir = j - 1;
			}
			else
			{
				istack[jstack] = j - 1;
				istack[jstack - 1] = l;
				l = i;
			}
		}
	}
	free_ivector(istack, 1, NSTACK);
}

void StreamPower::Indexx(int n, float* arr, int* indx)
{
	std::vector<float> arrv = ArrayToVector(arr, n);
	std::vector<int> indxv = Indexx(arrv);
	VectorToArray(indxv, indx);
}

std::vector<int> StreamPower::Indexx(std::vector<float>& arr)
{
	return SortFortranIndices(arr);
}

void StreamPower::tridag(float a[], float b[], float c[], float r[], float u[], unsigned long n)
{
	unsigned long j;
	float bet, *gam;

	gam = vector(1, n);
	u[1] = r[1] / (bet = b[1]);
	for (j = 2; j <= n; j++)
	{
		gam[j] = c[j - 1] / bet;
		bet = b[j] - a[j] * gam[j];
		u[j] = (r[j] - a[j] * u[j - 1]) / bet;
	}
	for (j = (n - 1); j >= 1; j--)
		u[j] -= gam[j + 1] * u[j + 1];
	free_vector(gam, 1, n);
}

void StreamPower::Tridag(float a[], float b[], float c[], float r[], float u[], unsigned long n)
{
	std::vector<float> av = ArrayToVector(a, n, true);
	std::vector<float> bv = ArrayToVector(b, n, true);
	std::vector<float> cv = ArrayToVector(c, n, true);
	std::vector<float> rv = ArrayToVector(r, n, true);
	std::vector<float> uv = ArrayToVector(u, n, true);

	Tridag(av, bv, cv, rv, uv, n);

	VectorToArray(av, a);
	VectorToArray(bv, b);
	VectorToArray(cv, c);
	VectorToArray(rv, r);
	VectorToArray(uv, u);
}


void StreamPower::Tridag(std::vector<float>& a, std::vector<float>& b, std::vector<float>& c, std::vector<float>& r, std::vector<float>& u, int n)
{
	unsigned long j;
	float bet;
	std::vector<float> gam;

	gam = Vector(1, n);
	u[1] = r[1] / (bet = b[1]);
	for (j = 2; j <= n; j++)
	{
		gam[j] = c[j - 1] / bet;
		bet = b[j] - a[j] * gam[j];
		u[j] = (r[j] - a[j] * u[j - 1]) / bet;
	}
	for (j = (n - 1); j >= 1; j--)
	{
		u[j] -= gam[j + 1] * u[j + 1];
	}
}

void StreamPower::setupgridneighbors()
{
	int i, j;

	_idown = ivector(1, lattice_size_x);
	_iup = ivector(1, lattice_size_x);
	_jup = ivector(1, lattice_size_y);
	_jdown = ivector(1, lattice_size_y);

	for (i = 1; i <= lattice_size_x; i++)
	{
		_idown[i] = i - 1;
		_iup[i] = i + 1;
	}
	_idown[1] = 1;
	_iup[lattice_size_x] = lattice_size_x;
	for (j = 1; j <= lattice_size_y; j++)
	{
		_jdown[j] = j - 1;
		_jup[j] = j + 1;
	}
	_jdown[1] = 1;
	_jup[lattice_size_y] = lattice_size_y;
}

void StreamPower::SetupGridNeighbors()
{
	int i, j;

	idown = IVector(1, lattice_size_x);
	iup = IVector(1, lattice_size_x);
	jup = IVector(1, lattice_size_y);
	jdown = IVector(1, lattice_size_y);

	for (i = 1; i <= lattice_size_x; i++)
	{
		idown[i] = i - 1;
		iup[i] = i + 1;
	}
	idown[1] = 1;
	iup[lattice_size_x] = lattice_size_x;
	for (j = 1; j <= lattice_size_y; j++)
	{
		jdown[j] = j - 1;
		jup[j] = j + 1;
	}
	jdown[1] = 1;
	jup[lattice_size_y] = lattice_size_y;
}

// Replace with TauDEM flood
/*
void fillinpitsandflats(int i, int j)
{
	float min;

	min = _topo[i][j];
	if (_topo[_iup[i]][j] < min) min = _topo[_iup[i]][j];
	if (_topo[_idown[i]][j] < min) min = _topo[_idown[i]][j];
	if (_topo[i][_jup[j]] < min) min = _topo[i][_jup[j]];
	if (_topo[i][_jdown[j]] < min) min = _topo[i][_jdown[j]];
	if (_topo[_iup[i]][_jup[j]] < min) min = _topo[_iup[i]][_jup[j]];
	if (_topo[_idown[i]][_jup[j]] < min) min = _topo[_idown[i]][_jup[j]];
	if (_topo[_idown[i]][_jdown[j]] < min) min = _topo[_idown[i]][_jdown[j]];
	if (_topo[_iup[i]][_jdown[j]] < min) min = _topo[_iup[i]][_jdown[j]];
	if ((_topo[i][j] <= min) && (i>1) && (j>1) && (i < lattice_size_x) && (j < lattice_size_y))
	{
		_topo[i][j] = min + fillincrement;
		fillinpitsandflats(i, j);
		fillinpitsandflats(_iup[i], j);
		fillinpitsandflats(_idown[i], j);
		fillinpitsandflats(i, _jup[j]);
		fillinpitsandflats(i, _jdown[j]);
		fillinpitsandflats(_iup[i], _jup[j]);
		fillinpitsandflats(_idown[i], _jup[j]);
		fillinpitsandflats(_idown[i], _jdown[j]);
		fillinpitsandflats(_iup[i], _jdown[j]);
	}
}

/* this is ok
void mfd_flowroute(int i, int j)

{
	float tot;

	tot = 0;
	if (_topo[i][j] > _topo[_iup[i]][j])
		tot += pow(_topo[i][j] - _topo[_iup[i]][j], 1.1f);
	if (_topo[i][j] > _topo[_idown[i]][j])
		tot += pow(_topo[i][j] - _topo[_idown[i]][j], 1.1f);
	if (_topo[i][j] > _topo[i][_jup[j]])
		tot += pow(_topo[i][j] - _topo[i][_jup[j]], 1.1f);
	if (_topo[i][j] > _topo[i][_jdown[j]])
		tot += pow(_topo[i][j] - _topo[i][_jdown[j]], 1.1f);
	if (_topo[i][j] > _topo[_iup[i]][_jup[j]])
		tot += pow((_topo[i][j] - _topo[_iup[i]][_jup[j]])*oneoversqrt2, 1.1f);
	if (_topo[i][j] > _topo[_iup[i]][_jdown[j]])
		tot += pow((_topo[i][j] - _topo[_iup[i]][_jdown[j]])*oneoversqrt2, 1.1f);
	if (_topo[i][j] > _topo[_idown[i]][_jup[j]])
		tot += pow((_topo[i][j] - _topo[_idown[i]][_jup[j]])*oneoversqrt2, 1.1f);
	if (_topo[i][j] > _topo[_idown[i]][_jdown[j]])
		tot += pow((_topo[i][j] - _topo[_idown[i]][_jdown[j]])*oneoversqrt2, 1.1f);
	if (_topo[i][j] > _topo[_iup[i]][j])
		_flow1[i][j] = pow(_topo[i][j] - _topo[_iup[i]][j], 1.1f) / tot;
	else _flow1[i][j] = 0;
	if (_topo[i][j] > _topo[_idown[i]][j])
		_flow2[i][j] = pow(_topo[i][j] - _topo[_idown[i]][j], 1.1f) / tot;
	else _flow2[i][j] = 0;
	if (_topo[i][j] > _topo[i][_jup[j]])
		_flow3[i][j] = pow(_topo[i][j] - _topo[i][_jup[j]], 1.1f) / tot;
	else _flow3[i][j] = 0;
	if (_topo[i][j] > _topo[i][_jdown[j]])
		_flow4[i][j] = pow(_topo[i][j] - _topo[i][_jdown[j]], 1.1f) / tot;
	else _flow4[i][j] = 0;
	if (_topo[i][j] > _topo[_iup[i]][_jup[j]])
		_flow5[i][j] = pow((_topo[i][j] - _topo[_iup[i]][_jup[j]])*oneoversqrt2, 1.1f) / tot;
	else _flow5[i][j] = 0;
	if (_topo[i][j] > _topo[_iup[i]][_jdown[j]])
		_flow6[i][j] = pow((_topo[i][j] - _topo[_iup[i]][_jdown[j]])*oneoversqrt2, 1.1f) / tot;
	else _flow6[i][j] = 0;
	if (_topo[i][j] > _topo[_idown[i]][_jup[j]])
		_flow7[i][j] = pow((_topo[i][j] - _topo[_idown[i]][_jup[j]])*oneoversqrt2, 1.1f) / tot;
	else _flow7[i][j] = 0;
	if (_topo[i][j] > _topo[_idown[i]][_jdown[j]])
		_flow8[i][j] = pow((_topo[i][j] - _topo[_idown[i]][_jdown[j]])*oneoversqrt2, 1.1f) / tot;
	else _flow8[i][j] = 0;
	_flow[_iup[i]][j] += _flow[i][j] * _flow1[i][j];
	_flow[_idown[i]][j] += _flow[i][j] * _flow2[i][j];
	_flow[i][_jup[j]] += _flow[i][j] * _flow3[i][j];
	_flow[i][_jdown[j]] += _flow[i][j] * _flow4[i][j];
	_flow[_iup[i]][_jup[j]] += _flow[i][j] * _flow5[i][j];
	_flow[_iup[i]][_jdown[j]] += _flow[i][j] * _flow6[i][j];
	_flow[_idown[i]][_jup[j]] += _flow[i][j] * _flow7[i][j];
	_flow[_idown[i]][_jdown[j]] += _flow[i][j] * _flow8[i][j];
}
*/

/* ok
void calculatealongchannel_slope(int i, int j)
{
	float down;

	down = 0;
	if (_topo[_iup[i]][j] - _topo[i][j] < down) down = _topo[_iup[i]][j] - _topo[i][j];
	if (_topo[_idown[i]][j] - _topo[i][j] < down) down = _topo[_idown[i]][j] - _topo[i][j];
	if (_topo[i][_jup[j]] - _topo[i][j] < down) down = _topo[i][_jup[j]] - _topo[i][j];
	if (_topo[i][_jdown[j]] - _topo[i][j] < down) down = _topo[i][_jdown[j]] - _topo[i][j];
	if ((_topo[_iup[i]][_jup[j]] - _topo[i][j])*oneoversqrt2 < down)
		down = (_topo[_iup[i]][_jup[j]] - _topo[i][j])*oneoversqrt2;
	if ((_topo[_idown[i]][_jup[j]] - _topo[i][j])*oneoversqrt2 < down)
		down = (_topo[_idown[i]][_jup[j]] - _topo[i][j])*oneoversqrt2;
	if ((_topo[_iup[i]][_jdown[j]] - _topo[i][j])*oneoversqrt2 < down)
		down = (_topo[_iup[i]][_jdown[j]] - _topo[i][j])*oneoversqrt2;
	if ((_topo[_idown[i]][_jdown[j]] - _topo[i][j])*oneoversqrt2 < down)
		down = (_topo[_idown[i]][_jdown[j]] - _topo[i][j])*oneoversqrt2;
	_slope[i][j] = fabs(down) / delt_ax;
}
*/


void StreamPower::hill_slopediffusioninit()
{
	int i, j, count;
	float term1;

	_ax = vector(1, lattice_size_x);
	_ay = vector(1, lattice_size_y);
	_bx = vector(1, lattice_size_x);
	_by = vector(1, lattice_size_y);
	_cx = vector(1, lattice_size_x);
	_cy = vector(1, lattice_size_y);
	_ux = vector(1, lattice_size_x);
	_uy = vector(1, lattice_size_y);
	_rx = vector(1, lattice_size_x);
	_ry = vector(1, lattice_size_y);
	D = 10000000.0;
	count = 0;
	term1 = D / (deltax*deltax);
	for (i = 1; i <= lattice_size_x; i++)
	{
		_ax[i] = -term1;
		_cx[i] = -term1;
		_bx[i] = 4 * term1 + 1;
		if (i == 1)
		{
			_bx[i] = 1;
			_cx[i] = 0;
		}
		if (i == lattice_size_x)
		{
			_bx[i] = 1;
			_ax[i] = 0;
		}
	}
	for (j = 1; j <= lattice_size_y; j++)
	{
		_ay[j] = -term1;
		_cy[j] = -term1;
		_by[j] = 4 * term1 + 1;
		if (j == 1)
		{
			_by[j] = 1;
			_cy[j] = 0;
		}
		if (j == lattice_size_y)
		{
			_by[j] = 1;
			_ay[j] = 0;
		}
	}
	while (count < 5)
	{
		count++;
		for (i = 1; i <= lattice_size_x; i++)
			for (j = 1; j <= lattice_size_y; j++)
				_topo2[i][j] = _topo[i][j];
		for (i = 1; i <= lattice_size_x; i++)
		{
			for (j = 1; j <= lattice_size_y; j++)
			{
				_ry[j] = term1*(_topo[_iup[i]][j] + _topo[_idown[i]][j]) + _topoold[i][j];
				if (j == 1)
					_ry[j] = _topoold[i][j];
				if (j == lattice_size_y)
					_ry[j] = _topoold[i][j];
			}
			tridag(_ay, _by, _cy, _ry, _uy, lattice_size_y);
			for (j = 1; j <= lattice_size_y; j++)
				_topo[i][j] = _uy[j];
		}
		for (i = 1; i <= lattice_size_x; i++)
			for (j = 1; j <= lattice_size_y; j++)
				_topo2[i][j] = _topo[i][j];
		for (j = 1; j <= lattice_size_y; j++)
		{
			for (i = 1; i <= lattice_size_x; i++)
			{
				_rx[i] = term1*(_topo[i][_jup[j]] + _topo[i][_jdown[j]]) + _topoold[i][j];
				if (i == 1)
					_rx[i] = _topoold[i][j];
				if (i == lattice_size_x)
					_rx[i] = _topoold[i][j];
			}
			tridag(_ax, _bx, _cx, _rx, _ux, lattice_size_x);
			for (i = 1; i <= lattice_size_x; i++)
				_topo[i][j] = _ux[i];
		}
	}
}

void StreamPower::HillSlopeDiffusionInit()
{
	int i, j, count;
	float term1;

	ax = Vector(1, lattice_size_x);
	ay = Vector(1, lattice_size_y);
	bx = Vector(1, lattice_size_x);
	by = Vector(1, lattice_size_y);
	cx = Vector(1, lattice_size_x);
	cy = Vector(1, lattice_size_y);
	ux = Vector(1, lattice_size_x);
	uy = Vector(1, lattice_size_y);
	rx = Vector(1, lattice_size_x);
	ry = Vector(1, lattice_size_y);
	D = 10000000.0;
	count = 0;
	term1 = D / (deltax*deltax);
	for (i = 1; i <= lattice_size_x; i++)
	{
		ax[i] = -term1;
		cx[i] = -term1;
		bx[i] = 4 * term1 + 1;
		if (i == 1)
		{
			bx[i] = 1;
			cx[i] = 0;
		}
		if (i == lattice_size_x)
		{
			bx[i] = 1;
			ax[i] = 0;
		}
	}
	for (j = 1; j <= lattice_size_y; j++)
	{
		ay[j] = -term1;
		cy[j] = -term1;
		by[j] = 4 * term1 + 1;
		if (j == 1)
		{
			by[j] = 1;
			cy[j] = 0;
		}
		if (j == lattice_size_y)
		{
			by[j] = 1;
			ay[j] = 0;
		}
	}
	while (count < 5)
	{
		count++;
		for (i = 1; i <= lattice_size_x; i++)
		{
			for (j = 1; j <= lattice_size_y; j++)
			{
				topo2[i][j] = topo[i][j];
			}
		}
		for (i = 1; i <= lattice_size_x; i++)
		{
			for (j = 1; j <= lattice_size_y; j++)
			{
				ry[j] = term1*(topo[iup[i]][j] + topo[idown[i]][j]) + topoold[i][j];
				if (j == 1)
					ry[j] = topoold[i][j];
				if (j == lattice_size_y)
					ry[j] = topoold[i][j];
			}
			Tridag(ay, by, cy, ry, uy, lattice_size_y);
			for (j = 1; j <= lattice_size_y; j++)
			{
				topo[i][j] = uy[j];
			}
		}
		for (i = 1; i <= lattice_size_x; i++)
		{
			for (j = 1; j <= lattice_size_y; j++)
			{
				topo2[i][j] = topo[i][j];
			}
		}
		for (j = 1; j <= lattice_size_y; j++)
		{
			for (i = 1; i <= lattice_size_x; i++)
			{
				rx[i] = term1*(topo[i][jup[j]] + topo[i][jdown[j]]) + topoold[i][j];
				if (i == 1)
				{
					rx[i] = topoold[i][j];
				}
				if (i == lattice_size_x)
				{
					rx[i] = topoold[i][j];
				}
			}
			Tridag(ax, bx, cx, rx, ux, lattice_size_x);
			for (i = 1; i <= lattice_size_x; i++)
			{
				topo[i][j] = ux[i];
			}
		}
	}
}


void StreamPower::Avalanche(int i, int j)
{
	if (topo[iup[i]][j] - _topo[i][j] > thresh)
		topo[iup[i]][j] = _topo[i][j] + thresh;
	if (topo[idown[i]][j] - _topo[i][j] > thresh)
		topo[idown[i]][j] = _topo[i][j] + thresh;
	if (topo[i][jup[j]] - _topo[i][j] > thresh)
		topo[i][jup[j]] = _topo[i][j] + thresh;
	if (topo[i][jdown[j]] - _topo[i][j] > thresh)
		topo[i][jdown[j]] = _topo[i][j] + thresh;
	if (topo[iup[i]][jup[j]] - _topo[i][j] > (thresh*sqrt2))
		topo[iup[i]][jup[j]] = _topo[i][j] + thresh*sqrt2;
	if (topo[iup[i]][jdown[j]] - _topo[i][j] > (thresh*sqrt2))
		topo[iup[i]][jdown[j]] = _topo[i][j] + thresh*sqrt2;
	if (topo[idown[i]][jup[j]] - _topo[i][j] > (thresh*sqrt2))
		topo[idown[i]][jup[j]] = _topo[i][j] + thresh*sqrt2;
	if (topo[idown[i]][jdown[j]] - _topo[i][j] > (thresh*sqrt2))
		topo[idown[i]][jdown[j]] = _topo[i][j] + thresh*sqrt2;
}

void StreamPower::avalanche(int i, int j)
{
	if (_topo[_iup[i]][j] - _topo[i][j] > thresh)
		_topo[_iup[i]][j] = _topo[i][j] + thresh;
	if (_topo[_idown[i]][j] - _topo[i][j] > thresh)
		_topo[_idown[i]][j] = _topo[i][j] + thresh;
	if (_topo[i][_jup[j]] - _topo[i][j] > thresh)
		_topo[i][_jup[j]] = _topo[i][j] + thresh;
	if (_topo[i][_jdown[j]] - _topo[i][j] > thresh)
		_topo[i][_jdown[j]] = _topo[i][j] + thresh;
	if (_topo[_iup[i]][_jup[j]] - _topo[i][j] > (thresh*sqrt2))
		_topo[_iup[i]][_jup[j]] = _topo[i][j] + thresh*sqrt2;
	if (_topo[_iup[i]][_jdown[j]] - _topo[i][j] > (thresh*sqrt2))
		_topo[_iup[i]][_jdown[j]] = _topo[i][j] + thresh*sqrt2;
	if (_topo[_idown[i]][_jup[j]] - _topo[i][j] > (thresh*sqrt2))
		_topo[_idown[i]][_jup[j]] = _topo[i][j] + thresh*sqrt2;
	if (_topo[_idown[i]][_jdown[j]] - _topo[i][j] > (thresh*sqrt2))
		_topo[_idown[i]][_jdown[j]] = _topo[i][j] + thresh*sqrt2;
}


/*
int start()
{
	FILE *fp1;
	float deltah, time, m_ax, duration;
	int printinterval, idum, i, j, t, step;
	fp1 = fopen("streampower_topo", "w");
	lattice_size_x = 250;
	lattice_size_y = 250;
	idum = -678;
	U = 1;                // m/kyr
	K = 0.05;             // kyr^-1
	printinterval = 100;
	delt_ax = 200.0;       // m
	thresh = 0.58*delt_ax; // 30 deg
	timestep = 1;         // kyr
	duration = 100;
	setupgridneighbors();
	_topo = matrix(1, lattice_size_x, 1, lattice_size_y);
	_topo2 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_topoold = matrix(1, lattice_size_x, 1, lattice_size_y);
	_slope = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow1 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow2 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow3 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow4 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow5 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow6 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow7 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow8 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_topovec = vector(1, lattice_size_x*lattice_size_y);
	_topovecind = ivector(1, lattice_size_x*lattice_size_y);
	for (i = 1; i <= lattice_size_x; i++)
		for (j = 1; j <= lattice_size_y; j++)
		{
			_topo[i][j] = 0.5*gasdev(&idum);
			_topoold[i][j] = _topo[i][j];
			_flow[i][j] = 1;
		}
	//construct diffusional landscape for initial _flow routing
	for (step = 1; step <= 10; step++)
	{
		hill_slopediffusioninit();
		for (i = 2; i <= lattice_size_x - 1; i++)
			for (j = 2; j <= lattice_size_y - 1; j++)
			{
				_topo[i][j] += 0.1;
				_topoold[i][j] += 0.1;
			}
	}
	time = 0;
	while (time < duration)
	{//perform landsliding
		for (j = 1; j <= lattice_size_y; j++)
			for (i = 1; i <= lattice_size_x; i++)
				_topovec[(j - 1)*lattice_size_x + i] = _topo[i][j];
		indexx(lattice_size_x*lattice_size_y, _topovec, _topovecind);
		t = 0;
		while (t < lattice_size_x*lattice_size_y)
		{
			t++;
			i = (_topovecind[t]) % lattice_size_x;
			if (i == 0) i = lattice_size_x;
			j = (_topovecind[t]) / lattice_size_x + 1;
			if (i == lattice_size_x) j--;
			avalanche(i, j);
		}
		for (j = 1; j <= lattice_size_y; j++)
			for (i = 1; i <= lattice_size_x; i++)
				_topoold[i][j] = _topo[i][j];
		for (j = 1; j <= lattice_size_y; j++)
			for (i = 1; i <= lattice_size_x; i++)
				fillinpitsandflats(i, j);
		for (j = 1; j <= lattice_size_y; j++)
			for (i = 1; i <= lattice_size_x; i++)
			{
				_flow[i][j] = 1;
				_topovec[(j - 1)*lattice_size_x + i] = _topo[i][j];
			}
		indexx(lattice_size_x*lattice_size_y, _topovec, _topovecind);
		t = lattice_size_x*lattice_size_y + 1;
		while (t > 1)
		{
			t--;
			i = (_topovecind[t]) % lattice_size_x;
			if (i == 0) i = lattice_size_x;
			j = (_topovecind[t]) / lattice_size_x + 1;
			if (i == lattice_size_x) j--;
			mfd_flowroute(i, j);
		}
		for (i = 2; i <= lattice_size_x - 1; i++)
			for (j = 2; j <= lattice_size_y - 1; j++)
			{
				_topo[i][j] += U*timestep;
				_topoold[i][j] += U*timestep;
			}
		//perform upwind erosion
		m_ax = 0;
		for (i = 2; i <= lattice_size_x - 1; i++)
			for (j = 2; j <= lattice_size_y - 1; j++)
			{
				calculatealongchannel_slope(i, j);
				deltah = timestep*K*sqrt(_flow[i][j])*delt_ax*_slope[i][j];
				_topo[i][j] -= deltah;
				if (_topo[i][j]<0) _topo[i][j] = 0;
				if (K*sqrt(_flow[i][j])*delt_ax>m_ax) m_ax = K*sqrt(_flow[i][j])*delt_ax;
			}
		time += timestep;
		if (m_ax > 0.3*delt_ax / timestep)
		{
			time -= timestep;
			timestep /= 2.0;
			for (i = 2; i <= lattice_size_x - 1; i++)
				for (j = 2; j <= lattice_size_y - 1; j++)
					_topo[i][j] = _topoold[i][j] - U*timestep;
		}
		else
		{
			if (m_ax < 0.03*delt_ax / timestep) timestep *= 1.2;
			for (j = 1; j <= lattice_size_y; j++)
				for (i = 1; i <= lattice_size_x; i++)
					_topoold[i][j] = _topo[i][j];
		}
		if (time > printinterval)
		{
			printinterval += 100;
			for (i = 1; i <= lattice_size_x; i++)
				for (j = 1; j <= lattice_size_y; j++)
				{
					fprintf(fp1, "%f\n", _topo[i][j]);
				}
		}
	}
	fclose(fp1);
	return 0;
}
*/

float** StreamPower::create_random_field()
{
	float** mat = matrix(1, lattice_size_x, 1, lattice_size_y);
	int idum = -678;
	for (int i = 1; i <= lattice_size_x; i++)
	{
		for (int j = 1; j <= lattice_size_y; j++)
		{
			mat[i][j] = 0.5*gasdev(&idum);
		}
	}
	return mat;
}

std::vector<std::vector<float>> StreamPower::CreateRandomField()
{
	std::vector<std::vector<float>> mat = Matrix(1, lattice_size_x, 1, lattice_size_y);
	std::default_random_engine generator;
	std::normal_distribution<float> distribution(0.0f, 1.0f);
	for (int i = 1; i <= lattice_size_x; i++)
	{
		for (int j = 1; j <= lattice_size_y; j++)
		{
			mat[i][j] = 0.5*Gasdev(generator, distribution);
		}
	}
	return mat;
}

void StreamPower::set_topo(float** t)
{
	for (int i = 1; i <= lattice_size_x; i++)
	{
		for (int j = 1; j <= lattice_size_y; j++)
		{
			_topo[i][j] = t[i][j];
			_topoold[i][j] = _topo[i][j];
			_flow[i][j] = 1;
		}
	}
	init_diffusion();
}

void StreamPower::SetTopo(std::vector<std::vector<float>> t)
{
	for (int i = 1; i <= lattice_size_x; i++)
	{
		for (int j = 1; j <= lattice_size_y; j++)
		{
			topo[i][j] = t[i][j];
			topoold[i][j] = topo[i][j];
			flow[i][j] = 1;
		}
	}
	InitDiffusion();
}

void StreamPower::init_diffusion()
{
	//construct diffusional landscape for initial flow routing
	for (int step = 1; step <= 10; step++)
	{
		hill_slopediffusioninit();
		for (int i = 2; i <= lattice_size_x - 1; i++)
		{
			for (int j = 2; j <= lattice_size_y - 1; j++)
			{
				_topo[i][j] += 0.1;
				_topoold[i][j] += 0.1;
			}
		}
	}
	
}

void StreamPower::InitDiffusion()
{
	//construct diffusional landscape for initial flow routing
	for (int step = 1; step <= 10; step++)
	{
		HillSlopeDiffusionInit();
		for (int i = 2; i <= lattice_size_x - 1; i++)
		{
			for (int j = 2; j <= lattice_size_y - 1; j++)
			{
				topo[i][j] += 0.1;
				topoold[i][j] += 0.1;
			}
		}
	}
}

void StreamPower::init()
{

	U = 1;                // m/kyr
	K = 0.05;             // kyr^-1
	int printinterval = 100;
	deltax = 200.0;       // m
	thresh = 0.58*deltax; // 30 deg
	timestep = 1;         // kyr
	int duration = 100;
	setupgridneighbors();
	_topo = matrix(1, lattice_size_x, 1, lattice_size_y);
	_topo2 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_topoold = matrix(1, lattice_size_x, 1, lattice_size_y);
	_slope = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow1 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow2 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow3 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow4 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow5 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow6 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow7 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_flow8 = matrix(1, lattice_size_x, 1, lattice_size_y);
	_topovec = vector(1, lattice_size_x*lattice_size_y);
	_topovecind = ivector(1, lattice_size_x*lattice_size_y);

}

void StreamPower::Init()
{

	U = 1;                // m/kyr
	K = 0.05;             // kyr^-1
	int printinterval = 100;
	deltax = 200.0;       // m
	thresh = 0.58*deltax; // 30 deg
	timestep = 1;         // kyr
	int duration = 100;
	SetupGridNeighbors();
	topo = Matrix(1, lattice_size_x, 1, lattice_size_y);
	topo2 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	topoold = Matrix(1, lattice_size_x, 1, lattice_size_y);
	slope = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow1 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow2 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow3 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow4 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow5 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow6 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow7 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow8 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	topovec = Vector(1, lattice_size_x*lattice_size_y);
	topovecind = IVector(1, lattice_size_x*lattice_size_y);
}

StreamPower::StreamPower() : StreamPower(50, 50) {}

StreamPower::StreamPower(int nx, int ny) : lattice_size_x(nx), lattice_size_y(ny)
{

}

StreamPower::~StreamPower() {}

float** StreamPower::get_topo()
{
	return _topo;
}

std::vector<std::vector<float>> StreamPower::GetTopo()
{
	return topo;
}