#pragma once

#include <vector>
#include <random>
#include <numeric>
#include <algorithm>

#define NR_END 1
#define FREE_ARG char*

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 100000

#define sqrt2 1.414213562373f
#define oneoversqrt2 0.707106781186f
#define fillincrement 0.01f

class StreamPower
{
public:

	int lattice_size_x, lattice_size_y;
	float U, K, D, timestep, delt_ax, thresh, thresholdarea;

	// old vars
	float **_flow1, **_flow2, **_flow3, **_flow4, **_flow5, **_flow6, **_flow7, **_flow8, **_flow;
	float **_topo, **_topoold, **_topo2, **_slope, *_ax, *_ay, *_bx, *_by, *_cx, *_cy, *_ux, *_uy;
	float *_rx, *_ry, *_topovec;
	int *_topovecind, *_iup, *_idown, *_jup, *_jdown;

	// new vars
	std::vector<int> iup, idown, jup, jdown;

	StreamPower();
	StreamPower(int nx, int ny);
	~StreamPower();

	static void free_ivector(int *v, long nl, long nh);
	static void free_vector(float *v, long nl, long nh);

	static float* vector(int nl, int nh);				// old implementation
	static std::vector<float> Vector(int nl, int nh);	// new implementation

	static int* ivector(int nl, int nh);				// old implementation
	static std::vector<int> IVector(int nl, int nh);	// new implementation

	static float** matrix(long nrl, long nrh, long ncl, long nch);							// old implementation
	static std::vector<std::vector<float>> Matrix(int nrl, int nrh, int ncl, int nch);		// new implementation

	static int** imatrix(long nrl, long nrh, long ncl, long nch);							// old implementation
	static std::vector<std::vector<int>> IMatrix(int nrl, int nrh, int ncl, int nch);		// new implementation

	static float ran3(int* idum);	// old implementation
	static float Ran3(std::default_random_engine& generator, std::uniform_real_distribution<float>& distribution);	// new implementation

	static float gasdev(int* idnum);	// old implementation
	static float Gasdev(std::default_random_engine& generator, std::normal_distribution<float>& distribution);	// new implementation

	static void indexx(int n, float* arr, int* indx);	// old implementation
	static void Indexx(int n, float* arr, int* indx);	// interface from old to new implementation
	static std::vector<int> Indexx(std::vector<float>& arr);	// new implementation

	static void tridag(float a[], float b[], float c[], float r[], float u[], unsigned long n); // old implementation
	static void Tridag(float a[], float b[], float c[], float r[], float u[], unsigned long n); // interface from old to new implementation
	static void Tridag(std::vector<float>& a, std::vector<float>& b, std::vector<float>& c, std::vector<float>& r, std::vector<float>& u, int n); // interface from old to new implementation

	void setupgridneighbors(); // old implementation
	void SetupGridNeighbors(); // new implementation
};

template <typename T> std::vector<T> Arr_ayToVector(T* a, int size)
{
	std::vector<T> v = std::vector<T>(size);
	for (int i = 0; i < size; i++)
	{
		v[i] = a[i];
	}
	return v;
}

template <typename T> std::vector<T> Arr_ayToVector(T* a, int size, bool fortranIndexing)
{
	std::vector<T> v;
	if (fortranIndexing)
	{
		v = std::vector<T>(size + 1);
	}
	else 
	{
		v = std::vector<T>(size);
	}

	for (int i = 0; i < size; i++)
	{
		v[i] = a[i];
	}
	return v;
}

template <typename T> void VectorToArr_ay(std::vector<T>& v, T* a)
{
	for (int i = 0; i < v.size(); i++)
	{
		a[i] = v[i];
	}
}

// http://stackover_flow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T> std::vector<int> SortIndices(const std::vector<T>& v)
{

	// initialize original index locations
	std::vector<int> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] < v[i2]; });

	return idx;
}

template <typename T> std::vector<int> SortFortranIndices(const std::vector<T>& v)
{

	// initialize original index locations
	std::vector<int> idx(v.size());
	std::iota(idx.begin()+1, idx.end(), 1);

	// sort indexes based on comparing values in v
	std::sort(idx.begin()+1, idx.end(), [&v](int i1, int i2) {return v[i1] < v[i2]; });

	return idx;
}
