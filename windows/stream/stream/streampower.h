#pragma once

#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include "Array2D.hpp"

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

	int lattice_size_x, lattice_size_y, duration, printinterval, printstep;
	float U, K, D, timestep, deltax, thresh, thresholdarea;

	// new vars
	float xllcorner, yllcorner, nodata;
	std::vector<int> iup, idown, jup, jdown, topovecind;
	std::vector<float> ax, ay, bx, by, cx, cy, ux, uy, rx, ry, topovec;
	std::vector<std::vector<float>> topo, topoold, topo2, slope, flow, flow1, flow2, flow3, flow4, flow5, flow6, flow7, flow8;
	Array2D<float> elevation;

	static std::vector<float> Vector(int nl, int nh);
	static std::vector<int> IVector(int nl, int nh);
	static std::vector<std::vector<float>> Matrix(int nrl, int nrh, int ncl, int nch);			
	static std::vector<std::vector<int>> IMatrix(int nrl, int nrh, int ncl, int nch);	

	static float Ran3(std::default_random_engine& generator, std::uniform_real_distribution<float>& distribution);
	static float Gasdev(std::default_random_engine& generator, std::normal_distribution<float>& distribution);

	static void Indexx(int n, float* arr, int* indx);	// interface from old to new implementation
	static std::vector<int> Indexx(std::vector<float>& arr);	// new implementation

	static std::vector<int> Indexx_C(std::vector<float>& arr);	// done

	static void Tridag(float a[], float b[], float c[], float r[], float u[], unsigned long n); // interface from old to new implementation
	static void Tridag(std::vector<float>& a, std::vector<float>& b, std::vector<float>& c, std::vector<float>& r, std::vector<float>& u, int n); // new implementation

	static void Tridag_C(std::vector<float>& a, std::vector<float>& b, std::vector<float>& c, std::vector<float>& r, std::vector<float>& u, int n); // done

	StreamPower(int nx, int ny);
	~StreamPower();

	std::vector<std::vector<float>> CreateRandomField();
	std::vector<std::vector<float>> CreateRandomField_C(); //done

	std::vector<std::vector<float>> ReadArcInfoASCIIGrid(char* fname);
	std::vector<std::vector<float>> ReadArcInfoASCIIGrid_C(char* fname); // done

	std::vector<std::vector<float>> GetTopo(); //done


	void Init(); //done

	void AssignVariables_C();
	void SetInitialValues_C(std::vector<std::vector<float>> t);

	void SetTopo(std::vector<std::vector<float>> t);
	void SetTopo_C(std::vector<std::vector<float>> t); // done

	void SetupGridNeighbors(); 
	void SetupGridNeighbors_C(); // done

	void HillSlopeDiffusionInit();
	void HillSlopeDiffusionInit_C(); //done
	
	void InitDiffusion();
	void InitDiffusion_C(); //done
	
	void Avalanche(int i, int j); //done
	
	void CalculateAlongChannelSlope(int i, int j); //done

	void MFDFlowRoute(int i, int j); //done

	void Flood(); // Barnes pit filling
	void Flood_C(); //done

	void Start(); 
	void Start_C(); // done

	void PrintState(char* fname);
	void PrintState_C(char* fname); // done
};

template <typename T> std::vector<T> ArrayToVector(T* a, int size)
{
	std::vector<T> v = std::vector<T>(size);
	for (int i = 0; i < size; i++)
	{
		v[i] = a[i];
	}
	return v;
}

template <typename T> std::vector<T> ArrayToVector(T* a, int size, bool fortranIndexing)
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

template <typename T> void VectorToArray(std::vector<T>& v, T* a)
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
