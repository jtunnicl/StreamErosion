#pragma once

#ifdef STREAMDLL_EXPORTS
#define STREAMDLL_API __declspec(dllexport) 
#else
#define STREAMDLL_API __declspec(dllimport) 
#endif

#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include "Array2D.hpp"

#define sqrt2 1.414213562373f
#define oneoversqrt2 0.707106781186f
#define fillincrement 0.01f

typedef struct StreamErosionModelParameters
{
	int timestep;
	int duration;
	double K;		// Diffusion kyr^-1
} Parameters;

class StreamPower
{
public:

	// parameters
	int lattice_size_x, lattice_size_y, duration, printinterval, printstep;
	float K, D, timestep, deltax, thresh, thresholdarea, time;
	float xllcorner, yllcorner, nodata;

	// internal variables
	std::vector<int> iup, idown, jup, jdown, topovecind;
	std::vector<float> ax, ay, bx, by, cx, cy, ux, uy, rx, ry, topovec;
	std::vector<std::vector<float>> topo, topoold, topo2, slope, flow, flow1, flow2, flow3, flow4, flow5, flow6, flow7, flow8, U;
	Array2D<float> elevation;

	static float Ran3(std::default_random_engine& generator, std::uniform_real_distribution<float>& distribution);
	static float Gasdev(std::default_random_engine& generator, std::normal_distribution<float>& distribution);
	static std::vector<int> Indexx(std::vector<float>& arr);	// done
	static void Tridag(std::vector<float>& a, std::vector<float>& b, std::vector<float>& c, std::vector<float>& r, std::vector<float>& u, int n); // done

	STREAMDLL_API StreamPower(Parameters p);
	STREAMDLL_API ~StreamPower();

	STREAMDLL_API std::vector<std::vector<float>> CreateRandomField();
	STREAMDLL_API std::vector<std::vector<float>> ReadArcInfoASCIIGrid(char* fname);
	STREAMDLL_API std::vector<std::vector<float>> GetTopo();

	STREAMDLL_API void SetTopo(std::vector<std::vector<float>> t);
	STREAMDLL_API void SetU(std::vector<std::vector<float>> u);
	STREAMDLL_API void SetU(float u);
	STREAMDLL_API void Start();
	STREAMDLL_API void Step();
	STREAMDLL_API void PrintState(char* fname);
	STREAMDLL_API void Init(int nx, int ny, float xllcorner, float yllcorner, float deltax, float nodata);
	void AssignVariables();
	void SetInitialValues(std::vector<std::vector<float>> t);
	void SetupGridNeighbors();
	void HillSlopeDiffusion();
	void InitDiffusion();
	void Avalanche(int i, int j);
	void CalculateAlongChannelSlope(int i, int j);
	void MFDFlowRoute(int i, int j);
	void Flood();

};

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
