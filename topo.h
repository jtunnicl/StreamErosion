#ifndef TOPO_H
#define TOPO_H

#include <vector>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;


class Topo
{

public:

	int nx, ny;                                   // Size of grid (elements)
	int xllcorner, yllcorner;
	int nodata;
    int elems;                                    // x * y
    double deltax;                                // Grid spacing (meters)

    int year, month, day;
    int hour, minute, second;
    int lat, longit;                               // Coordinates for solar calcs

	vector<vector<vector<double>>> strat;			  // Stratigraphy vector for each topo pixel

    vector<vector<double>> topo;              // Topographic grid
    vector<vector<double>> topo2;             // Computation results raster
    vector<vector<double>> topoold;           // Value holder, when updating
    vector<vector<double>> slope;             // Topographic gradient  
	vector<vector<double>> flow, flow1;       // Raster that holds flow magnitude within each cell
    vector<vector<double>> flow2, flow3;
    vector<vector<double>> flow4, flow5;
    vector<vector<double>> flow6, flow7, flow8;
    vector<vector<double>> topoVec;				// Sorting array for pit-filling
    vector<vector<double>> vegState;			// Veg growth raster
    vector<vector<double>> exposureAge;			// For ice-rich ground - ablation time
   
	vector<vector<int>> fixedPts;				// Points, such as bedrock, not subject to sediment movement
	vector<int> topoVecInd;				// Index of sorted values in topoVec

    vector<int> iup, idown;                  // Neighbour pixels
    vector<int> jup, jdown;                  // Neighbour pixels
    vector<int> ax, bx, cx;                  // Neighbour pixels
    vector<int> ay, by, cy;                  // Neighbour pixels
    vector<int> ux, rx;                      // Neighbour pixels
    vector<int> uy, ry;                      // Neighbour pixels	


	Topo();
	Topo(const char* dem);
	~Topo();


	void LoadDEM(ifstream& in, const string& delim);
	void indexX(int n, vector <double>& arr, vector<double>& indx);
	void IndexX();
	static void PrintMatrix(const vector<vector<double>>& mat, int nx, int ny);
	static void PrintMatrixRow(const vector<vector<double>>& mat, int nx, int row);

	//void triDag(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& r, vector<double>& u);
	void triDag();

protected:
	void Init();
	void SetupGridNeighbors();


};

#endif // TOPO_H
