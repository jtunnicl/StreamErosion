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

    Topo();

    int nx, ny;                                   // Size of grid (elements)
    int elems;                                    // x * y
    double deltax;                                // Grid spacing (meters)
    int year, month, day;
    int hour, minute, second;
    int lat, longit;                               // Coordinates for solar calcs

    vector <vector <double> > topo;              // Topographic grid
    vector <vector <double> > topo2;             // Computation results raster
    vector <vector <double> > topoold;           // Value holder, when updating
    vector <vector <double> > slope;             // Topographic gradient
    vector <vector <int> > fixedPts;             // Points, such as bedrock, not subject to sediment movement
    vector <vector <double> > flow, flow1;       // Raster that holds flow magnitude within each cell
    vector <vector <double> > flow2, flow3;
    vector <vector <double> > flow4, flow5;
    vector <vector <double> > flow6, flow7, flow8;
    vector <vector <double> > topoVec;           // Sorting array for pit-filling
    vector <vector <double> > topoVecInd;        // Index of sorted values in topoVec
    vector <vector <double> > vegState;          // Veg growth raster
    vector <vector <double> > exposureAge;       // For ice-rich ground - ablation time
    vector <vector <vector <double> > > strat;   // Stratigraphy vector for each topo pixel
    vector <double> iup, idown;                  // Neighbour pixels
    vector <double> jup, jdown;                  // Neighbour pixels
    vector <double> ax, bx, cx;                  // Neighbour pixels
    vector <double> ay, by, cy;                  // Neighbour pixels
    vector <double> ux, rx;                      // Neighbour pixels
    vector <double> uy, ry;                      // Neighbour pixels


    void loadMatrix(ifstream *is, const string& delim);

    void setupgridneighbors();
	
	void indexX(int n, vector <double>& arr, vector<double>& indx);
	void indexX();

	//void triDag(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& r, vector<double>& u);
	void triDag();

};

#endif // TOPO_H
