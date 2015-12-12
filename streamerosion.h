#ifndef STREAMEROSION_H
#define STREAMEROSION_H

//#include <QMainWindow>
//#include <QTimer>
#include <vector>
#include "topo.h"
//include "qcustomplot.h"
#include <cmath>

using namespace std;

typedef struct StreamErosionModelParameters {
	int timestep;
	int duration;
	double U;		// Uplift, m/kyr
	double K;		// Diffusion kyr^-1
	double deltax;	// m
	double deltah;
	double thresh;	// 30 deg
} Parameters;

class StreamErosion
{

protected:

    Topo* q;
	Parameters p;
    int idum;
	int time;
	int printinternval;

    float ran3(int *idum);
    float gasdev(int *idum);
    void avalanche( Topo *q, double thresh, int i, int j );
    void avalanche( );
    void hillslopediffusioninit(Topo *q, int i, int j );

public:

    StreamErosion();
    StreamErosion(Topo* q, Parameters& p);
    ~StreamErosion();
	void Start();

};


#endif // STREAMEROSION_H

