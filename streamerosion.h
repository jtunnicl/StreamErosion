#ifndef STREAMEROSION_H
#define STREAMEROSION_H

//#include <QMainWindow>
//#include <QTimer>
#include <vector>
#include "topo.h"
//include "qcustomplot.h"
#include <cmath>

using namespace std;

class StreamErosion
{

public:

    Topo q;
    int idum;
    int timestep, printinternval, duration;
    double U, K, deltax, thresh;

    float ran3(int *idum);
    float gasdev(int *idum);
    void avalanche( Topo *q, double thresh, int i, int j );
    void avalanche( );
    void hillslopediffusioninit(Topo *q, int i, int j );

    StreamErosion();
    StreamErosion(Topo* q);
    ~StreamErosion();

};


#endif // STREAMEROSION_H

