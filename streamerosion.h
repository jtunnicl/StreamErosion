#ifndef STREAMEROSION_H
#define STREAMEROSION_H

#include <QMainWindow>
#include <QTimer>
#include <vector>
#include "topo.h"
#include "qcustomplot.h"
#include <cmath>

using namespace std;



class StreamErosion
{

public:

    int idum;
    Topo q;
    int timestep, printinternval, duration;
    double U, K, deltax, thresh;

    float ran3(int *idum);

    float gasdev(int *idum);

    void avalanche( Topo *q, double thresh, int i, int j );

    void hillslopediffusioninit(Topo *q, int i, int j );

    explicit StreamErosion();
    ~StreamErosion();

}


#endif // STREAMEROSION_H

