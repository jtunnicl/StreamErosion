#include "flowrouting.h"
#include <cstdlib>
#include <iostream>

#define sqrt2 1.414213562373
#define oneoversqrt2 0.707106781186
#define fillincrement 0.01

FlowRouting::FlowRouting()
{
    // Constructor stuff, here.

}

void FlowRouting::fillinpitsandflats()
{
	
	cout << "Flowrouting::fillinpits and flats not implemented" << endl;
	exit(1);

}
void FlowRouting::fillinpitsandflats( Topo *q, int i, int j )
{
    double min;

    min = q->topo[i][j];

    if (q->topo[q->iup[i]][j] < min)
        min = q->topo[q->iup[i]][j];
    if (q->topo[q->idown[i]][j] < min)
        min = q->topo[q->idown[i]][j];
    if (q->topo[i][q->jup[j]] < min)
        min = q->topo[i][q->jup[j]];
    if (q->topo[i][q->jdown[j]] < min)
        min = q->topo[i][q->jdown[j]];
    if (q->topo[q->iup[i]][q->jup[j]] < min)
        min = q->topo[q->iup[i]][q->jup[j]];
    if (q->topo[q->idown[i]][q->jup[j]] < min)
        min = q->topo[q->idown[i]][q->jup[j]];
    if (q->topo[q->idown[i]][q->jdown[j]] < min)
        min = q->topo[q->idown[i]][q->jdown[j]];
    if (q->topo[q->iup[i]][q->jdown[j]] < min)
        min = q->topo[q->iup[i]][q->jdown[j]];
    if ((q->topo[i][j] <= min) && (i>1)&&(j>1) &&
            (i < q->nx) && (j < q->ny))
    {
        q->topo[i][j] = min + fillincrement;
        fillinpitsandflats(q, i, j);
        fillinpitsandflats(q, q->iup[i], j);
        fillinpitsandflats(q, q->idown[i], j);
        fillinpitsandflats(q, i, q->jup[j]);
        fillinpitsandflats(q, i, q->jdown[j]);
        fillinpitsandflats(q, q->iup[i], q->jup[j]);
        fillinpitsandflats(q, q->idown[i], q->jup[j]);
        fillinpitsandflats(q, q->idown[i], q->jdown[j]);
        fillinpitsandflats(q, q->iup[i], q->jdown[j]);
    }
}

void FlowRouting::mfdflowroute()
{
	
	cout << "Flowrouting::mfdflowroute and flats not implemented" << endl;
	exit(1);

}
void FlowRouting::mfdflowroute( Topo *q, int i, int j )
{
    float tot;

    tot=0;

    if (q->topo[i][j] > q->topo[q->iup[i]][j])
        tot += pow(q->topo[i][j]-q->topo[q->iup[i]][j],1.1);
    if (q->topo[i][j] > q->topo[q->idown[i]][j])
        tot += pow(q->topo[i][j]-q->topo[q->idown[i]][j],1.1);
    if (q->topo[i][j] > q->topo[i][q->jup[j]])
        tot += pow(q->topo[i][j]-q->topo[i][q->jup[j]],1.1);
    if (q->topo[i][j] > q->topo[i][q->jdown[j]])
        tot += pow(q->topo[i][j]-q->topo[i][q->jdown[j]],1.1);
    if (q->topo[i][j] > q->topo[q->iup[i]][q->jup[j]])
        tot += pow((q->topo[i][j]-q->topo[q->iup[i]][q->jup[j]]) * oneoversqrt2,1.1);
    if (q->topo[i][j] > q->topo[q->iup[i]][q->jdown[j]])
        tot += pow((q->topo[i][j]-q->topo[q->iup[i]][q->jdown[j]]) * oneoversqrt2,1.1);
    if (q->topo[i][j] > q->topo[q->idown[i]][q->jup[j]])
        tot += pow((q->topo[i][j]-q->topo[q->idown[i]][q->jup[j]]) * oneoversqrt2,1.1);
    if (q->topo[i][j] > q->topo[q->idown[i]][q->jdown[j]])
        tot += pow((q->topo[i][j]-q->topo[q->idown[i]][q->jdown[j]]) * oneoversqrt2,1.1);
    if (q->topo[i][j] > q->topo[q->iup[i]][j])
        q->flow1[i][j] = pow(q->topo[i][j]-q->topo[q->iup[i]][j],1.1) / tot;
    else
        q->flow1[i][j] = 0;

    if (q->topo[i][j] > q->topo[q->idown[i]][j])
        q->flow2[i][j] = pow(q->topo[i][j]-q->topo[q->idown[i]][j],1.1)/tot;
    else
        q->flow2[i][j] = 0;

    if (q->topo[i][j] > q->topo[i][q->jup[j]])
        q->flow3[i][j] = pow(q->topo[i][j]-q->topo[i][q->jup[j]],1.1)/tot;
    else
        q->flow3[i][j] = 0;

    if (q->topo[i][j] > q->topo[i][q->jdown[j]])
        q->flow4[i][j] = pow(q->topo[i][j]-q->topo[i][q->jdown[j]],1.1)/tot;
    else
        q->flow4[i][j] = 0;

    if (q->topo[i][j] > q->topo[q->iup[i]][q->jup[j]])
        q->flow5[i][j] = pow((q->topo[i][j]-q->topo[q->iup[i]][q->jup[j]]) * oneoversqrt2,1.1)/tot;
    else
        q->flow5[i][j] = 0;

    if (q->topo[i][j] > q->topo[q->iup[i]][q->jdown[j]])
        q->flow6[i][j] = pow((q->topo[i][j]-q->topo[q->iup[i]][q->jdown[j]]) * oneoversqrt2,1.1)/tot;
    else
        q->flow6[i][j]=0;

    if (q->topo[i][j] > q->topo[q->idown[i]][q->jup[j]])
        q->flow7[i][j] = pow((q->topo[i][j]-q->topo[q->idown[i]][q->jup[j]]) * oneoversqrt2,1.1)/tot;
    else
        q->flow7[i][j] = 0;

    if (q->topo[i][j] > q->topo[q->idown[i]][q->jdown[j]])
        q->flow8[i][j] = pow((q->topo[i][j]-q->topo[q->idown[i]][q->jdown[j]]) * oneoversqrt2,1.1)/tot;
    else
        q->flow8[i][j] = 0;

    q->flow[q->iup[i]][j] += q->flow[i][j] * q->flow1[i][j];
    q->flow[q->idown[i]][j] += q->flow[i][j] * q->flow2[i][j];
    q->flow[i][q->jup[j]] += q->flow[i][j] * q->flow3[i][j];
    q->flow[i][q->jdown[j]] += q->flow[i][j] * q->flow4[i][j];
    q->flow[q->iup[i]][q->jup[j]] += q->flow[i][j] * q->flow5[i][j];
    q->flow[q->iup[i]][q->jdown[j]] += q->flow[i][j] * q->flow6[i][j];
    q->flow[q->idown[i]][q->jup[j]] += q->flow[i][j] * q->flow7[i][j];
    q->flow[q->idown[i]][q->jdown[j]] += q->flow[i][j] * q->flow8[i][j];
}

void FlowRouting::calculatealongchannelslope()
{


}


void FlowRouting::calculatealongchannelslope( Topo *q, int i, int j )
{
    float down;

    down=0;

    if (q->topo[q->iup[i]][j] - q->topo[i][j] < down)
        down=q->topo[q->iup[i]][j] - q->topo[i][j];
    if (q->topo[q->idown[i]][j] - q->topo[i][j] < down)
        down=q->topo[q->idown[i]][j] - q->topo[i][j];
    if (q->topo[i][q->jup[j]] - q->topo[i][j] < down)
        down=q->topo[i][q->jup[j]] - q->topo[i][j];
    if (q->topo[i][q->jdown[j]] - q->topo[i][j] < down)
        down=q->topo[i][q->jdown[j]] - q->topo[i][j];
    if ((q->topo[q->iup[i]][q->jup[j]] - q->topo[i][j]) * oneoversqrt2 < down)
        down=(q->topo[q->iup[i]][q->jup[j]] - q->topo[i][j]) * oneoversqrt2;
    if ((q->topo[q->idown[i]][q->jup[j]]  -q->topo[i][j]) * oneoversqrt2 < down)
        down=(q->topo[q->idown[i]][q->jup[j]] - q->topo[i][j]) * oneoversqrt2;
    if ((q->topo[q->iup[i]][q->jdown[j]] - q->topo[i][j]) * oneoversqrt2 < down)
        down=(q->topo[q->iup[i]][q->jdown[j]] - q->topo[i][j]) * oneoversqrt2;
    if ((q->topo[q->idown[i]][q->jdown[j]] - q->topo[i][j]) * oneoversqrt2 < down)
        down=(q->topo[q->idown[i]][q->jdown[j]] - q->topo[i][j]) * oneoversqrt2;

    q->slope[i][j] = fabs(down) / q->deltax;
}

