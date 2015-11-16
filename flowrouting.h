#ifndef FLOWROUTING_H
#define FLOWROUTING_H

#include "topo.h"
#include "streamerosion.h"

using namespace std;

class FlowRouting
{

public:

    explicit FlowRouting();
    ~FlowRouting();

    void fillinpitsandflats(Topo *q, int i, int j );

    void mfdflowroute(Topo *q, int i, int j );

    void calculatealongchannelslope( Topo *q, int i, int j );

};

#endif // FLOWROUTING_H
