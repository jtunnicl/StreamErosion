#ifndef FLOWROUTING_H
#define FLOWROUTING_H

#include "topo.h"
//#include "streamerosion.h"

using namespace std;

class FlowRouting
{

// make these static
public:

    explicit FlowRouting();
    ~FlowRouting();

    static void fillinpitsandflats(Topo *q, int i, int j );
    static void fillinpitsandflats();

    static void mfdflowroute(Topo *q, int i, int j );
    static void mfdflowroute();

    static void calculatealongchannelslope( Topo *q, int i, int j );
    static void calculatealongchannelslope();

};

#endif // FLOWROUTING_H
