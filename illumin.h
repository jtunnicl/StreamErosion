#ifndef ILLUMIN_H
#define ILLUMIN_H

#include "topo.h"

class Illumin
{
public:
    Illumin();

    void solar ()             // Object that holds solar info
    {
        public:

            double dec;
            double sha;       // Solar Hour Angle is 0° at solar noon. Since the Earth rotates 15° per hour,
                              //      each hour away from solar noon corresponds to an angular motion of the sun in the sky of 15°.
                              //      In the morning the hour angle is negative, in the afternoon the hour angle is positive.
            double AST;       // Apparent Solar Time
            double alt;       // altitude
            double azm;       // azimuth
    }

    void solarRad( Topo *q )
    {
        int Year, Mnth, Day;
        int h, Minute, Lat, longit

    }

    void hillshape ( Topo *q )
    {

    }



};

#endif // ILLUMIN_H
