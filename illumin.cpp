#include "illumin.h"

#define PI 3.1415927

Illumin::Illumin()
{
    // Constructor stuff, here.

    // values for 'solar' need to be initialized, somewhere..

}

void Illumin::solar(double GHA, double alpha)       // Memory object for solar attributes
{
    double m, n;

    double dec = atan(tan(eps * DegRad) * sin(alpha * DegRad)) / DegRad;

    double sha = GHA - longit;

    double AST = 12.0 + (sha/15.0);

    m = sin(Lat * DegRad) * sin(dec * DegRad);
    n = cos(Lat * DegRad) * cos(dec * DegRad) * cos((sha+180.0) * DegRad);
    alt = asin(m-n)/DegRad;

    m = cos(Lat * DegRad) * tan(dec * DegRad);
    n = sin(Lat * DegRad) * cos((sha+180.0) * DegRad);
    azm = acos((m + n) * cos(dec * DegRad) / cos(alt * DegRad)) / DegRad;

    if ( sha > 0 )
    {
      azm = 360 - azm;
    }
}

void Illum::solarRad( Topo *q, solar *s )           // Update solar object (above), based on new time/date
{

    double m, n, t, eps, G;
    double Year1, Mnth1;
    double YStdMed = q->stdmed;
    double DegRad = PI / 180;
    double C, L, alpha;
    double EOT, GHA;

  // Calculate Greenwich Hour Angle
    double LocalTime = q->hour + ( q->minute / 60 );
    double UT = LocalTime + ( YStdMed / 15.0 );

    if (q->month > 2)
    {
        Year1 = q->year;
        Mnth1 = q->month - 3;
    }
    else
    {
        Year1 = q->year - 1;
        Mnth1 = q->month + 9;
    }

    m = 30.6 * Mnth1 + 0.5;
    n = 365.25 * (Year1 - 1976);
    t = ((UT/24) + q->day + m + n - 8707.5) / 36525;   // Eqn 1.4.2
    eps = 23.4393 - 0.013 * t;                      // obliquity of the ecliptic
    G = 357.528 + 35999.05 * t;                     // mean anomaly
    if (abs(G) > 360)
    {
        G = G % 360.0;
    }

    C = 1.915 * sin(G * DegRad)+ 0.020 * sin(2 * G * DegRad);     // Correction to centre
    L = 280.460 + 36000.770 * t + C;                // mean longitude

    if (abs(L) > 360)
    {
        L = L % 360.0;
    }

    alpha = L - 2.466 * sin(2 * L * DegRad) + 0.053 * sin(4 * L * DegRad);
                                                    // Right ascension of sun (apparent) in degrees
    EOT = (L - C - alpha) / 15;                     // Equation of Time Ranges from -15 - +15 -ish ?
    GHA = 15 * UT - 180.0 - C + L - alpha;          // Greenwich Hour Angle is the angle of Aries at Universal Time, measured from
                                                    // 0° to 360°, increasing towards the West.
    if (abs(GHA) > 360)
    {
        GHA = GHA % 360;
    }

    solar(GHA);                                     // Update solar model

}

void Illum::hillshape ( Topo *q, solar *s )
{

    // Whole bunch of MATLAB operation need to be C++ified, here

    double az = s->azm;
    double alt = s->alt;
    double zfac = 0.5;

    // convert light azimuth to polar angle in radians
    az = 360.0 - az + 90;
    az(az>=360) = az - 360;     // !!! MATLAB !!!
    az = az * (pi/180);

    // convert lighting altitude to zenith angle in radians
    alt = (90-alt) * (pi/180);

    // calculate slope and aspect in radians
    [asp,grad] = cart2pol(q.slope);                 //  !!! MATLAB !!!  cartesion to polar
    grad = atan( zfac * grad);                      //  "",  scale slope by z-factor

    // convert aspect angle
    asp( asp < pi ) = asp( asp < pi ) + ( pi / 2 ); //  !!! MATLAB !!!
    asp( asp < 0 ) = asp( asp < 0 ) + ( 2 * pi );

    // hillshade calculation
    topo_hs = 255.0*( (cos(alt).*cos(grad) ) + ( sin(alt).*sin(grad).*cos(az-asp)) ); // ESRIs algorithm
    topo_hs(topo_hs<0) = 0; // set hillshade values to min of 0.
    %topo_hs -= 180;
    %topo_hs(abs(topo_hs)>5) = 0;

    // zeros border
    topo_b = zeros(size(topo_hs));
    topo_b(2:end-1,2:end-1) += topo_hs(2:end-1,2:end-1);
    topo_hs = topo_b;
    imagesc(topo_hs);
    colormap gray;
}
