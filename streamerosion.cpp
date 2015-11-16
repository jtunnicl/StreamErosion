/***************************************************************************
**                                                                        **
**  Landscape Erosion Model                                               **
**  Based on code by Jon Pelletier                                        **
**                                                                        **
**  For use with SAGA project                                             **
**                                                                        **
****************************************************************************/

#include "streamerosion.h"
#include "ui_mainwindow.h"
#include <QDebug>
#include <QDesktopWidget>
#include <QScreen>
#include <QMessageBox>
#include <QMetaEnum>
#include <vector>
#include <fstream>
#include <iostream>

#define sqrt2 1.414213562373
#define oneoversqrt2 0.707106781186
#define fillincrement 0.01

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float StreamErosion::ran3(int *idum)                      // Random Number generator; Numerical Recipes (original C code)
{
    static int inext,inextp;
    static long ma[56];
    static int iff=0;
    long mj,mk;
    int i,ii,k;

    if (*idum < 0 || iff == 0)
    {
        iff = 1;
        mj = MSEED - (*idum < 0 ? - *idum: *idum);
        mj %= MBIG;
        ma[55] = mj;
        mk = 1;
        for (i=1;i<=54;i++)
        {
            ii=(21*i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ) mk += MBIG;
            mj = ma[ii];
        }
        for (k=1; k<=4; k++)
            for (i=1;i<=55;i++)
            {
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }
        inext = 0;
        inextp = 31;
        *idum = 1;
    }
    if (++inext == 56) inext = 1;
    if (++inextp == 56) inextp = 1;
    mj = ma[inext] - ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext]=mj;

    return mj * FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

float StreamErosion::gasdev(int *idum)                         // draws a standard normal random variable (NR) //
{
    static int iset = 0;
    static float gset;
    float fac, r, v1, v2;

    if  (iset == 0)
    {
        {
            v1 = 2.0 * ran3(idum) - 1.0;
            v2 = 2.0 * ran3(idum) - 1.0;
            r = v1 * v1 + v2 * v2;
        }
        while (r >= 1.0);
        fac = sqrt( -2.0 * log(r) / r);
        gset = v1 * fac;
        iset = 1;
        return v2 * fac;
    } else {
            iset = 0;
            return gset;
    }
}

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 100000

void StreamErosion::avalanche( Topo *q, double thresh, int i, int j)
{
    if (q->topo[iup[i]][j] - q->topo[i][j]>thresh)
        q->topo[iup[i]][j] = q->topo[i][j] + thresh;
    if (q->topo[idown[i]][j] - q->topo[i][j]>thresh)
        q->topo[idown[i]][j] = q->topo[i][j] + thresh;
    if (q->topo[i][jup[j]] - q->topo[i][j]>thresh)
        q->topo[i][jup[j]] = q->topo[i][j] + thresh;
    if (q->topo[i][jdown[j]] - q->topo[i][j]>thresh)
        q->topo[i][jdown[j]] = q->topo[i][j] + thresh;
    if (q->topo[iup[i]][jup[j]] - q->topo[i][j]>(thresh*sqrt2))
        q->topo[iup[i]][jup[j]] = q->topo[i][j] + thresh*sqrt2;
    if (q->topo[iup[i]][jdown[j]] - q->topo[i][j]>(thresh*sqrt2))
        q->topo[iup[i]][jdown[j]] = q->topo[i][j] + thresh*sqrt2;
    if (q->topo[idown[i]][jup[j]] - q->topo[i][j]>(thresh*sqrt2))
        q->topo[idown[i]][jup[j]]  =  q->topo[i][j] + thresh*sqrt2;
    if (q->topo[idown[i]][jdown[j]] - q->topo[i][j]>(thresh*sqrt2))
        q->topo[idown[i]][jdown[j]]  =  q->topo[i][j] + thresh*sqrt2;
}

void StreamErosion::hillslopediffusioninit(Topo *q, int i, int j)
{

    float D = 10000000.0;
    int count = 0;
    float term1 = D / (deltax * deltax);


    for ( i = 1; i <= q->nx; i++)
    {
        q->ax[i] = -term1;
        q->cx[i] = -term1;
        q->bx[i] = 4 * term1 + 1;
        if (i == 1)
        {
            q->bx[i] = 1;
            q->cx[i] = 0;
        }
        if (i==q->nx)
        {
            q->bx[i] = 1;
            q->ax[i] = 0;
        }
    }

    for ( j = 1; j <= q->ny; j++ )
    {
        q->ay[j] = -term1;
        q->cy[j] = -term1;
        q->by[j] = 4 * term1 + 1;
        if (j == 1)
        {
            q->by[j] = 1;
            q->cy[j] = 0;
        }
        if (j == q->ny)
        {
            q->by[j] = 1;
            q->ay[j] = 0;
        }
    }

    while (count < 5)
    {
        count++;
        for ( i = 1; i <= q->nx; i++ )
            for (j = 1; j <= q->ny; j++ )
                q->topo2[i][j] = q->topo[i][j];
        for ( i = 1; i <= q->nx; i++ )
        {
            for ( j = 1; j <= q->ny; j++ )
            {
                q->ry[j] = term1 * ( q->topo[q->iup[i]][j] + q->topo[q->idown[i]][j] ) + q->topoold[i][j];
                if ( j == 1 )
                    q->ry[j] = q->topoold[i][j];
                if ( j==q->ny )
                    q->ry[j] = q->topoold[i][j];
            }

            tridag( &q->ay, &q->by, &q->cy, &ry, &uy );

            for ( j = 1; j <= q->ny; j++ )
                q->topo[i][j] = q->uy[j];
        }

        for ( i = 1; i <= q->nx; i++ )
            for ( j = 1; j <= q->ny; j++ )
                q->topo2[i][j] = q->topo[i][j];

        for ( j = 1; j <= q->ny; j++ )
        {
            for ( i = 1; i <= q->nx; i++)
            {
                q->rx[i] = term1 *( q->topo[i][q->jup[j]] + q->topo[i][q->jdown[j]] ) + q->topoold[i][j];
                if ( i==1 )
                    q->rx[i] = q->topoold[i][j];
                if ( i==q->nx )
                    q->rx[i] = q->topoold[i][j];
            }

            tridag( &q->ax, &q->bx, &q->cx, &q->rx, &q->ux );

            for ( i = 1; i <= q->nx; i++ )
                q->topo[i][j] = q->ux[i];
        }
    }
}

StreamErosion::StreamErosion(Topo *q)
{

    int idum = -678;
    int i, j, max, t;
    int time = 1;
    int printinterval = 100;
    int timestep = 1;                  // kyr
    int duration = 1000;

    double U = 1;                      // Uplift, m/kyr
    double K = 0.05;                   // Diffusion kyr^-1
    double deltax = 200.0;             // m
    double deltah = 0.0;
    double thresh = 0.58 * deltax;     // 30 deg

    // Topographic matrix
    vector < vector < double > > topo;           // Model Domain

    // ASCII raster file
    ifstream is("H:\\putauaki_8m.txt");

    load_matrix(&is, &topo, " ");
    q->ny = static_cast<int>(topo.size());                 // size of y
    q->nx = static_cast<int>(topo[0].size());              // size of x

    // Vectors to represent neighbour cells
    vector <int> idown (q->nx + 1, 0);
    vector <int> iup (q->nx + 1, 0);
    vector <int> jdown (q->ny + 1, 0);
    vector <int> jup (q->ny + 1 ,0);

    // Topo matrices
    //vector < vector < double > > topo (q->ny, vector<double>(q->nx));
    vector < vector < double > > topo2 (q->ny, vector<double>(q->nx));
    vector < vector < double > > topoold (q->ny, vector<double>(q->nx));

    vector < vector < double > > slope (q->ny, vector<double>(q->nx));
    vector < vector < double > > flow (q->ny, vector<double>(q->nx));
    vector < vector < double > > flow1 (q->ny, vector<double>(q->nx));
    vector < vector < double > > flow2 (q->ny, vector<double>(q->nx));
    vector < vector < double > > flow3 (q->ny, vector<double>(q->nx));
    vector < vector < double > > flow4 (q->ny, vector<double>(q->nx));
    vector < vector < double > > flow5 (q->ny, vector<double>(q->nx));
    vector < vector < double > > flow6 (q->ny, vector<double>(q->nx));
    vector < vector < double > > flow7 (q->ny, vector<double>(q->nx));
    vector < vector < double > > flow8 (q->ny, vector<double>(q->nx));

    vector < double > topovec ( q->nx * q->ny, 0 );
    vector < int > topovecind ( q->nx * q->ny, 0 );

    setupColorMapDemo(ui->customPlot, &matrix);
    setWindowTitle("Stream Erosion");
    statusBar()->clearMessage();

    setupgridneighbors( q->nx, q->ny, &idown, &iup, &jdown, &jup );

  //   ******  Start Run  ******   //

    while ( time < duration )
    {                    
        for ( j = 1; j <= q->ny ; j++)                         // perform landsliding //
          for (i = 1; i <= q->nx; i++)
              topovec[ (j-1) * q->nx+i ] = topo[i][j];

        Indexx( topovec );

        t = 0;
        while ( t < q->nx * q->ny )
        {
            t++;
            i = ( topovecind[t] ) % q->nx;
            if ( i==0 ) i = q->nx;
            j = ( topovecind[t] ) / q->nx+1;
            if ( i==q->nx ) j--;
            avalanche( &topo, &idown, &iup, &jdown, &jup, thresh, i, j );                                 // Call to avalance subroutine
        }

        for ( j=1; j<=q->ny; j++ )
            for ( i=1; i<=q->nx; i++ )
                topoold[i][j]=topo[i][j];                      //  Save topo as topoold

        for (j=1; j<=q->ny; j++ )
            for (i=1; i<=q->nx; i++ )
                fillinpitsandflats( &topo, &idown, &iup, &jdown, &jup, i, j );                     //  Call to pits and flats routine

        for ( j=1; j<=q->ny; j++ )
            for ( i=1; i <= q->nx; i++ )
            {
                flow[i][j] = 1;
                topovec[ (j-1) * q->nx + i ] = topo[i][j];
            }

        Indexx( topovec );

        t = q->nx * q->ny+1;

        while ( t > 1 )
        {
            t--;
            i=(topovecind[t])%q->nx;
            if ( i==0 )
                i=q->nx;
            j=(topovecind[t])/q->nx+1;
            if ( i==q->nx)
                j--;
            mfdflowroute( &topo, &flow, &flow1, &flow2, &flow3, &flow4,
                          &flow5, &flow6, &flow7, &flow8, &iup, &idown,
                          &jdown, &jup, i,j);
        }

        for (i=2;i<=q->nx-1;i++)
            for (j=2;j<=q->ny-1;j++)
            {
                topo[i][j] += U*timestep;
                topoold[i][j] += U*timestep;
            }

        //   perform upwind erosion   //
        max = 0;
        for ( i=2; i<=q->nx-1; i++ )
            for (j=2; j<=q->ny-1; j++)
            {
                calculatealongchannelslope(&topo, &slope, &iup, &idown, &jdown, &jup, i,j);
                deltah = timestep * K * sqrt(flow[i][j]) * deltax * slope[i][j];
                topo[i][j] -= deltah;
                if (topo[i][j]<0)
                    topo[i][j] = 0;
                if (K*sqrt(flow[i][j])*deltax>max)
                    max = K * sqrt(flow[i][j]) * deltax;
            }

        time += timestep;
        if ( max > 0.3 * deltax/timestep )
        {
            time -= timestep;
            timestep /= 2.0;
            for ( i=2; i<=q->nx-1; i++ )
                for ( j=2; j<=q->ny-1; j++ )
                    topo[i][j] = topoold[i][j] - U * timestep;
        }
        else
        {
            if ( max < 0.03 *deltax / timestep)
                timestep *= 1.2;
            for ( j=1; j<=q->ny; j++ )
                for ( i=1; i<=q->nx; i++ )
                    topoold[i][j] = topo[i][j];
        }

//        if ( time > printinterval )
//        {
//            printinterval += 100;
//                for ( i=1; i<=q->nx; i++)
//                    for (j=1; j<=q->ny; j++)
//                    {
//                        fprintf( fp1, "%f\n", topo[i][j] );
//                    }
//        }
    }
}

