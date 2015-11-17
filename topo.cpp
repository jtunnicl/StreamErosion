#include "topo.h"
#include <iostream>
#include <cstdlib>


#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 100000


Topo::Topo()
{

    nx = 100;
    ny = 100;
    elems = nx * ny;
    deltax = 0;
    ifstream in("putauaki.txt"); // Load matrix into memory

    loadMatrix(&in, " ");

}

void Topo::loadMatrix(ifstream *in, const string& delim = " ")
{
    string line;
    string strnum;
    double number;

    // clear topo first
    topo.clear();

    // parse line by line
    while (getline(*in, line))
    {
        topo.push_back(vector<double>());

        for (string::const_iterator i = line.begin(); i != line.end(); ++ i)
        {
            // If i is not a delim, then append it to strnum
            if (delim.find(*i) == string::npos)
            {
                strnum += *i;
                if (i + 1 != line.end()) // If it's the last char, do not continue
                    continue;
            }

            // if strnum is still empty, it means the previous char is also a
            // delim (several delims appear together). Ignore this char.
            if (strnum.empty())
                continue;

            // If we reach here, we got a number. Convert it to double.
            number = stod(strnum);
            topo.back().push_back(number);

            strnum.clear();
        }
    }

    nx = static_cast<int>(topo.size());                // Assign dimension variables
    ny = static_cast<int>(topo[0].size());
    elems = nx * ny;

}

void Topo::setupgridneighbors()
{
    int i, j;

    idown = vector <double>(nx);
    iup = vector <double>(nx);
    jup = vector <double>(ny);
    jdown = vector <double>(ny);

    for (i = 0; i <= nx; i++)
    {
        idown[i] = i - 1;
        iup[i] = i + 1;
    }

    idown[0] = 0;
    iup[nx] = nx;

    for (j = 0; j <= ny; j++)
    {
        jdown[j] = j - 1;
        jup[j] = j + 1;
    }
    jdown[0] = 0;
    jup[ny] = ny;
}

void Topo::triDag()  //ax,bx,cx,rx,ux,lattice_size_x
{
    int j;
    int n = static_cast<int>(topo.size());
    double bet;
    vector <double> gam(n);

    if (bx[0] == 0.0)
        throw("Error 1 in tridag");

    ux[0] = rx[0] / ( bet = bx[0] );

    for (j = 1; j < n; j++)
    {
        gam[j] = cx[j-1] / bet;
        bet = bx[j] - ax[j] * gam[j];
        if (bet == 0.0)
            throw("Error 2 in tridag");

        ux[j]= (rx[j] - ax[j] * ux[j-1] ) / bet;
    }

    for ( j = (n-2); j >= 0 ; j-- )
        ux[j] -= gam[j+1] * ux[j+1];
}

void Topo::indexX()
{

	cout << "Topo::indexX not implemented" << endl;
	exit(1);
}

void Topo::indexX(int n, vector<double> &arr, vector<double> &indx)
{


#undef M
#undef NSTACK
#undef SWAP
}














