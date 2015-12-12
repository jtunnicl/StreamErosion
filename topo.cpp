#include "topo.h"
#include "util.h"
#include <iostream>
#include <cstdlib>

Topo::Topo() : nx(0), ny(0), elems(0), deltax(0)
{

}

void Topo::Init()
{

	// setup required data structures
	// strat = todo
	topo2 = vector<vector<double>>(ny, vector<double>(nx));
	topoold = vector<vector<double>>(ny, vector<double>(nx));
	slope = vector<vector<double>>(ny, vector<double>(nx));
	flow1 = vector<vector<double>>(ny, vector<double>(nx));
	flow2 = vector<vector<double>>(ny, vector<double>(nx));
	flow3 = vector<vector<double>>(ny, vector<double>(nx));
	flow4 = vector<vector<double>>(ny, vector<double>(nx));
	flow5 = vector<vector<double>>(ny, vector<double>(nx));
	flow6 = vector<vector<double>>(ny, vector<double>(nx));
	flow7 = vector<vector<double>>(ny, vector<double>(nx));
	flow8 = vector<vector<double>>(ny, vector<double>(nx));
	topoVec = vector<vector<double>>(ny, vector<double>(nx));
	vegState = vector<vector<double>>(ny, vector<double>(nx));
	exposureAge = vector<vector<double>>(ny, vector<double>(nx));

	topoVecInd = vector<int>(1, elems);

	// Vectors to represent neighbour cells
	vector<int> idown (nx + 1, 0);
    vector<int> iup (nx + 1, 0);
    vector<int> jdown (ny + 1, 0);
    vector<int> jup (ny + 1 ,0);

	SetupGridNeighbors();

}

Topo::Topo(const char* dem)
{
	// try to open file
	ifstream in(dem);
	if (in)
	{
		LoadDEM(in, " ");
		Init();
	}
	else
	{
		string msg = "Unable to open file: '" + string(dem) + "'";
		Util::Error(msg, 1);
	}

}

Topo::~Topo()
{

}

void Topo::PrintMatrix(const vector<vector<double>>& mat, int nrows, int ncols)
{
	for (int r = 0; r < nrows; r++)
	{
		for (int c = 0; c < ncols; c++)
		{
			cout << mat[r][c] << " ";
		} 
		cout << endl;
	}
}

void Topo::PrintMatrixRow(const vector<vector<double>>& mat, int row, int ncols)
{
	for (int c = 0; c < ncols; c++)
	{
		cout << mat[row][c] << " ";
	} 
	cout << endl;
}

void Topo::LoadDEM(ifstream& in, const string& delim)
{
	Util::Warning("Reading DEM without any checks or guarantees ...");
	string line;

	// clear topo
	topo.clear();

	// read 6 lines of metadata
	string key;
	in >> key; in >> nx; // ncols
	in >> key; in >> ny; // nrows
	in >> key; in >> xllcorner;
	in >> key; in >> yllcorner;
	in >> key; in >> deltax;
	in >> key; in >> nodata;

	// init topo
	topo = vector<vector<double>>(ny, vector<double>(nx));
	elems = nx * ny;

	// read data
	for (int y = 0; y < ny; y++)
	{
		for (int x = 0; x < nx; x++)
		{
			in >> topo[y][x];
		}
	}
	//PrintMatrx(topo, ny, nx);
	PrintMatrixRow(topo, ny-1, nx);
}

void Topo::SetupGridNeighbors()
{

    for (int i = 0; i <= nx; i++)
    {
        idown[i] = i - 1;
        iup[i] = i + 1;
    }
    idown[0] = 0;
    iup[nx] = nx;

    for (int j = 0; j <= ny; j++)
    {
        jdown[j] = j - 1;
        jup[j] = j + 1;
    }
    jdown[0] = 0;
    jup[ny] = ny;
}

/*
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
*/

void Topo::triDag()
{
	Util::Error("Topo::triDag not implemented", 1);
}
void Topo::IndexX()
{
	//Util::Error("Topo::indexX not implemented", 1);
	// sort topo and place sorted indices in topoVecInd

}


/*
DEPRECATED
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
*/












