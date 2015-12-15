#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "streampower.h"
#include "utility.h"
#include "priority_flood.hpp"



/* allocate a float vector with subscript range v[nl..nh] */
std::vector<float> StreamPower::Vector(int nl, int nh)
{
	int size = nh - nl + 1 + NR_END;
	return std::vector<float>(size);
}

/* allocate a float vector with subscript range v[nl..nh] */
std::vector<int> StreamPower::IVector(int nl, int nh)
{
	int size = nh - nl + 1 + NR_END;
	return std::vector<int>(size);
}

std::vector<std::vector<int>> StreamPower::IMatrix(int nrl, int nrh, int ncl, int nch)
{
	int rsize = nrh - nrl + 1 + NR_END;
	int csize = nch - ncl + 1 + NR_END;
	return std::vector<std::vector<int>>(rsize, std::vector<int>(csize));
}

std::vector<std::vector<float>> StreamPower::Matrix(int nrl, int nrh, int ncl, int nch)
{
	int rsize = nrh - nrl + 1 + NR_END;
	int csize = nch - ncl + 1 + NR_END;
	return std::vector<std::vector<float>>(rsize, std::vector<float>(csize));
}

float StreamPower::Ran3(std::default_random_engine& generator, std::uniform_real_distribution<float>& distribution)
{
	return distribution(generator);
}

float StreamPower::Gasdev(std::default_random_engine& generator, std::normal_distribution<float>& distribution)
{
	/*
		Assuming this is the same code from here: http://www.stat.berkeley.edu/~paciorek/diss/code/regression.binomial/gasdev.C
		We need to return a standard, normally distributed gaussian random number
	*/

	return distribution(generator);

}

void StreamPower::Indexx(int n, float* arr, int* indx)
{
	std::vector<float> arrv = ArrayToVector(arr, n);
	std::vector<int> indxv = Indexx(arrv);
	VectorToArray(indxv, indx);
}

std::vector<int> StreamPower::Indexx(std::vector<float>& arr)
{
	return SortFortranIndices(arr);
}

void StreamPower::Tridag(float a[], float b[], float c[], float r[], float u[], unsigned long n)
{
	std::vector<float> av = ArrayToVector(a, n, true);
	std::vector<float> bv = ArrayToVector(b, n, true);
	std::vector<float> cv = ArrayToVector(c, n, true);
	std::vector<float> rv = ArrayToVector(r, n, true);
	std::vector<float> uv = ArrayToVector(u, n, true);

	Tridag(av, bv, cv, rv, uv, n);

	VectorToArray(av, a);
	VectorToArray(bv, b);
	VectorToArray(cv, c);
	VectorToArray(rv, r);
	VectorToArray(uv, u);
}

void StreamPower::Tridag(std::vector<float>& a, std::vector<float>& b, std::vector<float>& c, std::vector<float>& r, std::vector<float>& u, int n)
{
	unsigned long j;
	float bet;
	std::vector<float> gam;

	gam = Vector(1, n);
	u[1] = r[1] / (bet = b[1]);
	for (j = 2; j <= n; j++)
	{
		gam[j] = c[j - 1] / bet;
		bet = b[j] - a[j] * gam[j];
		u[j] = (r[j] - a[j] * u[j - 1]) / bet;
	}
	for (j = (n - 1); j >= 1; j--)
	{
		u[j] -= gam[j + 1] * u[j + 1];
	}
}

void StreamPower::SetupGridNeighbors()
{
	int i, j;

	idown = IVector(1, lattice_size_x);
	iup = IVector(1, lattice_size_x);
	jup = IVector(1, lattice_size_y);
	jdown = IVector(1, lattice_size_y);

	for (i = 1; i <= lattice_size_x; i++)
	{
		idown[i] = i - 1;
		iup[i] = i + 1;
	}
	idown[1] = 1;
	iup[lattice_size_x] = lattice_size_x;
	for (j = 1; j <= lattice_size_y; j++)
	{
		jdown[j] = j - 1;
		jup[j] = j + 1;
	}
	jdown[1] = 1;
	jup[lattice_size_y] = lattice_size_y;
}

void StreamPower::Flood()
{
	// update elev
	for (int i = 1; i <= lattice_size_x; i++)
	{
		for (int j = 1; j <= lattice_size_y; j++)
		{
			elevation(i - 1, j - 1) = topo[i][j];
		}
	}

	// perform flooding
	original_priority_flood(elevation);

	// update topo
	for (int i = 0; i < lattice_size_x; i++)
	{
		for (int j = 0; j < lattice_size_y; j++)
		{
			topo[i + 1][j + 1] = elevation(i, j);
		}
	}

}

void StreamPower::MFDFlowRoute(int i, int j)
{
	float tot;

	tot = 0;
	if (topo[i][j] > topo[iup[i]][j])
		tot += pow(topo[i][j] - topo[iup[i]][j], 1.1f);
	if (topo[i][j] > topo[idown[i]][j])
		tot += pow(topo[i][j] - topo[idown[i]][j], 1.1f);
	if (topo[i][j] > topo[i][jup[j]])
		tot += pow(topo[i][j] - topo[i][jup[j]], 1.1f);
	if (topo[i][j] > topo[i][jdown[j]])
		tot += pow(topo[i][j] - topo[i][jdown[j]], 1.1f);
	if (topo[i][j] > topo[iup[i]][jup[j]])
		tot += pow((topo[i][j] - topo[iup[i]][jup[j]])*oneoversqrt2, 1.1f);
	if (topo[i][j] > topo[iup[i]][jdown[j]])
		tot += pow((topo[i][j] - topo[iup[i]][jdown[j]])*oneoversqrt2, 1.1f);
	if (topo[i][j] > topo[idown[i]][jup[j]])
		tot += pow((topo[i][j] - topo[idown[i]][jup[j]])*oneoversqrt2, 1.1f);
	if (topo[i][j] > topo[idown[i]][jdown[j]])
		tot += pow((topo[i][j] - topo[idown[i]][jdown[j]])*oneoversqrt2, 1.1f);
	if (topo[i][j] > topo[iup[i]][j])
		flow1[i][j] = pow(topo[i][j] - topo[iup[i]][j], 1.1f) / tot;
	else flow1[i][j] = 0;
	if (topo[i][j] > topo[idown[i]][j])
		flow2[i][j] = pow(topo[i][j] - topo[idown[i]][j], 1.1f) / tot;
	else flow2[i][j] = 0;
	if (topo[i][j] > topo[i][jup[j]])
		flow3[i][j] = pow(topo[i][j] - topo[i][jup[j]], 1.1f) / tot;
	else flow3[i][j] = 0;
	if (topo[i][j] > topo[i][jdown[j]])
		flow4[i][j] = pow(topo[i][j] - topo[i][jdown[j]], 1.1f) / tot;
	else flow4[i][j] = 0;
	if (topo[i][j] > topo[iup[i]][jup[j]])
		flow5[i][j] = pow((topo[i][j] - topo[iup[i]][jup[j]])*oneoversqrt2, 1.1f) / tot;
	else flow5[i][j] = 0;
	if (topo[i][j] > topo[iup[i]][jdown[j]])
		flow6[i][j] = pow((topo[i][j] - topo[iup[i]][jdown[j]])*oneoversqrt2, 1.1f) / tot;
	else flow6[i][j] = 0;
	if (topo[i][j] > topo[idown[i]][jup[j]])
		flow7[i][j] = pow((topo[i][j] - topo[idown[i]][jup[j]])*oneoversqrt2, 1.1f) / tot;
	else flow7[i][j] = 0;
	if (topo[i][j] > topo[idown[i]][jdown[j]])
		flow8[i][j] = pow((topo[i][j] - topo[idown[i]][jdown[j]])*oneoversqrt2, 1.1f) / tot;
	else flow8[i][j] = 0;
	flow[iup[i]][j] += flow[i][j] * flow1[i][j];
	flow[idown[i]][j] += flow[i][j] * flow2[i][j];
	flow[i][jup[j]] += flow[i][j] * flow3[i][j];
	flow[i][jdown[j]] += flow[i][j] * flow4[i][j];
	flow[iup[i]][jup[j]] += flow[i][j] * flow5[i][j];
	flow[iup[i]][jdown[j]] += flow[i][j] * flow6[i][j];
	flow[idown[i]][jup[j]] += flow[i][j] * flow7[i][j];
	flow[idown[i]][jdown[j]] += flow[i][j] * flow8[i][j];
}

void StreamPower::CalculateAlongChannelSlope(int i, int j)
{
	float down;

	down = 0;
	if (topo[iup[i]][j] - topo[i][j] < down) down = topo[iup[i]][j] - topo[i][j];
	if (topo[idown[i]][j] - topo[i][j] < down) down = topo[idown[i]][j] - topo[i][j];
	if (topo[i][jup[j]] - topo[i][j] < down) down = topo[i][jup[j]] - topo[i][j];
	if (topo[i][jdown[j]] - topo[i][j] < down) down = topo[i][jdown[j]] - topo[i][j];
	if ((topo[iup[i]][jup[j]] - topo[i][j])*oneoversqrt2 < down)
		down = (topo[iup[i]][jup[j]] - topo[i][j])*oneoversqrt2;
	if ((topo[idown[i]][jup[j]] - topo[i][j])*oneoversqrt2 < down)
		down = (topo[idown[i]][jup[j]] - topo[i][j])*oneoversqrt2;
	if ((topo[iup[i]][jdown[j]] - topo[i][j])*oneoversqrt2 < down)
		down = (topo[iup[i]][jdown[j]] - topo[i][j])*oneoversqrt2;
	if ((topo[idown[i]][jdown[j]] - topo[i][j])*oneoversqrt2 < down)
		down = (topo[idown[i]][jdown[j]] - topo[i][j])*oneoversqrt2;
	slope[i][j] = fabs(down) / deltax;
}

void StreamPower::HillSlopeDiffusionInit()
{
	int i, j, count;
	float term1;

	ax = Vector(1, lattice_size_x);
	ay = Vector(1, lattice_size_y);
	bx = Vector(1, lattice_size_x);
	by = Vector(1, lattice_size_y);
	cx = Vector(1, lattice_size_x);
	cy = Vector(1, lattice_size_y);
	ux = Vector(1, lattice_size_x);
	uy = Vector(1, lattice_size_y);
	rx = Vector(1, lattice_size_x);
	ry = Vector(1, lattice_size_y);
	D = 10000000.0;
	count = 0;
	term1 = D / (deltax*deltax);
	for (i = 1; i <= lattice_size_x; i++)
	{
		ax[i] = -term1;
		cx[i] = -term1;
		bx[i] = 4 * term1 + 1;
		if (i == 1)
		{
			bx[i] = 1;
			cx[i] = 0;
		}
		if (i == lattice_size_x)
		{
			bx[i] = 1;
			ax[i] = 0;
		}
	}
	for (j = 1; j <= lattice_size_y; j++)
	{
		ay[j] = -term1;
		cy[j] = -term1;
		by[j] = 4 * term1 + 1;
		if (j == 1)
		{
			by[j] = 1;
			cy[j] = 0;
		}
		if (j == lattice_size_y)
		{
			by[j] = 1;
			ay[j] = 0;
		}
	}
	while (count < 5)
	{
		count++;
		for (i = 1; i <= lattice_size_x; i++)
		{
			for (j = 1; j <= lattice_size_y; j++)
			{
				topo2[i][j] = topo[i][j];
			}
		}
		for (i = 1; i <= lattice_size_x; i++)
		{
			for (j = 1; j <= lattice_size_y; j++)
			{
				ry[j] = term1*(topo[iup[i]][j] + topo[idown[i]][j]) + topoold[i][j];
				if (j == 1)
					ry[j] = topoold[i][j];
				if (j == lattice_size_y)
					ry[j] = topoold[i][j];
			}
			Tridag(ay, by, cy, ry, uy, lattice_size_y);
			for (j = 1; j <= lattice_size_y; j++)
			{
				topo[i][j] = uy[j];
			}
		}
		for (i = 1; i <= lattice_size_x; i++)
		{
			for (j = 1; j <= lattice_size_y; j++)
			{
				topo2[i][j] = topo[i][j];
			}
		}
		for (j = 1; j <= lattice_size_y; j++)
		{
			for (i = 1; i <= lattice_size_x; i++)
			{
				rx[i] = term1*(topo[i][jup[j]] + topo[i][jdown[j]]) + topoold[i][j];
				if (i == 1)
				{
					rx[i] = topoold[i][j];
				}
				if (i == lattice_size_x)
				{
					rx[i] = topoold[i][j];
				}
			}
			Tridag(ax, bx, cx, rx, ux, lattice_size_x);
			for (i = 1; i <= lattice_size_x; i++)
			{
				topo[i][j] = ux[i];
			}
		}
	}
}

void StreamPower::Avalanche(int i, int j)
{
	if (topo[iup[i]][j] - topo[i][j] > thresh)
		topo[iup[i]][j] = topo[i][j] + thresh;
	if (topo[idown[i]][j] - topo[i][j] > thresh)
		topo[idown[i]][j] = topo[i][j] + thresh;
	if (topo[i][jup[j]] - topo[i][j] > thresh)
		topo[i][jup[j]] = topo[i][j] + thresh;
	if (topo[i][jdown[j]] - topo[i][j] > thresh)
		topo[i][jdown[j]] = topo[i][j] + thresh;
	if (topo[iup[i]][jup[j]] - topo[i][j] > (thresh*sqrt2))
		topo[iup[i]][jup[j]] = topo[i][j] + thresh*sqrt2;
	if (topo[iup[i]][jdown[j]] - topo[i][j] > (thresh*sqrt2))
		topo[iup[i]][jdown[j]] = topo[i][j] + thresh*sqrt2;
	if (topo[idown[i]][jup[j]] - topo[i][j] > (thresh*sqrt2))
		topo[idown[i]][jup[j]] = topo[i][j] + thresh*sqrt2;
	if (topo[idown[i]][jdown[j]] - topo[i][j] > (thresh*sqrt2))
		topo[idown[i]][jdown[j]] = topo[i][j] + thresh*sqrt2;
}

void StreamPower::Start()
{
	float deltah, time, max;
	int idum, i, j, t, step;
	time = 0;

	char fname[100];
	sprintf(fname, "erosion_%d.txt", 0);
	PrintState(fname);

	while (time < duration)
	{
		//perform landsliding
		for (j = 1; j <= lattice_size_y; j++)
		{
			for (i = 1; i <= lattice_size_x; i++)
			{
				topovec[(j - 1)*lattice_size_x + i] = topo[i][j];
			}
		}
		//indexx(lattice_size_x*lattice_size_y, _topovec, _topovecind);
		topovecind = Indexx(topovec);

		t = 0;
		while (t < lattice_size_x*lattice_size_y)
		{
			t++;
			i = (topovecind[t]) % lattice_size_x;
			if (i == 0)
			{
				i = lattice_size_x;
			}
			j = (topovecind[t]) / lattice_size_x + 1;
			if (i == lattice_size_x)
			{
				j--;
			}
			//avalanche(i, j);
			Avalanche(i, j);
		}
		for (j = 1; j <= lattice_size_y; j++)
		{
			for (i = 1; i <= lattice_size_x; i++)
			{
				topoold[i][j] = topo[i][j];
			}
		}

		// pit filling
		/*
		for (j = 1; j <= lattice_size_y; j++)
			for (i = 1; i <= lattice_size_x; i++)
				fillinpitsandflats(i, j);
		*/
		Flood();

		for (j = 1; j <= lattice_size_y; j++)
		{
			for (i = 1; i <= lattice_size_x; i++)
			{
				flow[i][j] = 1;
				topovec[(j - 1)*lattice_size_x + i] = topo[i][j];
			}
		}
		//indexx(lattice_size_x*lattice_size_y, _topovec, _topovecind);
		topovecind = Indexx(topovec);

		t = lattice_size_x*lattice_size_y + 1;
		while (t > 1)
		{
			t--;
			i = (topovecind[t]) % lattice_size_x;
			if (i == 0)
			{
				i = lattice_size_x;
			}
			j = (topovecind[t]) / lattice_size_x + 1;
			if (i == lattice_size_x)
			{
				j--;
			}
			//mfd_flowroute(i, j);
			MFDFlowRoute(i, j);
		}

		// perform uplift
		for (i = 2; i <= lattice_size_x - 1; i++)
		{
			for (j = 2; j <= lattice_size_y - 1; j++)
			{
				topo[i][j] += U*timestep;
				topoold[i][j] += U*timestep;
			}
		}

		//perform upwind erosion
		max = 0;
		for (i = 2; i <= lattice_size_x - 1; i++)
		{
			for (j = 2; j <= lattice_size_y - 1; j++)
			{
				//calculatealongchannel_slope(i, j);
				CalculateAlongChannelSlope(i, j);
				deltah = timestep*K*sqrt(flow[i][j])*deltax*slope[i][j];
				topo[i][j] -= deltah;
				if (topo[i][j]<0)
				{
					topo[i][j] = 0;
				}
				if (K*sqrt(flow[i][j])*deltax>max)
				{
					max = K*sqrt(flow[i][j])*deltax;
				}
			}
		}
		time += timestep;
		if (max > 0.3*deltax / timestep)
		{
			time -= timestep;
			timestep /= 2.0;
			for (i = 2; i <= lattice_size_x - 1; i++)
			{
				for (j = 2; j <= lattice_size_y - 1; j++)
				{
					topo[i][j] = topoold[i][j] - U*timestep;
				}
			}
		}
		else
		{
			if (max < 0.03*deltax / timestep)
			{
				timestep *= 1.2;
			}
			for (j = 1; j <= lattice_size_y; j++)
			{
				for (i = 1; i <= lattice_size_x; i++)
				{
					topoold[i][j] = topo[i][j];
				}

			}

		}
		//if (time > printinterval)
		//{
		char fname[100];
		//sprintf(fname, "erosion_%d.txt", printinterval);
		sprintf(fname, "erosion_%f.txt", time);
		PrintState(fname);
		//printinterval += printstep;
	//}
		std::cout << "Time: " << time << std::endl;

	}

}

void StreamPower::PrintState(char* fname)
{
	std::ofstream file;
	file.open(fname);
	// write arcgrid format
	file << "ncols " << lattice_size_y << std::endl;
	file << "nrows " << lattice_size_x << std::endl;
	file << "xllcorner " << xllcorner << std::endl;
	file << "yllcorner " << yllcorner << std::endl;
	file << "cellsize " << deltax << std::endl;
	file << "NODATA_value " << nodata << std::endl;
	for (int i = 1; i <= lattice_size_x; i++)
	{
		for (int j = 1; j <= lattice_size_y; j++)
		{
			file << topo[i][j] << " ";
		}
		file << std::endl;
	}
	file.close();
}

std::vector<std::vector<float>> StreamPower::CreateRandomField()
{
	std::vector<std::vector<float>> mat = Matrix(1, lattice_size_x, 1, lattice_size_y);
	std::default_random_engine generator;
	std::normal_distribution<float> distribution(0.0f, 1.0f);
	for (int i = 1; i <= lattice_size_x; i++)
	{
		for (int j = 1; j <= lattice_size_y; j++)
		{
			mat[i][j] = 0.5*Gasdev(generator, distribution);
		}
	}
	return mat;
}

void StreamPower::SetTopo(std::vector<std::vector<float>> t)
{

	topo = Matrix(1, lattice_size_x, 1, lattice_size_y);
	topo2 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	topoold = Matrix(1, lattice_size_x, 1, lattice_size_y);
	slope = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow1 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow2 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow3 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow4 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow5 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow6 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow7 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	flow8 = Matrix(1, lattice_size_x, 1, lattice_size_y);
	topovec = Vector(1, lattice_size_x*lattice_size_y);
	topovecind = IVector(1, lattice_size_x*lattice_size_y);
	elevation = Array2D<float>(lattice_size_x, lattice_size_y, -9999.0f);

	SetupGridNeighbors();

	for (int i = 1; i <= lattice_size_x; i++)
	{
		for (int j = 1; j <= lattice_size_y; j++)
		{
			topo[i][j] = t[i][j];
			topoold[i][j] = topo[i][j];
			flow[i][j] = 1;
		}
	}
	InitDiffusion();
}

std::vector<std::vector<float>> StreamPower::ReadArcInfoASCIIGrid(char* fname)
{
	std::ifstream in(fname);
	std::vector<std::vector<float>> t;
	std::string line;

	Util::Warning("Reading DEM without any checks or guarantees ...");

	// read 6 lines of metadata
	std::string key;
	in >> key; in >> lattice_size_y; // ncols //NOTE: Pelltier's code was originally written for [x][y] indexing; Saga uses [y][x].
	in >> key; in >> lattice_size_x; // nrows
	in >> key; in >> xllcorner;
	in >> key; in >> yllcorner;
	in >> key; in >> deltax;
	in >> key; in >> nodata;

	t = Matrix(1, lattice_size_x, 1, lattice_size_y);

	// read data
	for (int x = 1; x <= lattice_size_x; x++)
	{
		for (int y = 1; y <= lattice_size_y; y++)
		{
			in >> t[x][y];
		}
	}

	thresh = 0.58*deltax; // 30 deg // This may have to be adjusted for variable deltax (deltax was originally 200)

	Util::Info("Done reading DEM");

	return t;
}

void StreamPower::InitDiffusion()
{
	//construct diffusional landscape for initial flow routing
	for (int step = 1; step <= 10; step++)
	{
		HillSlopeDiffusionInit();
		for (int i = 2; i <= lattice_size_x - 1; i++)
		{
			for (int j = 2; j <= lattice_size_y - 1; j++)
			{
				topo[i][j] += 0.1;
				topoold[i][j] += 0.1;
			}
		}
	}
}

void StreamPower::Init()
{

	U = 1;                // m/kyr
	K = 0.05;             // kyr^-1
	printstep = 1;
	printinterval = printstep;
	deltax = 200.0;       // m
	nodata = -9999.0;
	xllcorner = 0;
	yllcorner = 0;
	thresh = 0.58*deltax; // 30 deg // This may have to be adjusted for variable deltax (deltax was originally 200)
	timestep = 1;         // kyr
	duration = 100;


}

StreamPower::StreamPower(int nx, int ny) : lattice_size_x(nx), lattice_size_y(ny)
{

}

StreamPower::~StreamPower() {}

std::vector<std::vector<float>> StreamPower::GetTopo()
{
	return topo;
}