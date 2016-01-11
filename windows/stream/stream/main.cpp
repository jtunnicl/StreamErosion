#include "streampower.h"

int main(int argc, char** argv)
{
	int nx = 100;
	int ny = 100;
	int duration = 1;

	char* fname = argv[1];

	Parameters p;
	p.timestep = 1;
	p.duration = 1;
	p.K = 0.05;

	StreamPower sp = StreamPower(p);
	sp.SetTopo(sp.ReadArcInfoASCIIGrid(fname));
	sp.SetU(1.0f);
	// shoulde also work as: sp.SetU(sp.ReadArcInfoASCIIGrid(ufname));  // WARNING: should have the same metadata as the topo DEM e.g. deltax, nrows, ncols.

	sp.Start();





}