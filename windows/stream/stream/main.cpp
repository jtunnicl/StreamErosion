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
	p.U = 1;
	p.K = 0.05;

	StreamPower sp = StreamPower(p);
	sp.SetTopo(sp.ReadArcInfoASCIIGrid(fname));
	sp.Start();





}