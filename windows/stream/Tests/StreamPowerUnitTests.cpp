#include "stdafx.h"
#include "streampower.h"
#include <cmath>

using namespace System;
using namespace System::Text;
using namespace System::Collections::Generic;
using namespace	Microsoft::VisualStudio::TestTools::UnitTesting;

namespace Tests
{
	[TestClass]
	public ref class StreamPowerUnitTests
	{
	public:

		[TestMethod]
		void TestReference()
		{

			double tolerance = 1e-3;

			char* finput = "test_input.asc";
			char* foutput = "test_output_double_win64.asc"; // after duration of 1ky

			Parameters p;
			p.timestep = 1;
			p.duration = 1;
			p.K = 0.05;

			StreamPower sp = StreamPower(p);
			sp.SetTopo(sp.ReadArcInfoASCIIGrid(finput));
			sp.SetU(1.0f);
			sp.Start();

			std::vector<std::vector<double>> tsim = sp.GetTopo();
			std::vector<std::vector<double>> ttest = sp.ReadArcInfoASCIIGrid(foutput);

			for (int i = 0; i < sp.lattice_size_x; i++)
			{
				for (int j = 0; j < sp.lattice_size_y; j++)
				{
					Assert::AreEqual(ttest[i][j], tsim[i][j], tolerance);
				}
			}

		}

	};
}