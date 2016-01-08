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
		void TestArrayToVector()
		{
			int size = 5;
			float farr[5] = { 1, 2, 3, 4, 5 };
			std::vector<float> fv = ArrayToVector(farr, size);
			Assert::AreEqual(farr[0], fv[0]);
			Assert::AreEqual(farr[4], fv[4]);
		}

		[TestMethod]
		void TestVectorToArray()
		{
			int size = 5;
			float farr[5] = { 0, 0, 0, 0, 0 };
			std::vector<float> fv = std::vector<float>();
			for (int i = 0; i < size; i++)
			{
				fv.push_back(i);
			}
			VectorToArray(fv, farr);
			Assert::AreEqual(fv[0], farr[0]);
			Assert::AreEqual(fv[4], farr[4]);
		}

		[TestMethod]
		void TestRandomFieldPitFilling()
		{
			int duration = 1;
			StreamPower sp = StreamPower(100, 100);
			sp.Init(duration);
			std::vector<std::vector<float>> rtopo = sp.CreateRandomField();
			sp.SetTopo(rtopo);
			sp.Flood();
			Assert::IsTrue(true);

		}

		[TestMethod]
		void TestSetupGridNeighbours_C()
		{
			int nx = 10;
			int ny = 10;
			int duration = 1;

			// fortran indexing
			StreamPower sp_f = StreamPower(nx, ny);
			sp_f.Init(duration);
			std::vector<std::vector<float>> rt_f = sp_f.CreateRandomField();
			sp_f.SetTopo(rt_f);

			// create matrix using c indexing
			std::vector<std::vector<float>> rt_c(nx, std::vector<float>(ny));
			for (int i = 1; i <= nx; i++)
			{
				for (int j = 1; j <= ny; j++)
				{
					rt_c[i - 1][j - 1] = rt_f[i][j];
				}
			}

			// c indexing
			StreamPower sp_c = StreamPower(nx, ny);
			sp_c.AssignVariables_C();
			sp_c.SetupGridNeighbors_C();

			std::vector<int> iup_f, idown_f, jup_f, jdown_f;
			std::vector<int> iup_c, idown_c, jup_c, jdown_c;

			iup_f = sp_f.iup;	iup_c = sp_c.iup;
			jup_f = sp_f.jup;	jup_c = sp_c.jup;
			idown_f = sp_f.idown;	idown_c = sp_c.idown;
			jdown_f = sp_f.jdown;	jdown_c = sp_c.jdown;

			for (int i = 1; i <= nx; i++)
			{
				Assert::AreEqual(iup_f[i] - 1, iup_c[i - 1]);
				Assert::AreEqual(idown_f[i] - 1, idown_c[i - 1]);
			}

			for (int j = 1; j <= ny; j++)
			{
				Assert::AreEqual(jup_f[j] - 1, jup_c[j - 1]);
				Assert::AreEqual(jdown_f[j] - 1, jdown_c[j - 1]);
			}
		}



		[TestMethod]
		void TestSetTopo_C()
		{
			int nx = 10;
			int ny = 10;
			int duration = 1;
			float tolerance = 1e-6f;

			// fortran indexing
			StreamPower sp_f = StreamPower(nx, ny);
			sp_f.Init(duration);
			std::vector<std::vector<float>> rt_f = sp_f.CreateRandomField();
			sp_f.SetTopo(rt_f);

			// create matrix using c indexing
			std::vector<std::vector<float>> rt_c(nx, std::vector<float>(ny));
			for (int i = 1; i <= nx; i++)
			{
				for (int j = 1; j <= ny; j++)
				{
					rt_c[i - 1][j - 1] = rt_f[i][j];
				}
			}

			// c indexing
			StreamPower sp_c = StreamPower(nx, ny);
			sp_c.Init(duration);
			sp_c.SetTopo_C(rt_c);

			// check that set topo is working correctly
			std::vector<std::vector<float>> t_f = sp_f.GetTopo();
			std::vector<std::vector<float>> t_c = sp_c.GetTopo();
			for (int i = 1; i <= nx; i++)
			{
				for (int j = 1; j <= ny; j++)
				{
					Assert::AreEqual(t_f[i][j], t_c[i - 1][j - 1], tolerance);
				}
			}
		}


		[TestMethod]
		void TestCIndexing()
		{

			int nx = 10;
			int ny = 10;
			int duration = 10;
			float tolerance = 1e-6f;

			// fortran indexing
			StreamPower sp_f = StreamPower(nx, ny);
			sp_f.Init(duration);
			std::vector<std::vector<float>> rt_f = sp_f.CreateRandomField();
			sp_f.SetTopo(rt_f);

			// create matrix using c indexing
			std::vector<std::vector<float>> rt_c(nx, std::vector<float>(ny));
			for (int i = 1; i <= nx; i++)
			{
				for (int j = 1; j <= ny; j++)
				{
					rt_c[i - 1][j - 1] = rt_f[i][j];
				}
			}

			// c indexing
			StreamPower sp_c = StreamPower(nx, ny);
			sp_c.Init(duration);
			sp_c.SetTopo_C(rt_c);

			// run main loop
			sp_f.Start();
			sp_c.Start_C();
			std::vector<std::vector<float>> t_f = sp_f.GetTopo();
			std::vector<std::vector<float>> t_c = sp_c.GetTopo();

			for (int i = 1; i <= nx; i++)
			{
				for (int j = 1; j <= ny; j++)
				{
					Assert::AreEqual(t_f[i][j], t_c[i - 1][j - 1], tolerance);
				}
			}
		}

		[TestMethod]
		void TestReference()
		{

			float tolerance = 1e-3f;
			int duration = 1;

			char* finput = "test_input.asc";
			char* foutput = "test_output.asc"; // after duration of 1ky ~ 13 timesteps
			StreamPower sp = StreamPower();
			sp.Init(duration);
			sp.SetTopo(sp.ReadArcInfoASCIIGrid(finput));
			sp.Start();

			std::vector<std::vector<float>> tsim = sp.GetTopo();
			std::vector<std::vector<float>> ttest = sp.ReadArcInfoASCIIGrid(foutput);

			for (int i = 1; i <= sp.lattice_size_x; i++)
			{
				for (int j = 1; j <= sp.lattice_size_y; j++)
				{
					Assert::AreEqual(ttest[i][j], tsim[i][j], tolerance);
				}
			}

		}


		[TestMethod]
		void TestReference_C()
		{

			float tolerance = 1e-3f;
			int duration = 1;

			char* finput = "test_input.asc";
			char* foutput = "test_output.asc"; // after duration of 1ky ~ 13 timesteps
			StreamPower sp = StreamPower();
			sp.Init(duration);
			sp.SetTopo_C(sp.ReadArcInfoASCIIGrid_C(finput));
			sp.Start_C();

			std::vector<std::vector<float>> tsim = sp.GetTopo();
			std::vector<std::vector<float>> ttest = sp.ReadArcInfoASCIIGrid_C(foutput);

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