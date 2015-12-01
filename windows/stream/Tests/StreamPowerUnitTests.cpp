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
		void TestVector()
		{
			int nl = 1;
			int nh = 10;
			float* v_old = StreamPower::vector(nl, nh);
			std::vector<float> v_new = StreamPower::Vector(nl, nh);
			
			// check first and last index
			Assert::AreEqual(v_old[nl], v_new[nl]);
			Assert::AreEqual(v_new[nh], v_new[nh]);;
		}

		[TestMethod]
		void TestIVector()
		{
			int nl = 1;
			int nh = 10;
			int* v_old = StreamPower::ivector(nl, nh);
			std::vector<int> v_new = StreamPower::IVector(nl, nh);

			// check first and last index
			Assert::AreEqual(v_old[nl], v_new[nl]);
			Assert::AreEqual(v_old[nh], v_new[nh]);;
		}

		[TestMethod]
		void TestMatrix()
		{
			// 10 x 10 mati_rx
			int nrl = 1;
			int nrh = 10;
			int ncl = 1;
			int nch = 10;

			float** m_old = StreamPower::matrix(nrl, nrh, ncl, nch);
			std::vector<std::vector<float>> m_new = StreamPower::Matrix(nrl, nrh, ncl, nch);

			// check first and last indicies
			Assert::AreEqual(m_old[nrl][ncl], m_new[nrl][ncl]);		// top left
			Assert::AreEqual(m_old[nrl][nch], m_new[nrl][nch]);		// top right
			Assert::AreEqual(m_old[nrh][ncl], m_new[nrh][ncl]);		// bottom left
			Assert::AreEqual(m_old[nrh][nch], m_new[nrh][nch]);		// bottom right

		}

		[TestMethod]
		void TestIMatrix()
		{
			// 10 x 10 mati_rx
			int nrl = 1;
			int nrh = 10;
			int ncl = 1;
			int nch = 10;

			int** m_old = StreamPower::imatrix(nrl, nrh, ncl, nch);
			std::vector<std::vector<int>> m_new= StreamPower::IMatrix(nrl, nrh, ncl, nch);

			// check first and last indicies
			Assert::AreEqual(m_old[nrl][ncl], m_new[nrl][ncl]);		// top left
			Assert::AreEqual(m_old[nrl][nch], m_new[nrl][nch]);		// top right
			Assert::AreEqual(m_old[nrh][ncl], m_new[nrh][ncl]);		// bottom left
			Assert::AreEqual(m_old[nrh][nch], m_new[nrh][nch]);		// bottom right
			
		}

		[TestMethod]
		void TestRan3()
		{
			unsigned int trials = 1e6;
			float tolerance = 1e-3;
			std::default_random_engine generator;
			std::uniform_real_distribution<float> distribution(0.0f, 1.0f);

			std::vector<float> v_old = std::vector<float>(trials);
			std::vector<float> v_new = std::vector<float>(trials);

			int idum;
			for (int i = 0; i < trials; i++)
			{
				v_old[i] = StreamPower::ran3(&idum);
				v_new[i] = StreamPower::Ran3(generator, distribution);
			}

			Assert::AreEqual(trials, v_old.size());
			Assert::AreEqual(trials, v_new.size());

			// mean
			float mean_old = std::accumulate(v_old.begin(), v_old.end(), 0.0f) / v_old.size();
			float mean_new = std::accumulate(v_new.begin(), v_new.end(), 0.0f) / v_new.size();

			// variance
			std::vector<float> diff_old = std::vector<float>(trials);
			std::vector<float> diff_new = std::vector<float>(trials);
			float d_old, d_new;
			for (int i = 0; i < trials; i++)
			{
				d_old = v_old[i] - mean_old;
				d_new = v_new[i] - mean_new;
				diff_old[i] = d_old * d_old;
				diff_new[i] = d_new * d_new;
			}
			float stdev_old = sqrt(	std::accumulate(diff_old.begin(), diff_old.end(), 0.0f) / diff_old.size() );
			float stdev_new = sqrt(	std::accumulate(diff_new.begin(), diff_new.end(), 0.0f) / diff_new.size() );

			Assert::AreEqual(mean_old, mean_new, tolerance);
			Assert::AreEqual(stdev_old, stdev_new, tolerance);
		}

		[TestMethod]
		void TestGasdev()
		{
			unsigned int trials = 1e6;
			float tolerance = 1e-2;
			std::default_random_engine generator;
			std::normal_distribution<float> distribution(0.0f, 1.0f);

			std::vector<float> v_old = std::vector<float>(trials);
			std::vector<float> v_new = std::vector<float>(trials);

			int idum;
			for (int i = 0; i < trials; i++)
			{
				v_old[i] = StreamPower::gasdev(&idum);
				v_new[i] = StreamPower::Gasdev(generator, distribution);
			}

			Assert::AreEqual(trials, v_old.size());
			Assert::AreEqual(trials, v_new.size());

			// mean
			float mean_old = std::accumulate(v_old.begin(), v_old.end(), 0.0f) / v_old.size();
			float mean_new = std::accumulate(v_new.begin(), v_new.end(), 0.0f) / v_new.size();

			// variance
			std::vector<float> diff_old = std::vector<float>(trials);
			std::vector<float> diff_new = std::vector<float>(trials);
			float d_old, d_new;
			for (int i = 0; i < trials; i++)
			{
				d_old = v_old[i] - mean_old;
				d_new = v_new[i] - mean_new;
				diff_old[i] = d_old * d_old;
				diff_new[i] = d_new * d_new;
			}
			float stdev_old = sqrt(std::accumulate(diff_old.begin(), diff_old.end(), 0.0f) / diff_old.size());
			float stdev_new = sqrt(std::accumulate(diff_new.begin(), diff_new.end(), 0.0f) / diff_new.size());

			Assert::AreEqual(mean_old, mean_new, tolerance);
			Assert::AreEqual(stdev_old, stdev_new, tolerance);
		}

		[TestMethod]
		void TestArr_ayToVector()
		{
			int size = 5;
			float farr[5] = {1, 2, 3, 4, 5};
			std::vector<float> fv = Arr_ayToVector(farr, size);
			Assert::AreEqual(farr[0], fv[0]);
			Assert::AreEqual(farr[4], fv[4]);
		}

		[TestMethod]
		void TestVectorToArray()
		{
			int size = 5;
			float farr[5] = {0, 0, 0, 0, 0};
			std::vector<float> fv = std::vector<float>();
			for (int i = 0; i < size; i++)
			{
				fv.push_back(i);
			}
			VectorToArr_ay(fv, farr);
			Assert::AreEqual(fv[0], farr[0]);
			Assert::AreEqual(fv[4], farr[4]);
		}

		[TestMethod]
		void TestIndexx()
		{
			int nl = 1;
			int nh = 5;
			int size = nh - nl + 1 + NR_END;
			float* farr = StreamPower::vector(nl, nh);
			for (int i = 1; i <= nh; i++)
			{
				farr[i] = i;
			}

			int* indx = StreamPower::ivector(nl , nh);
			StreamPower::indexx(nh, farr, indx);

			std::vector<float> fv = Arr_ayToVector(farr, size);
			std::vector<int> indxv = StreamPower::Indexx(fv);
			for (int i = 1; i <= nh; i++)
			{
				Assert::AreEqual(indx[i], indxv[i]);
			}
		}

		[TestMethod]
		void TestTridag()
		{
			int size = 5;
			float a[6] = {0.0f, 0.1f, 0.3f, -0.2f, 1.5f, 0.0f};
			float b[6] = {0.0f, 0.2f, 0.4f, -0.3f, 2.5f, 0.0f};
			float c[6] = {0.0f, 0.3f, 0.5f, -0.4f, 3.5f, 0.0f};
			float r[6] = {0.0f, 0.4f, 0.6f, -0.5f, 4.5f, 0.0f};
			float u[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
			float uv[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

			StreamPower::tridag(a, b, c, r, u, size);
			StreamPower::Tridag(a, b, c, r, uv, size);

			for (int i = 1; i <= size; i++)
			{
				Assert::AreEqual(u[i], uv[i]);
			}
			
		}

		[TestMethod]
		void TestSetupGridNeighbours()
		{
			int nx = 10;
			int ny = 10;
			
			StreamPower sp_old = StreamPower(nx, ny);
			sp_old.setupgridneighbors();

			StreamPower sp_new = StreamPower(nx, ny);
			sp_new.SetupGridNeighbors();

			
			for (int i = 1; i <= nx; i++)
			{
				Assert::AreEqual(sp_old._iup[i], sp_new.iup[i]);
				Assert::AreEqual(sp_old._idown[i], sp_new.idown[i]);
			}

			for (int j = 1; j <= ny; j++)
			{
				Assert::AreEqual(sp_old._jup[j], sp_new.jup[j]);
				Assert::AreEqual(sp_old._jdown[j], sp_new.jdown[j]);
			}
			
		}
	};
}
