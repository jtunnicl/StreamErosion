#include "stdafx.h"
#include "streampower.h"
#include <numeric>
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
		void ShouldFail()
		{
			Assert::AreEqual(3, StreamPower::Count());
		}

		[TestMethod]
		void ShouldPass()
		{
			Assert::AreEqual(2, StreamPower::Count());
		}

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
			// 10 x 10 matirx
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
			// 10 x 10 matirx
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
	};
}
