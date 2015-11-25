#include "stdafx.h"
#include "streampower.h"

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
			Assert::AreEqual(m_old[nrl][ncl], m_new[nrl][ncl]);
			Assert::AreEqual(m_old[nrl][nch], m_new[nrl][nch]);
			Assert::AreEqual(m_old[nrh][ncl], m_new[nrh][ncl]);
			Assert::AreEqual(m_old[nrh][nch], m_new[nrh][nch]);
			


		}
	};
}
