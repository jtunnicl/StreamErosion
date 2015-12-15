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
			float farr[5] = {1, 2, 3, 4, 5};
			std::vector<float> fv = ArrayToVector(farr, size);
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
			VectorToArray(fv, farr);
			Assert::AreEqual(fv[0], farr[0]);
			Assert::AreEqual(fv[4], farr[4]);
		}

		[TestMethod]
		void TestRandomFieldPitFilling()
		{
			StreamPower sp = StreamPower(100, 100);
			sp.Init();
			std::vector<std::vector<float>> rtopo = sp.CreateRandomField();
			sp.SetTopo(rtopo);
			sp.Flood();
			Assert::IsTrue(true);

		}


	};
}
