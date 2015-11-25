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
		void TestFloatVectorCreation()
		{
			int nl = 1;
			int nh = 10;
			float* v_old = StreamPower::vector(nl, nh);
			std::vector<float> v_new = StreamPower::Vector(nl, nh);
			
			// check first and last index
			Assert::AreEqual(v_old[nl], v_new[nl]);
			Assert::AreEqual(v_new[nh], v_new[nh]);;
		}
	};
}
