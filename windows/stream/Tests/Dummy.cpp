#include "stdafx.h"
#include "streampower.h"

using namespace System;
using namespace System::Text;
using namespace System::Collections::Generic;
using namespace	Microsoft::VisualStudio::TestTools::UnitTesting;

namespace Tests
{
	[TestClass]
	public ref class Dummy
	{
	public: 
		[TestMethod]
		void ShouldFail()
		{
			Assert::AreEqual(StreamPower::Count(), 3);
		}

		[TestMethod]
		void ShouldPass()
		{
			Assert::AreEqual(StreamPower::Count(), 2);
		}
	};
}
