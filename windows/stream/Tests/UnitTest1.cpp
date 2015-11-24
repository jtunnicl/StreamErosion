#include "stdafx.h"
#include "streampower.h"

using namespace System;
using namespace System::Text;
using namespace System::Collections::Generic;
using namespace	Microsoft::VisualStudio::TestTools::UnitTesting;

namespace Tests
{
	[TestClass]
	public ref class UnitTest1
	{
	public: 
		[TestMethod]
		void TestMethod1()
		{
			Assert::AreEqual(StreamPower::Count(), 3);
		}
	};
}
