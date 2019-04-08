using PerformanceTracer.DataModel;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PerformanceTracerTests
{
    class TestThreadContainer: ThreadContainer
    {
        static public bool TestParseLog(string str, ref Log log)
        {
            return ParseLog(str, ref log);
        }

        public void TestAddLog(Log log)
        {
            AddLog(log);
        }
    }
}
