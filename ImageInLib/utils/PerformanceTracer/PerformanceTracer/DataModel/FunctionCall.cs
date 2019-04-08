using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PerformanceTracer.DataModel
{
    public class FunctionCall
    {
        public double Start { get; set; }
        public double Stop { get; set; }
        public string Description { get; set; }
        public uint Level { get; set; }
        public uint ThreadID { get; set; }
        public FunctionCall()
        {
            Start = -1;
            Stop = -1;
            Description = "";
            Level = 0;
            ThreadID = 0;
        }
    }
}

