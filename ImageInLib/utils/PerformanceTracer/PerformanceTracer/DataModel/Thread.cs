using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;

namespace PerformanceTracer.DataModel
{
    public class Thread
    {
        public Thread(uint ithr_id) 
        {
            ThreadId = ithr_id;
            FunctionCalls = new Dictionary<uint, FunctionCall>();
        }
        public Dictionary<uint, FunctionCall> FunctionCalls { get; set; }
        public uint ThreadId { get; set; }
    }
}
