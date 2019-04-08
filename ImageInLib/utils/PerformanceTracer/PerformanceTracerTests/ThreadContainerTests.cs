using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace PerformanceTracerTests
{
    [TestClass]
    public class ThreadContainerTests
    {
        [TestMethod]
        public void ProcessLogFile()
        {
            PerformanceTracer.DataModel.ThreadContainer thr_cont = new PerformanceTracer.DataModel.ThreadContainer();

            Assert.IsNotNull(thr_cont);
        }

        [TestMethod]
        public void ParseLog()
        {
            TestThreadContainer thr_cont = new TestThreadContainer();
            Assert.IsNotNull(thr_cont);

            TestThreadContainer.Log log = new TestThreadContainer.Log();
            Assert.IsNotNull(log);

            //test null string
            Assert.IsFalse(TestThreadContainer.TestParseLog(null, ref log));

            //test different decimal separator
            //and parsed values
            Assert.IsTrue(TestThreadContainer.TestParseLog(String.Format("0.0\t1\t2\tstart\tFUN0()"), ref log));
            Assert.IsTrue(log.time == 0);
            Assert.IsTrue(log.thread_id == 1);
            Assert.IsTrue(log.call_id == 2);
            Assert.IsTrue(log.log_action == "start");
            Assert.IsTrue(log.description == "FUN0()");

            Assert.IsTrue(TestThreadContainer.TestParseLog(String.Format("1,0\t3\t4\tend\tFUN0()"), ref log));
            Assert.IsTrue(log.time == 1.0);
            Assert.IsTrue(log.thread_id == 3);
            Assert.IsTrue(log.call_id == 4);
            Assert.IsTrue(log.log_action == "end");
            Assert.IsTrue(log.description == "FUN0()");

            //test uncorrect format
            Assert.IsFalse(TestThreadContainer.TestParseLog(String.Format("-0.1\t0\t0\tstart\tFUN0()"), ref log));
            Assert.IsFalse(TestThreadContainer.TestParseLog(String.Format("\t0\t0\tstart\tFUN0()"), ref log));
            Assert.IsFalse(TestThreadContainer.TestParseLog(String.Format("str1\tstr2\tstr3\tstart\tFUN0()"), ref log));
        }

        [TestMethod]
        public void AddLog()
        {
            TestThreadContainer thr_cont = new TestThreadContainer();
            Assert.IsNotNull(thr_cont);

            TestThreadContainer.Log log = new TestThreadContainer.Log();
            
            log.time = 0.2;
            log.thread_id = 0;
            log.call_id = 0;
            log.description = "FUN1()";
            log.log_action = "start";
            thr_cont.TestAddLog(log);

            //is the 1st thred added
            Assert.IsTrue(thr_cont.Threads.ContainsKey(log.thread_id));

            log.time = 0.3;
            log.call_id = 1;
            log.description = "FUN2()";
            thr_cont.TestAddLog(log);
            Assert.IsTrue(thr_cont.Threads.ContainsKey(log.thread_id));
            Assert.IsTrue(thr_cont.Threads.Count == 1);

            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls.Count == 2);
            //test thread 0 - call 0
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls.ContainsKey(0));
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls[0].Start == 0.2);
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls[0].Stop == -1);
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls[0].Level == 0);
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls[0].Description == "FUN1()");
            //test thread 0 - call 1
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls.ContainsKey(1));
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls[1].Start == 0.3);
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls[1].Stop == -1);
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls[1].Level == 1);
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls[1].Description == "FUN2()");

            log.thread_id = 1;
            log.time = 0.32;
            log.call_id = 0;
            log.description = "FUN1()";
            thr_cont.TestAddLog(log);
            Assert.IsTrue(thr_cont.Threads.ContainsKey(log.thread_id));
            Assert.IsTrue(thr_cont.Threads.Count == 2);
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls.Count == 1);
            //test thread 1 - call 0
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls.ContainsKey(0));
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls[0].Start == 0.32);
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls[0].Stop == -1);
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls[0].Level == 0);
            Assert.IsTrue(thr_cont.Threads[log.thread_id].FunctionCalls[0].Description == "FUN1()");


            log.time = 0.4;
            log.log_action = "end";
            log.thread_id = 0;
            log.call_id = 1;
            thr_cont.TestAddLog(log);

            log.time = 0.42;
            log.thread_id = 1;
            log.call_id = 0;
            thr_cont.TestAddLog(log);

            log.time = 0.45;
            log.thread_id = 0;
            thr_cont.TestAddLog(log);

            Assert.IsTrue(thr_cont.Threads.Count == 2);
            Assert.IsTrue(thr_cont.Threads[0].FunctionCalls.Count == 2);
            //test thread 0 - call 0
            Assert.IsTrue(thr_cont.Threads[0].FunctionCalls.ContainsKey(0));
            Assert.IsTrue(thr_cont.Threads[0].FunctionCalls[0].Start == 0.2);
            Assert.IsTrue(thr_cont.Threads[0].FunctionCalls[0].Stop == 0.45);
            Assert.IsTrue(thr_cont.Threads[0].FunctionCalls[0].Level == 0);
            Assert.IsTrue(thr_cont.Threads[0].FunctionCalls[0].Description == "FUN1()");
            //test thread 0 - call 1
            Assert.IsTrue(thr_cont.Threads[0].FunctionCalls.ContainsKey(1));
            Assert.IsTrue(thr_cont.Threads[0].FunctionCalls[1].Start == 0.3);
            Assert.IsTrue(thr_cont.Threads[0].FunctionCalls[1].Stop == 0.4);
            Assert.IsTrue(thr_cont.Threads[0].FunctionCalls[1].Level == 1);
            Assert.IsTrue(thr_cont.Threads[0].FunctionCalls[1].Description == "FUN2()");

            Assert.IsTrue(thr_cont.Threads.ContainsKey(log.thread_id));
            Assert.IsTrue(thr_cont.Threads[1].FunctionCalls.Count == 1);
            //test thread 1 - call 0
            Assert.IsTrue(thr_cont.Threads[1].FunctionCalls.ContainsKey(0));
            Assert.IsTrue(thr_cont.Threads[1].FunctionCalls[0].Start == 0.32);
            Assert.IsTrue(thr_cont.Threads[1].FunctionCalls[0].Stop == 0.42);
            Assert.IsTrue(thr_cont.Threads[1].FunctionCalls[0].Level == 0);
            Assert.IsTrue(thr_cont.Threads[1].FunctionCalls[0].Description == "FUN1()");

            //test incorrect log action
            log.log_action = "helmut";
            try
            {
                thr_cont.TestAddLog(log);
            }
            catch
            {
                Assert.IsTrue(true);
                return;
            }

            Assert.IsTrue(false);
        }
    }
}
