using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace PerformanceTracer.DataModel
{
    public class ThreadContainer
    {
        
        public struct Log
        {
            public double time;
            public uint thread_id;
            public uint call_id;
            public string log_action;
            public string description;
        }

        public Dictionary<uint, Thread> Threads { get; set; }

        public ThreadContainer() { Threads = new Dictionary<uint, Thread>(); }

        //reads log file (with format defined in CPerformanceLog) and fill container Threads
        public void ProcessLogFile(string FileName)
        {
            StreamReader sr = null;

            try
            {
                sr = new StreamReader(FileName);
            }
            catch (Exception ex)
            {
                string str_ex = ex.Message;
                MessageBox.Show(str_ex);
                return;
            }

            string str = null;
            Log log = new Log();

            do
            {
                try
                {
                    str = sr.ReadLine();
                }
                catch
                {
                    break;
                }

                if (str == null)
                    break;

                

                if(!ParseLog(str, ref log))
                {
                    MessageBox.Show("Incorrect file format");
                    Threads.Clear();
                    return;
                }

                AddLog(log);

            } while (str != null);

            sr.Close();
        }

        static protected bool ParseLog(string str, ref Log log)
        {
            if (str == null)
            {
                return false;
            }

            string[] log_items = str.Split(new char[] { '\t' });
            //decimal separator will be used in order to solve problem with current culture setting and double conversion
            char dec_sep = Convert.ToChar(CultureInfo.CurrentCulture.NumberFormat.NumberDecimalSeparator);
            //assure usage of correct decimal separator
            if (dec_sep == '.')
                log_items[0] = log_items[0].Replace(',', dec_sep);
            else
                log_items[0] = log_items[0].Replace('.', dec_sep);

            try
            {
                log.time = double.Parse(log_items[0]);

                if (log.time < 0)
                    return false;

                log.thread_id = uint.Parse(log_items[1]);
                log.call_id = uint.Parse(log_items[2]);
            }
            catch
            {
                return false;
            }

            log.log_action = log_items[3];
            log.description = log_items[4];

            return true;
        }

        //Adds log to container Threads
        protected void AddLog(Log log)
        {
            bool badd_new_thread = false;
            bool badd_new_call = false;
            bool bnew_thread_added = false;

            //identify wether new thread/function call should be added to dontainer
            //and add new thread if it is indicated
            try
            {
                badd_new_thread = !Threads.ContainsKey(log.thread_id);
                badd_new_call = badd_new_thread || !Threads[log.thread_id].FunctionCalls.ContainsKey(log.call_id);

                if (badd_new_thread)
                {
                    Threads.Add(log.thread_id, new Thread(log.thread_id));
                    bnew_thread_added = true;
                }
            }
            catch
            {
                return;
            }

            try
            {
                if (badd_new_call)
                {
                    Threads[log.thread_id].FunctionCalls.Add(log.call_id, new FunctionCall());
                    Threads[log.thread_id].FunctionCalls[log.call_id].ThreadID = log.thread_id;                   
                }
            }
            catch 
            {
                if (bnew_thread_added)
                {
                    Threads.Remove(log.thread_id);
                }

                return;
            }

            uint level = 1000;
            uint tmp_call_id = log.call_id - 1;

            //compute level of function call
            do
            {
                try
                {
                    if (!Threads[log.thread_id].FunctionCalls.ContainsKey(tmp_call_id))
                        level = 0;
                    else if (Threads[log.thread_id].FunctionCalls[tmp_call_id].Stop == -1)
                    {
                        level = Threads[log.thread_id].FunctionCalls[tmp_call_id].Level + 1;
                        break;
                    }                                      
                }
                catch
                { 
                }

                tmp_call_id--;
            } while (level > 0);
           

            if (log.log_action == "start")
            {
                Threads[log.thread_id].FunctionCalls[log.call_id].Start = log.time;
                Threads[log.thread_id].FunctionCalls[log.call_id].Level = level;
                Threads[log.thread_id].FunctionCalls[log.call_id].Description = log.description;
            }
            else if (log.log_action == "end")
                Threads[log.thread_id].FunctionCalls[log.call_id].Stop = log.time;
            else
                throw new NotImplementedException();
        }
    }
}
