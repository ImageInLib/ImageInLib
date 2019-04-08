using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;

namespace PerformanceTracer.ViewModel
{
    class Thread:INotifyPropertyChanged
    {
        DataModel.Thread Thread_dm;
        protected double TimeOffset = 0;

        public Thread(DataModel.Thread ithread, Brush Color)
        {
            Thread_dm = ithread;
            MaxTime = 0.0;
            MinTime = double.MaxValue;

            Levels = 0;
            ThreadObjects = new List<object>();
            

            foreach (KeyValuePair<uint, DataModel.FunctionCall> entry in ithread.FunctionCalls)
            {
                if (Levels < entry.Value.Level)
                    Levels = entry.Value.Level;

                if (MaxTime < entry.Value.Stop)
                    MaxTime = entry.Value.Stop;

                if (MinTime > entry.Value.Start)
                    MinTime = entry.Value.Start;

                FunctionCall new_call = new FunctionCall(entry.Value);
                new_call.Color = Color;
                ThreadObjects.Add(new_call);
            }

            //ThreadObjects.Add(new ThreadLabel(this));

            Levels++;
        }

        public List<object> ThreadObjects { get; set; }
       
        public double Width 
        { 
            get 
            { return (MaxTime - StartTimeOffset)  * Zoom * FunctionCall.TimeScale; 
            }
            
        }
        public uint Height { get { return Levels * FunctionCall.LevelHeight; } }

       
        public double MinTime { get; protected set; }
        public double MaxTime { get; protected set; }
        public uint Levels    {get; protected set;}
        

        public uint ThreadID { get { return Thread_dm.ThreadId; } }

        double CurrenZoom = 1.0;

        internal double Zoom
        {
            get
            {
                return CurrenZoom;
            }

            set
            {
                if (CurrenZoom == value)
                    return;

                foreach (var entry in ThreadObjects)
                {
                    FunctionCall fun_call = entry as FunctionCall;

                    if (fun_call == null)
                        continue;

                    fun_call.Zoom = value;
                }

                CurrenZoom = value;
                
                if (PropertyChanged != null)
                {
                    PropertyChanged(this, new PropertyChangedEventArgs("Width"));
                    //not needed
                    //PropertyChanged(this, new PropertyChangedEventArgs("Zoom"));
                }
            }
        }

        public double StartTimeOffset 
        { 
            get
            {
                return TimeOffset;
            }
            set
            {
                if(TimeOffset == value)
                    return;

                TimeOffset = value;

                foreach (var entry in ThreadObjects)
                {
                    FunctionCall fun_call = entry as FunctionCall;

                    if (fun_call == null)
                        continue;

                    fun_call.StartTimeOffset = TimeOffset;
                }

               }
        }

        public event PropertyChangedEventHandler PropertyChanged;
    }
}
