using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;

namespace PerformanceTracer.ViewModel
{
    class ThreadLine : INotifyPropertyChanged
    {
        public SortedDictionary<double, ViewModel.Thread> m_Threads;
        public List<object> LineObjects { get; set; }
        private double MaxTime { get; set; }
        public double MinTime { get; set; }
        private uint Levels { get; set; }
        protected double TimeOffset = 0;
        public Brush Color { get; set; }

        public ThreadLine(ViewModel.Thread iTread_vm)
        {
            MaxTime = 0.0;
            MinTime = double.MaxValue;
            Levels = 0;
            LineObjects = new List<object>();
            m_Threads = new SortedDictionary<double, ViewModel.Thread>();
            AddThread(iTread_vm);
        }
        public bool AddThread(ViewModel.Thread iTread_vm)
        {
            bool add = false;
            double start_time = 0;
            double end_time = double.MaxValue;
           
            //there is no threads in line just add new one
            if (m_Threads.Count() == 0)
            {
                add = true;
            }
            else
            {
                //check if is not space prior existing thread
                KeyValuePair<double, ViewModel.Thread> element = m_Threads.First();
                end_time = element.Value.MinTime;
            }

            int num_of_treads = m_Threads.Count();
            for (int i = 0; i < num_of_treads; i++)
            {                
                if( start_time <= iTread_vm.MinTime &&
                    iTread_vm.MaxTime <= end_time)
                {
                    add = true;
                    break;
                }
                start_time = m_Threads.ElementAt(i).Value.MaxTime;
                if(i+1 >= num_of_treads )
                    end_time = double.MaxValue;
                else
                    end_time = m_Threads.ElementAt(i+1).Value.MinTime;
            }
            //if it possible past last thread 
            if (start_time <= iTread_vm.MinTime &&
                   iTread_vm.MaxTime <= end_time)
            {
                add = true;
            }


            if (!add)
                return false;

            m_Threads.Add(iTread_vm.MinTime,iTread_vm);

            MaxTime = iTread_vm.MaxTime;
            if (MinTime > iTread_vm.MinTime)
                MinTime = iTread_vm.MinTime;

            if (Levels < iTread_vm.Levels)
                Levels = iTread_vm.Levels;

            foreach (var item in iTread_vm.ThreadObjects)
            {
                FunctionCall fun_call = item as FunctionCall;
                if (fun_call == null)
                    continue;
                LineObjects.Add(fun_call);
            }
            return true;
        }

        public double Width
        {
            get
            {
                return (MaxTime - StartTimeOffset) * Zoom * FunctionCall.TimeScale;
            }
        }
        public uint Height { get { return Levels * FunctionCall.LevelHeight; } }

        public double StartTimeOffset
        {
            get
            {
                return TimeOffset;
            }
            set
            {
                if (TimeOffset == value)
                    return;

                TimeOffset = value;

                foreach (var entry in LineObjects)
                {
                    FunctionCall fun_call = entry as FunctionCall;

                    if (fun_call == null)
                        continue;

                    fun_call.StartTimeOffset = TimeOffset;
                }
            }
        }
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

                foreach (var entry in LineObjects)
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
        public event PropertyChangedEventHandler PropertyChanged;
    }
}
