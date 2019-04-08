using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Reflection;
using System.Windows.Media;
using PerformanceTracer.DataModel;

namespace PerformanceTracer.ViewModel
{
    class ThreadLines : INotifyPropertyChanged
    {
        private List<ViewModel.Thread> m_ThreadContainer = null;
        public List<ThreadLine> Lines { get; protected set; }
        protected double CurrentShift = 0;
        protected double CurrentWindowWidth = 0;
        public TimeLine Ruler { get; set; }
        public double Zoom = 1.0;
        public double MinTime { get; set; }
        public double Shift
        {
            get
            {
                return CurrentShift;
            }
            set
            {
                if (CurrentShift == value)
                    return;

                CurrentShift = value;

                UpdateWindowSize();
            }
        }
        protected double CurrentWindowMax = 0;
        protected double CurrentWindowMin = 0;

        private void UpdateWindowSize()
        {
            CurrentWindowMax = (MinTime + (CurrentShift + WindowWidth) / (Zoom * FunctionCall.TimeScale));
            CurrentWindowMin = (MinTime + (CurrentShift) / (Zoom * FunctionCall.TimeScale));

            if (PropertyChanged != null)
            {
                PropertyChanged(this, new PropertyChangedEventArgs("WindowMax"));
                PropertyChanged(this, new PropertyChangedEventArgs("WindowMin"));
            }

            if (Ruler != null)
                Ruler.UpdateMinMax(WindowMin,WindowMax,Zoom,MinTime, Shift);
        }
        public double WindowWidth
        {
            get
            {
                return CurrentWindowWidth;
            }
            set
            {
                if (CurrentWindowWidth == value)
                    return;

                CurrentWindowWidth = value;

                UpdateWindowSize();
            }
        }
        public double WindowMax
        {
            get
            {
                return CurrentWindowMax;
            }
            set
            {
                CurrentWindowMax = value;
            }
        }
        public double WindowMin
        {
            get
            {
                return CurrentWindowMin;
            }
            set
            {
                CurrentWindowMin = value;
            }
        }
        public ThreadLines(List<ViewModel.Thread> iThreadsList_vm)
        {
            Lines = new List<ThreadLine>();
            m_ThreadContainer = iThreadsList_vm;
            MinTime = double.MaxValue;

            foreach (ViewModel.Thread thread in m_ThreadContainer)
            {
                bool new_line_needed = true; 
               
                foreach(ThreadLine line in  Lines)
                {
                    if (line.AddThread(thread))
                    {
                        new_line_needed = false;
                        break;
                    }
                }
                if(new_line_needed)
                    Lines.Add(new ThreadLine(thread));

                if (MinTime > thread.MinTime)
                    MinTime = thread.MinTime;
            }

            WindowWidth = 0;
            UpdateWindowSize();
            Ruler = new TimeLine(WindowMin, WindowMax, Zoom, MinTime, Shift);

        }
        internal void SetZoom(double CurrZoom)
        {
            foreach (ThreadLine entry in Lines)
            {
                entry.Zoom = CurrZoom;
            }

            Shift = Shift * CurrZoom / Zoom;
            Zoom = CurrZoom;

            UpdateWindowSize();

            if (PropertyChanged != null)
            {
                PropertyChanged(this, new PropertyChangedEventArgs("Shift"));
            }
        }
        public event PropertyChangedEventHandler PropertyChanged;
    }
}
