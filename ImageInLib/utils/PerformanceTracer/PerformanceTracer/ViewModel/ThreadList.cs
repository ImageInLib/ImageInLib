using PerformanceTracer.DataModel;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;

namespace PerformanceTracer.ViewModel
{
    class ThreadList : INotifyPropertyChanged
    {
        private ThreadContainer ThreadCont = null;
        public static uint MarginX = 0;
        private uint colorIndex = 0;
        protected double CurrentShift = 0;
        protected double CurrentWindowWidth = 0;

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

        private void UpdateWindowSize()
        {
            CurrentWindowMax = (MinTime + (CurrentShift + WindowWidth) / (Zoom* FunctionCall.TimeScale));
            CurrentWindowMin = (MinTime + (CurrentShift) / (Zoom * FunctionCall.TimeScale));

            if (Ruler != null)
                Ruler.UpdateMinMax(CurrentWindowMin, CurrentWindowMax, Zoom, MinTime, Shift);

            if (PropertyChanged != null)
            {
                PropertyChanged(this, new PropertyChangedEventArgs("WindowMax"));
                PropertyChanged(this, new PropertyChangedEventArgs("WindowMin"));
            }


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

        protected double CurrentWindowMax = 0;
        protected double CurrentWindowMin = 0;

        public string WindowMax 
        { 
            get
            {
                return CurrentWindowMax.ToString();
            }
            set
            {
                CurrentWindowMax = double.Parse(value);
            }
        }

        public string WindowMin
        {
            get
            {
                return CurrentWindowMin.ToString();
            }
            set
            {
                CurrentWindowMin = double.Parse(value);
            }
        }

        public TimeLine Ruler { get; set; }
        public double MinTime { get; set; }

        public double Zoom = 1.0;

        public ThreadList(ThreadContainer thread_cont)
        {
            Threads = new List<Thread>();
            ThreadCont = thread_cont;

            foreach (KeyValuePair<uint, DataModel.Thread> entry in ThreadCont.Threads)
            {
                Thread thr = new Thread(entry.Value, GetNextBrush());
                Threads.Add(thr);
            }

            MinTime = double.MaxValue;
            foreach(Thread thr in Threads)
            {

                if (MinTime > thr.MinTime)
                    MinTime = thr.MinTime;
            }

            foreach(Thread thr in Threads)
            {
                thr.StartTimeOffset = MinTime;
            }

            WindowWidth = 0;
            UpdateWindowSize();
            Ruler = new TimeLine(double.Parse(WindowMin), double.Parse(WindowMax), Zoom, MinTime, Shift);
        }

        public List<Thread> Threads { get; protected set; }

        //get next brush
        private Brush GetNextBrush()
        {
            Type brushesType = typeof(System.Windows.Media.Brushes);

            // Get all static properties
            var properties = brushesType.GetProperties(BindingFlags.Static | BindingFlags.Public);

            if (colorIndex >= properties.Length)
                colorIndex = 0;

            return (Brush)properties[colorIndex++].GetValue(null, null);
        }


        internal void SetZoom(double CurrZoom)
        {
            foreach (Thread entry in Threads)
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
