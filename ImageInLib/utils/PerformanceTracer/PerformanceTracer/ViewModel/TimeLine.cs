using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PerformanceTracer.ViewModel
{
    class TimeLine : INotifyPropertyChanged
    {
        protected double CurrentMin = 0;
        protected double CurrentMax = 0;
        protected double CurrentZoom = 1.0;

        private ObservableCollection<Tick> m_Ticks = null;
        public ObservableCollection<Tick> Ticks 
        {
            get {
                return m_Ticks; 
            }
            protected set
            {
                m_Ticks = value;
            }
        }

        public TimeLine(double iMin, double iMax, double iZoom, double iTimeOffset, double iShift)
        {
            UpdateMinMax(iMin, iMax, iZoom, iTimeOffset, iShift);

            m_Ticks = new ObservableCollection<Tick>();
            UpdateTicks();

        }

        private void UpdateTicks()
        {
            double one = FunctionCall.TimeScale * CurrentZoom;
            double min = 0.25*FunctionCall.TimeScale;
            double c = one/min;
            double exp = Math.Log(c, 10);
            double step_base = Math.Pow(10, -Math.Floor(exp));
            double step = 0.5 * step_base;
            
            if (m_Ticks == null)
                return;

            m_Ticks.Clear();
            
            double nmax = Max;
            double tick_size = 1.0;
            double eps = 0.00001;
            string str = "";
            double x_coord = 0;
            bool bshow_text = true;
            double last_text_coord = -1;
            const double text_step = 50;
            
            for (double x = 0; x <= nmax; x += step)
            {
                if (x < Min)
                    continue;

                x_coord = (x - TimeOffset) * FunctionCall.TimeScale * CurrentZoom - Shift;

                if (Math.Abs(x % step_base ) < eps || Math.Abs(x % step_base - step_base) < eps)
                {
                    tick_size = 1.0;
                    double diff = (x * one - last_text_coord);

                    if (last_text_coord < 0 || last_text_coord >= 0 && diff - eps >= text_step)
                        {
                            last_text_coord = x * one;
                            bshow_text = true;
                            tick_size = 1.0;
                        }
                        else
                            bshow_text = false;
                }
                else 
                {
                    tick_size = 0.5;
                    bshow_text = false;
                }

                if (bshow_text)
                {
                    str = Math.Round(x, 7).ToString();
                }
                else
                    str = "";

                m_Ticks.Add(new Tick(x_coord, str, tick_size, bshow_text));
            }
        }

        public void UpdateMinMax(double iMin,double iMax,double iZoom,double iTimeOffset, double iShift)
        {
            Min = iMin;
            Max = iMax;
            CurrentZoom = iZoom;
            TimeOffset = iTimeOffset;
            Shift = iShift;

            UpdateTicks();
        }

        public double Min
        {
            get
            {
                return CurrentMin;
            }
            set
            {
                if (Min == value)
                    return;

                CurrentMin = value;

                if (PropertyChanged != null)
                {
                    PropertyChanged(this, new PropertyChangedEventArgs("Min"));

                }
            }
        }

        public double Max
        {
            get
            {
                return CurrentMax;
            }
            set
            {
                if (Max == value)
                    return;

                CurrentMax = value;

                if (PropertyChanged != null)
                {
                    PropertyChanged(this, new PropertyChangedEventArgs("Max"));

                }
            }
        }

        public event PropertyChangedEventHandler PropertyChanged;
        public double TimeOffset { get; set; }

        protected double Shift { get; set; }
        public static double Height { get { return 50; } }
    }
}
