using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PerformanceTracer.ViewModel
{
    class ThreadLabel:INotifyPropertyChanged
    {
        protected Thread Thr = null;
        private double OffsetX = 0;
        private double OffsetY = 0;

        public ThreadLabel(Thread ithread)
        {
            Thr = ithread;

        }

        public double Top { get { return OffsetY; } }
        public double Bottom { get { return Top+Thr.Height; } }
        public double Left { get { return OffsetX; } }
        public double Right { get { return Left; } }
        public string Label { get { return Thr.ThreadID.ToString(); } }
        public double Width { get { return Right - Left; } }
        public double Height { get { return ThreadList.MarginX + Bottom - Top; } }

        internal void SetViewOffset(double offX, double offY)
        {
            OffsetX = offX;
            OffsetY = offY;

            if (PropertyChanged != null)
            {
                PropertyChanged(this, new PropertyChangedEventArgs("Top"));
                PropertyChanged(this, new PropertyChangedEventArgs("Bottom"));
                PropertyChanged(this, new PropertyChangedEventArgs("Left"));
                PropertyChanged(this, new PropertyChangedEventArgs("Label"));
            }
        }

        public event PropertyChangedEventHandler PropertyChanged;
    }
}
