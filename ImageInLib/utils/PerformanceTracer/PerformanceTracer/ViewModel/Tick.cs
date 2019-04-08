using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PerformanceTracer.ViewModel
{
    class Tick
    {
        protected double Length = 10;
        protected double YOffset = 0;
        public double X {get; set;}
        public double HorizMargin { protected get; set; }
        public string Description { get; protected set; }

        public Tick(double XCoord, string iDescription, double iTickSize, bool iShowText)
        {
            HorizMargin = 30;
            X = XCoord;

            if (iShowText)
            {
                Description = iDescription;
                Thickness = 2;
            }
            else
                Thickness = 1;

            Length *= iTickSize;

        }

        public double Top 
        {
            get
            {
                return 0;
            }
        }
        
        public double Bottom 
        {
            get
            {
                return Length;
            }
            
        }

        public double BorderTop { get { return BorderBottom - Length - 20; } }
        public double BorderBottom { 
            get 
            { 
                return TimeLine.Height; 
            } 
        }
        public double BorderLeft 
        { 
            get 
            { 
                return X - HorizMargin; 
            } 
        }
        public double BorderRight { get { return X + HorizMargin; } }
        public double BorderWidth
        {
            get
            {
                return BorderRight - BorderLeft;
            }
        }
        public double BorderHeight
        {
            get
            {
                return BorderBottom - BorderTop;
            }
        }

        public int Thickness { get; protected set; }
    }
}
