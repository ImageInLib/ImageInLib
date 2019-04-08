using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;

namespace PerformanceTracer.ViewModel
{
    class FunctionCall: INotifyPropertyChanged
    {
        public event PropertyChangedEventHandler PropertyChanged;
        private uint Level = 0;
        private double CurrenZoom = 1.0;
        DataModel.FunctionCall FunCall_dm;
         public Brush Color { get; set; }

        public FunctionCall(DataModel.FunctionCall ifun_call)
        {
            FunCall_dm = ifun_call;
            Level = ifun_call.Level;

            string str = String.Format("Thread: {0}\n {1}\nStart time: {2} ms\nDuration: {3:0.######} ms", FunCall_dm.ThreadID.ToString(), FunCall_dm.Description, FunCall_dm.Start, FunCall_dm.Stop - FunCall_dm.Start);
            Description = str;
            Zoom = 1;
            Color = Brushes.Red;
        }

        public static uint TimeScale = 25;
        public static uint LevelHeight = 50;

        public double Top { get { return Level * LevelHeight; } }
        public double Bottom { get { return (Level + 1) * LevelHeight; } }
        public double Left { get { return (FunCall_dm.Start-StartTimeOffset) * TimeScale * Zoom; } }
        public double Right { get { return (FunCall_dm.Stop-StartTimeOffset)* TimeScale * Zoom; } }
        public double Width { get {
            double width = Right - Left;
            //minimal call width
            const double thr_x = 5;

            if (width < thr_x)
                return thr_x;
            return width; 
        } }
        public double Height { get { return Bottom - Top; } 
        }
        public string Description { get; private set; }
        public string BasicDescription { get { return FunCall_dm.Description;}}
        public double Zoom 
        { 
            get
            {
                return CurrenZoom;
            }
            set
            {
                if (CurrenZoom == value)
                    return;

                CurrenZoom = value;

                if (PropertyChanged != null)
                {
                    PropertyChanged(this, new PropertyChangedEventArgs("Zoom"));
                    PropertyChanged(this, new PropertyChangedEventArgs("Left"));
                    PropertyChanged(this, new PropertyChangedEventArgs("Right"));
                    PropertyChanged(this, new PropertyChangedEventArgs("Width"));
                }
            }
        }

        public double StartTimeOffset { protected get; set; }
    }
}
