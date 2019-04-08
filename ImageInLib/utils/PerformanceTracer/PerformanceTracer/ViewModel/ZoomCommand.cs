using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Input;

namespace PerformanceTracer.ViewModel
{
    class ZoomCommand : ICommand
    {
        View view = null;

        public ZoomCommand(View v)
        {
            view = v;
            view.PropertyChanged += view_PropertyChanged;
        }

        void view_PropertyChanged(object sender, System.ComponentModel.PropertyChangedEventArgs e)
        {
            if (e.PropertyName == "Threads")
            {
                if (CanExecuteChanged != null)
                    CanExecuteChanged(this, null);
            }
        }

        public event EventHandler CanExecuteChanged;

        public void Execute(object parameter)
        {
            double zoom_change = (double)parameter;

            if (view.CurrentZoom >= 1.0)
                view.CurrentZoom = view.CurrentZoom + zoom_change;
            else
            {
                if(zoom_change < 0)
                    view.CurrentZoom = view.CurrentZoom * Math.Abs(zoom_change);
                else
                    view.CurrentZoom = view.CurrentZoom / Math.Abs(zoom_change);
            }
            return;
        }

        public bool CanExecute(object parameter)
        {
            if (view == null || view.Threads == null)
                return false;

            if (CanExecuteChanged == null)
                return true;

            return true;
        }
    }
}
