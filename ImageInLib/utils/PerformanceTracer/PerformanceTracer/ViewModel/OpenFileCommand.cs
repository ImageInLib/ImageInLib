using Microsoft.Win32;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Input;

namespace PerformanceTracer.ViewModel
{
    class OpenFileCommand : ICommand
    {
        View view = null;
        public OpenFileCommand(View v)
        {
            view = v;
        }

        public bool CanExecute(object parameter)
        {
            if (CanExecuteChanged == null)
                return true;

            return true;
        }

        public event EventHandler CanExecuteChanged;

        public void Execute(object parameter)
        {
            string file_name = null;

            OpenFileDialog openFileDialog = new OpenFileDialog();
            if (openFileDialog.ShowDialog() == true)
            {
                file_name = openFileDialog.FileName;
                DataModel.ThreadContainer thread_cont = new DataModel.ThreadContainer();
                thread_cont.ProcessLogFile(file_name);

                view.Threads = new ThreadList(thread_cont);
                view.ResetZoom();                
                view.Lines = new ThreadLines(view.Threads.Threads);
            }
            //throw new NotImplementedException();
        }
    }
}
