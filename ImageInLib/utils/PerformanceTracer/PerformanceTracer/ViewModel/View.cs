using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Input;

namespace PerformanceTracer.ViewModel
{
    class View : INotifyPropertyChanged
    {
        public enum ViewMode
        {
            ALL_THREADS,
            ACTIVE_THREADS_ONLY
        };

        public View()
        {
            OpenFile = new OpenFileCommand(this);
            Zoom = new ZoomCommand(this);
            CurrentViewMode = ViewMode.ALL_THREADS;
            this.PropertyChanged += OnMyViewPropertyChanged;
        }

        void OnMyViewPropertyChanged(object sender, PropertyChangedEventArgs e)
        {
            if (e.PropertyName == "Threads")
                PropertyChanged(this, new PropertyChangedEventArgs("ActiveView"));
            if (e.PropertyName == "CurrentViewMode")
                PropertyChanged(this, new PropertyChangedEventArgs("ActiveView"));
            if (e.PropertyName == "CurrentZoom")
                PropertyChanged(this, new PropertyChangedEventArgs("ZoomLabel"));
        }

        public ICommand OpenFile { get; protected set; }
        public ICommand Zoom { get; protected set; }

        ThreadList m_List = null;
        ThreadLines m_Lines = null;

        private const double DefaultZoom = 1.0;
        protected double CurrZoom = DefaultZoom;

        public ThreadList Threads
        {
            get { return m_List; }
            set
            {
                if (m_List == value)
                    return;
                m_List = value;
                if (PropertyChanged != null)
                    PropertyChanged(this, new PropertyChangedEventArgs("Threads"));
            }
        }

        public ThreadLines Lines
        {
            get { return m_Lines; }
            set
            {
                if (m_Lines == value)
                    return;
                m_Lines = value;
                if (PropertyChanged != null)
                    PropertyChanged(this, new PropertyChangedEventArgs("Lines"));
            }
        }

        public event PropertyChangedEventHandler PropertyChanged;

        protected ViewMode CurrentVM = ViewMode.ALL_THREADS;
        protected ViewMode CurrentViewMode
        {
            get
            {
                return CurrentVM;
            }
            set
            {
                if (CurrentVM == value)
                    return;

                CurrentVM = value;
                if (PropertyChanged != null)
                    PropertyChanged(this, new PropertyChangedEventArgs("CurrentViewMode"));
            }
        }

        public bool IsAllThreadsMode
        {
            get
            {
                return CurrentViewMode == ViewMode.ALL_THREADS;
            }
            set
            {
                if (value)
                {
                    CurrentViewMode = ViewMode.ALL_THREADS;
                }
            }
        }
        public bool IsActiveThreadsOnlyMode
        {
            get
            {
                return CurrentViewMode == ViewMode.ACTIVE_THREADS_ONLY;
            }
            set
            {
                if (value)
                {
                    CurrentViewMode = ViewMode.ACTIVE_THREADS_ONLY;

                }
            }
        }

        public object ActiveView
        {
            get
            {
                if (CurrentViewMode == ViewMode.ALL_THREADS)
                    return Threads;
                else
                    return Lines;
            }
        }

        public double CurrentZoom
        {
            get
            {
                return CurrZoom;
            }
            set
            {
                if (CurrZoom == value || Threads == null)
                    return;

                CurrZoom = value;
                Threads.SetZoom(CurrZoom);
                Lines.SetZoom(CurrZoom);

                if (PropertyChanged != null)
                    PropertyChanged(this, new PropertyChangedEventArgs("CurrentZoom"));
            }
        }
        public string ZoomLabel
        {
            get
            {
                return String.Format("Zoom: {0} %", (CurrentZoom * 100).ToString());
            }
        }

        public void ResetZoom()
        {
            CurrentZoom = DefaultZoom;
        }

        public double ZoomInStep { get { return 0.5; } }
        public double ZoomOutStep { get { return -ZoomInStep; } }
    }

}
