using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Controls.Primitives;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace ThreadsGraphControl
{
    [TemplatePart(Name = ThreadsGraph.PART_ItemsControl, Type = typeof(ItemsPresenter))]
    [TemplatePart(Name = ThreadsGraph.PART_Ruler, Type = typeof(ContentPresenter))]

    public class ThreadsGraph : ListBox, IScrollInfo, System.ComponentModel.INotifyPropertyChanged
    {
        public const string PART_ItemsControl = "PART_ItemsControl";
        public const string PART_Ruler = "PART_Ruler";

        static ThreadsGraph()
        {
            DefaultStyleKeyProperty.OverrideMetadata(typeof(ThreadsGraph), new FrameworkPropertyMetadata(typeof(ThreadsGraph)));
        }

        public static readonly DependencyProperty ItemHeaderTemplateProperty = DependencyProperty.Register(
            "ItemHeaderTemplate",
            typeof(DataTemplate),
            typeof(ThreadsGraph),
            new PropertyMetadata((s, e) => ((ThreadsGraph)s).ApplyItemHeaderTemplate())
            );

        public DataTemplate ItemHeaderTemplate
        {
            get { return (DataTemplate)GetValue(ItemHeaderTemplateProperty); }
            set { SetValue(ItemHeaderTemplateProperty, value); }
        }

        public void ApplyItemHeaderTemplate()
        {
            InvalidateMeasure();
            InvalidateArrange();
        }

        public static readonly DependencyProperty RulerTemplateProperty = DependencyProperty.Register(
            "RulerTemplate",
            typeof(DataTemplate),
            typeof(ThreadsGraph),
            new PropertyMetadata((s, e) => ((ThreadsGraph)s).ApplyItemHeaderTemplate())
            );

        public DataTemplate RulerTemplate
        {
            get { return (DataTemplate)GetValue(RulerTemplateProperty); }
            set { SetValue(RulerTemplateProperty, value); }
        }

        public static readonly DependencyProperty XOffsetProperty = DependencyProperty.Register(
            "XOffset",
            typeof(double),
            typeof(ThreadsGraph),
            new PropertyMetadata((s, e) => ((ThreadsGraph)s).OnXOffsetChanged())
            );

        public void OnXOffsetChanged()
        {
            SetHorizontalOffset(XOffset);
            if (ScrollOwner != null)
                ScrollOwner.InvalidateScrollInfo();
            SetVisibleXMax();

            if (PropertyChanged != null)
                PropertyChanged(this, new System.ComponentModel.PropertyChangedEventArgs("XOffset"));
        }

        private void SetVisibleXMax()
        {
            VisibleXMax = _viewport.Width + _offset.X;
        }

        public double XOffset
        {
            get { return (double)GetValue(XOffsetProperty); }
            set { SetValue(XOffsetProperty, value); }
        }

        public static readonly DependencyProperty VisibleXMaxProperty = DependencyProperty.Register(
    "VisibleXMax",
    typeof(double),
    typeof(ThreadsGraph),
    new PropertyMetadata((s, e) => ((ThreadsGraph)s).OnVisibleXMaxChanged())
    );

        public void OnVisibleXMaxChanged()
        {
        }

        public double VisibleXMax
        {
            get { return (double)GetValue(VisibleXMaxProperty); }
            set { SetValue(VisibleXMaxProperty, value); }
        }

        public ThreadsGraph()
        {
        }

        protected override DependencyObject GetContainerForItemOverride()
        {
            return new ThreadGraphListItem(this);
        }


        protected override bool IsItemItsOwnContainerOverride(object item)
        {
            return (item as ThreadGraphListItem != null);
        }

        public bool CanHorizontallyScroll
        {
            get;
            set;
        }

        public bool CanVerticallyScroll
        {
            get;
            set;
        }

        public double ExtentHeight
        {
            get
            {
                return _extent.Height;
            }
        }

        public double ExtentWidth
        {
            get { return _extent.Width; }
        }

        public double HorizontalOffset
        {
            get { return _offset.X; }
        }

        public void LineDown()
        {
            _offset.Y += 1;
            ValidateOffset();
            InvalidateArrange();
        }

        public void LineLeft()
        {
            _offset.X -= 1;
            ValidateOffset();
            XOffset = _offset.X;
            InvalidateArrange();
        }

        public void LineRight()
        {
            _offset.X += 1;
            ValidateOffset();
            XOffset = _offset.X;
            InvalidateArrange();
        }

        public void LineUp()
        {
            _offset.Y -= 1;
            ValidateOffset();
            InvalidateArrange();
        }

        public Rect MakeVisible(Visual visual, Rect rectangle)
        {
            return rectangle;
            //throw new NotImplementedException();
        }

        public void MouseWheelDown()
        {
            LineDown();
        }

        public void MouseWheelLeft()
        {
            LineLeft();
        }

        public void MouseWheelRight()
        {
            LineRight();
        }

        public void MouseWheelUp()
        {
            LineUp();
        }

        public void PageDown()
        {
            _offset.Y += _viewport.Height;
            ValidateOffset();
            InvalidateArrange();
        }

        public void PageLeft()
        {
            _offset.Y -= _viewport.Width;
            ValidateOffset();
            InvalidateArrange();
        }

        public void PageRight()
        {
            _offset.Y += _viewport.Width;
            ValidateOffset();
            InvalidateArrange();
        }

        public void PageUp()
        {
            _offset.Y -= _viewport.Height;
            ValidateOffset();
            InvalidateArrange();
        }

        public ScrollViewer ScrollOwner
        {
            get;
            set;
        }

        public void SetHorizontalOffset(double offset)
        {
            _offset.X = offset;
            ValidateOffset();
            XOffset = _offset.X;
            InvalidateArrange();
        }

        public void SetVerticalOffset(double offset)
        {
            _offset.Y = offset;
            ValidateOffset();
            InvalidateArrange();
        }

        public double VerticalOffset
        {
            get { return _offset.Y; }
        }

        public double ViewportHeight
        {
            get { return _viewport.Height; }
        }

        public double ViewportWidth
        {
            get { return _viewport.Width; }
        }

        private Point _offset = new Point(0, 0);
        private void ValidateOffset()
        {
            if (_offset.X < 0)
                _offset.X = 0;

            if (_offset.Y < 0)
                _offset.Y = 0;

            double max_x = _extent.Width - _viewport.Width;
            double max_y = _extent.Height - _viewport.Height;

            if (max_x < 0)
                max_x = 0;
            if (max_y < 0)
                max_y = 0;

            if (_offset.X >= max_x)
                _offset.X = max_x;

            if (_offset.Y >= max_y)
                _offset.Y = max_y;
        }

        private Size _extent = new Size(0, 0);
        private Size _viewport = new Size(0, 0);

        private static Size InfiniteSize = new Size(double.PositiveInfinity, double.PositiveInfinity);

        protected override Size MeasureOverride(Size availableSize)
        {
            FrameworkElement ItemsControl = GetTemplateChild(PART_ItemsControl) as FrameworkElement;
            FrameworkElement ruler = GetTemplateChild(PART_Ruler) as FrameworkElement;
            IScrollInfo info = ItemsControl as IScrollInfo;

            ruler.Measure(availableSize);
            ItemsControl.Measure(InfiniteSize);

            Size itemsControlSize = ItemsControl.DesiredSize;
            Size rulerSize = ruler.DesiredSize;
            Size extent = new Size(0, 0);
            extent.Width = itemsControlSize.Width;
            if (rulerSize.Width > extent.Width)
                extent.Width = rulerSize.Width;

            extent.Height = itemsControlSize.Height + rulerSize.Height;
            if (extent != _extent)
            {
                _extent = extent;
                ScrollOwner.InvalidateScrollInfo();
            }

            if (availableSize != _viewport)
            {
                _viewport = availableSize;

                SetVisibleXMax();

                ScrollOwner.InvalidateScrollInfo();
            }
            return availableSize;
        }

        protected override Size ArrangeOverride(Size arrangeBounds)
        {
            Size retVal = base.ArrangeOverride(arrangeBounds);

            FrameworkElement ItemsControl = GetTemplateChild(PART_ItemsControl) as FrameworkElement;
            FrameworkElement ruler = GetTemplateChild(PART_Ruler) as FrameworkElement;

            double rullerHeight = ruler.DesiredSize.Height;

            Rect rulerRect = new Rect(arrangeBounds);
            rulerRect.Height = rullerHeight;
            rulerRect.Y = 0;
            ruler.Arrange(rulerRect);

            Rect itemsRect = new Rect(arrangeBounds);
            itemsRect.Height = ItemsControl.DesiredSize.Height;
            itemsRect.Y -= _offset.Y - rullerHeight;
            ItemsControl.InvalidateArrange();
            ItemsControl.Arrange(itemsRect);

            Rect itemsClip = new Rect(arrangeBounds);
            itemsClip.Height -= rullerHeight;
            itemsClip.Y += _offset.Y;
            ItemsControl.Clip = new RectangleGeometry(itemsClip);



            return retVal;
        }

        public event System.ComponentModel.PropertyChangedEventHandler PropertyChanged;
    }
}
