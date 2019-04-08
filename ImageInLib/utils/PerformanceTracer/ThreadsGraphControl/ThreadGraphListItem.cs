using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace ThreadsGraphControl
{
    [TemplatePart(Name = ThreadGraphListItem.PART_Header, Type = typeof(ContentPresenter))]
    [TemplatePart(Name = ThreadGraphListItem.PART_Content, Type = typeof(ContentPresenter))]

    public class ThreadGraphListItem : ListBoxItem
    {
        public const string PART_Header = "PART_Header";
        public const string PART_Content = "PART_Content";

        static ThreadGraphListItem()
        {
            DefaultStyleKeyProperty.OverrideMetadata(typeof(ThreadGraphListItem), new FrameworkPropertyMetadata(typeof(ThreadGraphListItem)));
        }

        public static readonly DependencyProperty ItemHeaderTemplateProperty = DependencyProperty.Register(
            "ItemHeaderTemplate",
            typeof(DataTemplate),
            typeof(ThreadGraphListItem),
            new PropertyMetadata((s, e) => ((ThreadGraphListItem)s).ApplyItemHeaderTemplate())
            );

        public DataTemplate ItemHeaderTemplate
        {
            get { return (DataTemplate)GetValue(ItemHeaderTemplateProperty); }
            set { SetValue(ItemHeaderTemplateProperty, value); }
        }

        public void ApplyItemHeaderTemplate()
        {
        }

        ThreadsGraph _graph = null;
        public ThreadGraphListItem(ThreadsGraph graph)
        {
            _graph = graph;
            _graph.PropertyChanged += (sender, parameter) =>
                {
                    if (parameter.PropertyName == "XOffset")
                        InvalidateArrange();
                };

            ItemHeaderTemplate = _graph.ItemHeaderTemplate;
        }

        protected override Size MeasureOverride(Size availableSize)
        {
            FrameworkElement Header = GetTemplateChild(PART_Header) as FrameworkElement;
            FrameworkElement Content = GetTemplateChild(PART_Content) as FrameworkElement;

            Header.Measure(availableSize);
            Content.Measure(availableSize);

            Size HeaderSize = Header.DesiredSize;
            Size ContentSize = Content.DesiredSize;

            Size result = new Size(HeaderSize.Width+ContentSize.Width, HeaderSize.Height);

            if(result.Height < ContentSize.Height)
                result.Height = ContentSize.Height;
            return result;
        }

        protected override Size ArrangeOverride(Size arrangeBounds)
        {
            Size retVal = base.ArrangeOverride(arrangeBounds);

            FrameworkElement Header = GetTemplateChild(PART_Header) as FrameworkElement;
            FrameworkElement Content = GetTemplateChild(PART_Content) as FrameworkElement;

            Rect headerRect = new Rect(arrangeBounds);
            headerRect.Width = Header.DesiredSize.Width;
            Header.Arrange(headerRect);

            Rect contentRect = new Rect(arrangeBounds);
            contentRect.Width = Content.DesiredSize.Width;
            contentRect.X += headerRect.Width;
            contentRect.X -= _graph.HorizontalOffset;

            Rect contentClip = new Rect(arrangeBounds);
            contentClip.Width -= headerRect.Width;
            contentClip.X += _graph.HorizontalOffset;
            Content.Clip = new RectangleGeometry(contentClip);
            
            Content.Arrange(contentRect);

            return retVal;
        }

    }
}
