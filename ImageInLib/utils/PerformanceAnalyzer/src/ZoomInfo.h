#pragma once

static double _100ns = 1e-7;

namespace DataModel {
	class CZoomInfo
	{
	public:
		CZoomInfo();
		~CZoomInfo();

	protected:
		double m_dOriginalZoomCoef;
		double m_dActualZoomCoef;
		double m_dZoomRangeStart;
		double m_dZoomRangeWidth;
		double m_dZoomRangeMinWidth;
		double m_dDataWidth;
		double m_dDataShift;

	public:
		void Reset();

		double DataWidth() const { return m_dDataWidth; }
		double DataShift() const { return m_dDataShift; }
		bool SetDataShift(double dDataShift);

		double ZoomRangeStart() const { return m_dZoomRangeStart + m_dDataShift; }
		double ZoomRange() const { return m_dZoomRangeWidth; }
		double ZoomMinRange() const { return m_dZoomRangeMinWidth; }
		double SetZoomRange(double dRangeStart, double dRangeWidth);
		double SetMinZoomRange(double dMinRange);

		void UpdateZoomCoefs(double dDataWndWidth, double dDataWidth);

		double ShiftZoomRangeToPosition(double dStartPosition)
		{
			return SetZoomRange(dStartPosition, ZoomRange());
		}


		double TimePositionToPixelPosition(double dTimePosition) const
		{
			return dTimePosition * m_dActualZoomCoef;
		}

		double PixelPositionToTimePosition(double dPixelPosition) const
		{
			return dPixelPosition / m_dActualZoomCoef;
		}

	protected:
		bool UpdateActualZoom();

	};
}