#ifndef __RULER_H
#define __RULER_H

class Ruler
{
public:
	Ruler(CDC & dc, double minimum, double range, const CRect & rulerRect, int leftOffset);
	~Ruler();
	void Draw();

private:
	void DrawSegments(bool bBig, double step, double leftPoint, double rightPoint, bool bInverse);

private:
	CDC & m_dc;
	double m_min;
	double m_range;
	CRect m_rulerRect;
	int m_leftOffset;

	static CBrush m_rulerBrush;
	static CPen m_rulerPen;
};

#endif
