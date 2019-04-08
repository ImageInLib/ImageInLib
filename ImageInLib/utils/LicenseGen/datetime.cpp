/*
 * Calendar time to Julian date conversions.
 * Julian date is commonly used in astronomical applications,
 *	since it is numerically accurate and computationally simple.
 * The algorithms here will accurately convert between Julian day
 *	and calendar date for all non-negative Julian days
 *	(i.e. from Nov 23, -4713 on).
 *
 * Ref: Explanatory Supplement to the Astronomical Almanac, 1992.
 *	University Science Books, 20 Edgehill Rd. Mill Valley CA 94941.
 *
 * Use the algorithm by Henry Fliegel.
 */
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#ifdef _WIN32
#include "windows.h"
#endif
#include "time.h"

int date2j(int y, int m, int d)
{
	int m12 = (m - 14) / 12;
	return ((1461 * (y + 4800 + m12)) / 4 + (367 * (m - 2 - 12 * (m12))) / 12
			- (3 * ((y + 4900 + m12) / 100)) / 4 + d - 32075);
}

void j2date(int jd, int *year, int *month, int *day)
{
	int	j, y, m, d;
	int i, l, n;

	l = jd + 68569;
	n = (4 * l) / 146097;
	l -= (146097 * n + 3) / 4;
	i = (4000 * (l + 1)) / 1461001;
	l += 31 - (1461 * i) / 4;
	j = (80 * l) / 2447;
	d = l - (2447 * j) / 80;
	l = j / 11;
	m = (j + 2) - (12 * l);
	y = 100 * (n - 49) + i + l;

	*year = y;
	*month = m;
	*day = d;
}

int j2day(int date)
{
	return (date + 1) % 7;
}


long ConvertStringToDate(const char *str)
{
	char year[5];
	char month[3];
	char day[3];
	if (*str == 0)
		return -1;
	if (strspn(str, "0123456789.") != strlen(str))
		return -1;
	if (strlen(str) == 8 && strchr(str, '.') == NULL ) {
		year[4] = 0; memmove(year, str, 4);
		month[2] = 0; memmove(month, str + 4, 2);
		day[2] = 0; memmove(day, str + 6, 2);
	} else if (strlen(str) == 10 && strchr(str, '.') == str+4 && strchr(str, '.') == str+6) {
		year[4] = 0; memmove(year, str, 4);
		month[2] = 0; memmove(month, str + 5, 2);
		day[2] = 0; memmove(day, str + 8, 2);
	} else {
		return -1;
	}
	if (strchr(year, '.') != NULL || strchr(month, '.') != NULL || strchr(day, '.') != NULL)
		return -1;
	return date2j(atoi(year), atoi(month), atoi(day));
}

long ConvertStringToTime(const char *str)
{
	char hour[3];
	char min[3];
	char sec[3];
	const char *frac = "";
	if (*str == 0)
		return -1;
	if (strspn(str, "0123456789.:") != strlen(str))
		return -1;
	if (strchr(str,':') != NULL) {
		hour[0] = str[0]; hour[1] = str[1]; hour[2] = 0;
		if (str[2] != 0) {
			if (str[2] != ':')
				return 0;
			min[0] = str[3]; min[1] = str[4]; min[2] = 0;
			if (str[5] != 0) {
				if (str[5] != ':')
					return 0;
				sec[0] = str[6]; sec[1] = str[7]; sec[2] = 0;
				if (str[8] == '.')
					frac = str+6;
			} else {
				*sec = 0;
			}
		} else {
			*min = 0;
			*sec = 0;
		}
	} else {
		hour[0] = str[0]; hour[1] = str[1]; hour[2] = 0;
		min[0] = str[2]; min[1] = str[3]; min[2] = 0;
		sec[0] = str[4]; sec[1] = str[5]; sec[2] = 0;
		if (str[6] == '.')
			frac = str+6;
	}
	return long((atoi(hour) * 3600 + atoi(min) * 60 + atoi(sec) + atof(frac)) * 10000);
}


double ConvertStringToDateTime(const char *str)
{
	char year[5] = "";
	char month[3] = "";
	char day[3] = "";
	char hour[3] = "";
	char min[3] = "";
	char sec[3] = "";
	const char *frac = "";
	int len;
	double dateTime;

	if (strspn(str, "0123456789.+-") != strlen(str))
		return -1;
	len = strlen(str);

	if (len < 4)
		return -1;
	year[4] = 0; memmove(year, str, 4);
	if (strspn(year, "0123456789") != 4)
		return -1;
	str += 4; len -= 4;
	if (len && *str != '+' && *str != '-') {
		if (len < 2)
			return -1;
		month[0] = str[0]; month[1] = str[1]; month[2] = 0;
		if (strspn(month, "0123456789") != 2)
			return -1;
		str += 2; len -= 2;
		if (len && *str != '+' && *str != '-') {
			if (len < 2)
				return -1;
			day[0] = str[0]; day[1] = str[1]; day[2] = 0;
			if (strspn(day, "0123456789") != 2)
				return -1;
			str += 2; len -= 2;
			if (len && *str != '+' && *str != '-') {
				if (len < 2)
					return -1;
				hour[0] = str[0]; hour[1] = str[1]; hour[2] = 0;
				if (strspn(hour, "0123456789") != 2)
					return -1;
				str += 2; len -= 2;
				if (len && *str != '+' && *str != '-') {
					if (len < 2)
						return -1;
					min[0] = str[0]; min[1] = str[1]; min[2] = 0;
					if (strspn(min, "0123456789") != 2)
						return -1;
					str += 2; len -= 2;
					if (len && *str != '+' && *str != '-') {
						if (len < 2)
							return -1;
						sec[0] = str[0]; sec[1] = str[1]; sec[2] = 0;
						if (strspn(sec, "0123456789") != 2)
							return -1;
						str += 2; len -= 2;
						if (*str == '.')
							frac = str;
					}
				}
			}
		}
	}
	dateTime = double(date2j(atoi(year), atoi(month), atoi(day))) * 86400.0 +
		atoi(hour) * 3600 + atoi(min) * 60 + atoi(sec) + atof(frac);
	if (len == 5 && (*str == '+' || *str == '-')) {
		long offset;
		offset = str[1] * 36000 + str[2] * 3600 + str[3] * 600 + str[4] * 60;
		if (*str == '+')
			dateTime += offset;
		else 
			dateTime -= offset;
	}
	return dateTime;
}

void ConvertDateToString(long date, char *str)
{
	int year, month, day;
	if (date == -1) {
		*str = 0;
		return;
	}
	j2date(date, &year, &month, &day);
	sprintf(str, "%04d%02d%02d", year, month, day);
}

void ConvertTimeToString(long time, char *str)
{
	int hour, min, sec, frac;
	if (time == -1) {
		*str = 0;
		return;
	}
	frac = time % 10000;
	time /= 10000;
	hour = time / 3600;
	min = (time % 3600) / 60;
	sec = time % 60;
	sprintf(str, "%02d%02d%02d.%04d", hour, min, sec, frac);
}

void ConvertDateTimeToString(double dateTime, char *str)
{
	int year, month, day, hour, min, sec, frac;
	if (dateTime == -1) {
		*str = 0;
		return;
	}
	j2date(int(dateTime / 86400), &year, &month, &day);
	double time = fmod(dateTime, 86400);
	hour = int(time / 3600);
	min = int(fmod(time, 3600) / 60);
	sec = int(fmod(time, 60));
	frac = int((time - floor(time)) * 1000000);
	sprintf(str, "%04d%02d%02d%02d%02d%02d.%06d", year, month, day, hour, min, sec, frac);
}


double GetCurrentDateTime(void)
{
#ifdef _WIN32
	SYSTEMTIME t;
	GetSystemTime(&t);
	return double(date2j(t.wYear, t.wMonth, t.wDay)) * 86400 +
		t.wHour * 3600 + t.wMinute * 60 + t.wSecond + t.wMilliseconds/1000.0;
#else
	time_t t = time(NULL);
	struct tm *utc = gmtime(&t);
	return double(date2j(utc->tm_year + 1900, utc->tm_mon + 1, utc->tm_mday)) * 86400 +
		utc->tm_hour * 3600 + utc->tm_min * 60 + utc->tm_sec;
#endif
}

double GetLocalDateTime(void)
{
#ifdef _WIN32
	SYSTEMTIME t;
	GetLocalTime(&t);
	return double(date2j(t.wYear, t.wMonth, t.wDay)) * 86400 +
		t.wHour * 3600 + t.wMinute * 60 + t.wSecond + t.wMilliseconds/1000.0;
#else
	time_t t = time(NULL);
	struct tm *utc = localtime(&t);
	return double(date2j(utc->tm_year + 1900, utc->tm_mon + 1, utc->tm_mday)) * 86400 +
		utc->tm_hour * 3600 + utc->tm_min * 60 + utc->tm_sec;
#endif
}
