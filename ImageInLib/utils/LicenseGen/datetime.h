#ifndef __datetime_h__
#define __datetime_h__

long ConvertStringToDate(const char *str);
long ConvertStringToTime(const char *str);
double ConvertStringToDateTime(const char *str);
void ConvertDateToString(long date, char *str);
void ConvertTimeToString(long time, char *str);
void ConvertDateTimeToString(double dateTime, char *str);
double GetCurrentDateTime(void);
double GetLocalDateTime(void);

#endif
