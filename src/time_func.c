#include "biort.h"

// Adapted from https://www.geeksforgeeks.org/find-number-of-days-between-two-given-dates/

// To store number of days in all months from January to Dec.
const int month_days[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

// Count number of leap years before the given date (from AD 0)
int CountLeapYears(int y, int m)
{
    int             years = y;

    // Check if the current year needs to be considered for the count of leap years or not
    if (m <= 2)
    {
        years--;
    }

    // A year is a leap year if it is a multiple of 4, multiple of 400 and not a multiple of 100.
    return years / 4 - years / 100 + years / 400;
}

// This function returns number of days between two given dates
int GetDifference(int ymd1, int ymd2)
{
    int             y1, m1, d1;
    int             y2, m2, d2;
    int             n1, n2;
    int             km;

    y1 = (int)(ymd1 / 10000);
    m1 = (int)((ymd1 - y1 * 10000) / 100);
    d1 = ymd1 - y1 * 10000 - m1 * 100;

    y2 = (int)(ymd2 / 10000);
    m2 = (int)((ymd2 - y2 * 10000) / 100);
    d2 = ymd2 - y2 * 10000 - m2 * 100;

    // Count total number of days before first date
    // Initialize count using years and day
    n1 = y1 * 365 + d1;

    // Add days for months in given date
    for (km = 0; km < m1 - 1; km++)
    {
        n1 += month_days[km];
    }

    // Since every leap year is of 366 days, add a day for every leap year
    n1 += CountLeapYears(y1, m1);

    // Similarly, count total number of days before second date
    n2 = y2 * 365 + d2;

    for (km = 0; km < m2 - 1; km++)
    {
        n2 += month_days[km];
    }
    n2 += CountLeapYears(y2, m2);

    // Return difference between two counts
    return n2 - n1;
}
