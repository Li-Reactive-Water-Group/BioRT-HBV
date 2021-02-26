#include "biort.h"

int ReadParam(const char buffer[], const char keyword[], char type, const char fn[], int lno, void *value)
{
    char            optstr[MAXSTRING];
    int             success = 1;
    int             match;

    if (NULL == value)
    {
        match = sscanf(buffer, "%s", optstr);
        if (strcasecmp(keyword, optstr) != 0)
        {
            biort_printf(VL_ERROR, "Expected keyword \"%s\", detected keyword \"%s\".\n", keyword, optstr);
            success = 0;
        }
        else if (match != 1)
        {
            success = 0;
        }
    }
    else
    {
        switch (type)
        {
            case 'd':
                match = sscanf(buffer, "%s %lf", optstr, (double *)value);
                if (strcasecmp(keyword, optstr) != 0)
                {
                    biort_printf(VL_ERROR, "Expected keyword \"%s\", detected keyword \"%s\".\n", keyword, optstr);
                    success = 0;
                }
                else if (match != 2)
                {
                    success = 0;
                }
                break;
            case 'i':
                match = sscanf(buffer, "%s %d", optstr, (int *)value);
                if (strcasecmp(keyword, optstr) != 0)
                {
                    biort_printf(VL_ERROR, "Expected keyword \"%s\", detected keyword \"%s\".\n", keyword, optstr);
                    success = 0;
                }
                else if (match != 2)
                {
                    success = 0;
                }
                break;
            case 's':
                match = sscanf(buffer, "%s %[^\n]", optstr, (char *)value);
                if (strcasecmp(keyword, optstr) != 0)
                {
                    biort_printf(VL_ERROR, "Expected keyword \"%s\", detected keyword \"%s\".\n", keyword, optstr);
                    success = 0;
                }
                else if (match != 2)
                {
                    success = 0;
                }
                break;
            case 'w':
                match = sscanf(buffer, "%s %s", optstr, (char *)value);
                if (strcasecmp(keyword, optstr) != 0)
                {
                    biort_printf(VL_ERROR, "Expected keyword \"%s\", detected keyword \"%s\".\n", keyword, optstr);
                    success = 0;
                }
                else if (match != 2)
                {
                    success = 0;
                }
                break;
            default:
                biort_printf(VL_ERROR, "Error: Keyword type \'%c\' is not defined.\n", type);
                exit(EXIT_FAILURE);
        }
    }

    if (!success)
    {
        biort_printf(VL_ERROR, "File %s format error at Line %d.\n", fn, lno);
        exit(EXIT_FAILURE);
    }

    return success;
}