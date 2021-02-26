#include "biort.h"

int roundi(double x)
{
    return (int)((x < 0.0) ? x - 0.5 : x + 0.5);
}

void Wrap(char *str)
{
    char            word[MAXSTRING];

    sprintf(word, "'%s'", str);
    strcpy(str, word);
}

void Unwrap(const char wrapped_str[], char str[])
{
    int             i, j = 0;

    for (i = 0; i < (int)strlen(wrapped_str); i++)
    {
        if (wrapped_str[i] != '\'')
        {
            str[j] = wrapped_str[i];
            j++;
        }
    }

    str[j] = '\0';
}

void FreeStruct(int nsub, int *steps[], subcatch_struct subcatch[])
{
    int             ksub;

    free(*steps);

    for (ksub = 0; ksub < nsub; ksub++)
    {
        free(subcatch[ksub].ws[ksub]);
        free(subcatch[ksub].q[ksub]);
    }

    free(subcatch[ksub].ws);
    free(subcatch[ksub].q);
}
