#include "biort.h"

void PrintHeader(FILE *fp, const rttbl_struct *rttbl, const chemtbl_struct chemtbl[])
{
    int             kspc;
    char            chemn[MAXSTRING];
    char            tempstr[MAXSTRING];

    // Soil concentration file header
    fprintf(fp, "%-15s",  "TIME");

    for (kspc = 0; kspc < rttbl->num_spc; kspc++)
    {
        Unwrap(chemtbl[kspc].name, chemn);

        sprintf(tempstr, "%s_inf", chemn);
        fprintf(fp, "\t%-15s", tempstr);

    }
    for (kspc = 0; kspc < rttbl->num_stc + rttbl->num_ssc; kspc++)
    {
        Unwrap(chemtbl[kspc].name, chemn);

        sprintf(tempstr, "%s_UZ", chemn);
        fprintf(fp, "\t%-15s", tempstr);

    }
    for (kspc = 0; kspc < rttbl->num_stc + rttbl->num_ssc; kspc++)
    {
        Unwrap(chemtbl[kspc].name, chemn);

        sprintf(tempstr, "%s_LZ", chemn);
        fprintf(fp, "\t%-15s", tempstr);

    }
    for (kspc = 0; kspc < rttbl->num_stc + rttbl->num_ssc; kspc++)
    {
        Unwrap(chemtbl[kspc].name, chemn);

        sprintf(tempstr, "%s_riv", chemn);
        fprintf(fp, "\t%-15s", tempstr);

    }
    fprintf(fp, "\n");

    fflush(fp);
}

void PrintDailyResults(FILE *fp, int step, int nsub, const rttbl_struct *rttbl, const subcatch_struct subcatch[])
{
    int             kspc;
    int             ksub;

    // Soil concentration file header
    fprintf(fp, "%-15d",  step);

    for (ksub = 0; ksub < nsub; ksub++)
    {
        for (kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch[ksub].chms[SNSM].prim_conc[kspc]);
        }
        for (kspc = 0; kspc < rttbl->num_stc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch[ksub].chms[UZ].prim_conc[kspc]);
        }
        for (kspc = 0; kspc < rttbl->num_ssc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch[ksub].chms[UZ].sec_conc[kspc]);
        }
        for (kspc = 0; kspc < rttbl->num_stc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch[ksub].chms[LZ].prim_conc[kspc]);
        }
        for (kspc = 0; kspc < rttbl->num_ssc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch[ksub].chms[LZ].sec_conc[kspc]);
        }
        for (kspc = 0; kspc < rttbl->num_stc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch[ksub].chms[STREAM].prim_conc[kspc]);
        }
        for (kspc = 0; kspc < rttbl->num_ssc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch[ksub].chms[STREAM].sec_conc[kspc]);
        }
    }
    fprintf(fp, "\n");

    fflush(fp);
}