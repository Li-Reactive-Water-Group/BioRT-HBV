#include "biort.h"

void PrintHeader(FILE *fp, int transpt, const rttbl_struct *rttbl, const chemtbl_struct chemtbl[])
{
    int             kspc;
    char            chemn[MAXSTRING];
    char            tempstr[MAXSTRING];

    // Soil concentration file header
    fprintf(fp, "%-15s",  "TIME");

    if (transpt == KIN_REACTION)
    {
        for (kspc = 0; kspc < rttbl->num_spc; kspc++)  // 2021-06-29
        {
            Unwrap(chemtbl[kspc].name, chemn);

            sprintf(tempstr, "%s_SNOW", chemn);
            fprintf(fp, "\t%-15s", tempstr);

        }

        for (kspc = 0; kspc < rttbl->num_spc; kspc++)  // 2021-06-29
        {
            Unwrap(chemtbl[kspc].name, chemn);

            sprintf(tempstr, "%s_SURFACE", chemn);
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

        for (kspc = 0; kspc < rttbl->num_min; kspc++)  // 2021-05-14
        {
            Unwrap(chemtbl[kspc + rttbl->num_stc - rttbl->num_min].name, chemn);

            sprintf(tempstr, "%s_rate_UZ", chemn);
            fprintf(fp, "\t%-23s", tempstr);
        }

        //for (kspc = 0; kspc < rttbl->num_min; kspc++)
        //{
            //Unwrap(chemtbl[kspc + rttbl->num_stc - rttbl->num_min].name, chemn);

            //sprintf(tempstr, "%s_rate_UZ", chemn);
            //fprintf(fp, "\t%-23s", tempstr);
        //}

        for (kspc = 0; kspc < rttbl->num_min; kspc++)
        {
            Unwrap(chemtbl[kspc + rttbl->num_stc - rttbl->num_min].name, chemn);

            sprintf(tempstr, "%s_rate_LZ", chemn);
            fprintf(fp, "\t%-23s", tempstr);
        }
    }
    else    // In transport mode, only print primary species
    {

        for (kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            Unwrap(chemtbl[kspc].name, chemn);

            sprintf(tempstr, "%s_SNOW", chemn);
            fprintf(fp, "\t%-15s", tempstr);

        }
        for (kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            Unwrap(chemtbl[kspc].name, chemn);

            sprintf(tempstr, "%s_SURFACE", chemn);
            fprintf(fp, "\t%-15s", tempstr);

        }
        for (kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            Unwrap(chemtbl[kspc].name, chemn);

            sprintf(tempstr, "%s_UZ", chemn);
            fprintf(fp, "\t%-15s", tempstr);

        }
        for (kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            Unwrap(chemtbl[kspc].name, chemn);

            sprintf(tempstr, "%s_LZ", chemn);
            fprintf(fp, "\t%-15s", tempstr);

        }
        for (kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            Unwrap(chemtbl[kspc].name, chemn);

            sprintf(tempstr, "%s_riv", chemn);
            fprintf(fp, "\t%-15s", tempstr);

        }
    }
    fprintf(fp, "\n");

    // UNITS
    fprintf(fp, "%-15s",  "YYYYMMDD");

    if (transpt == KIN_REACTION)
    {
        for (kspc = 0; kspc < 2 * (rttbl->num_spc); kspc++)   // 2021-09-27, SNOW +SURFACE
        {
            fprintf(fp, "\t%-15s", "mol/L");
        }
        for (kspc = 0; kspc < 3 * (rttbl->num_stc + rttbl->num_ssc); kspc++)   // UZ + LZ + STREAM
        {
            fprintf(fp, "\t%-15s", "mol/L");
        }
        for (kspc = 0; kspc < 2 * rttbl->num_min; kspc++)  // UZ + LZ
        {
            fprintf(fp, "\t%-23s", "mol/m2/day");
        }
    }
    else    // In transport mode, only print primary species
    {
        for (kspc = 0; kspc < 5 * rttbl->num_spc; kspc++)  // 2021-06-29
        {
            fprintf(fp, "\t%-15s", "mol/L");
        }
    }
    fprintf(fp, "\n");

    fflush(fp);
}

void PrintDailyResults(FILE *fp, int transpt, int step, int nsub, const rttbl_struct *rttbl,
    const subcatch_struct subcatch[])
{
    int             kspc;
    int             ksub;

    // Soil concentration file header
    fprintf(fp, "%-15d",  step);

    for (ksub = 0; ksub < nsub; ksub++)
    {

        if (transpt == KIN_REACTION)
        {

            for (kspc = 0; kspc < rttbl->num_spc; kspc++)   // 2021-06-29
            {
                fprintf(fp, "\t%-15lg", subcatch[ksub].chms[SNOW].prim_conc[kspc]);
            }
            for (kspc = 0; kspc < rttbl->num_spc; kspc++)   // 2021-06-29
            {
                fprintf(fp, "\t%-15lg", subcatch[ksub].chms[SURFACE].prim_conc[kspc]);
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
            //for (kspc = 0; kspc < rttbl->num_min; kspc++)  // 2021-05-14
            //{
            //    fprintf(fp, "\t%-23lg", subcatch[ksub].react_rate[SURFACE][kspc]);
            //}
            for (kspc = 0; kspc < rttbl->num_min; kspc++)
            {
                fprintf(fp, "\t%-23lg", subcatch[ksub].react_rate[UZ][kspc]);
            }
            for (kspc = 0; kspc < rttbl->num_min; kspc++)
            {
                fprintf(fp, "\t%-23lg", subcatch[ksub].react_rate[LZ][kspc]);
            }
        }
        else    // In transport mode, only print primary species
        {
            for (kspc = 0; kspc < rttbl->num_spc; kspc++)  // 2021-06-29
            {
                fprintf(fp, "\t%-15lg", subcatch[ksub].chms[SNOW].prim_conc[kspc]);
            }
            for (kspc = 0; kspc < rttbl->num_spc; kspc++)  // 2021-06-29
            {
                fprintf(fp, "\t%-15lg", subcatch[ksub].chms[SURFACE].prim_conc[kspc]);
            }
            for (kspc = 0; kspc < rttbl->num_spc; kspc++)
            {
                fprintf(fp, "\t%-15lg", subcatch[ksub].chms[UZ].prim_conc[kspc]);
            }
            for (kspc = 0; kspc < rttbl->num_spc; kspc++)
            {
                fprintf(fp, "\t%-15lg", subcatch[ksub].chms[LZ].prim_conc[kspc]);
            }
            for (kspc = 0; kspc < rttbl->num_spc; kspc++)
            {
                fprintf(fp, "\t%-15lg", subcatch[ksub].chms[STREAM].prim_conc[kspc]);
            }
        }
    }
    fprintf(fp, "\n");

    fflush(fp);
}
