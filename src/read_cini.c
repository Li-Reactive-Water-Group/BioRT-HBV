#include "biort.h"

void ReadCini(const char dir[], int nsub, const chemtbl_struct *chemtbl, rttbl_struct *rttbl,
    subcatch_struct subcatch[])
{
    char            fn[MAXSTRING];
    char            cmdstr[MAXSTRING];
    FILE           *fp;
    int             lno = 0;
    double          dummy[MAXSPS];

    sprintf(fn, "input/%s/cini.txt", dir);
    fp = fopen(fn, "r");

    // Read precipitation concentration
    FindLine(fp, "PRECIPITATION", &lno, cmdstr);
    ReadConc(fp, rttbl->num_stc, chemtbl, &lno, subcatch[0].prcp_conc, dummy);

    // Read upper zone concentration
    FindLine(fp, "UZ", &lno, cmdstr);
    ReadConc(fp, rttbl->num_stc, chemtbl, &lno, subcatch[0].chms[UZ].tot_conc, subcatch[0].chms[UZ].ssa);

    // Read lower zone concentration
    FindLine(fp, "LZ", &lno, cmdstr);
    ReadConc(fp, rttbl->num_stc, chemtbl, &lno, subcatch[0].chms[LZ].tot_conc, subcatch[0].chms[LZ].ssa);

    fclose(fp);
}

void ReadConc(FILE *fp, int num_stc, const chemtbl_struct chemtbl[], int *lno, double tot_conc[], double ssa[])
{
    char            cmdstr[MAXSTRING];
    char            temp_str[MAXSTRING];
    int             ind;
    int             convert;
    int             kspc;

    for (kspc = 0; kspc < num_stc; kspc++)
    {
        NextLine(fp, cmdstr, lno);
        sscanf(cmdstr, "%s", temp_str);
        if (strcmp("pH", temp_str) == 0)
        {
            convert = 1;
        }

        ind = FindChem(temp_str, num_stc, chemtbl);
        if (ind < 0)
        {
            biort_printf(VL_ERROR, "Error finding chemical %s.\n", temp_str);
            exit(EXIT_FAILURE);
        }

        if (chemtbl[ind].itype == MINERAL)
        {
            if (sscanf(cmdstr, "%*s %lf %*s %lf", &tot_conc[ind], &ssa[ind]) !=2)
            {
                biort_printf(VL_ERROR, "Error reading initial condition in %s at Line %d.\n", "cini.txt", *lno);
            }
        }
        else
        {
            if (sscanf(cmdstr, "%*s %lf", &tot_conc[ind]) != 1)
            {
                biort_printf(VL_ERROR, "Error reading initial condition in %s at Line %d.\n", "cini.txt", *lno);
            }
        }

        tot_conc[ind] = (strcmp(chemtbl[ind].name, "pH") == 0 && convert == 1) ?
            pow(10, -tot_conc[ind]) : tot_conc[ind];
    }
}
