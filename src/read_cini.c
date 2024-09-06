#include "biort.h"

void ReadCini(const char dir[], int nsub, const chemtbl_struct *chemtbl, rttbl_struct *rttbl,
    subcatch_struct subcatch[])
{
    char            fn[MAXSTRING];
    char            cmdstr[MAXSTRING];
    FILE           *fp;
    int             lno = 0;
    int             i;
    double          dummy[MAXSPS];

    sprintf(fn, "input/%s/cini.txt", dir);
    fp = fopen(fn, "r");

    // Read precipitation concentration
    FindLine(fp, "PRECIPITATION", &lno, cmdstr);
    ReadConc(fp, rttbl->num_stc, chemtbl, &lno, subcatch[0].prcp_conc, dummy, dummy, dummy, dummy, dummy, dummy);

    // Read surface concentration  2021-05-07
    FindLine(fp, "SF", &lno, cmdstr);
    ReadConc(fp, rttbl->num_stc, chemtbl, &lno, subcatch[0].chms[SURFACE].tot_conc, subcatch[0].chms[SURFACE].ssa,subcatch[0].chms[SURFACE].k_cini, subcatch[0].chms[SURFACE].q10, subcatch[0].chms[SURFACE].sw_thld, subcatch[0].chms[SURFACE].sw_exp, subcatch[0].chms[SURFACE].n_alpha);

    // Read upper zone concentration
    FindLine(fp, "SZ", &lno, cmdstr);
    ReadConc(fp, rttbl->num_stc, chemtbl, &lno, subcatch[0].chms[UZ].tot_conc, subcatch[0].chms[UZ].ssa, subcatch[0].chms[UZ].k_cini, subcatch[0].chms[UZ].q10, subcatch[0].chms[UZ].sw_thld, subcatch[0].chms[UZ].sw_exp, subcatch[0].chms[UZ].n_alpha);

    // Read lower zone concentration
    FindLine(fp, "DZ", &lno, cmdstr);
    ReadConc(fp, rttbl->num_stc, chemtbl, &lno, subcatch[0].chms[LZ].tot_conc, subcatch[0].chms[LZ].ssa, subcatch[0].chms[LZ].k_cini, subcatch[0].chms[LZ].q10, subcatch[0].chms[LZ].sw_thld, subcatch[0].chms[LZ].sw_exp, subcatch[0].chms[LZ].n_alpha);

    //To check if UZ and LZ has same k or rate for each mineral
    for (i = 0; i < rttbl->num_min; i++){
        if (subcatch[0].chms[LZ].k_cini[i + rttbl->num_stc - rttbl->num_min]!=subcatch[0].chms[UZ].k_cini[i + rttbl->num_stc - rttbl->num_min]){
            biort_printf(VL_ERROR, "k for minerals in cini should be same in both UZ and LZ \n");
            exit(EXIT_FAILURE);
        }
    }

    fclose(fp);
}

void ReadConc(FILE *fp, int num_stc, const chemtbl_struct chemtbl[], int *lno, double tot_conc[], double ssa[], double k_cini[], double q10[], double sw_thld[], double sw_exp[], double n_alpha[])
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
            if (sscanf(cmdstr, "%*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf ", &tot_conc[ind], &ssa[ind], &k_cini[ind], &q10[ind], &sw_thld[ind], &sw_exp[ind], &n_alpha[ind]) !=7)
            {
                biort_printf(VL_ERROR, "Error reading initial condition in %s at Line %d.\n", "cini.txt", *lno);
                exit(EXIT_FAILURE);
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
