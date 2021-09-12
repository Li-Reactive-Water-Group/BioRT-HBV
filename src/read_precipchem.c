#include "biort.h"

void ReadPrecipChem(const char dir[], int nsub, int *nsteps, int *steps[], subcatch_struct subcatch[], int num_stc, const chemtbl_struct chemtbl[],int mode)
{
    FILE           *fp;
    char            fn[MAXSTRING];
    char            cmdstr[MAXSTRING];
    char            temp_str[MAXSTRING];
    int             ksub;
    int             kstep;
    int             kspc;
    int             pH_index = 0;
    int             pH_convert = 0;
    int             ind;
    int             ntime;

    ntime = *nsteps;

    if (mode == 0)
    {
        sprintf(fn, "input/%s/precipchem.txt", dir);
        biort_printf(VL_NORMAL, "\nREADING TIME-SERIES PRECIPITATION CHEMISTRY from \"precipchem.txt\"\n");
    }

    else if (mode == 1)

    {
        sprintf(fn, "input/%s/Numexp_precipchem.txt", dir);
        biort_printf(VL_NORMAL, "\nREADING TIME-SERIES PRECIPITATION CHEMISTRY from \"Numexp_precipchem.txt\"\n");
    }

    fp = fopen(fn, "r");

    *nsteps = CountLines(fp, cmdstr, 0) - 1;

    rewind(fp);

    if (ntime != *nsteps & mode == 1){
        biort_printf(VL_ERROR,"\nNumber of time steps in \"Numexp_precipchem.txt\" should be same as in \"Numexp_Results.txt\" file.\n");
        exit(EXIT_FAILURE);
    }
    *steps = (int *)malloc(*nsteps * sizeof(int));

    for (ksub = 0; ksub < nsub; ksub++)
    {
        subcatch[ksub].prcp_conc_time = (double **)malloc(*nsteps * sizeof(double *));

        for (kstep = 0; kstep < *nsteps; kstep++)
        {
            subcatch[ksub].prcp_conc_time[kstep] = (double *)malloc(num_stc * sizeof(double));
        }
    }

    for (ksub = 0; ksub < nsub; ksub++)
    {

        // read header to locate pH position
        for (kspc = 0; kspc < num_stc + 1; kspc++)  // add one more column of date
        {
            fscanf(fp, "%s", temp_str);

            if (strcmp("pH", temp_str) == 0)
            {
                pH_convert = 1;
                pH_index = kspc - 1;
            }

            // also check chemical species, 0629
            if (kspc > 0)
            {
                ind = FindChem(temp_str, num_stc, chemtbl);
                if (ind < 0)
                {
                    biort_printf(VL_ERROR, "Error finding chemical %s.\n", temp_str);
                    exit(EXIT_FAILURE);
                }
            }

        }

        for (kstep = 0; kstep < *nsteps; kstep++)
        {

            if (ksub == 0)
            {
                fscanf(fp, "%d", &((*steps)[kstep]));    // Read model steps
            }
            else
            {
                fscanf(fp, "%*d");
            }

            for (kspc = 0; kspc < num_stc; kspc++)  // Read precipitation chemistry
            {
                if (kspc == pH_index && pH_convert == 1)
                {
                    fscanf(fp, "%lf", &subcatch[ksub].prcp_conc_time[kstep][kspc]);
                    //printf("  step = %d, converting time-series precipitation pH (%lf) to ", kstep, subcatch[ksub].prcp_conc_time[kstep][kspc]);
                    subcatch[ksub].prcp_conc_time[kstep][kspc] = pow(10, -subcatch[ksub].prcp_conc_time[kstep][kspc]);
                    //printf("H+ concentration (%lf) \n", subcatch[ksub].prcp_conc_time[kstep][kspc]);
                }
                else
                {
                    fscanf(fp, "%lf", &subcatch[ksub].prcp_conc_time[kstep][kspc]);
                }
            }

        }

        fclose(fp);
    }
}
