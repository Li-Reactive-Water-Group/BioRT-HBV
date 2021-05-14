#include "biort.h"

void ReadHbvResults(const char dir[], int nsub, int *nsteps, int *steps[], subcatch_struct subcatch[])
{
    FILE           *fp;
    char            fn[MAXSTRING];
    char            cmdstr[MAXSTRING];
    int             ksub;
    int             kstep;
    int             lno = 0;

    sprintf(fn, "input/%s/Results.txt", dir);
    fp = fopen(fn, "r");

    *nsteps = CountLines(fp, cmdstr, 0) - 1;

    rewind(fp);

    *steps = (int *)malloc(*nsteps * sizeof(int));

    for (ksub = 0; ksub < nsub; ksub++)
    {
        subcatch[ksub].ws = (double **)malloc(*nsteps * sizeof(double *));
        subcatch[ksub].q = (double **)malloc(*nsteps * sizeof(double *));
        subcatch[ksub].tmp = (double *)malloc(*nsteps * sizeof(double));

        for (kstep = 0; kstep < *nsteps; kstep++)
        {
            subcatch[ksub].ws[kstep] = (double *)malloc(NWS * sizeof(double));
            subcatch[ksub].q[kstep] = (double *)malloc(NQ * sizeof(double));
        }
    }

    for (ksub = 0; ksub < nsub; ksub++)
    {
        NextLine(fp, cmdstr, &lno);     // Skip header line

        for (kstep = 0; kstep < *nsteps; kstep++)
        {
            double          snow, sm;

            if (ksub == 0)
            {
                fscanf(fp, "%d", &((*steps)[kstep]));    // Read model steps
            }
            else
            {
                fscanf(fp, "%*d");
            }

            fscanf(fp, "%*f %*f");  // Skip "Qsim" and "Qobs"
            fscanf(fp, "%lf %lf", &subcatch[ksub].q[kstep][PRECIP], &subcatch[ksub].tmp[kstep]);    // Read precip. and
                                                                                                    // air temperature
            fscanf(fp, "%*lf %*lf");   // Skip "AET" and "PET"
            fscanf(fp, "%lf %*lf %lf", &snow, &sm);     // Read snow and soil moisture
            subcatch[ksub].ws[kstep][SURFACE] = snow + sm;  // 2021-05-14
            fscanf(fp, "%lf", &subcatch[ksub].q[kstep][RECHG]);  // Read recharge
            fscanf(fp, "%*lf %*lf");   // Skip upper and lower zone storages
            fscanf(fp, "%lf %lf %lf", &subcatch[ksub].q[kstep][Q0],
                &subcatch[ksub].q[kstep][Q1], &subcatch[ksub].q[kstep][Q2]);  // Read Q0, Q1, and Q2
            fscanf(fp, "%*lf %*lf");   // Skip "Qsim_rain" and "Qsim_snow"

            // In the upper zone, HBV first calculates percolation to the lower zone.
            //   percolation = MIN(perc0, SUZ0 + recharge).                         (1)
            // If there is water left in the upper zone, calculate Q0 and Q1:
            //   Q0 = MAX(SUZ0 + recharge - percoloation - UZL, 0) * K0,            (2)
            //   Q1 = (SUZ0 + recharge - percoloation) * K1,                        (3)
            //   SUZ = SUZ0 + recharge - percolation - Q0 - Q1,                     (4)
            // i.e.,
            //   SUZ = Q1 / K1 - Q0 - Q1.                                           (5)
            // In the lower zone, percoloation is first added to the lower zone storage. Then Q2 is calculated:
            //   Q2 = (SLZ0 + percolation) * K2,                                    (6)
            //   SLZ = SLZ0 + percololation - Q2,                                   (7)
            // i.e.,
            //   SLZ = Q2 / k2 - Q2.                                                (8)
            //
            // HBV Light outputs SUZ and SLZ with low precision (one decimal digit), and does not output percolation
            // rate. Therefore, SUZ, SLZ, and percolation rate are calculated using Equations (1), (4), and (8).

            subcatch[ksub].ws[kstep][UZ] = subcatch[ksub].q[kstep][Q1] / subcatch[ksub].k1 -
                (subcatch[ksub].q[kstep][Q0] + subcatch[ksub].q[kstep][Q1]);
            subcatch[ksub].ws[kstep][LZ] = subcatch[ksub].q[kstep][Q2] / subcatch[ksub].k2 -
                subcatch[ksub].q[kstep][Q2];
            subcatch[ksub].q[kstep][PERC] = (kstep == 0) ?
                0.0 : MIN(subcatch[ksub].perc, subcatch[ksub].ws[kstep - 1][UZ] + subcatch[ksub].q[kstep][RECHG]);
        }

        // Add residual moisture
        for (kstep = 0; kstep < *nsteps; kstep++)
        {
            subcatch[ksub].ws[kstep][SURFACE] += STORAGE_MIN;   // 2021-05-14
            subcatch[ksub].ws[kstep][UZ] += subcatch[ksub].res_uz;
            subcatch[ksub].ws[kstep][LZ] += subcatch[ksub].res_lz;
        }

        fclose(fp);
    }
}

void ReadHbvParam(const char dir[], int nsub, subcatch_struct subcatch[])
{
    FILE           *fp;
    char            fn[MAXSTRING];
    char            cmdstr[MAXSTRING];
    char            tag[MAXSTRING];
    int             lno = 0;
    int             ksub;
    double          value;

    sprintf(fn, "input/%s/Parameter.xml", dir);
    fp = fopen(fn, "r");

    biort_printf(VL_NORMAL, "\nHBV MODEL PARAMETERS\n");

    for (ksub = 0; ksub < nsub; ksub++)
    {
        FindLine(fp, "<SubCatchmentParameters>", &lno, cmdstr);

        while (!feof(fp))
        {
            NextLine(fp, cmdstr, &lno);
            ParseLine(cmdstr, tag, &value);
            if (strcmp(tag, "PERC") == 0)
            {
                subcatch[ksub].perc = value;
                biort_printf(VL_NORMAL, "  Percolation rate is %.2lf mm day-1\n", subcatch[ksub].perc);
            }
            else if (strcmp(tag, "K1") == 0)
            {
                subcatch[ksub].k1 = value;
                biort_printf(VL_NORMAL, "  K1 is %.2lf day-1\n", subcatch[ksub].k1);
            }
            else if (strcmp(tag, "K2") == 0)
            {
                subcatch[ksub].k2 = value;
                biort_printf(VL_NORMAL, "  K2 is %.2lf day-1\n", subcatch[ksub].k2);
            }
            else if (strcmp(tag, "MAXBAS") == 0)
            {
                subcatch[ksub].maxbas = value;
                biort_printf(VL_NORMAL, "  Routing parameter is %.2lf\n", subcatch[ksub].maxbas);
                break;
            }
        }
    }

    fclose(fp);
}

void ParseLine(const char cmdstr[], char tag[], double *value)
{
    char           *temp;
    int             tag_index1, tag_index2, tag_index3;
    int             k;

    temp = strchr(cmdstr, '<');
    tag_index1 = (int)(temp - cmdstr);

    temp = strchr(cmdstr, '>');
    tag_index2 = (int)(temp - cmdstr);

    temp = strchr(temp, '<');
    tag_index3 = (int)(temp - cmdstr);

    for (k = tag_index1 + 1; k < tag_index2; k++)
    {
        tag[k - tag_index1 - 1] = cmdstr[k];
    }
    tag[tag_index2 - tag_index1 - 1] = '\0';

    for (k = tag_index2 + 1; k < tag_index3; k++)
    {
        temp[k - tag_index2 - 1] = cmdstr[k];
    }
    temp[tag_index3 - tag_index2 - 1] = '\0';
    sscanf(temp, "%lf", value);
}
