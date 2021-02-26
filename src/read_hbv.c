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
            fscanf(fp, "%lf", &subcatch[ksub].q[kstep][PRECIP]);     // Read precipitation
            fscanf(fp, "%*lf %*lf %*lf");   // Skip "Temperature", "AET", and "PET"
            fscanf(fp, "%lf %*lf %lf", &snow, &sm);     // Read snow and soil moisture
            subcatch[ksub].ws[kstep][SNSM] = snow + sm;
            fscanf(fp, "%lf", &subcatch[ksub].q[kstep][RECHG]);  // Read recharge
            fscanf(fp, "%lf %*lf", &subcatch[ksub].ws[kstep][UZ]);   // Read upper zone storages
            fscanf(fp, "%lf %lf %lf", &subcatch[ksub].q[kstep][Q0],
                &subcatch[ksub].q[kstep][Q1], &subcatch[ksub].q[kstep][Q2]);  // Read Q0, Q1, and Q2
            fscanf(fp, "%*lf %*lf");   // Skip "Qsim_rain" and "Qsim_snow"

            // HBV output of storage terms have low precisions (one decimal digit). So we calculate lower zone storage
            // from Q2 and K2
            subcatch[ksub].ws[kstep][LZ] = subcatch[ksub].q[kstep][Q2] * (1.0 / subcatch[ksub].k2 - 1.0);

            // Add residual moisture
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
                printf("  Percolation rate is %.2lf mm day-1\n", subcatch[ksub].perc);
            }
            else if (strcmp(tag, "K2") == 0)
            {
                subcatch[ksub].k2 = value;
                printf("  K2 is %.2lf day-1\n", subcatch[ksub].k2);
            }
            else if (strcmp(tag, "MAXBAS") == 0)
            {
                subcatch[ksub].maxbas = value;
                printf("  Routing parameter is %.2lf\n", subcatch[ksub].maxbas);
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
