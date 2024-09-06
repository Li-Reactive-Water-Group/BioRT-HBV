#include "biort.h"

void ReadHbvResults(const char dir[], int nsub,  double stepsize, int *nsteps, double *steps[], subcatch_struct subcatch[], int mode)//double stepsize,
{
    FILE           *fp;
    char            fn[MAXSTRING];
    char            cmdstr[MAXSTRING];
    int             ksub;
    int             kstep;
    int             lno = 0;
    int             i;
    int             len_numexp=1;
    int             numexp_file_flag=0;
    int             date;
    int             hour;
    int             min;
    int             sec;

    if (mode==0)
    {
        sprintf(fn, "input/%s/Results.txt", dir);
        fp = fopen(fn, "r");
    } else if (mode == 1)
    {
        sprintf(fn, "input/%s/Numexp_Results.txt", dir);
        if (access(fn , F_OK ) != -1)
        {
            fp = fopen(fn, "r");
            numexp_file_flag = 1;
            biort_printf(VL_NORMAL, "\nHydrology in %s will be used for numerical experiment.\n", fn);
        } else
        {
            sprintf(fn, "input/%s/Numexp_precipchem.txt", dir);
            fp = fopen(fn, "r");
            len_numexp = CountLines(fp, cmdstr, 0) - 1;
            fclose(fp);
            sprintf(fn, "input/%s/Results.txt", dir);
            fp = fopen(fn, "r");
            numexp_file_flag = 0;
            biort_printf(VL_NORMAL, "\nHydrology in %s will be used for numerical experiment.\n", fn);
        }
    }


    *nsteps = CountLines(fp, cmdstr, 0) - 1;

    rewind(fp);

    *steps = (double *)malloc(*nsteps * sizeof(double));

    if (mode == 1 & numexp_file_flag == 0)
    {
      if (len_numexp % *nsteps != 0){
          biort_printf(VL_ERROR,"\nNumber of time steps in \"Numexp_precipchem.txt\" should be a multiple of time steps in \"Results.txt\" file, if \"Numexp_precipchem.txt\" file is not provided.\n");
          exit(EXIT_FAILURE);
    }
        len_numexp /= *nsteps;
    }

    for (ksub = 0; ksub < nsub; ksub++)
    {
        subcatch[ksub].ws = (double **)malloc(len_numexp * *nsteps * sizeof(double *));
        subcatch[ksub].q = (double **)malloc(len_numexp * *nsteps * sizeof(double *));
        subcatch[ksub].tmp = (double *)malloc(len_numexp * *nsteps * sizeof(double));

        for (kstep = 0; kstep < len_numexp * *nsteps; kstep++)
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
                if (stepsize < 86400)
                {
                    fscanf(fp, "%d %d:%d:%d",&date, &hour, &min, &sec);
                    ((*steps)[kstep])=date+((hour*3600+min*60+sec)/86400.0);
                }
                else 
                {
                    fscanf(fp, "%lf", &((*steps)[kstep]));
                    //((*steps)[kstep])=date;// Read model steps
                }
                
            }
            else
            {
                fscanf(fp, "%*d");
            }

            fscanf(fp, "%*f %*f");  // Skip "Qsim" and "Qobs"
            fscanf(fp, "%lf %lf", &subcatch[ksub].q[kstep][PRECIP], &subcatch[ksub].tmp[kstep]);    // Read precip. and
                                                                                                    // air temperature
            fscanf(fp, "%*f %*f");   // Skip "AET" and PET"
            fscanf(fp, "%lf %*lf %lf", &snow, &sm);     // Read snow and soil moisture
            subcatch[ksub].ws[kstep][SNOW] = snow;  // 2021-05-14
            subcatch[ksub].ws[kstep][SM] = sm;
            fscanf(fp, "%lf", &subcatch[ksub].q[kstep][RECHG]);  // Read recharge
            fscanf(fp, "%*lf %*lf");   // Skip upper and lower zone storages
            fscanf(fp, "%lf %lf %lf", &subcatch[ksub].q[kstep][Q0],
                &subcatch[ksub].q[kstep][Q1], &subcatch[ksub].q[kstep][Q2]);  // Read Q0, Q1, and Q2
            fscanf(fp, "%*lf %*lf");   // Skip "Qsim_rain" and "Qsim_snow"

            //Incoming precipitation might become snow or enter the soil zone directly.
            //  Psnow = PRECIP * SFCF (if temp < TT)
            //  Prain = PRECIP (if temp >=TT)
            //  Psnow + snow0 = snowmelt + snow
            //
            // Solving the above equations,
            //
            subcatch[ksub].q[kstep][Prain] = (subcatch[ksub].tmp[kstep] < subcatch[ksub].tt) ?
                0.0 : subcatch[ksub].q[kstep][PRECIP];
            
            subcatch[ksub].q[kstep][Psnow] = (subcatch[ksub].tmp[kstep] < subcatch[ksub].tt) ?
                subcatch[ksub].q[kstep][PRECIP] * subcatch[ksub].sfcf : 0.0;
            
            subcatch[ksub].q[kstep][snowmelt] = (kstep == 0) ?
                0.0 : MAX(0, subcatch[ksub].q[kstep][Psnow] + subcatch[ksub].ws[kstep - 1][SNOW] - subcatch[ksub].ws[kstep][SNOW]);
            //biort_printf(VL_NORMAL, "  snowmelt %.2f m3 m-3.\n", subcatch[ksub].q[kstep][snowmelt]);
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
            subcatch[ksub].ws[kstep][SURFACE] = subcatch[ksub].q[kstep][Q0];
            subcatch[ksub].ws[kstep][UZ] = subcatch[ksub].q[kstep][Q1] / subcatch[ksub].k1 -
                (subcatch[ksub].q[kstep][Q0] + subcatch[ksub].q[kstep][Q1]);
            subcatch[ksub].ws[kstep][LZ] = subcatch[ksub].q[kstep][Q2] / subcatch[ksub].k2 -
                subcatch[ksub].q[kstep][Q2];
            subcatch[ksub].q[kstep][PERC] = (kstep == 0) ?
                0.0 : MIN(subcatch[ksub].perc, subcatch[ksub].ws[kstep - 1][UZ] + subcatch[ksub].q[kstep][RECHG]);

        }
        
        //Change PERC value for 1st time step using water storage for last timestep
        
        subcatch[ksub].q[0][PERC] = MIN(subcatch[ksub].perc, subcatch[ksub].ws[kstep - 1][UZ] + subcatch[ksub].q[0][RECHG]);

        // Add 1. residual moisture to LZ & UZ and 2. SM to UZ
        for (kstep = 0; kstep < *nsteps; kstep++)
        {
            //subcatch[ksub].ws[kstep][SURFACE] += subcatch[ksub].res_surface;   // 2021-05-14
            subcatch[ksub].ws[kstep][UZ] += subcatch[ksub].res_uz;
            subcatch[ksub].ws[kstep][LZ] += subcatch[ksub].res_lz;
            subcatch[ksub].ws[kstep][UZ] += subcatch[ksub].ws[kstep][SM];

        }
        for (kstep = 0; kstep < *nsteps; kstep++)
        {
            //if (subcatch[ksub].ws[kstep][SURFACE] > (subcatch[ksub].d_surface * subcatch[ksub].porosity_surface))
            //{
                //biort_printf(VL_NORMAL, "\nWater storage in SURFACE exceeds maximum water storage capacity at line %d.\n", kstep);
            //}
            if (subcatch[ksub].ws[kstep][UZ] > (subcatch[ksub].d_uz * subcatch[ksub].porosity_uz))
            {
                biort_printf(VL_ERROR, "\nWater storage in UZ exceeds maximum water storage capacity at line %d.\n", kstep);
                exit(EXIT_FAILURE);
            }
            if (subcatch[ksub].ws[kstep][LZ] > (subcatch[ksub].d_lz * subcatch[ksub].porosity_lz))
            {
                biort_printf(VL_ERROR, "\nWater storage in LZ exceeds maximum water storage capacity at line %d.\n", kstep);
                exit(EXIT_FAILURE);
            }
        }

        fclose(fp);
    }

    if (mode == 1 & numexp_file_flag == 0)
    {
        if (len_numexp > 1)
        {
            for (i = 1; i < len_numexp; i++)
            {
                for (ksub = 0; ksub < nsub; ksub++)
                {
                    for (kstep = 0; kstep < *nsteps; kstep++)
                    {
                        subcatch[ksub].ws[(i * *nsteps)+kstep][SNOW] = subcatch[ksub].ws[kstep][SNOW];
                        //subcatch[ksub].ws[(i * *nsteps)+kstep][SURFACE] = subcatch[ksub].ws[kstep][SURFACE];
                        subcatch[ksub].ws[(i * *nsteps)+kstep][UZ] = subcatch[ksub].ws[kstep][UZ];
                        subcatch[ksub].ws[(i * *nsteps)+kstep][LZ] = subcatch[ksub].ws[kstep][LZ];
                        subcatch[ksub].q[(i * *nsteps)+kstep][PRECIP] = subcatch[ksub].q[kstep][PRECIP];
                        subcatch[ksub].q[(i * *nsteps)+kstep][RECHG] = subcatch[ksub].q[kstep][RECHG];
                        subcatch[ksub].q[(i * *nsteps)+kstep][PERC] = subcatch[ksub].q[kstep][PERC];
                        subcatch[ksub].q[(i * *nsteps)+kstep][Q0] = subcatch[ksub].q[kstep][Q0];
                        subcatch[ksub].q[(i * *nsteps)+kstep][Q1] = subcatch[ksub].q[kstep][Q1];
                        subcatch[ksub].q[(i * *nsteps)+kstep][Q2] = subcatch[ksub].q[kstep][Q2];
                        subcatch[ksub].tmp[(i * *nsteps)+kstep] = subcatch[ksub].tmp[kstep];
                        //biort_printf(VL_NORMAL, "\n%d kstep: PERC%d\n", kstep, subcatch[ksub].q[kstep][PERC]);
                    }
                }
            }
        }
    }
    *nsteps *= len_numexp;
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
        FindLine(fp, "<SubCatchmentVegetationZoneParameters>", &lno, cmdstr);

        while (!feof(fp))
        {
            NextLine(fp, cmdstr, &lno);
            ParseLine(cmdstr, tag, &value);
            if (strcmp(tag,"TT") == 0)
            {
                subcatch[ksub].tt = value;
                biort_printf(VL_NORMAL, "  TT is %.2lf\n", subcatch[ksub].tt);
            }
            else if (strcmp(tag,"SFCF") == 0)
            {
                subcatch[ksub].sfcf = value;
                biort_printf(VL_NORMAL, "  SFCF is %.2lf\n", subcatch[ksub].sfcf);
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
