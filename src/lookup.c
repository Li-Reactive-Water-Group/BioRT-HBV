#include "biort.h"

void Lookup(FILE *fp, const calib_struct *calib, chemtbl_struct chemtbl[], kintbl_struct kintbl[], rttbl_struct *rttbl, subcatch_struct subcatch[])
{
    char            cmdstr[MAXSTRING];
    int             i, j, k;
    int             ind;
    int             keq_position = BADVAL;
    int             total_temp_points;
    int             lno = 0;
    double          pot_dep[MAXSPS][MAXSPS];// dependency of all possible kinetic species on primary species
    double          keq_kin_all[MAXSPS];    // Keq's of all possible kinetic species

    biort_printf(VL_NORMAL, "\nDATABASE\n");

    // Find temperature point in database
    NextLine(fp, cmdstr, &lno);
    while(MatchWrappedKey(cmdstr, "'temperature points'") != 0)
    {
        NextLine(fp, cmdstr, &lno);
    }
    ReadTempPoints(cmdstr, rttbl->tmp, &total_temp_points, &keq_position);
    if (keq_position == BADVAL)
    {
        biort_printf(VL_ERROR, "Error reading temperature points in %s near Line %d", ".cdbs", lno);
        exit(EXIT_FAILURE);
    }

    // Read Debye Huckel parameters from database
    rttbl->adh = BADVAL;
    rttbl->bdh = BADVAL;
    rttbl->bdt = BADVAL;

    while (MatchWrappedKey(cmdstr, "'Debye-Huckel adh'") != 0)
    {
        NextLine(fp, cmdstr, &lno);
    }
    ReadDHParam(cmdstr, keq_position, &rttbl->adh);
    if (roundi(rttbl->adh) == BADVAL)
    {
        biort_printf(VL_ERROR, "Error reading Debye Huckel parameters in %s near Line %d", ".cdbs", lno);
        exit(EXIT_FAILURE);
    }

    while (MatchWrappedKey(cmdstr, "'Debye-Huckel bdh'") != 0)
    {
        NextLine(fp, cmdstr, &lno);
    }
    ReadDHParam(cmdstr, keq_position, &rttbl->bdh);
    if (roundi(rttbl->bdh) == BADVAL)
    {
        biort_printf(VL_ERROR, "Error reading Debye Huckel parameters in %s near Line %d", ".cdbs", lno);
        exit(EXIT_FAILURE);
    }

    while (MatchWrappedKey(cmdstr, "'Debye-Huckel bdt'") != 0)
    {
        NextLine(fp, cmdstr, &lno);
    }
    ReadDHParam(cmdstr, keq_position, &rttbl->bdt);
    if (roundi(rttbl->bdt) == BADVAL)
    {
        biort_printf(VL_ERROR, "Error reading Debye Huckel parameters in %s near Line %d", ".cdbs", lno);
        exit(EXIT_FAILURE);
    }

    biort_printf(VL_NORMAL, " Debye-Huckel Parameters set to A=%6.4f; B=%6.4f; b=%6.4f\n\n",
        rttbl->adh, rttbl->bdh, rttbl->bdt);


    for (i = 0; i < MAXSPS; i++)
    {
        for (j = 0; j < MAXSPS; j++)
        {
            pot_dep[i][j] = 0.0;    // num_min x num_stc
            rttbl->dep_mtx[i][j] = 0.0;     // num_ssc x num_sdc
            rttbl->dep_kin[i][j] = 0.0;     // (um_mkr + num_akr) x num_stc
            rttbl->conc_contrib[i][j] = 0.0;    // num_stc x (num_stc + num_ssc)
#if NOT_YET_IMPLEMENTED
            rttbl->Totalconck[i][j] = 0.0;  // num_stc x (num_stc + num_ssc)
#endif
        }
        // Keqs of equilibrium: kinetic and kinetic all
        keq_kin_all[i] = 0.0;   // num_min
        rttbl->keq[i] = 0.0;    // num_ssc
        rttbl->keq_kin[i] = 0.0;    // num_mkr + num_akr
    }

    // Update species parameters
    for (i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        if (strcmp(chemtbl[i].name, "pH") == 0)
        {
            strcpy(chemtbl[i].name, "H+");
        }
        Wrap(chemtbl[i].name);

        chemtbl[i].charge = 0.0;
        chemtbl[i].size_fac = 1.0;
        chemtbl[i].molar_mass = BADVAL;
    }

    // Read parameters for primary species
    while (MatchWrappedKey(cmdstr, "'End of primary'") != 0)
    {
        NextLine(fp, cmdstr, &lno);
        ReadPrimary(cmdstr, rttbl->num_stc, chemtbl);
    }

    // Read parameters for secondary species
    while (MatchWrappedKey(cmdstr, "'End of secondary'") != 0)
    {
        NextLine(fp, cmdstr, &lno);
        ReadSecondary(cmdstr, total_temp_points, keq_position, chemtbl, rttbl);
    }

    // Read parameters for minerals
    while (MatchWrappedKey(cmdstr, "'End of minerals'") != 0)
    {
        NextLine(fp, cmdstr, &lno);
        ReadMinerals(cmdstr, total_temp_points, keq_position, pot_dep, keq_kin_all, chemtbl, rttbl);
    }

    // Build dependencies
    for (i = 0; i < rttbl->num_mkr + rttbl->num_akr; i++)
    {
        ind = kintbl[i].position - rttbl->num_stc + rttbl->num_min;
        biort_printf(VL_NORMAL, " Selecting the kinetic species %s from all possible species.\n\n",
            chemtbl[kintbl[i].position].name);
        rttbl->keq_kin[i] = keq_kin_all[ind];
        for (k = 0; k < rttbl->num_stc; k++)
        {
            rttbl->dep_kin[i][k] = pot_dep[ind][k];
        }
    }

    while (strcmp(cmdstr, "End of surface complexation\r\n") != 0 &&
        strcmp(cmdstr, "End of surface complexation\n") != 0)
    {
        NextLine(fp, cmdstr, &lno);
        ReadAdsorption(cmdstr, total_temp_points, keq_position, chemtbl, rttbl);
    }

    while (strcmp(cmdstr, "End of exchange\r\n") != 0 && strcmp(cmdstr, "End of exchange\n") != 0)
    {
        NextLine(fp, cmdstr, &lno);
        ReadCationEchg(cmdstr, calib->xsorption, chemtbl, rttbl);
    }

    for (i = 0; i < MAXSPS; i++)
    {
        for (j = 0; j < MAXDEP; j++)
        {
            kintbl[i].dep_index[j] = 0;
            kintbl[i].monod_index[j] = 0;
            kintbl[i].inhib_index[j] = 0;
        }
    }

    for (i = 0; i < rttbl->num_mkr; i++)
    {
        // Initialize kinetic reaction type
        kintbl[i].type = BADVAL;

        FindLine(fp, "BOF", &lno, ".cdbs");
        NextLine(fp, cmdstr, &lno);
        while (strcmp(cmdstr, "Begin mineral kinetics\r\n") != 0 && strcmp(cmdstr, "Begin mineral kinetics\n") != 0)
        {
            NextLine(fp, cmdstr, &lno);
        }

        while (strcmp(cmdstr, "End of mineral kinetics\r\n") != 0 && strcmp(cmdstr, "End of mineral kinetics\n") != 0)
        {
            ReadMinKin(fp, rttbl->num_stc, subcatch[0].chms[UZ].k_cini[i + rttbl->num_stc - rttbl->num_min], &lno, cmdstr, chemtbl, &kintbl[i]);
            NextLine(fp, cmdstr, &lno);
        }

        if (kintbl[i].type == BADVAL)
        {
            biort_printf(VL_ERROR, "Error finding mineral kinetic %s -label %s in the database.\n",
                chemtbl[kintbl[i].position].name, kintbl[i].label);
            exit(EXIT_FAILURE);
        }
    }

    for (i = 0; i < rttbl->num_stc; i++)
    {
        rttbl->conc_contrib[i][i] = 1.0;
    }

    for (i = rttbl->num_stc; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        for (j = 0; j < rttbl->num_sdc; j++)
        {
            rttbl->conc_contrib[j][i] += rttbl->dep_mtx[i - rttbl->num_stc][j];
        }
    }
    if (rttbl->num_ssc > 0)
    {
        biort_printf(VL_NORMAL, " \n Dependency Matrix!\n");

        biort_printf(VL_NORMAL, "%-15s", " ");
        for (i = 0; i < rttbl->num_sdc; i++)
        {
            biort_printf(VL_NORMAL, "%-14s", chemtbl[i].name);
        }
        biort_printf(VL_NORMAL, "\n");

        for (i = 0; i < rttbl->num_ssc; i++)
        {
            biort_printf(VL_NORMAL, " %-14s", chemtbl[i + rttbl->num_stc].name);
            for (j = 0; j < rttbl->num_sdc; j++)
            {
                biort_printf(VL_NORMAL, "%-14.2f", rttbl->dep_mtx[i][j]);
            }
            biort_printf(VL_NORMAL, " %6.2f\n", rttbl->keq[i]);
        }
    }

    biort_printf(VL_NORMAL, " \n Total Concentration Matrix!\n");
    biort_printf(VL_NORMAL, "%-18s", " ");
    for (i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        biort_printf(VL_NORMAL, "%-14s", chemtbl[i].name);
    }
    biort_printf(VL_NORMAL, "\n");
    for (i = 0; i < rttbl->num_stc; i++)
    {
        biort_printf(VL_NORMAL, " Sum%-14s", chemtbl[i].name);
        for (j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
        {
            biort_printf(VL_NORMAL, "%-14.2f", rttbl->conc_contrib[i][j]);
        }
        biort_printf(VL_NORMAL, "\n");
    }

    biort_printf(VL_NORMAL, " \n Kinetic Mass Matrix!\n");
    biort_printf(VL_NORMAL, "%-15s", " ");
    for (i = 0; i < rttbl->num_stc; i++)
    {
        biort_printf(VL_NORMAL, "%-14s", chemtbl[i].name);
    }
    biort_printf(VL_NORMAL, "\n");
    for (j = 0; j < rttbl->num_mkr + rttbl->num_akr; j++)
    {
        biort_printf(VL_NORMAL, " %-14s", chemtbl[kintbl[j].position].name);
        for (i = 0; i < rttbl->num_stc; i++)
        {
            biort_printf(VL_NORMAL, "%-14f", rttbl->dep_kin[j][i]);
        }
        biort_printf(VL_NORMAL, " Keq = %-6.2f\n", rttbl->keq_kin[j]);
    }

#if NOT_YET_IMPLEMENTED
    // Use calibration coefficient to produce new Keq values for 1) CO2, 2) other kinetic reaction
    double          Cal_PCO2 = 1.0;
    double          Cal_Keq = 1.0;
    for (i = 0; i < rttbl->num_akr + rttbl->num_mkr; i++)
    {
        rttbl->KeqKinect[i] += (!strcmp(chemtbl[i + rttbl->num_spc + rttbl->num_ads + rttbl->num_cex].name,
            "'CO2(*g)'")) ? log10(Cal_PCO2) : log10(Cal_Keq);
    }

    printf("\n Kinetic Mass Matrix (calibrated Keq)! \n");
    printf("%-15s", " ");
    for (i = 0; i < rttbl->num_stc; i++)
    {
        printf("%-14s", chemtbl[i].name);
    }
    printf("\n");
    for (j = 0; j < rttbl->num_mkr + rttbl->num_akr; j++)
    {
        printf(" %-14s", chemtbl[j + rttbl->num_spc + rttbl->num_ads + rttbl->num_cex].name);
        for (i = 0; i < rttbl->num_stc; i++)
        {
            printf("%-14.2f", rttbl->Dep_kinetic[j][i]);
        }
        printf(" Keq = %-6.2f\n", rttbl->KeqKinect[j]);
    }
    printf("\n");
#endif

    biort_printf(VL_NORMAL, " \n Mass action species type determination (0: immobile, 1: mobile, 2: Mixed) \n");
    for (i = 0; i < rttbl->num_spc; i++)
    {
        chemtbl[i].mtype = (chemtbl[i].itype == AQUEOUS) ? MOBILE_MA : IMMOBILE_MA;

        for (j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
        {
            chemtbl[i].mtype = (rttbl->conc_contrib[i][j] != 0 && chemtbl[j].itype != chemtbl[i].mtype) ?
                MIXED_MA : chemtbl[i].mtype;
        }
        biort_printf(VL_NORMAL, " %12s    %10d\n", chemtbl[i].name, chemtbl[i].mtype);
    }

    biort_printf(VL_NORMAL,
        " \n Individual species type determination (1: aqueous, 2: adsorption, 3: ion exchange, 4: solid)\n");
    for (i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        biort_printf(VL_NORMAL, " %12s    %10d\n", chemtbl[i].name, chemtbl[i].itype);
    }
}

void ReadTempPoints(const char cmdstr[], double tmp, int *total_points, int *keq_position)
{
    int             bytes_now;
    int             bytes_consumed = 0;
    int             i;
    double          val;

    if (sscanf(cmdstr + bytes_consumed, "'%*[^']'%n", &bytes_now) != 0)
    {
        return;
    }
    bytes_consumed += bytes_now;

    if (sscanf(cmdstr + bytes_consumed, "%d%n", total_points, &bytes_now) != 1)
    {
        return;
    }
    bytes_consumed += bytes_now;

    for (i = 0; i < *total_points; i++)
    {
        if (sscanf(cmdstr + bytes_consumed, "%lf%n", &val, &bytes_now) != 1)
        {
            return;
        }
        bytes_consumed += bytes_now;

        if (fabs(val - tmp) < 1.0E-5)
        {
            *keq_position = i + 1;
            biort_printf(VL_NORMAL, "\n Temperature point %.2f found in database (Pos. %d).\n\n", val, *keq_position);
            return;
        }
    }
}

void ReadDHParam(const char cmdstr[], int tmp_position, double *param)
{
    int             bytes_now;
    int             bytes_consumed = 0;
    int             i;

    if (sscanf(cmdstr + bytes_consumed, "'%*[^']'%n", &bytes_now) != 0)
    {
        return;
    }
    bytes_consumed += bytes_now;

    for (i = 0; i < tmp_position; i++)
    {
        if (sscanf(cmdstr + bytes_consumed, "%lf%n", param, &bytes_now) != 1)
        {
            return;
        }
        bytes_consumed += bytes_now;
    }
}

void ReadPrimary(const char cmdstr[], int num_stc, chemtbl_struct chemtbl[])
{
    int             i;

    for (i = 0; i < num_stc; i++)
    {
        if (MatchWrappedKey(cmdstr, chemtbl[i].name) == 0)
        {
            if (sscanf(cmdstr, "'%*[^']' %lf %lf %lf",
                &chemtbl[i].size_fac, &chemtbl[i].charge, &chemtbl[i].molar_mass) != 3)
            {
                biort_printf(VL_ERROR, "Error reading primary species parameters for %s\n", chemtbl[i].name);
                exit(EXIT_FAILURE);
            }
            biort_printf(VL_NORMAL, " Primary species %s found in database.\n molar_mass = %6.4f\n\n",
                chemtbl[i].name, chemtbl[i].molar_mass);
            break;
        }
    }
}

void ReadSecondary(const char cmdstr[], int npoints, int keq_position, chemtbl_struct chemtbl[], rttbl_struct *rttbl)
{
    int             bytes_now;
    int             bytes_consumed = 0;
    int             i, j, k;
    int             ind;
    int             ndep;
    double          dep;
    char            chemn[MAXSTRING];

    sscanf(cmdstr + bytes_consumed, "'%[^']' %d%n", chemn, &ndep, &bytes_now);
    bytes_consumed += bytes_now;
    Wrap(chemn);

    for (i = 0; i < rttbl->num_ssc; i++)
    {
        ind = i + rttbl->num_stc;

        if (MatchWrappedKey(chemn, chemtbl[ind].name) == 0)
        {
            biort_printf(VL_NORMAL, " Secondary species %s found in database!\n", chemtbl[ind].name);
            biort_printf(VL_NORMAL, " %s", cmdstr);
            chemtbl[ind].itype = AQUEOUS;

            for (j = 0; j < ndep; j++)
            {
                sscanf(cmdstr + bytes_consumed, "%lf '%[^']'%n", &dep, chemn, &bytes_now);
                bytes_consumed += bytes_now;
                Wrap(chemn);

                for (k = 0; k < rttbl->num_sdc; k++)
                {
                    if (MatchWrappedKey(chemtbl[k].name, chemn) == 0)
                    {
                        rttbl->dep_mtx[i][k] = dep;
                    }
                }
            }

            for (j = 0; j < npoints; j++)
            {
                if (j + 1 == keq_position)
                {
                    sscanf(cmdstr + bytes_consumed, "%lf%n", &rttbl->keq[i], &bytes_now);
                    bytes_consumed += bytes_now;
                    biort_printf(VL_NORMAL, " Keq = %6.4f\n", rttbl->keq[i]);
                }
                else
                {
                    sscanf(cmdstr + bytes_consumed, "%*f%n", &bytes_now);
                    bytes_consumed += bytes_now;
                }
            }

            sscanf(cmdstr + bytes_consumed, "%lf %lf %lf",
                &chemtbl[ind].size_fac, &chemtbl[ind].charge, &chemtbl[ind].molar_mass);
            biort_printf(VL_NORMAL, " molar_mass = %6.4f, Charge = %6.4f, SizeFactor = %6.4f\n\n",
                chemtbl[ind].molar_mass, chemtbl[ind].charge, chemtbl[ind].size_fac);

            break;
        }
    }
}

void ReadMinerals(const char cmdstr[], int npoints, int keq_position, double pot_dep[MAXSPS][MAXSPS],
    double keq_kin_all[], chemtbl_struct chemtbl[], rttbl_struct *rttbl)
{
    int             bytes_now;
    int             bytes_consumed = 0;
    int             i, j, k, l;
    int             ind;
    int             ndep;
    double          dep;
    double          mvol;
    double          keq_sec = 0.0;   // fix by Kayal, 2021-05-26
    char            chemn[MAXSTRING];

    sscanf(cmdstr + bytes_consumed, "'%[^']' %lf %d%n", chemn, &mvol, &ndep, &bytes_now);
    bytes_consumed += bytes_now;
    Wrap(chemn);

    for (i = 0; i < rttbl->num_min; i++)
    {
        ind = i + rttbl->num_spc + rttbl->num_ads + rttbl->num_cex;

        if (MatchWrappedKey(chemn, chemtbl[ind].name) == 0)
        {
            biort_printf(VL_NORMAL, " Mineral %s found in database!\n", chemtbl[ind].name);
            biort_printf(VL_NORMAL, " %s", cmdstr);

            chemtbl[ind].molar_vol = mvol;
            chemtbl[ind].itype = MINERAL;

            for (j = 0; j < ndep; j++)
            {
                sscanf(cmdstr + bytes_consumed, "%lf '%[^']'%n", &dep, chemn, &bytes_now);
                bytes_consumed += bytes_now;
                Wrap(chemn);

                for (k = 0; k < rttbl->num_stc + rttbl->num_ssc; k++)
                {
                    if (MatchWrappedKey(chemtbl[k].name, chemn) == 0)
                    {
                        if (k < rttbl->num_stc)
                        {
                            pot_dep[i][k] = dep;
                        }
                        else
                        {
                            for (l = 0; l < rttbl->num_spc; l++)
                            {
                                pot_dep[i][l] += dep * rttbl->dep_mtx[k - rttbl->num_stc][l];
                            }
                            //keq_kin_all[i] += dep * rttbl->keq[k - rttbl->num_stc];
                            keq_sec += dep * rttbl->keq[k - rttbl->num_stc];  // fix by for adjusting Keq when reaction is written in terms of secondary species in database, 2021-05-26
                        }

                        break;
                    }
                }
            }

            for (j = 0; j < npoints; j++)
            {
                if (j + 1 == keq_position)
                {
                    sscanf(cmdstr + bytes_consumed, "%lf%n", &keq_kin_all[i], &bytes_now);
                    bytes_consumed += bytes_now;
                    if (keq_sec != 0)   // fix for adjusting Keq when reaction is written in terms of secondary species in database, 2021-05-26
                    {
                        keq_kin_all[i] += keq_sec;
                    }
                    keq_sec = 0;
                }
                else
                {
                    sscanf(cmdstr + bytes_consumed, "%*f%n", &bytes_now);
                    bytes_consumed += bytes_now;
                }
            }

            sscanf(cmdstr + bytes_consumed, "%lf", &chemtbl[ind].molar_mass);

            pot_dep[i][ind] = -1.0;

            biort_printf(VL_NORMAL, " Keq = %6.4f\n", keq_kin_all[i]);
            chemtbl[ind].charge = 0;
            biort_printf(VL_NORMAL, " molar_mass = %6.4f, molar_vol = %6.4f\n\n",
                chemtbl[ind].molar_mass, chemtbl[ind].molar_vol);

            break;
        }
    }
}

void ReadAdsorption(const char cmdstr[], int npoints, int keq_position, chemtbl_struct chemtbl[], rttbl_struct *rttbl)
{
    int             bytes_now;
    int             bytes_consumed = 0;
    int             i, j, k;
    int             ind;
    int             ndep;
    double          dep;
    char            chemn[MAXSTRING];

    if (sscanf(cmdstr + bytes_consumed, "'%[^']' %d%n", chemn, &ndep, &bytes_now) != 2)
    {
        return;
    }

    bytes_consumed += bytes_now;
    Wrap(chemn);

    for (i = 0; i < rttbl->num_ssc; i++)
    {
        ind = i + rttbl->num_stc;

        if (strcmp(chemtbl[ind].name, chemn) == 0)
        {
            biort_printf(VL_NORMAL, " Secondary surface complexation %s found in database!\n", chemtbl[ind].name);
            biort_printf(VL_NORMAL, " %s", cmdstr);
            chemtbl[ind].itype = ADSORPTION;
            for (j = 0; j < ndep; j++)
            {
                sscanf(cmdstr + bytes_consumed, "%lf '%[^']'%n", &dep, chemn, &bytes_now);
                bytes_consumed += bytes_now;
                Wrap(chemn);

                for (k = 0; k < rttbl->num_sdc; k++)
                {
                    if (strcmp(chemtbl[k].name, chemn) == 0)
                    {
                        rttbl->dep_mtx[i][k] = dep;
                        break;
                    }
                }
            }

            for (j = 0; j < npoints; j++)
            {
                if (j + 1 == keq_position)
                {
                    sscanf(cmdstr + bytes_consumed, "%lf%n", &rttbl->keq[i], &bytes_now);
                    bytes_consumed += bytes_now;
                    biort_printf(VL_NORMAL, " Keq = %6.4f\n", rttbl->keq[i]);
                }
                else
                {
                    sscanf(cmdstr + bytes_consumed, "%*f%n", &bytes_now);
                    bytes_consumed += bytes_now;
                }
            }
        }
    }
}

void ReadCationEchg(const char cmdstr[], double calval, chemtbl_struct chemtbl[], rttbl_struct *rttbl)
{
    int             bytes_now;
    int             bytes_consumed = 0;
    int             i, j, k;
    int             ind;
    int             ndep;
    double          dep;
    char            chemn[MAXSTRING];

    if (sscanf(cmdstr + bytes_consumed, "'%[^']' %d%n", chemn, &ndep, &bytes_now) != 2)
    {
        return;
    }
    bytes_consumed += bytes_now;
    Wrap(chemn);

    for (i = 0; i < rttbl->num_ssc; i++)
    {
        ind = i + rttbl->num_stc;

        if (strcmp(chemtbl[ind].name, chemn) == 0)
        {
            biort_printf(VL_NORMAL, " Secondary ion exchange %s found in database!\n", chemtbl[ind].name);
            biort_printf(VL_NORMAL, " %s", cmdstr);
            chemtbl[ind].itype = CATION_ECHG;
            for (j = 0; j < ndep; j++)
            {
                sscanf(cmdstr + bytes_consumed, "%lf '%[^']'%n", &dep, chemn, &bytes_now);
                bytes_consumed += bytes_now;
                Wrap(chemn);

                for (k = 0; k < rttbl->num_sdc; k++)
                {
                    if (strcmp(chemtbl[k].name, chemn) == 0)
                    {
                        rttbl->dep_mtx[i][k] = dep;
                        break;
                    }
                }
            }
            sscanf(cmdstr + bytes_consumed, "%lf", &rttbl->keq[i]);
            biort_printf(VL_NORMAL, " Keq = %6.4f \n", rttbl->keq[i]);

            rttbl->keq[i] += calval;
            biort_printf(VL_NORMAL, " After calibration: Keq = %6.4f \n", rttbl->keq[i]);
        }
    }
}

void ReadMinKin(FILE *fp, int num_stc, double calval, int *lno, char cmdstr[], chemtbl_struct chemtbl[],
    kintbl_struct *kintbl)
{
    int             bytes_now;
    int             bytes_consumed = 0;
    double          temp;
    char            chemn[MAXSTRING];
    char            optstr[MAXSTRING];
    char            label[MAXSTRING];

    sscanf(cmdstr, "%s", chemn);
    Wrap(chemn);

    if (strcmp(chemtbl[kintbl->position].name, chemn) == 0)
    {
        NextLine(fp, cmdstr, lno);
        sscanf(cmdstr, "%*s = %s", label);

        if (strcmp(kintbl->label, label) == 0)
        {
            biort_printf(VL_NORMAL, " \n Mineral kinetics %s %s found in database!\n",
                chemtbl[kintbl->position].name, kintbl->label);

            // For mineral kinetics, all species have the following lines:
            //   label, type, rate(25C), and activation.
            // Optional lines are:
            //   dependence, monod_terms, inhibition, and biomass
            // Read type
            NextLine(fp, cmdstr, lno);
            sscanf(cmdstr, "%s = %s", optstr, label);
            if (strcmp(optstr, "type") == 0)
            {
                if (strcmp(label, "tst") == 0)
                {
                    kintbl->type = TST;
                }
                else if (strcmp(label, "PrecipitationOnly") == 0)
                {
                    kintbl->type = PRCP_ONLY;
                }
                else if (strcmp(label, "DissolutionOnly") == 0)
                {
                    kintbl->type = DISS_ONLY;
                }
                else if (strcmp(label, "monod") == 0)
                {
                    kintbl->type = MONOD;
                }
            }
            else
            {
                biort_printf(VL_ERROR, "Error: Cannot find mineral kinetics type in .cdbs near Line %d.\n", *lno);
                exit(EXIT_FAILURE);
            }

            // Read rate
            NextLine(fp, cmdstr, lno);
            sscanf(cmdstr, "%s = %lf", optstr, &kintbl->rate);
            //if (strcmp(optstr, "rate(25C)") == 0)
            //{
            //    biort_printf(VL_NORMAL, " Rate is %f\n", kintbl->rate);

            //    kintbl->rate += calval;
            //    biort_printf(VL_NORMAL, " After calibration: Rate is %f, calib->Rate = %f \n", kintbl->rate, calval);
            //}
            //else
            //{
            //    biort_printf(VL_ERROR, "Error: Cannot find mineral kinetics rate in .cdbs near Line %d.\n", *lno);
            //    exit(EXIT_FAILURE);
            //}
            
            //Read rate from cini file
            kintbl->rate = calval;
            biort_printf(VL_NORMAL, " Rate is %f\n", kintbl->rate);

            // Read activation
            NextLine(fp, cmdstr, lno);
            sscanf(cmdstr, "%s = %lf", optstr, &kintbl->actv);
            if (strcmp(optstr, "activation") == 0)
            {
                biort_printf(VL_NORMAL, " Activation is %f\n", kintbl->actv);
            }
            else
            {
                biort_printf(VL_ERROR, "Error: Cannot find mineral kinetics activation in .cdbs near Line %d.\n", *lno);
                exit(EXIT_FAILURE);
            }

            kintbl->ndep = 0;
            kintbl->nmonod = 0;
            kintbl->ninhib = 0;

            NextLine(fp, cmdstr, lno);
            while (cmdstr[0] != '+')
            {
                bytes_consumed = 0;

                sscanf(cmdstr, "%s :%n", optstr, &bytes_now);
                bytes_consumed += bytes_now;

                if (strcmp(optstr, "dependence") == 0)
                {
                    // Read dependence
                    // Assume that all mineral kinetics only depend on one species
                    if (sscanf(cmdstr + bytes_consumed, "%s %lf", chemn, &temp) == 2)
                    {
                        Wrap(chemn);
                        kintbl->dep_index[0] = FindChem(chemn, num_stc, chemtbl);
                        if (kintbl->dep_index[0] < 0)
                        {
                            biort_printf(VL_ERROR, "Error finding dependence in species table.\n");
                            exit(EXIT_FAILURE);
                        }
                        kintbl->ndep = 1;
                        kintbl->dep_power[0] = temp;
                        biort_printf(VL_NORMAL, " Dependency: %s %f\n", chemn, kintbl->dep_power[0]);
                    }
                    else
                    {
                        biort_printf(VL_NORMAL, " No dependency.\n");
                    }
                }
                else if (strcmp(optstr, "biomass") == 0)
                {
                    // Biomass term
                    sscanf(cmdstr, "%*s : %s", chemn);
                    Wrap(chemn);
                    biort_printf(VL_NORMAL, " Biomass species: %s \n", chemn);
                    kintbl->biomass_index = FindChem(chemn, num_stc, chemtbl);
                    biort_printf(VL_NORMAL, " Biomass species position: %d \n", kintbl->biomass_index);
                }
                else if (strcmp(optstr, "monod_terms") == 0)
                {
                    // Monod term
                    while (sscanf(cmdstr + bytes_consumed, "%s %lf%n", chemn, &temp, &bytes_now) == 2)
                    {
                        bytes_consumed += bytes_now;

                        Wrap(chemn);
                        kintbl->monod_index[kintbl->nmonod] = FindChem(chemn, num_stc, chemtbl);
                        if (kintbl->monod_index[kintbl->nmonod] < 0)
                        {
                            biort_printf(VL_ERROR, "Error finding monod_terms in species table.\n");
                            exit(EXIT_FAILURE);
                        }
                        kintbl->monod_para[kintbl->nmonod] = temp;
                        biort_printf(VL_NORMAL, " Monod term: %s %f\n", chemn, kintbl->monod_para[kintbl->nmonod]);
                        kintbl->nmonod++;
                    }
                }
                else if (strcmp(optstr, "inhibition") == 0)
                {
                    // Inhibition term
                    while (sscanf(cmdstr + bytes_consumed, "%s %lf%n", chemn, &temp, &bytes_now) == 2)
                    {
                        bytes_consumed += bytes_now;

                        Wrap(chemn);
                        kintbl->inhib_index[kintbl->ninhib] = FindChem(chemn, num_stc, chemtbl);
                        if (kintbl->inhib_index[kintbl->ninhib] < 0)
                        {
                            biort_printf(VL_ERROR, "Error finding inhibition term in species table.\n");
                            exit(EXIT_FAILURE);
                        }
                        kintbl->inhib_para[kintbl->ninhib] = temp;
                        biort_printf(VL_NORMAL, " Inhibition term: %s %f\n", chemn, kintbl->inhib_para[kintbl->ninhib]);
                        kintbl->ninhib++;
                    }
                }

                NextLine(fp, cmdstr, lno);
            }
        }
    }
}
