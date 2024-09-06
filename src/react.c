#include "biort.h"


void SurfaceReaction(int kstep, int nsub,int ctrl_sfreaction, double stepsize, const double steps[], const chemtbl_struct chemtbl[],
    const kintbl_struct kintbl[], const rttbl_struct *rttbl, subcatch_struct subcatch[])
{
    int             ksub;
    int             kspc;
    double          satn;
    double          temp;
    double          substep;
    double          depth;
    double          porosity;
    double          Zw;
    
    
    for (ksub = 0; ksub < nsub; ksub++)
    {
        if (subcatch[ksub].ws[kstep][SURFACE] > 0.0 & ctrl_sfreaction == 1)
        {
            
            for (kspc = 0; kspc < rttbl->num_spc; kspc++)
            {
                subcatch[ksub].chms[STREAM].tot_mol[kspc] -= subcatch[ksub].chms[SURFACE].tot_conc[kspc] * subcatch[ksub].q[kstep][Q0];
                subcatch[ksub].chms[SURFACE].tot_mol[kspc] = subcatch[ksub].chms[SURFACE].tot_conc[kspc] * subcatch[ksub].q[kstep][Q0];
                
            }
            
            temp = subcatch[ksub].tmp[kstep];
            porosity=subcatch[ksub].porosity_surface;
            depth=subcatch[ksub].d_surface;//subcatch[ksub].ws[kstep][SURFACE];
            satn=1;
            Zw=0;
            
            for (kspc = 0; kspc < MAXSPS; kspc++)
            {
                subcatch[ksub].react_rate[SURFACE][kspc] = BADVAL;   // Set reaction rate to -999
            }
            
            substep = ReactControl(chemtbl, kintbl, rttbl, stepsize, porosity, depth, satn, temp, Zw,
                    subcatch[ksub].react_rate[SURFACE], &subcatch[ksub].chms[SURFACE]);
            
            if (substep < 0.0)
            {
               biort_printf(VL_NORMAL, "%f %s zone reaction failed with a substep of %.1lf s.\n",
                       steps[kstep], "SURFACE", -substep);
            }
            if (substep > 0.0)
            {
                biort_printf(VL_VERBOSE, "%f %s zone reaction passed with a minimum step of %.1lf s.\n",
                        steps[kstep], "SURFACE", substep);
            }
            
            for (kspc = 0; kspc < rttbl->num_spc; kspc++)
            {
                subcatch[ksub].chms[STREAM].tot_mol[kspc] += subcatch[ksub].chms[SURFACE].tot_conc[kspc] * subcatch[ksub].q[kstep][Q0];
                subcatch[ksub].chms[SURFACE].tot_mol[kspc] = 0.0;
                subcatch[ksub].chms[STREAM].tot_conc[kspc] =
                    (subcatch[ksub].q[kstep][Q0] + subcatch[ksub].q[kstep][Q1] + subcatch[ksub].q[kstep][Q2] > 0.0) ?
                    subcatch[ksub].chms[STREAM].tot_mol[kspc] /
                    (subcatch[ksub].q[kstep][Q0] + subcatch[ksub].q[kstep][Q1] + subcatch[ksub].q[kstep][Q2]) : ZERO_CONC;
            }
        } else
        {
            for (kspc = 0; kspc < MAXSPS; kspc++)
            {
                subcatch[ksub].react_rate[SURFACE][kspc] = ZERO_CONC;   // Set reaction rate to 0
            }
        }
    }
        

}

void Reaction(int kstep, int nsub, double stepsize, const double steps[], const chemtbl_struct chemtbl[],
    const kintbl_struct kintbl[], const rttbl_struct *rttbl, subcatch_struct subcatch[])
{
    int             ksub;
    int             kzone;
    int             kspc;
    double          satn;
    double          temp;
    double          substep;
    double          depth;
    double          porosity;
    double          Zw;
    const int       NZONES = 2;   // 2021-05-14

    for (ksub = 0; ksub < nsub; ksub++)
    {
        //ftemp = SoilTempFactor(rttbl->q10, subcatch[ksub].tmp[kstep]);
        temp = subcatch[ksub].tmp[kstep];

        for (kzone = UZ; kzone < UZ + NZONES; kzone++)   // 2021-05-14
        {
            switch (kzone)
            {
                //case SURFACE:   // 2021-05-14
                //    depth = subcatch[ksub].d_surface;
                //    porosity = subcatch[ksub].porosity_surface;
                //    break;
                case UZ:
                    depth = subcatch[ksub].d_uz;
                    porosity = subcatch[ksub].porosity_uz;
                    Zw = depth - (subcatch[ksub].ws[kstep][UZ]/porosity);
                    break;
                case LZ:
                    depth = subcatch[ksub].d_lz;
                    porosity = subcatch[ksub].porosity_lz;
                    Zw = depth - (subcatch[ksub].ws[kstep][LZ]/porosity);
                    break;
            }

            for (kspc = 0; kspc < MAXSPS; kspc++)
            {
                subcatch[ksub].react_rate[kzone][kspc] = BADVAL;   // Set reaction rate to -999
            }


            satn = subcatch[ksub].ws[kstep][kzone] / (depth * porosity);  // add porosity for saturation calculation

            satn = MIN(satn, 1.0);

            //biort_printf(VL_NORMAL, "%d %d %s zone reaction has saturation of %.1lf s.\n",
            //              steps[kstep],kzone, satn);

            if (satn > 1.0E-2)
            {
                substep = ReactControl(chemtbl, kintbl, rttbl, stepsize, porosity, depth, satn, temp, Zw,
                    subcatch[ksub].react_rate[kzone], &subcatch[ksub].chms[kzone]);

                if (substep < 0.0)
                {
                   if (kzone == SURFACE)   // 2021-05-14
                    {
                        biort_printf(VL_NORMAL, "%d %s zone reaction failed with a substep of %.1lf s.\n",
                          steps[kstep], "SURFACE", -substep);
                    } else {
                        biort_printf(VL_NORMAL, "%d %s zone reaction failed with a substep of %.1lf s.\n",
                          steps[kstep], (kzone == UZ) ? "Upper" : "Lower", -substep);
                    }
                }
                if (substep > 0.0)
                {
                    if (kzone == SURFACE)   // 2021-05-14
                    {
                      biort_printf(VL_VERBOSE, "%d %s zone reaction passed with a minimum step of %.1lf s.\n",
                        steps[kstep], "SURFACE", substep);
                    } else {
                      biort_printf(VL_VERBOSE, "%d %s zone reaction passed with a minimum step of %.1lf s.\n",
                        steps[kstep], (kzone == UZ) ? "Upper" : "Lower", substep);
                    }
                }
            }
        }
    }
}

int SolveReact(double stepsize, const chemtbl_struct chemtbl[], const kintbl_struct kintbl[], const rttbl_struct *rttbl,
    double satn, double temp, double porosity, double Zw, chmstate_struct *chms)
{
    int             i, j, k;
    int             kmonod, kinhib;
    int             control;
    int             min_pos;
    int             pivot_flg;
    double          monodterm = 1.0, inhibterm = 1.0;
    double          residue[MAXSPS];
    double          residue_t[MAXSPS];
    double          tmpconc[MAXSPS];
    double          tot_conc[MAXSPS];
    double          area[MAXSPS];
    double          ftemp[MAXSPS];
    double          fsw[MAXSPS];
    double          fzw[MAXSPS];
    double          gamma[MAXSPS];
    double          rate_pre[MAXSPS];
    double          iap[MAXSPS];
    double          dependency[MAXSPS];
    double          rate_spe[MAXSPS];
    double          rate_spet[MAXSPS];
    double          tmpval;
    double          inv_sat;
    double          imat;
    double          iroot;
    double          temp_keq;
    double          adh;
    double          bdh;
    double          bdt;
    double          max_error;
    double          tot_cec;
    realtype      **jcb;
    sunindextype    p[MAXSPS];
    realtype        x[MAXSPS];
    const int       SUFEFF = 1;
    const double    TMPPRB = 1.0E-2;
    const double    TMPPRB_INV = 1.0 / TMPPRB;

    inv_sat = 1.0 / satn;//inv_sat(L of porous space/L of water)= 1/sat(L of water/L of porous space)


    for (i = 0; i < rttbl->num_min; i++)
    {
        //area(m2/L of porous media) = ssa(m2/g of min) * mineral_conc(mol of min/L of porous media) * molar_mass(g of min/mol of min)
        area[i] = chms->ssa[i + rttbl->num_stc - rttbl->num_min] *
            chms->prim_conc[i + rttbl->num_stc - rttbl->num_min] *
            chemtbl[i + rttbl->num_stc - rttbl->num_min].molar_mass;
        ftemp[i] = SoilTempFactor(chms->q10[i + rttbl->num_stc - rttbl->num_min], temp);
        fzw[i] = WTDepthFactor(Zw, chms->n_alpha[i + rttbl->num_stc - rttbl->num_min]);
    }

    if (SUFEFF)
    {
        if (satn < 1.0)
        {
            for (i = 0; i < rttbl->num_min; i++)
            {
                fsw[i] = SoilMoistFactor(satn, chms->sw_thld[i + rttbl->num_stc - rttbl->num_min],chms->sw_exp[i + rttbl->num_stc - rttbl->num_min]);
            }
        }
        if (satn == 1.0)
        {
            for (i = 0; i < rttbl->num_min; i++)
            {
                fsw[i]= 1.0;
            }
        }
    }   // Lichtner's 2 third law if SUF_EFF is turned on

    for (i = 0; i < rttbl->num_stc; i++)
    {
        rate_spe[i] = 0.0;
    }

    for (i = 0; i < rttbl->num_mkr + rttbl->num_akr; i++)
    {
        min_pos = kintbl[i].position - rttbl->num_stc + rttbl->num_min;

        if (kintbl[i].type == TST)
        {
            iap[i] = 0.0;
            for (j = 0; j < rttbl->num_stc; j++)
            {
                iap[i] += log10(chms->prim_actv[j]) * rttbl->dep_kin[i][j];
            }
            iap[i] = pow(10, iap[i]);

            temp_keq = pow(10, rttbl->keq_kin[i]);

            dependency[i] = 1.0;
            for (k = 0; k < kintbl[i].ndep; k++)
            {
                dependency[i] *= pow(chms->prim_actv[kintbl[i].dep_index[k]], kintbl[i].dep_power[k]);
            }

            // Calculate the predicted rate depending on the type of rate law
            //   rate_pre: rate per reaction (mol L-1 porous media s-1)
            //   area: m2 of min/ L of porous media
            //   rate: mol/m2 of min/s
            //   dependency: dimensionless
            rate_pre[i] = area[min_pos] * pow(10, kintbl[i].rate) * dependency[i] * (1.0 - iap[i] / temp_keq) * ftemp[min_pos] * fsw[min_pos] * fzw[min_pos];
        }
        else if (kintbl[i].type == MONOD)
        {
            monodterm = 1.0;    // re-set for new species
            inhibterm = 1.0;    //re-set for new species

            // Calculate rate
            for (kmonod = 0; kmonod < kintbl[i].nmonod; kmonod++)
            {
                monodterm *= chms->prim_conc[kintbl[i].monod_index[kmonod]] /
                    (chms->prim_conc[kintbl[i].monod_index[kmonod]] + kintbl[i].monod_para[kmonod]);
            }

            for (kinhib = 0; kinhib < kintbl[i].ninhib; kinhib++)
            {
                inhibterm *= kintbl[i].inhib_para[kinhib] /
                    (kintbl[i].inhib_para[kinhib] + chms->prim_conc[kintbl[i].inhib_index[kinhib]]);
            }

            // Based on CrunchTope
            rate_pre[i] = area[min_pos] * pow(10, kintbl[i].rate) * monodterm * inhibterm * ftemp[min_pos] * fsw[min_pos] * fzw[min_pos];
        }

        for (j = 0; j < rttbl->num_stc; j++)
        {
            rate_spe[j] += rate_pre[i] * rttbl->dep_kin[i][j];
        }
    }

    for (i = 0; i < rttbl->num_mkr + rttbl->num_akr; i++)
    {
        min_pos = kintbl[i].position - rttbl->num_stc + rttbl->num_min;

        if (rate_pre[i] < 0.0)
        {
            // Mineral cutoff when mineral is disappearing
            area[min_pos] = (chms->prim_conc[kintbl[i].position] < 1.0E-8) ? 0.0 : area[min_pos];
        }
    }

    for (i = 0; i < rttbl->num_spc; i++)
    {
        // Aqueous species, saturation term for aqueous volume
        // rate(mol/m2 water/s)= rate(mol/m2 pm/s)*inv_sat(L of porous space/L of water)/porosity(L of porous space/L of pm)
        rate_spe[i] *= (chemtbl[i].itype == AQUEOUS) ? (inv_sat/porosity) : 1.0;
    }

    jcb = newDenseMat(rttbl->num_stc - rttbl->num_min, rttbl->num_stc - rttbl->num_min);

    adh = rttbl->adh;
    bdh = rttbl->bdh;
    bdt = rttbl->bdt;

    tot_cec = 0.0;
    imat = 0;
    for (i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        tmpconc[i] = (i < rttbl->num_stc) ? log10(chms->prim_conc[i]) : log10(chms->sec_conc[i - rttbl->num_stc]);

        tot_cec += (chemtbl[i].itype == CATION_ECHG) ? pow(10, tmpconc[i]) : 0.0;

        imat += 0.5 * pow(10, tmpconc[i]) * chemtbl[i].charge * chemtbl[i].charge;
    }
    iroot = sqrt(imat);

    for (i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        switch (chemtbl[i].itype)
        {
            case AQUEOUS:
                gamma[i] = (-adh * chemtbl[i].charge * chemtbl[i].charge * iroot) /
                    (1.0 + bdh * chemtbl[i].size_fac * iroot) + bdt * imat;
                break;
            case ADSORPTION:
                gamma[i] = log10(satn);
                break;
            case CATION_ECHG:
                gamma[i] = -log10(tot_cec);
                break;
            case MINERAL:
                gamma[i] = -tmpconc[i];
                break;
        }
    }

    control = 0;
    max_error = 0.0;
    do
    {
        for (i = 0; i < rttbl->num_ssc; i++)
        {
            tmpval = 0.0;
            for (j = 0; j < rttbl->num_sdc; j++)
            {
                tmpval += (tmpconc[j] + gamma[j]) * rttbl->dep_mtx[i][j];
            }
            tmpval -= rttbl->keq[i] + gamma[i + rttbl->num_stc];
            tmpconc[i + rttbl->num_stc] = tmpval;
        }

        for (j = 0; j < rttbl->num_stc; j++)
        {
            rate_spet[j] = 0.0;
        }

        for (i = 0; i < rttbl->num_mkr + rttbl->num_akr; i++)
        {
            min_pos = kintbl[i].position - rttbl->num_stc + rttbl->num_min;

            if (kintbl[i].type == TST)
            {
                iap[i] = 0.0;
                for (j = 0; j < rttbl->num_stc; j++)
                {
                    iap[i] += (chemtbl[j].itype != MINERAL) ? (tmpconc[j] + gamma[j]) * rttbl->dep_kin[i][j] : 0.0;
                }
                iap[i] = pow(10, iap[i]);
                temp_keq = pow(10, rttbl->keq_kin[i]);

#if NOT_YET_IMPLEMENTED
                if (iap[i] < temp_keq)
                {
                    rct_drct[i] = 1.0;
                }
                else if (iap[i] > temp_keq)
                {
                    rct_drct[i] = -1.0;
                }
                else
                {
                    rct_drct[i] = 0.0;
                }
#endif

                dependency[i] = 0.0;
                for (k = 0; k < kintbl[i].ndep; k++)
                {
                    dependency[i] += (tmpconc[kintbl[i].dep_index[k]] +
                        gamma[kintbl[i].dep_index[k]]) * kintbl[i].dep_power[k];
                }
                dependency[i] = pow(10, dependency[i]);

                // Calculate predicted rate depending on type of rate law
                // rate_pre: in mol / L water / s
                // area: m2/L water
                // rate: mol/m2/s
                // dependency: dimensionless
                rate_pre[i] = area[min_pos] * pow(10, kintbl[i].rate) * dependency[i] * (1.0 - (iap[i] / temp_keq)) * ftemp[min_pos] * fsw[min_pos] * fzw[min_pos];
            }
            else if (kintbl[i].type == MONOD)
            {
                monodterm = 1.0;
                inhibterm = 1.0;

                // Calculate rate
                for (kmonod = 0; kmonod < kintbl[i].nmonod; kmonod++)
                {
                    monodterm *= chms->prim_conc[kintbl[i].monod_index[kmonod]] /
                        (chms->prim_conc[kintbl[i].monod_index[kmonod]] +
                        kintbl[i].monod_para[kmonod]);
                }

                for (kinhib = 0; kinhib < kintbl[i].ninhib; kinhib++)
                {
                    inhibterm *= kintbl[i].inhib_para[kinhib] / (kintbl[i].inhib_para[kinhib] +
                        chms->prim_conc[kintbl[i].inhib_index[kinhib]]);
                }

                // Based on CrunchTope
                rate_pre[i] = area[min_pos] * pow(10, kintbl[i].rate) * monodterm * inhibterm * ftemp[min_pos] * fsw[min_pos] * fzw[min_pos];
            }

            for (j = 0; j < rttbl->num_stc; j++)
            {
                rate_spet[j] += rate_pre[i] * rttbl->dep_kin[i][j];
            }
            // Adjust the unit of the calculated rate. Note that for mineral, the unit of rate and the unit of
            // concentration are mol/L porous media. For the aqueous species, the unit of the rate and the unit of the
            // concentration are mol/L pm and mol/L water respectively.
        }

        for (i = 0; i < rttbl->num_spc; i++)
        {
            // rate(mol/m2 water/s)= rate(mol/m2 pm/s)*inv_sat(L of porous space/L of water)/porosity(L of porous space/L of pm)
            rate_spet[i] *= (chemtbl[i].itype == AQUEOUS) ? (inv_sat/porosity) : 1.0;
        }

        for (i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
        {
            tmpval = 0.0;
            for (j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
            {
                tmpval += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
            }
            tot_conc[i] = tmpval;
            residue[i] = tmpval - (chms->tot_conc[i] + (rate_spe[i] + rate_spet[i]) * stepsize * 0.5);
        }

        if (control % SKIP_JACOB == 0)
        {
            for (k = 0; k < rttbl->num_stc - rttbl->num_min; k++)
            {
                tmpconc[k] += TMPPRB;
                for (i = 0; i < rttbl->num_ssc; i++)
                {
                    tmpval = 0.0;
                    for (j = 0; j < rttbl->num_sdc; j++)
                    {
                        tmpval += (tmpconc[j] + gamma[j]) * rttbl->dep_mtx[i][j];
                    }
                    tmpval -= rttbl->keq[i] + gamma[i + rttbl->num_stc];
                    tmpconc[i + rttbl->num_stc] = tmpval;
                }
                for (i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
                {
                    tmpval = 0.0;
                    for (j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
                    {
                        tmpval += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
                    }
                    residue_t[i] = tmpval - (chms->tot_conc[i] + (rate_spe[i] + rate_spet[i]) * stepsize * 0.5);
                    jcb[k][i] = (residue_t[i] - residue[i]) * TMPPRB_INV;
                }
                tmpconc[k] -= TMPPRB;
            }
        }
        for (i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
        {
            x[i] = -residue[i];
        }

        pivot_flg = denseGETRF(jcb, rttbl->num_stc - rttbl->num_min, rttbl->num_stc - rttbl->num_min, p);
        if (pivot_flg != 0)
        {
            destroyMat(jcb);
            return 1;
        }

        denseGETRS(jcb, rttbl->num_stc - rttbl->num_min, p, x);

        max_error = 0.0;
        for (i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
        {
            if (fabs(x[i]) < 0.3)
            {
                tmpconc[i] += x[i];
            }
            else
            {
                tmpconc[i] += (x[i] < 0) ? -0.3 : 0.3;
            }

            max_error = MAX(fabs(residue[i] / tot_conc[i]), max_error);
        }

        control++;
        if (control > 10)
        {
            destroyMat(jcb);
            return 1;
        }
    } while (max_error > TOLERANCE);

    destroyMat(jcb);

    for (i = 0; i < rttbl->num_ssc; i++)
    {
        tmpval = 0.0;
        for (j = 0; j < rttbl->num_sdc; j++)
        {
            tmpval += (tmpconc[j] + gamma[j]) * rttbl->dep_mtx[i][j];
        }
        tmpval -= rttbl->keq[i] + gamma[i + rttbl->num_stc];
        tmpconc[i + rttbl->num_stc] = tmpval;
    }

    for (i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
    {
        tmpval = 0.0;
        for (j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
        {
            tmpval += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
        }
        tot_conc[i] = tmpval;
        residue[i] = tmpval - chms->tot_conc[i];
    }

    for (i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        if (i < rttbl->num_stc)
        {
            if (chemtbl[i].itype == MINERAL)
            {
                chms->tot_conc[i] += (rate_spe[i] + rate_spet[i]) * stepsize * 0.5;
                chms->prim_actv[i] = 1.0;
                chms->prim_conc[i] = chms->tot_conc[i];
            }
            else
            {
                chms->prim_conc[i] = pow(10, tmpconc[i]);
                chms->prim_actv[i] = pow(10, tmpconc[i] + gamma[i]);
                chms->tot_conc[i] = tot_conc[i];
            }
        }
        else
        {
            chms->sec_conc[i - rttbl->num_stc] = pow(10, tmpconc[i]);
#if TEMP_DISABLED
            chms->s_actv[i - rttbl->num_stc] = pow(10, tmpconc[i] + gamma[i]);
#endif
        }
    }

    return 0;
}

double ReactControl(const chemtbl_struct chemtbl[], const kintbl_struct kintbl[], const rttbl_struct *rttbl,
    double stepsize, double porosity, double depth, double satn, double temp, double Zw, double react_rate[],
    chmstate_struct *chms)
{
    int             flag;
    int             kspc;
    double          substep;
    double          step_counter = 0.0;
    double          conc0[MAXSPS];

    // Copy initial mineral concentration to array
    for (kspc = 0; kspc < rttbl->num_min; kspc++)
    {
        conc0[kspc] = chms->tot_conc[kspc + rttbl->num_stc - rttbl->num_min];
    }

    substep = stepsize;

    while (1.0 - step_counter / stepsize > 1.0E-10 && substep > 1.0E-30)
    {
        flag = SolveReact(substep, chemtbl, kintbl, rttbl, satn, temp, porosity, Zw, chms);

        if (flag == 0)
        {
            // Reaction passed with current step
            step_counter += substep;
        }
        else
        {
            substep *= 0.5;
        }
    }

    if (roundi(step_counter) != roundi(stepsize))   // Reactions fail
    {
        return -substep;
    }
    else    // Reactions succeed
    {
        for (kspc = 0; kspc < rttbl->num_min; kspc++)
        {
            // Calculate reaction rate (mole/m2-pm/day) = mol/L-pm/day * mm of depth   * (1m/1000mm) * (1000L/1m3)
            react_rate[kspc] =
                (chms->tot_conc[kspc + rttbl->num_stc - rttbl->num_min] - conc0[kspc]) * depth;

        }

        for (kspc = 0; kspc <rttbl->num_spc; kspc++)
        {
            chms->tot_mol[kspc] = chms->tot_conc[kspc] * porosity * satn * depth; // tot_mol (moles-mm of water/L water ) =tot_conc (mol/L water) * satn * porosity * depth
        }
        if (roundi(substep) != roundi(stepsize))
        {
            return substep;
        }
        else
        {
            return 0;
        }
    }
}

double SoilTempFactor(double q10, double stc)
{
    return pow(q10, (stc - 20.0) / 10.0);
}

double SoilMoistFactor(double satn, double sw_thld, double sw_exp)
{
    double fsw;
    fsw     =   (satn <= sw_thld) ?
                    pow(satn/sw_thld, sw_exp) : pow((1.0 - satn) / (1.0 - sw_thld), sw_exp);
    return fsw;
}

double WTDepthFactor(double Zw, double n_alpha){
    double fzw;
    fzw     =   (n_alpha==0) ? 1 : exp(-fabs(n_alpha)*pow(Zw,n_alpha/fabs(n_alpha)));
    return    fzw;
}
