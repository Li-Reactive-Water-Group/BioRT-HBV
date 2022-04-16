#include "biort.h"

void CopyConstSubcatchProp(int nsub, const subcatch_struct subcatch[], subcatch_struct subcatch_numexp[])
{
    subcatch_numexp[0].k1 = subcatch[0].k1;
    subcatch_numexp[0].k2 = subcatch[0].k2;
    subcatch_numexp[0].maxbas = subcatch[0].maxbas;
    subcatch_numexp[0].perc = subcatch[0].perc;
    subcatch_numexp[0].tt = subcatch[0].tt;
    subcatch_numexp[0].sfcf = subcatch[0].sfcf;
    //subcatch_numexp[0].porosity_surface = subcatch[0].porosity_surface;
    subcatch_numexp[0].porosity_uz = subcatch[0].porosity_uz;
    subcatch_numexp[0].porosity_lz = subcatch[0].porosity_lz;
    //subcatch_numexp[0].res_surface = subcatch[0].res_surface;
    subcatch_numexp[0].res_uz = subcatch[0].res_uz;
    subcatch_numexp[0].res_lz = subcatch[0].res_lz;
    //subcatch_numexp[0].d_surface = subcatch[0].d_surface;
    subcatch_numexp[0].d_uz = subcatch[0].d_uz;
    subcatch_numexp[0].d_lz = subcatch[0].d_lz;
}

void CopyInitChemSubcatch(int nsub, rttbl_struct *rttbl, const subcatch_struct subcatch[], subcatch_struct subcatch_numexp[])
{
    int ksub;
    int kspc;

    for (ksub = 0; ksub < nsub; ksub++)
    {
        for (kspc = 0; kspc < rttbl->num_stc; kspc++)
        {
            subcatch_numexp[ksub].chms[SNOW].ssa[kspc] = subcatch[ksub].chms[SNOW].ssa[kspc];
            subcatch_numexp[ksub].chms[UZ].ssa[kspc] = subcatch[ksub].chms[UZ].ssa[kspc];
            subcatch_numexp[ksub].chms[LZ].ssa[kspc] = subcatch[ksub].chms[LZ].ssa[kspc];
            
            subcatch_numexp[ksub].chms[UZ].q10[kspc] = subcatch[ksub].chms[UZ].q10[kspc];
            subcatch_numexp[ksub].chms[LZ].q10[kspc] = subcatch[ksub].chms[LZ].q10[kspc];
            
            subcatch_numexp[ksub].chms[UZ].sw_thld[kspc] = subcatch[ksub].chms[UZ].sw_thld[kspc];
            subcatch_numexp[ksub].chms[LZ].sw_thld[kspc] = subcatch[ksub].chms[LZ].sw_thld[kspc];
            
            subcatch_numexp[ksub].chms[UZ].sw_exp[kspc] = subcatch[ksub].chms[UZ].sw_exp[kspc];
            subcatch_numexp[ksub].chms[LZ].sw_exp[kspc] = subcatch[ksub].chms[LZ].sw_exp[kspc];
            
            subcatch_numexp[ksub].chms[UZ].n_alpha[kspc] = subcatch[ksub].chms[UZ].n_alpha[kspc];
            subcatch_numexp[ksub].chms[LZ].n_alpha[kspc] = subcatch[ksub].chms[LZ].n_alpha[kspc];

            subcatch_numexp[ksub].chms[SNOW].prim_conc[kspc] = subcatch[ksub].chms[SNOW].prim_conc[kspc];
            subcatch_numexp[ksub].chms[UZ].prim_conc[kspc] = subcatch[ksub].chms[UZ].prim_conc[kspc];
            subcatch_numexp[ksub].chms[LZ].prim_conc[kspc] = subcatch[ksub].chms[LZ].prim_conc[kspc];
            subcatch_numexp[ksub].chms[STREAM].prim_conc[kspc] = subcatch[ksub].chms[STREAM].prim_conc[kspc];

            subcatch_numexp[ksub].chms[SNOW].prim_actv[kspc] = subcatch[ksub].chms[SNOW].prim_actv[kspc];
            subcatch_numexp[ksub].chms[UZ].prim_actv[kspc] = subcatch[ksub].chms[UZ].prim_actv[kspc];
            subcatch_numexp[ksub].chms[LZ].prim_actv[kspc] = subcatch[ksub].chms[LZ].prim_actv[kspc];
            subcatch_numexp[ksub].chms[STREAM].prim_actv[kspc] = subcatch[ksub].chms[STREAM].prim_actv[kspc];

            subcatch_numexp[ksub].chms[SNOW].tot_conc[kspc] = subcatch[ksub].chms[SNOW].tot_conc[kspc];
            subcatch_numexp[ksub].chms[UZ].tot_conc[kspc] = subcatch[ksub].chms[UZ].tot_conc[kspc];
            subcatch_numexp[ksub].chms[LZ].tot_conc[kspc] = subcatch[ksub].chms[LZ].tot_conc[kspc];
            subcatch_numexp[ksub].chms[STREAM].tot_conc[kspc] = subcatch[ksub].chms[STREAM].tot_conc[kspc];

            subcatch_numexp[ksub].chms[SNOW].tot_mol[kspc] = subcatch[ksub].chms[SNOW].tot_mol[kspc];
            subcatch_numexp[ksub].chms[UZ].tot_mol[kspc] = subcatch[ksub].chms[UZ].tot_mol[kspc];
            subcatch_numexp[ksub].chms[LZ].tot_mol[kspc] = subcatch[ksub].chms[LZ].tot_mol[kspc];
            subcatch_numexp[ksub].chms[STREAM].tot_mol[kspc] = subcatch[ksub].chms[STREAM].tot_mol[kspc];
        }

        for (kspc = 0; kspc < rttbl->num_ssc; kspc++)
        {
            subcatch_numexp[ksub].chms[SNOW].sec_conc[kspc] = subcatch[ksub].chms[SNOW].sec_conc[kspc];
            subcatch_numexp[ksub].chms[UZ].sec_conc[kspc] = subcatch[ksub].chms[UZ].sec_conc[kspc];
            subcatch_numexp[ksub].chms[LZ].sec_conc[kspc] = subcatch[ksub].chms[LZ].sec_conc[kspc];
        }

    }

}


//void RunNumExp(int nsteps, int nsub, const chemtbl_struct chemtbl[], const kintbl_struct kintbl[], rttbl_struct *rttbl, const ctrl_struct *ctrl, subcatch_struct subcatch[], FILE *fp)
//{
  //  int kstep;
    //for (kstep = 0; kstep < nsteps; kstep++)
    //{


        //biort_printf(VL_NORMAL, "\n%d.\n",kstep);
        // Transport and routing
        //Transpt(kstep, nsub, &rttbl, &ctrl, subcatch); //; // 2021-05-20



    //}
//}
