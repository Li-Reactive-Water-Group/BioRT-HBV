#include "biort.h"

// Calculates mass in the snow + soil zone, and the transport of mass due to recharge
void Transpt(int step, int nsub, rttbl_struct *rttbl, subcatch_struct subcatch[])
{
    int             ksub;
    int             kspc;
    double          perc;
    double          conc_rechg;

    for (ksub = 0; ksub < nsub; ksub++)
    {
        if (step > 0)
        {
            // Percolation rate is not an output variable in HBV-Light. Calculate percolation rate first using mass
            //balance
            perc = subcatch[ksub].ws[step][LZ] - subcatch[ksub].ws[step - 1][LZ] + subcatch[ksub].q[step][Q2];

            perc = MIN(perc, subcatch[ksub].perc);
            perc = MAX(perc, 0.0);
        }

        // Transport for primary species
        for (kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            if (step > 0)
            {
                // Precipitation
                subcatch[ksub].chms[SNSM].tot_mol[kspc] +=
                    subcatch[ksub].prcp_conc[kspc] * subcatch[ksub].q[step][PRECIP];

                // Recharge
                // When there is water in snow/soil moisture storage, use snow/soil moisture concentration for recharge
                // When there is no water in snow/soil moisture storage, use precipitation concentration for recharge
                conc_rechg = (subcatch[ksub].ws[step - 1][SNSM] > 0.0) ?
                    subcatch[ksub].chms[SNSM].tot_conc[kspc] : subcatch[ksub].prcp_conc[kspc];

                subcatch[ksub].chms[SNSM].tot_mol[kspc] -= conc_rechg * subcatch[ksub].q[step][RECHG];
                subcatch[ksub].chms[UZ].tot_mol[kspc] += conc_rechg * subcatch[ksub].q[step][RECHG];

                // Q0
                // Note that Q0 should always use snow/soil concentration (or precipitation concentration if there is no
                // snow/soil water)
                subcatch[ksub].chms[UZ].tot_mol[kspc] -= conc_rechg * subcatch[ksub].q[step][Q0];

                // Q1
                subcatch[ksub].chms[UZ].tot_mol[kspc] -=
                    subcatch[ksub].chms[UZ].tot_conc[kspc] * subcatch[ksub].q[step][Q1];

                // Percolation
                subcatch[ksub].chms[UZ].tot_mol[kspc] -= subcatch[ksub].chms[UZ].tot_conc[kspc] * perc;
                subcatch[ksub].chms[LZ].tot_mol[kspc] += subcatch[ksub].chms[UZ].tot_conc[kspc] * perc;

                // Q2
                subcatch[ksub].chms[LZ].tot_mol[kspc] -=
                    subcatch[ksub].chms[LZ].tot_conc[kspc] * subcatch[ksub].q[step][Q2];
            }


            // ROUTING
            subcatch[ksub].chms[STREAM].tot_mol[kspc] = conc_rechg * subcatch[ksub].q[step][Q0];

            subcatch[ksub].chms[STREAM].tot_mol[kspc] +=
                subcatch[ksub].chms[UZ].tot_conc[kspc] * subcatch[ksub].q[step][Q1];

            subcatch[ksub].chms[STREAM].tot_mol[kspc] +=
                subcatch[ksub].chms[LZ].tot_conc[kspc] * subcatch[ksub].q[step][Q2];

            // UPDATE CONCENTRATIONS
            subcatch[ksub].chms[SNSM].tot_conc[kspc] = (subcatch[ksub].ws[step][SNSM] > 0.0) ?
                subcatch[ksub].chms[SNSM].tot_mol[kspc] / subcatch[ksub].ws[step][SNSM] : ZERO_CONC;
            subcatch[ksub].chms[UZ].tot_conc[kspc] =
                subcatch[ksub].chms[UZ].tot_mol[kspc] / subcatch[ksub].ws[step][UZ];
            subcatch[ksub].chms[LZ].tot_conc[kspc] =
                subcatch[ksub].chms[LZ].tot_mol[kspc] / subcatch[ksub].ws[step][LZ];

            subcatch[ksub].chms[STREAM].tot_conc[kspc] =
                (subcatch[ksub].q[step][Q0] + subcatch[ksub].q[step][Q1] + subcatch[ksub].q[step][Q2] > 0.0) ?
                subcatch[ksub].chms[STREAM].tot_mol[kspc] /
                (subcatch[ksub].q[step][Q0] + subcatch[ksub].q[step][Q1] + subcatch[ksub].q[step][Q2]) : ZERO_CONC;

        }
    }
}

void UpdatePrimConc(int nsub, const rttbl_struct *rttbl, const ctrl_struct *ctrl, subcatch_struct subcatch[])
{
    int             ksub;
    int             kspc;

    for (ksub = 0; ksub < nsub; ksub++)
    {
        for (kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            subcatch[ksub].chms[SNSM].prim_conc[kspc] = subcatch[ksub].chms[SNSM].tot_conc[kspc];

            if (ctrl->transpt == TRANSPORT_ONLY)
            {
                subcatch[ksub].chms[UZ].prim_conc[kspc] = subcatch[ksub].chms[UZ].tot_conc[kspc];
                subcatch[ksub].chms[LZ].prim_conc[kspc] = subcatch[ksub].chms[LZ].tot_conc[kspc];
            }

            subcatch[ksub].chms[STREAM].prim_conc[kspc] = subcatch[ksub].chms[STREAM].tot_conc[kspc];
        }
    }
}
