#include "biort.h"

// Calculates mass from transport
// To increase stability, transport is calculated using backward-in-time scheme, in which
//   d(CV) / dt = C_in * q_in - C * q_out + f
// i.e.,
//   C * V - C_0 * V_0 = C_in * Q_in - C * Q_out + F.
// Thus,
//   C = (C_in * Qin + F + C0 * V0) / (V + Q_out).
void Transpt(int step, int nsub, rttbl_struct *rttbl, subcatch_struct subcatch[])
{
    int             ksub;
    int             kspc;
    double          conc_temp;

    for (ksub = 0; ksub < nsub; ksub++)
    {
        // Transport for primary species
        for (kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            // Precipitation
            subcatch[ksub].chms[SNSM].tot_mol[kspc] +=
                subcatch[ksub].prcp_conc[kspc] * subcatch[ksub].q[step][PRECIP];

            // Temperary concentration in snow/soil
            conc_temp = subcatch[ksub].chms[SNSM].tot_mol[kspc] /
                (subcatch[ksub].ws[step][SNSM] + subcatch[ksub].q[step][RECHG]);

            // Recharge
            subcatch[ksub].chms[SNSM].tot_mol[kspc] -= conc_temp * subcatch[ksub].q[step][RECHG];
            subcatch[ksub].chms[UZ].tot_mol[kspc] += conc_temp * subcatch[ksub].q[step][RECHG];

            // Q0
            subcatch[ksub].chms[UZ].tot_mol[kspc] -= conc_temp * subcatch[ksub].q[step][Q0];
            subcatch[ksub].chms[STREAM].tot_mol[kspc] = conc_temp * subcatch[ksub].q[step][Q0];

            // Temperary concentration in UZ
            conc_temp = subcatch[ksub].chms[UZ].tot_mol[kspc] /
                (subcatch[ksub].ws[step][UZ] + subcatch[ksub].q[step][Q1] + subcatch[ksub].q[step][PERC]);

            // Q1
            subcatch[ksub].chms[UZ].tot_mol[kspc] -= conc_temp * subcatch[ksub].q[step][Q1];
            subcatch[ksub].chms[STREAM].tot_mol[kspc] += conc_temp * subcatch[ksub].q[step][Q1];

            // Percolation
            subcatch[ksub].chms[UZ].tot_mol[kspc] -= conc_temp * subcatch[ksub].q[step][PERC];
            subcatch[ksub].chms[LZ].tot_mol[kspc] += conc_temp * subcatch[ksub].q[step][PERC];

            // Temperary concentration in LZ
            conc_temp = subcatch[ksub].chms[LZ].tot_mol[kspc] /
                (subcatch[ksub].ws[step][LZ] + subcatch[ksub].q[step][Q2]);

            // Q2
            subcatch[ksub].chms[LZ].tot_mol[kspc] -= conc_temp * subcatch[ksub].q[step][Q2];
            subcatch[ksub].chms[STREAM].tot_mol[kspc] += conc_temp * subcatch[ksub].q[step][Q2];

            // UPDATE CONCENTRATIONS
            subcatch[ksub].chms[SNSM].tot_conc[kspc] =
                subcatch[ksub].chms[SNSM].tot_mol[kspc] / subcatch[ksub].ws[step][SNSM];
            subcatch[ksub].chms[UZ].tot_conc[kspc] =
                subcatch[ksub].chms[UZ].tot_mol[kspc] / subcatch[ksub].ws[step][UZ];
            subcatch[ksub].chms[LZ].tot_conc[kspc] =
                subcatch[ksub].chms[LZ].tot_mol[kspc] / subcatch[ksub].ws[step][LZ];

            //subcatch[ksub].chms[STREAM].tot_conc[kspc] =
                //(subcatch[ksub].q[step][Q0] + subcatch[ksub].q[step][Q1] + subcatch[ksub].q[step][Q2] > 0.0) ?
                //subcatch[ksub].chms[STREAM].tot_mol[kspc] /
                //(subcatch[ksub].q[step][Q0] + subcatch[ksub].q[step][Q1] + subcatch[ksub].q[step][Q2]) : ZERO_CONC;
                
            subcatch[ksub].chms[STREAM].tot_conc[kspc] = 
                (subcatch[ksub].q[step][Q0] + subcatch[ksub].q[step][Q1] + subcatch[ksub].q[step][Q2] > 0.0) ?
                (subcatch[ksub].chms[SNSM].prim_conc[kspc] * subcatch[ksub].q[step][Q0] +  subcatch[ksub].chms[UZ].prim_conc[kspc] * subcatch[ksub].q[step][Q1] + subcatch[ksub].chms[LZ].prim_conc[kspc] * subcatch[ksub].q[step][Q2]) / (subcatch[ksub].q[step][Q0] + subcatch[ksub].q[step][Q1] + subcatch[ksub].q[step][Q2]) : ZERO_CONC;

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
