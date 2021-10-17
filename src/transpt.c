#include "biort.h"

// Calculates mass from transport
// To increase stability, transport is calculated using backward-in-time scheme, in which
//   d(CV) / dt = C_in * q_in - C * q_out + f
// i.e.,
//   C * V - C_0 * V_0 = C_in * Q_in - C * Q_out + F.
// Thus,
//   C = (C_in * Qin + F + C0 * V0) / (V + Q_out).
void Transpt(int step, int nsub, rttbl_struct *rttbl, const ctrl_struct *ctrl, subcatch_struct subcatch[])
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
            if (ctrl->precipchem == 0)  // constant precipitation chemistry mode
            {
                subcatch[ksub].chms[SNOW].tot_mol[kspc] += subcatch[ksub].prcp_conc[kspc] * subcatch[ksub].q[step][Psnow];
            } else if (ctrl->precipchem == 1)   // time-series precipitation chemistry mode
            {
                subcatch[ksub].chms[SNOW].tot_mol[kspc] += subcatch[ksub].prcp_conc_time[step][kspc] * subcatch[ksub].q[step][Psnow];
            }

            // Temperary concentration in snow

            if (subcatch[ksub].ws[step][SNOW] + subcatch[ksub].q[step][snowmelt] == 0)
            {
                conc_temp = 0;
                subcatch[ksub].chms[SNOW].tot_mol[kspc] = 0;
            }
            else
            {
                conc_temp = subcatch[ksub].chms[SNOW].tot_mol[kspc] /
                    (subcatch[ksub].ws[step][SNOW] + subcatch[ksub].q[step][snowmelt]);

            // Recharge
                subcatch[ksub].chms[SNOW].tot_mol[kspc] -= conc_temp * subcatch[ksub].q[step][snowmelt];
            }

            if ((subcatch[ksub].q[step][snowmelt] + subcatch[ksub].q[step][Prain]) == 0)
            {
                conc_temp = ZERO_CONC;
            } else if (ctrl->precipchem == 0)  // constant precipitation chemistry mode
            {
                subcatch[ksub].chms[UZ].tot_mol[kspc] += ((conc_temp * subcatch[ksub].q[step][snowmelt]) + (subcatch[ksub].prcp_conc[kspc] * subcatch[ksub].q[step][Prain]));
                conc_temp = ((conc_temp * subcatch[ksub].q[step][snowmelt]) + (subcatch[ksub].prcp_conc[kspc] * subcatch[ksub].q[step][Prain])) /
                    (subcatch[ksub].q[step][snowmelt] + subcatch[ksub].q[step][Prain]);
            } else if (ctrl->precipchem == 1)   // time-series precipitation chemistry mode
            {
                subcatch[ksub].chms[UZ].tot_mol[kspc] += ((conc_temp * subcatch[ksub].q[step][snowmelt]) + (subcatch[ksub].prcp_conc_time[step][kspc] * subcatch[ksub].q[step][Prain]));
                conc_temp = ((conc_temp * subcatch[ksub].q[step][snowmelt]) + (subcatch[ksub].prcp_conc_time[step][kspc] * subcatch[ksub].q[step][Prain])) /
                    (subcatch[ksub].q[step][snowmelt] + subcatch[ksub].q[step][Prain]);
            }


            subcatch[ksub].chms[SURFACE].tot_conc[kspc] = conc_temp;

            // Q0
            subcatch[ksub].chms[UZ].tot_mol[kspc] -= conc_temp * subcatch[ksub].q[step][Q0];
            //subcatch[ksub].chms[STREAM].tot_mol[kspc] = conc_temp * subcatch[ksub].q[step][Q0];

            // Temperary concentration in UZ
            conc_temp = subcatch[ksub].chms[UZ].tot_mol[kspc] /
                (subcatch[ksub].ws[step][UZ] + subcatch[ksub].q[step][Q1] + subcatch[ksub].q[step][PERC]);

            // Q1
            subcatch[ksub].chms[UZ].tot_mol[kspc] -= conc_temp * subcatch[ksub].q[step][Q1];
            //subcatch[ksub].chms[STREAM].tot_mol[kspc] += conc_temp * subcatch[ksub].q[step][Q1];

            // Percolation
            subcatch[ksub].chms[UZ].tot_mol[kspc] -= conc_temp * subcatch[ksub].q[step][PERC];
            subcatch[ksub].chms[LZ].tot_mol[kspc] += conc_temp * subcatch[ksub].q[step][PERC];

            // Temperary concentration in LZ
            conc_temp = subcatch[ksub].chms[LZ].tot_mol[kspc] /
                (subcatch[ksub].ws[step][LZ] + subcatch[ksub].q[step][Q2]);

            // Q2
            subcatch[ksub].chms[LZ].tot_mol[kspc] -= conc_temp * subcatch[ksub].q[step][Q2];
            //subcatch[ksub].chms[STREAM].tot_mol[kspc] += conc_temp * subcatch[ksub].q[step][Q2];

            // UPDATE CONCENTRATIONS
            subcatch[ksub].chms[SNOW].tot_conc[kspc] =
                (subcatch[ksub].ws[step][SNOW] == 0) ? ZERO_CONC : subcatch[ksub].chms[SNOW].tot_mol[kspc] / subcatch[ksub].ws[step][SNOW];
            //subcatch[ksub].chms[SURFACE].tot_conc[kspc] =
            //    (subcatch[ksub].ws[step][SNOW] == 0) ? ZERO_CONC : subcatch[ksub].chms[SURFACE].tot_mol[kspc] / subcatch[ksub].q[step][SURFACE];
            subcatch[ksub].chms[UZ].tot_conc[kspc] =
                subcatch[ksub].chms[UZ].tot_mol[kspc] / subcatch[ksub].ws[step][UZ];
            subcatch[ksub].chms[LZ].tot_conc[kspc] =
                subcatch[ksub].chms[LZ].tot_mol[kspc] / subcatch[ksub].ws[step][LZ];

            //subcatch[ksub].chms[STREAM].tot_conc[kspc] =
            //    (subcatch[ksub].q[step][Q0] + subcatch[ksub].q[step][Q1] + subcatch[ksub].q[step][Q2] > 0.0) ?
            //    subcatch[ksub].chms[STREAM].tot_mol[kspc] /
            //    (subcatch[ksub].q[step][Q0] + subcatch[ksub].q[step][Q1] + subcatch[ksub].q[step][Q2]) : ZERO_CONC;
            
           
            
        }
    }
}

void UpdatePrimConc(int step, int nsub, const rttbl_struct *rttbl, const ctrl_struct *ctrl, subcatch_struct subcatch[])
{
    int             ksub;
    int             kspc;

    for (ksub = 0; ksub < nsub; ksub++)
    {
        for (kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            subcatch[ksub].chms[SNOW].prim_conc[kspc] = subcatch[ksub].chms[SNOW].tot_conc[kspc];
            subcatch[ksub].chms[SURFACE].prim_conc[kspc] = subcatch[ksub].chms[SURFACE].tot_conc[kspc];
            if (ctrl->transpt == TRANSPORT_ONLY)
            {
                subcatch[ksub].chms[UZ].prim_conc[kspc] = subcatch[ksub].chms[UZ].tot_conc[kspc];
                subcatch[ksub].chms[LZ].prim_conc[kspc] = subcatch[ksub].chms[LZ].tot_conc[kspc];
            }
            subcatch[ksub].chms[STREAM].tot_conc[kspc] =
                (subcatch[ksub].q[step][Q0] + subcatch[ksub].q[step][Q1] + subcatch[ksub].q[step][Q2] > 0.0) ?
                (subcatch[ksub].chms[SURFACE].prim_conc[kspc] * subcatch[ksub].q[step][Q0] +  subcatch[ksub].chms[UZ].prim_conc[kspc] * subcatch[ksub].q[step][Q1] + subcatch[ksub].chms[LZ].prim_conc[kspc] * subcatch[ksub].q[step][Q2]) / (subcatch[ksub].q[step][Q0] + subcatch[ksub].q[step][Q1] + subcatch[ksub].q[step][Q2]) : ZERO_CONC;

            
            subcatch[ksub].chms[STREAM].prim_conc[kspc] = subcatch[ksub].chms[STREAM].tot_conc[kspc];
        }
    }
}
