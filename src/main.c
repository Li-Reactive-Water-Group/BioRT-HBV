#include "biort.h"

int             verbose_mode;

int main(int argc, char *argv[])
{
    char            dir[MAXSTRING];         // name of input directory
    char            fn[MAXSTRING];          // name of output file
    char            timestr[MAXSTRING];     // time stamp
    int             nsub = 1;               // number of sub-catchments
    int             nsteps;                 // number of simulation steps
    int            *steps;                  // simulation steps
    int             nsteps_numexp;                 // number of simulation steps
    int            *steps_numexp;                  // simulation steps
    int             kstep;                  // loop index
    int             kcycle;                 // loop index
    time_t          rawtime;
    struct tm      *timestamp;
    rttbl_struct    rttbl;
    chemtbl_struct  chemtbl[MAXSPS];
    kintbl_struct   kintbl[MAXSPS];
    subcatch_struct *subcatch;
    subcatch_struct *subcatch_numexp;
    calib_struct    calib;
    ctrl_struct     ctrl;
    FILE           *fp;

    // Read command line arguments
    ParseCmdLineParam(argc, argv, dir);

    // Allocate
    subcatch = (subcatch_struct *)malloc(nsub * sizeof(subcatch_struct));
    subcatch_numexp = (subcatch_struct *)malloc(nsub * sizeof(subcatch_struct));

    calib.rate = 0.0;
    calib.xsorption = 0.0;
    calib.ssa = 1.0;

    // Read subcatchment soil parameters
    ReadSoil(dir, nsub, subcatch);

    // Read HBV simulation steps, water states and fluxes
    ReadHbvParam(dir, nsub, subcatch);
    ReadHbvResults(dir, nsub, &nsteps, &steps, subcatch, 0);

    // Read chemistry control file
    ReadChem(dir, &ctrl, &rttbl, chemtbl, kintbl);

    // Read chemistry initial conditions
    ReadCini(dir, nsub, chemtbl, &rttbl, subcatch);

    // Read time-series precipitation chemistry if defined in chem.txt  2021-05-20
    if (ctrl.precipchem == 1) {
        //printf("read in time-series precipitation %d", ctrl.precipchem);
        ReadPrecipChem(dir, nsub, &nsteps, &steps, subcatch, rttbl.num_stc, chemtbl, 0);
    }


    if (ctrl.precipchem_numexp == 1){
        ctrl.recycle = ctrl.recycle - 1;
        CopyConstSubcatchProp(nsub, subcatch, subcatch_numexp);
        ReadHbvResults(dir, nsub, &nsteps_numexp,  &steps_numexp, subcatch_numexp, 1);
        ReadPrecipChem(dir, nsub, &nsteps_numexp, &steps_numexp, subcatch_numexp, rttbl.num_stc, chemtbl, 1);
    }

    // Initialize RT structures
    InitChem(dir, nsub, &calib, &ctrl, chemtbl, kintbl, &rttbl, subcatch);

    // Create output directory when necessary
    mkdir("output");

    //
    time(&rawtime);
    timestamp = localtime(&rawtime);
    strftime(timestr, 11, "%y%m%d%H%M", timestamp);

    // Open output file
    sprintf(fn, "output/%s_results_%s.txt", dir, timestr);
    fp = fopen(fn, "w");
    PrintHeader(fp, ctrl.transpt, &rttbl, chemtbl);

    biort_printf(VL_NORMAL, "\nHBV-BioRT %s simulation started.\n", dir);

    // Loop through forcing cycles
    for (kcycle = 0; kcycle <= ctrl.recycle; kcycle++)
    {
        // Loop through model steps to calculate reactive transport
        for (kstep = 0; kstep < nsteps; kstep++)
        {
            // Transport and routing
            Transpt(kstep, nsub, chemtbl, &rttbl, &ctrl, subcatch); // 2021-05-20

            // Transport changes total concentrations. Primary concentrations needs to be updated using total
            // concentrations
            UpdatePrimConc(kstep, nsub, &rttbl, &ctrl, subcatch);
            
            StreamSpeciation(kstep, nsub, chemtbl, &ctrl, &rttbl, subcatch);

            PrintDailyResults(fp, ctrl.transpt, steps[kstep], nsub, &rttbl, subcatch);
            
            if (ctrl.transpt == KIN_REACTION)
            {
                // In reaction mode, simulate reaction for soil, and speciation for stream
                Reaction(kstep, nsub, 86400.0, steps, chemtbl, kintbl, &rttbl, subcatch);
            }
            else
            {
                Speciation(nsub, chemtbl, &ctrl, &rttbl, subcatch);
            }

            
        }
    }

    biort_printf(VL_NORMAL, "\nHBV-BioRT %s simulation succeeded.\n", dir);

    fclose(fp);

    if (ctrl.precipchem_numexp == 1){
        CopyInitChemSubcatch(nsub, &rttbl, subcatch, subcatch_numexp);

        //sprintf(fn, "output/%s_results_numexp.txt", dir);
        sprintf(fn, "output/%s_results_numexp_%s.txt", dir, timestr);
        fp = fopen(fn, "w");
        PrintHeader(fp, ctrl.transpt, &rttbl, chemtbl);

        biort_printf(VL_NORMAL, "\nHBV-BioRT %s numerical experiment started.\n", dir);

        for (kstep = 0; kstep < nsteps_numexp; kstep++){


            Transpt(kstep, nsub, chemtbl, &rttbl, &ctrl, subcatch_numexp); // 2021-05-20

            // Transport changes total concentrations. Primary concentrations needs to be updated using total
            // concentrations
            UpdatePrimConc(kstep, nsub, &rttbl, &ctrl, subcatch_numexp);
            
            StreamSpeciation(kstep, nsub, chemtbl, &ctrl, &rttbl, subcatch_numexp);

            PrintDailyResults(fp, ctrl.transpt, steps_numexp[kstep], nsub, &rttbl, subcatch_numexp);

            if (ctrl.transpt == KIN_REACTION)
            {
                // In reaction mode, simulate reaction for soil, and speciation for stream
                Reaction(kstep, nsub, 86400.0, steps, chemtbl, kintbl, &rttbl, subcatch_numexp);
            }
            else
            {
                Speciation(nsub, chemtbl, &ctrl, &rttbl, subcatch_numexp);
            }

        }

        biort_printf(VL_NORMAL, "\nHBV-BioRT %s numerical experiment succeeded.\n", dir);

        fclose(fp);
    }

    FreeStruct(nsub, nsteps, &steps, subcatch);
    FreeStruct(nsub, nsteps_numexp, &steps_numexp, subcatch_numexp);

    free(subcatch);
    free(subcatch_numexp);
}
