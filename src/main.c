#include "biort.h"

int main(int argc, char *argv[])
{
    char            dir[MAXSTRING];         // name of input directory
    char            fn[MAXSTRING];          // name of output file
    int             nsub = 1;               // number of sub-catchments
    int             nsteps;                 // number of simulation steps
    int            *steps;                  // simulation steps
    int             kstep;                  // loop index
    int             kcycle;                 // loop index
    rttbl_struct    rttbl;
    chemtbl_struct  chemtbl[MAXSPS];
    kintbl_struct   kintbl[MAXSPS];
    subcatch_struct *subcatch;
    calib_struct    calib;
    ctrl_struct     ctrl;
    FILE           *fp;

    // Get name of directory from command line
    if (argc != 2)
    {
        printf("Please specify name of directory for your simulation.\n");
        return EXIT_FAILURE;
    }
    strcpy(dir, argv[1]);

    // Allocate
    subcatch = (subcatch_struct *)malloc(nsub * sizeof(subcatch_struct));

    calib.rate = 0.0;
    calib.xsorption = 0.0;
    calib.ssa = 1.0;

    // Read subcatchment soil parameters
    ReadSoil(dir, nsub, subcatch);

    // Read HBV simulation steps, water states and fluxes
    ReadHbvParam(dir, nsub, subcatch);
    ReadHbvResults(dir, nsub, &nsteps, &steps, subcatch);

    // Read chemistry control file
    ReadChem(dir, &ctrl, &rttbl, chemtbl, kintbl);

    // Read chemistry initial conditions
    ReadCini(dir, nsub, chemtbl, &rttbl, subcatch);

    // Initialize RT structures
    InitChem(dir, nsub, &calib, &ctrl, chemtbl, kintbl, &rttbl, subcatch);

    // Open output file
    sprintf(fn, "%s_results.txt", dir);
    fp = fopen(fn, "w");
    PrintHeader(fp, &rttbl, chemtbl);

    // Loop through forcing cycles
    for (kcycle = 0; kcycle <= ctrl.recycle; kcycle++)
    {
        // Loop through model steps to calculate reactive transport
        for (kstep = 0; kstep < nsteps; kstep++)
        {
            // Transport and routing
            Transpt(kstep, nsub, &rttbl, subcatch);

            // Transport changes total concentrations. Primary concentrations needs to be updated using total
            // concentrations
            UpdatePrimConc(nsub, &rttbl, &ctrl, subcatch);

            if (ctrl.transpt == KIN_REACTION)
            {
                // In reaction mode, simulate reaction for soil, and speciation for stream
                Reaction(kstep, nsub, 86400.0, chemtbl, kintbl, &rttbl, subcatch);

                StreamSpeciation(nsub, chemtbl, &ctrl, &rttbl, subcatch);
            }

            PrintDailyResults(fp, steps[kstep], nsub, &rttbl, subcatch);
        }
    }

    fclose(fp);
    FreeStruct(nsub, &steps, subcatch);
    free(subcatch);
}
