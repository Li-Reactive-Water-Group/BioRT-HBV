#include "biort.h"

void ReadSoil(const char dir[], int nsub, subcatch_struct sub[])
{
    char            cmdstr[MAXSTRING];
    char            fn[MAXSTRING];
    int             lno = 0;
    FILE           *fp;

    // READ SOIL.TXT FILE
    sprintf(fn, "input/%s/soil.txt", dir);
    fp = fopen(fn, "r");

    biort_printf(VL_NORMAL, "SOIL PARAMTERS\n");

    //NextLine(fp, cmdstr, &lno);   // 2021-05-14
    //ReadParam(cmdstr, "POROSITY_SURFACE", 'd', fn, lno, &sub[0].porosity_surface);
    //biort_printf(VL_NORMAL, "  Surface porosity is %.2f m3 m-3.\n", sub[0].porosity_surface);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "POROSITY_UZ", 'd', fn, lno, &sub[0].porosity_uz);
    biort_printf(VL_NORMAL, "  Upper zone porosity is %.2f m3 m-3.\n", sub[0].porosity_uz);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "POROSITY_LZ", 'd', fn, lno, &sub[0].porosity_lz);
    biort_printf(VL_NORMAL, "  Lower zone porosity is %.2f m3 m-3.\n", sub[0].porosity_lz);

    //NextLine(fp, cmdstr, &lno);  // 2021-05-14
    //ReadParam(cmdstr, "WS_PASSIVE_SURFACE", 'd', fn, lno, &sub[0].res_surface);
    //sub[0].res_surface = MAX(sub[0].res_surface, STORAGE_MIN);
    //biort_printf(VL_NORMAL, "  Surface passive water storage is %.2f mm.\n", sub[0].res_surface);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "WS_PASSIVE_UZ", 'd', fn, lno, &sub[0].res_uz);
    sub[0].res_uz = MAX(sub[0].res_uz, STORAGE_MIN);
    biort_printf(VL_NORMAL, "  Upper zone passive water storage is %.2f mm.\n", sub[0].res_uz);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "WS_PASSIVE_LZ", 'd', fn, lno, &sub[0].res_lz);
    sub[0].res_lz = MAX(sub[0].res_lz, STORAGE_MIN);
    biort_printf(VL_NORMAL, "  Lower zone passive water storage is %.2f mm.\n", sub[0].res_lz);

    //NextLine(fp, cmdstr, &lno);  // 2021-05-14
    //ReadParam(cmdstr, "DEPTH_SURFACE", 'd', fn, lno, &sub[0].d_surface);
    //biort_printf(VL_NORMAL, "  Depth of surface zone is %.2f mm.\n", sub[0].d_surface);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "DEPTH_UZ", 'd', fn, lno, &sub[0].d_uz);
    biort_printf(VL_NORMAL, "  Depth of upper zone is %.2f mm.\n", sub[0].d_uz);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "DEPTH_LZ", 'd', fn, lno, &sub[0].d_lz);
    biort_printf(VL_NORMAL, "  Depth of lower zone is %.2f mm.\n", sub[0].d_lz);

    fclose(fp);
}
