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

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "POROSITY_UZ", 'd', fn, lno, &sub[0].porosity_uz);
    biort_printf(VL_NORMAL, "  Upper zone porosity is %.2f m3 m-3.\n", sub[0].porosity_uz);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "POROSITY_LZ", 'd', fn, lno, &sub[0].porosity_lz);
    biort_printf(VL_NORMAL, "  Lower zone porosity is %.2f m3 m-3.\n", sub[0].porosity_lz);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "WS_IMMOBILE_UZ", 'd', fn, lno, &sub[0].ws_immobile_uz);
    sub[0].ws_immobile_uz = MAX(sub[0].ws_immobile_uz, STORAGE_MIN);
    biort_printf(VL_NORMAL, "  Upper zone immobile water storage is %.2f mm.\n", sub[0].ws_immobile_uz);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "WS_IMMOBILE_LZ", 'd', fn, lno, &sub[0].ws_immobile_lz);
    sub[0].ws_immobile_lz = MAX(sub[0].ws_immobile_lz, STORAGE_MIN);
    biort_printf(VL_NORMAL, "  Lower zone immobile water storage is %.2f mm.\n", sub[0].ws_immobile_lz);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "WS_MAX_UZ", 'd', fn, lno, &sub[0].ws_max_uz);
    biort_printf(VL_NORMAL, "  Upper zone maximum water storage capacity is %.2f mm.\n", sub[0].ws_max_uz);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "WS_MAX_LZ", 'd', fn, lno, &sub[0].ws_max_lz);
    biort_printf(VL_NORMAL, "  Lower zone maximum water storage capacity is %.2f mm.\n", sub[0].ws_max_lz);

    fclose(fp);
}
