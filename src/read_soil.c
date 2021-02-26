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

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "POROSITY_UZ", 'd', fn, lno, &sub[0].porosity_uz);
    printf("  Upper zone porosity is %.2f m3 m-3.\n", sub[0].porosity_uz);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "POROSITY_LZ", 'd', fn, lno, &sub[0].porosity_lz);
    printf("  Lower zone porosity is %.2f m3 m-3.\n", sub[0].porosity_lz);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "RES_UZ", 'd', fn, lno, &sub[0].res_uz);
    sub[0].res_uz = MAX(sub[0].res_uz, 1.0);
    printf(" Upper zone residual moisture is %.2f mm.\n", sub[0].res_uz);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "RES_LZ", 'd', fn, lno, &sub[0].res_lz);
    sub[0].res_lz = MAX(sub[0].res_lz, 1.0);
    printf(" Lower zone residual moisture is %.2f mm.\n", sub[0].res_lz);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "D_UZ", 'd', fn, lno, &sub[0].d_uz);
    printf("  Upper zone depth is %.2f m.\n", sub[0].d_uz);

    NextLine(fp, cmdstr, &lno);
    ReadParam(cmdstr, "D_LZ", 'd', fn, lno, &sub[0].d_lz);
    printf("  Lower zone depth is %.2f m.\n", sub[0].d_lz);

    fclose(fp);
}
