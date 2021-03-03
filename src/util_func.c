#include "biort.h"
#include "optparse.h"

int roundi(double x)
{
    return (int)((x < 0.0) ? x - 0.5 : x + 0.5);
}

void Wrap(char *str)
{
    char            word[MAXSTRING];

    sprintf(word, "'%s'", str);
    strcpy(str, word);
}

void Unwrap(const char wrapped_str[], char str[])
{
    int             i, j = 0;

    for (i = 0; i < (int)strlen(wrapped_str); i++)
    {
        if (wrapped_str[i] != '\'')
        {
            str[j] = wrapped_str[i];
            j++;
        }
    }

    str[j] = '\0';
}

void FreeStruct(int nsub, int nsteps, int *steps[], subcatch_struct subcatch[])
{
    int             ksub;
    int             kstep;

    free(*steps);

    for (ksub = 0; ksub < nsub; ksub++)
    {
        for (kstep = 0; kstep < nsteps; kstep++)
        {
            free(subcatch[ksub].ws[kstep]);
            free(subcatch[ksub].q[kstep]);
        }

        free(subcatch[ksub].ws);
        free(subcatch[ksub].q);
        free(subcatch[ksub].tmp);
    }
}

void ParseCmdLineParam(int argc, char *argv[], char dir[])
{
    int             option;
    struct optparse options;
    struct optparse_long longopts[] = {
        {"brief",      'b', OPTPARSE_NONE},
        {"silent",     's', OPTPARSE_NONE},
        {"version",    'V', OPTPARSE_NONE},
        {"verbose",    'v', OPTPARSE_NONE},
        {0, 0, 0}
    };

    optparse_init(&options, argv);

    while ((option = optparse_long(&options, longopts, NULL)) != -1)
    {
        switch (option)
        {
            case 'v':
                // Verbose mode
                verbose_mode = VL_VERBOSE;
                break;
            case 'b':
                // Brief mode
                verbose_mode = VL_BRIEF;
                break;
            case 's':
                // Silent mode
                verbose_mode = VL_SILENT;
                break;
            case 'V':
                // Print version number
                printf("HBV-BioRT Version %s\n", VERSION);
                exit(EXIT_SUCCESS);
                break;
            case '?':
                biort_printf(VL_ERROR, "Option not recognizable %s\n", options.errmsg);
                exit(EXIT_FAILURE);
                break;
            default:
                break;
        }

        fflush(stdout);
    }

    if (options.optind >= argc)
    {
        biort_printf(VL_ERROR, "Error:You must specify the name of input directory!\n"
            "Usage: ./biort [-b] [-v] [-V]"
            " <project directory>\n"
            "    -b Brief mode\n"
            "    -V Version number\n"
            "    -v Verbose mode\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        // Parse remaining arguments
        strcpy(dir, optparse_arg(&options));
    }
}