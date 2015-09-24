/*****************************************************************************
FILE: lndsrbm.c
  
PURPOSE: Contains functions for handling the computation of the cloud-related
QA values using the surface reflectance and brightness temperature values.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
2/10/2014    Gail Schmidt     Original development (based on FORTRAN and ksh
                              code from original lndsrbm application)

NOTES:
  1. The XML metadata format read by this application follows the ESPA internal
     metadata format found in ESPA Raw Binary Format v1.0.doc.  The schema for
     the ESPA internal metadata format is available at
     http://espa.cr.usgs.gov/static/schema/espa_internal_metadata_v1_0.xsd.
*****************************************************************************/
#include <getopt.h>
#include <math.h>
#include <sys/stat.h>
#include "error_handler.h"
#include "espa_metadata.h"
#include "parse_metadata.h"
#include "raw_binary_io.h"

typedef signed short int16;
typedef unsigned char uint8;

/******************************************************************************
MODULE: usage

PURPOSE: Prints the usage information for this application.

RETURN VALUE:
Type = None

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
1/14/2014    Gail Schmidt     Original Development

NOTES:
******************************************************************************/
void usage ()
{
    printf ("lndsrbm updates the cloud-related QA values in the input "
            "product, using the already computed and available surface "
            "reflectance and brightness temperature values.\n\n");
    printf ("usage: lndsrbm "
            "--temp_center=scene_center_temperature_in_kelvin "
            "--dx=deltax --dy=deltay "
            "--xml=input_xml_filename\n");
    printf ("\nwhere the following parameters are required:\n");
    printf ("    -temp_center: temperature at the scene center (Kelvin)\n");
    printf ("    -dx: delta x (for northern adjustment)\n");
    printf ("    -dy: delta y (for northern adjustment)\n");
    printf ("    -xml: name of the input XML metadata file which follows "
            "the ESPA internal raw binary schema\n");
    printf ("\nExample: lndsrbm "
            "--center_temp=250.037186 --dx=-27.5052 --dy=91.7454 "
            "--xml=LE70230282011250EDC00.xml\n");
}


/******************************************************************************
MODULE:  get_args

PURPOSE:  Gets the command-line arguments and validates that the required
arguments were specified.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error getting the command-line arguments or a command-line
                argument and associated value were not specified
SUCCESS         No errors encountered

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
2/10/2014    Gail Schmidt     Original development
2/11/2014    Gail Schmidt     northern adjustment is now calculated from delta
                              x and y values; computation is borrowed from
                              compadjn.f

NOTES:
  1. Memory is allocated for the xml file.  This should be a character pointer
     set to NULL on input.  The caller is responsible for freeing the allocated
     memory upon successful return.
******************************************************************************/
short get_args
(
    int argc,             /* I: number of cmd-line args */
    char *argv[],         /* I: string of cmd-line args */
    float *center_temp,   /* O: address of the scene center temp (Kelvin) */
    float *north_adj,     /* O: addr of adjustment for true north (degrees) */
    char **xml_infile     /* O: address of input XML filename */
)
{
    int c;                           /* current argument index */
    int option_index;                /* index for the command-line option */
    float dx, dy;                    /* delta x and y values for northern
                                        adjustment computation */
    float fac;                       /* adjustment factor */
    char errmsg[STR_SIZE];           /* error message */
    char FUNC_NAME[] = "get_args";   /* function name */
    static struct option long_options[] =
    {
        {"xml", required_argument, 0, 'i'},
        {"center_temp", required_argument, 0, 't'},
        {"dx", required_argument, 0, 'x'},
        {"dy", required_argument, 0, 'y'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /* Loop through all the cmd-line options */
    dx = -9999.0;
    dy = -9999.0;
    opterr = 0;   /* turn off getopt_long error msgs as we'll print our own */
    while (1)
    {
        /* optstring in call to getopt_long is empty since we will only
           support the long options */
        c = getopt_long (argc, argv, "", long_options, &option_index);
        if (c == -1)
        {   /* Out of cmd-line options */
            break;
        }

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
     
            case 'h':  /* help */
                usage ();
                return (ERROR);
                break;

            case 'i':  /* XML infile */
                *xml_infile = strdup (optarg);
                break;
     
            case 't':  /* scene center temperature */
                *center_temp = atof (optarg);
                break;
     
            case 'x':  /* delta x */
                dx = atof (optarg);
                break;
     
            case 'y':  /* delta y */
                dy = atof (optarg);
                break;
     
            case '?':
            default:
                sprintf (errmsg, "Unknown option %s", argv[optind-1]);
                error_handler (true, FUNC_NAME, errmsg);
                usage ();
                return (ERROR);
                break;
        }
    }

    /* Make sure the XML file was specified */
    if (*xml_infile == NULL)
    {
        sprintf (errmsg, "XML input file is a required argument");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    /* Compute the northern adjustment */
    fac = atan (1.0) / 45.0;
    *north_adj = atan (dx / dy) / fac;


    return (SUCCESS);
}


/*****************************************************************************
MODULE: lndsrbm
  
PURPOSE: Reads in surface reflectance bands 1, 2, 3, 5, and 6 along with QA
bands for cloud, cloud shadow, adjacent cloud, snow, and fill.  Recomputes and
overwrites the cloud, cloud shadow, and adjacent cloud QA bands.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error computing the cloud QA
SUCCESS         Successfully computed the cloud QA

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
2/10/2014    Gail Schmidt     Rewrote the original cmrbv1.0.f in C code to
                              better work with the ESPA internal file format.
                              FORTRAN doesn't play well with raw binary files.
                              Modified to address the fact that TOA band 6
                              (temperature) data is now in Kelvin vs. degrees
                              Celsius

NOTES:
*****************************************************************************/
int main (int argc, char **argv)
{
    char errmsg[STR_SIZE];           /* error message */
    char FUNC_NAME[] = "lndsrbm";    /* function name */
    char *xml_infile = NULL; /* input XML filename */
    float center_temp;       /* scene center temp (Kelvin) */
    float tclear;            /* clear temperature (Celcius) */
    float t6;                /* band 6 temperature (Celcius) */
    float anom;              /* band 1 and 3 combination */
    float north_adj;         /* adjustment for true north (degrees) */
    float pclear;            /* percentage of pixels which are clear pixels */
    float cfac;              /* cloud factor */
    float dtr;               /* arctangent of 1.0 / 45. 0 */
    float facj;              /* cloud height factor in the line direction */
    float fack;              /* cloud height factor in the sample direction */
    float tcloud;   /* temperature of the current pixel (celsius) */
    float cldh;     /* cloud height (based on temperature of the cloud) */
    double mclear;           /* average/mean temp of the clear pixels */
    uint8 *cloud_qa = NULL;       /* cloud QA data */
    uint8 *cloud_shad_qa = NULL;  /* cloud shadow QA data */
    uint8 *cloud_adja_qa = NULL;  /* adjacent cloud QA data */
    uint8 *snow_qa = NULL;        /* snow QA data */
    uint8 *fill_qa = NULL;        /* fill QA data */
    uint8 *tmpbit_qa = NULL;      /* temporary bit for QA data */
    uint8 QA_OFF = 0;        /* value for QA turned off */
    uint8 QA_ON = 255;       /* value for QA turned on */
    int16 *band1 = NULL;     /* band 1 data */
    int16 *band2 = NULL;     /* band 2 data */
    int16 *band3 = NULL;     /* band 3 data */
    int16 *band5 = NULL;     /* band 5 data */
    int16 *band6 = NULL;     /* temperature (band6) data (Kelvin) */
    int mband5;         /* storage for the band 5 value */
    int mband5_pix;     /* storage for the pixel location of band 5 value */
    int rep_indx=-1;    /* band index in XML file for the current product */
    int cldhmin;        /* minimum bound of the cloud height */
    int cldhmax;        /* maximum bound of the cloud height */
    int icldh;          /* looping variable for cloud height */
    int ib;             /* looping variable for bands */
    int i;              /* looping variable for pixels */
    int il, is;         /* looping variables for lines and samples */
    int j, k;           /* looping variables for 3x3 and 5x5 blocks */
    int pix;            /* location of current pixel */
    long nbcloud;       /* count of the cloud pixels */
    long nbclear;       /* count of the clear (non-cloud) pixels */
    long nbval;         /* count of the non-fill pixels */
    FILE *cloud_fp = NULL;       /* cloud QA file */
    FILE *cloud_shad_fp = NULL;  /* cloud shadow QA file */
    FILE *cloud_adja_fp = NULL;  /* adjacent cloud QA file */
    FILE *snow_fp = NULL;        /* snow QA file */
    FILE *fill_fp = NULL;        /* fill QA file */
    FILE *band1_fp = NULL;       /* band 1 file */
    FILE *band2_fp = NULL;       /* band 2 file */
    FILE *band3_fp = NULL;       /* band 3 file */
    FILE *band5_fp = NULL;       /* band 5 file */
    FILE *band6_fp = NULL;       /* temperature (band6) file */
    Espa_internal_meta_t xml_metadata;  /* XML metadata structure to be
                                   populated by reading the XML metadata file */
    Espa_global_meta_t *gmeta = NULL;   /* pointer to global metadata */
    Espa_band_meta_t *bmeta = NULL;     /* pointer to the band metadata array
                                           within the output structure */

    /* Read the command-line arguments */
    if (get_args (argc, argv, &center_temp, &north_adj, &xml_infile) != SUCCESS)
    {   /* get_args already printed the error message */
        exit (ERROR);
    }
    printf ("north_adj: %f\n", north_adj);

    /* Validate the input metadata file */
    if (validate_xml_file (xml_infile) != SUCCESS)
    {  /* Error messages already written */
        exit (ERROR);
    }

    /* Initialize the metadata structure */
    init_metadata_struct (&xml_metadata);

    /* Parse the metadata file into our internal metadata structure; also
       allocates space as needed for various pointers in the global and band
       metadata */
    if (parse_metadata (xml_infile, &xml_metadata) != SUCCESS)
    {  /* Error messages already written */
        exit (ERROR);
    }
    gmeta = &xml_metadata.global;

    /* Look for band1 in the SR product and use for our representative band */
    for (ib = 0; ib < xml_metadata.nbands; ib++)
    {
        if (!strcmp (xml_metadata.band[ib].name, "sr_band1") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
        {
            /* this is the index we'll use for band info from the XML
               strcuture */
            rep_indx = ib;
            break;
        }
    }
    if (rep_indx == -1)
    {
        strcpy (errmsg, "Error finding sr_band1 band in the XML file");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    bmeta = &xml_metadata.band[rep_indx];

    /* Initialize the QA and band arrays to hold all the band data at a time */
    printf ("Allocating memory ...\n");
    cloud_qa = calloc (bmeta->nlines * bmeta->nsamps, sizeof (uint8));
    if (cloud_qa == NULL)
    {
        strcpy (errmsg, "Error allocating memory for cloud QA band.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    cloud_shad_qa = calloc (bmeta->nlines * bmeta->nsamps,
        sizeof (uint8));
    if (cloud_shad_qa == NULL)
    {
        strcpy (errmsg, "Error allocating memory for cloud shadow QA band.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    cloud_adja_qa = calloc (bmeta->nlines * bmeta->nsamps,
        sizeof (uint8));
    if (cloud_adja_qa == NULL)
    {
        strcpy (errmsg, "Error allocating memory for adjacent cloud QA band.");
        error_handler (true, FUNC_NAME, errmsg);
    }

    snow_qa = calloc (bmeta->nlines * bmeta->nsamps, sizeof (uint8));
    if (snow_qa == NULL)
    {
        strcpy (errmsg, "Error allocating memory for snow QA band.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    fill_qa = calloc (bmeta->nlines * bmeta->nsamps, sizeof (uint8));
    if (fill_qa == NULL)
    {
        strcpy (errmsg, "Error allocating memory for fill QA band.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    tmpbit_qa = calloc (bmeta->nlines * bmeta->nsamps, sizeof (uint8));
    if (tmpbit_qa == NULL)
    {
        strcpy (errmsg, "Error allocating memory for tmpbit QA band.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    band1 = calloc (bmeta->nlines * bmeta->nsamps, sizeof (int16));
    if (band1 == NULL)
    {
        strcpy (errmsg, "Error allocating memory for band 1 band.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    band2 = calloc (bmeta->nlines * bmeta->nsamps, sizeof (int16));
    if (band2 == NULL)
    {
        strcpy (errmsg, "Error allocating memory for band 2 band.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    band3 = calloc (bmeta->nlines * bmeta->nsamps, sizeof (int16));
    if (band3 == NULL)
    {
        strcpy (errmsg, "Error allocating memory for band 3 band.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    band5 = calloc (bmeta->nlines * bmeta->nsamps, sizeof (int16));
    if (band5 == NULL)
    {
        strcpy (errmsg, "Error allocating memory for band 5 band.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    band6 = calloc (bmeta->nlines * bmeta->nsamps, sizeof (int16));
    if (band6 == NULL)
    {
        strcpy (errmsg, "Error allocating memory for band 6 band.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Open the raw binary file for each of these QA, surface reflectance,
       and brightness temp bands.  The cloud-related files need to be open for
       reading and writing.  All other bands for read only. */
    printf ("Opening the input data ...\n");
    for (ib = 0; ib < xml_metadata.nbands; ib++)
    {
        if (!strcmp (xml_metadata.band[ib].name, "sr_band1") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            band1_fp = open_raw_binary (xml_metadata.band[ib].file_name, "rb");

        if (!strcmp (xml_metadata.band[ib].name, "sr_band2") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            band2_fp = open_raw_binary (xml_metadata.band[ib].file_name, "rb");

        if (!strcmp (xml_metadata.band[ib].name, "sr_band3") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            band3_fp = open_raw_binary (xml_metadata.band[ib].file_name, "rb");

        if (!strcmp (xml_metadata.band[ib].name, "sr_band5") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            band5_fp = open_raw_binary (xml_metadata.band[ib].file_name, "rb");

        if ((!strcmp (xml_metadata.band[ib].name, "toa_band6") ||
             !strcmp (xml_metadata.band[ib].name, "toa_band61")) &&
            !strcmp (xml_metadata.band[ib].product, "toa_bt"))
            band6_fp = open_raw_binary (xml_metadata.band[ib].file_name, "rb");

        if (!strcmp (xml_metadata.band[ib].name, "sr_cloud_qa") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            cloud_fp = open_raw_binary (xml_metadata.band[ib].file_name, "rb+");

        if (!strcmp (xml_metadata.band[ib].name, "sr_cloud_shadow_qa") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            cloud_shad_fp = open_raw_binary (xml_metadata.band[ib].file_name,
                "rb+");

        if (!strcmp (xml_metadata.band[ib].name, "sr_adjacent_cloud_qa") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            cloud_adja_fp = open_raw_binary (xml_metadata.band[ib].file_name,
                "rb+");

        if (!strcmp (xml_metadata.band[ib].name, "sr_snow_qa") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            snow_fp = open_raw_binary (xml_metadata.band[ib].file_name, "rb");

        if (!strcmp (xml_metadata.band[ib].name, "sr_fill_qa") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            fill_fp = open_raw_binary (xml_metadata.band[ib].file_name, "rb");
    }

    /* Make sure all the band and QA files are open */
    if (band1_fp == NULL)
    {
        strcpy (errmsg, "Error opening band 1 or obtaining filename from XML.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (band2_fp == NULL)
    {
        strcpy (errmsg, "Error opening band 2 or obtaining filename from XML.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (band3_fp == NULL)
    {
        strcpy (errmsg, "Error opening band 3 or obtaining filename from XML.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (band5_fp == NULL)
    {
        strcpy (errmsg, "Error opening band 5 or obtaining filename from XML.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (band6_fp == NULL)
    {
        strcpy (errmsg, "Error opening band 6 or obtaining filename from XML.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (cloud_fp == NULL)
    {
        strcpy (errmsg, "Error opening cloud QA or obtaining filename from "
            "XML.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (cloud_shad_fp == NULL)
    {
        strcpy (errmsg, "Error opening cloud shadow QA or obtaining filename "
            "from XML.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (cloud_adja_fp == NULL)
    {
        strcpy (errmsg, "Error opening cloud adjacent QA or obtaining filename "
            "from XML.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (snow_fp == NULL)
    {
        strcpy (errmsg, "Error opening snow QA or obtaining filename from "
            "XML.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (fill_fp == NULL)
    {
        strcpy (errmsg, "Error opening fill QA or obtaining filename from "
            "XML.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the band and QA files */
    printf ("Reading the input data ...\n");
    if (read_raw_binary (band1_fp, bmeta->nlines, bmeta->nsamps,
        sizeof (int16), band1) != SUCCESS)
    {
        strcpy (errmsg, "Reading band 1 image data.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (band2_fp, bmeta->nlines, bmeta->nsamps,
        sizeof (int16), band2) != SUCCESS)
    {
        strcpy (errmsg, "Reading band 2 image data.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (band3_fp, bmeta->nlines, bmeta->nsamps,
        sizeof (int16), band3) != SUCCESS)
    {
        strcpy (errmsg, "Reading band 3 image data.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (band5_fp, bmeta->nlines, bmeta->nsamps,
        sizeof (int16), band5) != SUCCESS)
    {
        strcpy (errmsg, "Reading band 5 image data.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (band6_fp, bmeta->nlines, bmeta->nsamps,
        sizeof (int16), band6) != SUCCESS)
    {
        strcpy (errmsg, "Reading band 6 image data.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (cloud_fp, bmeta->nlines, bmeta->nsamps,
        sizeof (uint8), cloud_qa) != SUCCESS)
    {
        strcpy (errmsg, "Reading cloud QA data.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (cloud_shad_fp, bmeta->nlines, bmeta->nsamps,
        sizeof (uint8), cloud_shad_qa) != SUCCESS)
    {
        strcpy (errmsg, "Reading cloud shadow QA data.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (cloud_adja_fp, bmeta->nlines, bmeta->nsamps,
        sizeof (uint8), cloud_adja_qa) != SUCCESS)
    {
        strcpy (errmsg, "Reading adjacent cloud QA data.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (snow_fp, bmeta->nlines, bmeta->nsamps,
        sizeof (uint8), snow_qa) != SUCCESS)
    {
        strcpy (errmsg, "Reading snow QA data.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (fill_fp, bmeta->nlines, bmeta->nsamps,
        sizeof (uint8), fill_qa) != SUCCESS)
    {
        strcpy (errmsg, "Reading fill QA data.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the non-cloud file pointers */
    close_raw_binary (snow_fp);
    close_raw_binary (fill_fp);
    close_raw_binary (band1_fp);
    close_raw_binary (band2_fp);
    close_raw_binary (band3_fp);
    close_raw_binary (band5_fp);
    close_raw_binary (band6_fp);

    /* Convert the center temp to celcius */
    tclear = center_temp - 273.15;

    /* Reset the cloud, cloud shadow, and adjacent cloud bits */
    for (i = 0; i < bmeta->nlines * bmeta-> nsamps; i++)
    {
       cloud_qa[i] = QA_OFF;
       cloud_adja_qa[i] = QA_OFF;
       cloud_shad_qa[i] = QA_OFF;
       tmpbit_qa[i] = QA_OFF;
    }

    /* Update the cloud mask */
    printf ("Updating cloud mask ...\n");
    nbclear = 0;
    mclear = 0.0;
    nbval = 0;
    nbcloud = 0;
      
    /* Compute the average temperature of the clear data (celsius) */
    for (i = 0; i < bmeta->nlines * bmeta-> nsamps; i++)
    {
        /* Only use the non-fill pixels for this average */
        if (fill_qa[i] == QA_OFF)
        {
            nbval++;
            anom = band1[i] - band3[i] * 0.5;
            t6 = band6[i] * 0.1 - 273.15;   /* convert to celsius */
            if (snow_qa[i] == QA_ON)
                continue;
            else
            {  /* not snow */
                if ((anom > 300.0) && (band5[i] > 300.0) && (t6 < tclear))
                {  /* cloud */
                    cloud_qa[i] = QA_ON;
                    nbcloud++;
                }
                else
                {
                    if ((band1[i] > 3000.0) && (t6 < tclear))
                    {  /* cloud */
                        cloud_qa[i] = QA_ON;
                        nbcloud++;
                    }
                    else
                    {  /* not cloud (clear) */
                        mclear = mclear + t6 * 0.0001;
                        nbclear++;
                    }
                }
            }
        }
    }

    /* Determine the average/mean temp of the clear pixels and the percentage */
    if (nbclear > 0)
        mclear = mclear * 10000.0 / nbclear;
    pclear = nbclear * 100.0 / nbval;

    /* Write the info for the user */
    printf ("number of non-fill pixels: %ld\n", nbval);
    printf ("number of cloud pixels: %ld\n", nbcloud);
    printf ("number of non-cloud pixels: %ld\n", nbclear);
    printf ("nbclear-nbval: %ld\n", nbclear-nbval);
    printf ("average clear temperature (celcius): %f\n", mclear);
    printf ("%% clear pixels: %f\n", pclear);

    /* Reset the clear temp if greater than 5% of the pixels are clear */
    if (pclear > 5.0)
        tclear = mclear;
             
    /* Update the adjacent cloud bit; only set for non-fill pixels */
    printf ("Updating adjacent cloud bit ...\n");
    for (il = 0; il < bmeta->nlines; il++)
    {
        for (is = 0; is < bmeta->nsamps; is++)
        {
            /* Look for adjacent clouds of cloud pixels */
            pix = il * bmeta->nsamps + is;
            if (cloud_qa[pix] == QA_ON)
            {
                /* Look at the adjacent 5x5 block for adjacent clouds */
                for (j = il-5; j <= il+5; j++)
                {
                    /* Make sure the current line is within the bounds of the
                       image */
                    if (j >= 0 && j < bmeta->nlines)
                    {
                        pix = j * bmeta->nsamps + (is-5);
                        for (k = is-5; k <= is+5; k++, pix++)
                        {
                            /* Make sure the current sample is within the
                               bounds of the image */
                            if (k >= 0 && k < bmeta->nsamps)
                            {
                                /* If this pixel is not cloud or fill then set
                                   it to adjacent cloud */
                                if (cloud_qa[pix] != QA_ON &&
                                    fill_qa[pix] != QA_ON)
                                    cloud_adja_qa[pix] = QA_ON;
                            }
                        }  /* end for k */
                    }
                }  /* end for j */
            }
        }  /* end for is */
    }  /* end for il */
       
    /* Compute the cloud shadow (using temp in degrees Celsius) */
    mband5 = 9999;
    cfac = 6.0;
    dtr = atan (1.0) / 45.0;
    printf ("Looking for cloud shadow ...\n");
    facj = cos (gmeta->solar_azimuth * dtr) * tan (gmeta->solar_zenith * dtr) /
        bmeta->pixel_size[0];
    fack = sin (gmeta->solar_azimuth * dtr) * tan (gmeta->solar_zenith * dtr) /
        bmeta->pixel_size[0];
    for (il = 0; il < bmeta->nlines; il++)
    {
        for (is = 0; is < bmeta->nsamps; is++)
        {
            pix = il * bmeta->nsamps + is;
            if (cloud_qa[pix] == QA_ON)
            {
                /* Convert the temperature to celsius and use that to compute
                   the cloud height */
                tcloud = band6[pix] * 0.1 - 273.15;
                cldh = (tclear - tcloud) * 1000.0 / cfac;
                if (cldh < 0.0)
                    cldh = 0.0;
                cldhmin = (int) (cldh - 1000.0);
                cldhmax = (int) (cldh + 1000.0);
                mband5 = 9999.0;
                mband5_pix = 9999;

                /* Loop through the min to max cloud height values and
                   determine the cloud shadows from the height and sun
                   angle factors */
                for (icldh = cldhmin * 0.1; icldh <= cldhmax * 0.1; icldh++)
                {
                    cldh = icldh * 10.0;
                    j = (int) (il + facj * cldh);
                    k = (int) (is - fack * cldh);

                    /* Make sure the current pixel is within the bounds
                       of the image */
                    if ((j >= 0) && (j < bmeta->nlines) &&
                        (k >= 0) && (k < bmeta->nsamps))
                    {
                        pix = j * bmeta->nsamps + k;
                        if ((band5[pix] < 800.0) &&
                            (band2[pix] - band3[pix] < 100.0))
                        {
                            if ((cloud_adja_qa[pix] != QA_ON) &&
                                (cloud_qa[pix] != QA_ON) &&
                                (cloud_shad_qa[pix] != QA_ON) &&
                                (fill_qa[pix] != QA_ON))
                            {
                                /* Store the value of band5 as well as the
                                   pixel location */
                                if (band5[pix] < mband5)
                                {
                                    mband5 = band5[pix];
                                    mband5_pix = pix;
                                }
                            }
                        }
                    }
                }  /* end for icldh */

                if (mband5 < 9999)
                {
                    pix = mband5_pix;
                    cloud_shad_qa[pix] = QA_ON;
                }
            }
        }  /* end for is */
    }  /* end for il */

    /* Dilate the cloud shadow */
    printf ("Dilating cloud shadow ...\n");
    for (il = 0; il < bmeta->nlines; il++)
    {
        for (is = 0; is < bmeta->nsamps; is++)
        {
            /* Look for cloud shadows of cloud pixels */
            pix = il * bmeta->nsamps + is;
            if (cloud_shad_qa[pix] == QA_ON)
            {
                /* Look at the adjacent 3x3 block for adjacent clouds */
                for (j = il-3; j < il+3; j++)
                {
                    /* Make sure the current line is within the bounds of the
                       image */
                    if (j >= 0 && j < bmeta->nlines)
                    {
                        pix = j * bmeta->nsamps + (is-3);
                        for (k = is-3; k < is+3; k++, pix++)
                        {
                            /* Make sure the current pixel is within the bounds
                               of the image */
                            if (k >= 0 && k < bmeta->nsamps)
                            {
                                /* If this pixel is not cloud, adjacent cloud,
                                   cloud shadow, or fill then set the tmpbit
                                   to on */
                                if ((cloud_adja_qa[pix] != QA_ON) &&
                                    (cloud_qa[pix] != QA_ON) &&
                                    (cloud_shad_qa[pix] != QA_ON) &&
                                    (fill_qa[pix] != QA_ON))
                                    tmpbit_qa[pix] = QA_ON;
                            }
                        }  /* end for k */
                    }
                }  /* end for j */
            }
        }  /* end for is */
    }  /* end for il */

    /* Update the cloud shadow and clear QA for fill pixels */
    printf ("Updating cloud shadow ...\n");
    for (i = 0; i < bmeta->nlines * bmeta-> nsamps; i++)
    {
        /* Update the cloud shadow */
        if (tmpbit_qa[i] == QA_ON)
        {
            cloud_shad_qa[i] = QA_ON;
            tmpbit_qa[i] = QA_OFF;
        }

        /* Clear the QA info for the fill pixels */
        if (fill_qa[i] == QA_ON)
        {
            cloud_qa[i] = QA_OFF;
            cloud_adja_qa[i] = QA_OFF;
            cloud_shad_qa[i] = QA_OFF;
        }
    }

    /* Reset the cloud file pointers so they are reading for updating with
       the new QA values */
    rewind (cloud_fp);
    rewind (cloud_shad_fp);
    rewind (cloud_adja_fp);

    /* Write the updated cloud, cloud shadow, and adjacent cloud QA values back
       to the file */
    if (write_raw_binary (cloud_fp, bmeta->nlines, bmeta->nsamps,
        sizeof (uint8), cloud_qa) != SUCCESS)
    {
        strcpy (errmsg, "Updating cloud QA file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (write_raw_binary (cloud_shad_fp, bmeta->nlines, bmeta->nsamps,
        sizeof (uint8), cloud_shad_qa) != SUCCESS)
    {
        strcpy (errmsg, "Updating cloud shadow QA file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (write_raw_binary (cloud_adja_fp, bmeta->nlines, bmeta->nsamps,
        sizeof (uint8), cloud_adja_qa) != SUCCESS)
    {
        strcpy (errmsg, "Updating adjacent cloud QA file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the cloud file pointers */
    close_raw_binary (cloud_fp);
    close_raw_binary (cloud_shad_fp);
    close_raw_binary (cloud_adja_fp);

    /* Free the data pointers */
    free (cloud_qa);
    free (cloud_shad_qa);
    free (cloud_adja_qa);
    free (snow_qa);
    free (fill_qa);
    free (tmpbit_qa);
    free (band1);
    free (band2);
    free (band3);
    free (band5);
    free (band6);

    /* Free the metadata structure */
    free_metadata (&xml_metadata);
    free (xml_infile);

    /* Successful completion */
    exit (SUCCESS);
}

