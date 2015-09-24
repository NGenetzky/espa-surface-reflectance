#include <stdio.h> 
#include <stdlib.h> 
#include <ctype.h>
#include "espa_metadata.h"
#include "parse_metadata.h"

int main(int argc, char **argv)
{
    char *filename;      /* XML filename */
    int ib;              /* band looping variable */
    int rep_indx=-1;     /* band index in XML file for the current product */
    Espa_internal_meta_t xml_metadata;  /* XML metadata structure */
    Espa_global_meta_t *gmeta = NULL;   /* pointer to global metadata */
    Espa_band_meta_t *bmeta = NULL;     /* pointer to the band metadata array
                                           within the output structure */

    if (argc < 2)
    {
        printf("usage: %s <XML file>\n", argv[0]);
        exit (1);
    }
    filename = argv[1];

    /* Validate the input metadata file */
    if (validate_xml_file (filename) != SUCCESS)
    {  /* Error messages already written */
        exit (1);
    }

    /* Initialize the metadata structure */
    init_metadata_struct (&xml_metadata);

    /* Parse the metadata file into our internal metadata structure; also
       allocates space as needed for various pointers in the global and band
       metadata */
    if (parse_metadata (filename, &xml_metadata) != SUCCESS)
    {  /* Error messages already written */
        printf("Error parsing XML file: %s", filename);
        exit (1);
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
        printf("Error finding sr_band1 band in the XML file");
        exit (1);
    }
    bmeta = &xml_metadata.band[rep_indx];

    printf (":WestBoundingCoordinate = %lf\n",
        gmeta->bounding_coords[ESPA_WEST]);
    printf (":EastBoundingCoordinate = %lf\n",
        gmeta->bounding_coords[ESPA_EAST]);
    printf (":NorthBoundingCoordinate = %lf\n",
        gmeta->bounding_coords[ESPA_NORTH]);
    printf (":SouthBoundingCoordinate = %lf\n",
        gmeta->bounding_coords[ESPA_SOUTH]);
    printf (":AcquisitionDate = %sT%s\n", gmeta->acquisition_date,
      gmeta->scene_center_time);
    printf (":XDim_Grid = %d\n", bmeta->nsamps);
    printf (":YDim_Grid = %d\n", bmeta->nlines);
    printf (":SolarZenith = %f\n", gmeta->solar_zenith);
    printf (":SolarAzimuth = %f\n", gmeta->solar_azimuth);
    printf (":PixelSize = %f\n", bmeta->pixel_size[0]);

    /* Look for rep_band in the SR product and write the filename to the
       file */
    for (ib = 0; ib < xml_metadata.nbands; ib++)
    {
        if (!strcmp (xml_metadata.band[ib].name, "sr_band1") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            printf (":Band1 = %s\n", xml_metadata.band[ib].file_name);

        if (!strcmp (xml_metadata.band[ib].name, "sr_band2") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            printf (":Band2 = %s\n", xml_metadata.band[ib].file_name);

        if (!strcmp (xml_metadata.band[ib].name, "sr_band3") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            printf (":Band3 = %s\n", xml_metadata.band[ib].file_name);

        if (!strcmp (xml_metadata.band[ib].name, "sr_band4") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            printf (":Band4 = %s\n", xml_metadata.band[ib].file_name);

        if (!strcmp (xml_metadata.band[ib].name, "sr_band5") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            printf (":Band5 = %s\n", xml_metadata.band[ib].file_name);

        if ((!strcmp (xml_metadata.band[ib].name, "toa_band6") ||
             !strcmp (xml_metadata.band[ib].name, "toa_band61")) &&
            !strcmp (xml_metadata.band[ib].product, "toa_bt"))
            printf (":Band6 = %s\n", xml_metadata.band[ib].file_name);

        if (!strcmp (xml_metadata.band[ib].name, "sr_cloud_qa") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            printf (":CloudQA = %s\n", xml_metadata.band[ib].file_name);

        if (!strcmp (xml_metadata.band[ib].name, "sr_cloud_shadow_qa") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            printf (":CloudShadowQA = %s\n", xml_metadata.band[ib].file_name);

        if (!strcmp (xml_metadata.band[ib].name, "sr_adjacent_cloud_qa") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            printf (":CloudAdjQA = %s\n", xml_metadata.band[ib].file_name);

        if (!strcmp (xml_metadata.band[ib].name, "sr_snow_qa") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            printf (":SnowQA = %s\n", xml_metadata.band[ib].file_name);

        if (!strcmp (xml_metadata.band[ib].name, "sr_land_water_qa") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            printf (":LandWaterQA = %s\n", xml_metadata.band[ib].file_name);

        if (!strcmp (xml_metadata.band[ib].name, "sr_fill_qa") &&
            !strcmp (xml_metadata.band[ib].product, "sr_refl"))
            printf (":FillQA = %s\n", xml_metadata.band[ib].file_name);
    }

    exit (0);
}

