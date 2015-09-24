/*
!C****************************************************************************

!File: output.c
  
!Description: Functions creating and writting data to the product output file.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

 Revision 1.1 2012/09/27
 Gail Schmidt - Updated metadata to output information for each of the
 individual QA bands vs. a single packed bit QA band.

 Revision 1.2 2012/12/27
 Gail Schmidt - Modified to write the original LPGS L1G/L1T metadata filename
 for future reference by the user

 Revision 1.3 2013/01/22
 Gail Schmidt - Modified to utilize only one Version value - LEDAPSVersion

 Revision 1.4 2013/03/15
 Gail Schmidt - Removed NumberOfBands and BandNumbers from the lndsr output
 metadata.  Those only applied to the surface reflectance bands themselves
 and don't fully represent the final output bands for the lndsr product,
 therefore potentially leading to user confusion.

 Revision 1.5 2013/03/22
 Gail Schmidt, USGS EROS
 Modified to output the UL and LR corner lat/longs.  We are already writing
 the bounding coords, however for ascending scenes and scenes in the polar
 regions, the scenes are flipped upside down.  The bounding coords will be
 correct in North represents the northernmost latitude and South represents
 the southernmost latitude.  However, the UL corner in this case would be
 more south than the LR corner.  Comparing the UL and LR corners will allow
 the user to determine if the scene is flipped.

 Revision 1.5 2013/08/05
 Gail Schmidt, USGS EROS
 Modified to output the gain and bias values for the reflectance and
   brightness temperature bands.

 Revision 2.0 2014/02/03
 Gail Schmidt, USGS EROS
 Modified application to utilize the ESPA internal raw binary format.

!Team Unique Header:
  This software was developed by the MODIS Land Science Team Support 
  Group for the Laboratory for Terrestrial Physics (Code 922) at the 
  National Aeronautics and Space Administration, Goddard Space Flight 
  Center, under NASA Task 92-012-00.

 ! References and Credits:
  ! MODIS Science Team Member:
      Christopher O. Justice
      MODIS Land Science Team           University of Maryland
      justice@hermes.geog.umd.edu       Dept. of Geography
      phone: 301-405-1600               1113 LeFrak Hall
                                        College Park, MD, 20742

  ! Developers:
      Robert E. Wolfe (Code 922)
      MODIS Land Team Support Group     Raytheon ITSS
      robert.e.wolfe.1@gsfc.nasa.gov    4400 Forbes Blvd.
      phone: 301-614-5508               Lanham, MD 20770  
  
!END****************************************************************************
*/

#include <stdlib.h>
#include "output.h"
#include "input.h"
#include "error.h"


Output_t *OpenOutput(Espa_internal_meta_t *in_meta, Input_t *input,
  Param_t *param, Lut_t *lut)
/* 
!C******************************************************************************

!Description: 'OutputFile' sets up the 'output' data structure and opens the
 output file for write access.
 
!Input Parameters:
 in_meta        input XML metadata structure (band-related info)
 input          input structure with input image metadata (nband, iband, size)
 param          input paramter information (LEDAPS version)
 lut            lookup table (fill and saturation values)

!Output Parameters:
 (returns)      'output' data structure or NULL when an error occurs

!Revision:
revision on 1/27/2014 by Gail Schmidt, USGS/EROS
 - modified to utilize the ESPA internal raw binary structure

revision on 11/13/2014 by Gail Schmidt, USGS/EROS
 - modified to not write a fill_value for the output surface reflectance QA
   bands.  These are either on/off for that quality item, so technically there
   isn't a fill value.  The user will need to refer to the fill QA band.  Also
   ENVI doesn't include the fill value as part of the histogram, so you have
   to play weird games with the histograms in ENVI to be able to see both the
   0 and 255 when one of the values is flagged as fill.

revision on 1/15/2015 by Gail Schmidt, USGS/EROS
 - descriptions of the on/off QA values were flip-flopped. These have been
   fixed.

!END****************************************************************************
*/
{
  Output_t *this = NULL;       /* pointer to output structure */
  char *mychar = NULL;         /* pointer to '_' */
  char scene_name[STR_SIZE];   /* scene name for the current scene */
  int ib;             /* looping variables */
  int nband;          /* number of bands for this dataset */
  int nband_tot;      /* number of total bands with QA, for processing */
  int nband_out;      /* number of total bands with QA, for writing/output */
  int nband_out_extra; /* number of extra QA bands for writing/output */
  int rep_indx=-1;    /* band index in XML file for the current product */
  char production_date[MAX_DATE_LEN+1]; /* current date/time for production */
  time_t tp;          /* time structure */
  struct tm *tm;      /* time structure for UTC time */
  Espa_band_meta_t *bmeta = NULL;  /* pointer to the band metadata array
                         within the output structure */
  char *band_name_extra[NBAND_SR_EXTRA] = {"atmos_opacity", "fill_qa", "ddv_qa",
    "cloud_qa", "cloud_shadow_qa", "snow_qa", "land_water_qa",
    "adjacent_cloud_qa", "nb_dark_pixels", "avg_dark_sr_b7", "std_dark_sr_b7"};

  /* Determine the number of output bands. Don't plan to write the last 3 QA
     bands (nb_dark_pixels, avg_dark_sr_b7, or std_dark_sr_b7) */
  nband = input->nband;
  nband_tot = nband + NBAND_SR_EXTRA;
  nband_out = nband + NBAND_SR_EXTRA - 3;
  nband_out_extra = nband_out - nband;

  /* Check parameters */
  if (input->size.l < 1)
    RETURN_ERROR("invalid number of output lines", "OpenOutput", NULL);

  if (input->size.s < 1)
    RETURN_ERROR("invalid number of samples per output line", "OpenOutput",
      NULL);

  if (nband < 1 || nband > NBAND_REFL_MAX)
    RETURN_ERROR("invalid number of bands", "OpenOutput", NULL);

  /* Create the Output data structure */
  this = (Output_t *) malloc (sizeof(Output_t));
  if (this == NULL) 
    RETURN_ERROR("allocating Output data structure", "OpenOutput", NULL);

  /* Find the representative band for metadata information */
  for (ib = 0; ib < in_meta->nbands; ib++)
  {
    if (!strcmp (in_meta->band[ib].name, "toa_band1") &&
        !strcmp (in_meta->band[ib].product, "toa_refl"))
    {
      /* this is the index we'll use for band info from the XML strcuture */
      rep_indx = ib;
      break;
    }
  }
  if (rep_indx == -1)
    RETURN_ERROR("finding toa_band1 band in the XML file", "OpenOutput", NULL);

  /* Initialize the internal metadata for the output product. The global
     metadata won't be updated, however the band metadata will be updated
     and used later for appending to the original XML file. */
  init_metadata_struct (&this->metadata);

  /* Allocate memory for the total bands */
  if (allocate_band_metadata (&this->metadata, nband_out) != SUCCESS)
    RETURN_ERROR("allocating band metadata", "OpenOutput", NULL);
  bmeta = this->metadata.band;

  /* Determine the scene name */
  strcpy (scene_name, in_meta->band[rep_indx].file_name);
  mychar = strchr (scene_name, '_');
  if (mychar != NULL)
    *mychar = '\0';

  /* Get the current date/time (UTC) for the production date of each band */
  if (time (&tp) == -1)
    RETURN_ERROR("getting time", "OpenOutput", NULL);

  tm = gmtime (&tp);
  if (tm == NULL)
    RETURN_ERROR("converting time to UTC", "OpenOutput", NULL);

  if (strftime (production_date, MAX_DATE_LEN, "%Y-%m-%dT%H:%M:%SZ", tm) == 0)
    RETURN_ERROR("formating production date/time", "OpenOutput", NULL);

  /* Populate the data structure */
  this->open = false;
  this->nband_tot = nband_tot;
  this->nband_out = nband_out;
  this->size.l = input->size.l;
  this->size.s = input->size.s;
  for (ib = 0; ib < nband_out; ib++) {
    strncpy (bmeta[ib].short_name, in_meta->band[rep_indx].short_name, 3);
    bmeta[ib].short_name[3] = '\0';
    strcpy (bmeta[ib].product, "sr_refl");
    strcpy (bmeta[ib].source, "toa_refl");
    strcat (bmeta[ib].short_name, "SR");
    bmeta[ib].nlines = this->size.l;
    bmeta[ib].nsamps = this->size.s;
    bmeta[ib].pixel_size[0] = in_meta->band[rep_indx].pixel_size[0];
    bmeta[ib].pixel_size[1] = in_meta->band[rep_indx].pixel_size[1];
    strcpy (bmeta[ib].pixel_units, "meters");
    sprintf (bmeta[ib].app_version, "LEDAPS_%s", param->LEDAPSVersion);
    strcpy (bmeta[ib].production_date, production_date);

    if (ib < nband)  /* image band */
    {
      bmeta[ib].data_type = ESPA_INT16;
      bmeta[ib].fill_value = lut->output_fill;
      bmeta[ib].saturate_value = lut->output_satu;
      strcpy (bmeta[ib].category, "image");
      sprintf (bmeta[ib].name, "sr_band%d", input->meta.iband[ib]);
      bmeta[ib].scale_factor = lut->scale_factor;
      bmeta[ib].add_offset = lut->add_offset;
      sprintf (bmeta[ib].long_name, "band %d surface reflectance",
        input->meta.iband[ib]);
      strcpy (bmeta[ib].data_units, lut->units);
      bmeta[ib].valid_range[0] = lut->min_valid_sr;
      bmeta[ib].valid_range[1] = lut->max_valid_sr;
      bmeta[ib].calibrated_nt = lut->calibrated_nt;
    }
    else if (ib == nband)  /* atmospheric opacity */
    {
      bmeta[ib].data_type = ESPA_INT16;
      bmeta[ib].fill_value = lut->output_fill;
      strcpy (bmeta[ib].category, "image");
      sprintf (bmeta[ib].name, "sr_%s", band_name_extra[ib-nband]);
      bmeta[ib].scale_factor = lut->atmos_opacity_scale_factor;
      strcpy (bmeta[ib].long_name, band_name_extra[ib-nband]);
      strcpy (bmeta[ib].data_units, lut->units);
      bmeta[ib].valid_range[0] = lut->min_valid_sr;
      bmeta[ib].valid_range[1] = lut->max_valid_sr;
    }
    else  /* QA bands */
    {
      bmeta[ib].data_type = ESPA_UINT8;
      strcpy (bmeta[ib].category, "qa");
      sprintf (bmeta[ib].name, "sr_%s", band_name_extra[ib-nband]);
      strcpy (bmeta[ib].long_name, band_name_extra[ib-nband]);
      strcpy (bmeta[ib].data_units, "quality/feature classification");
      bmeta[ib].valid_range[0] = 0;
      bmeta[ib].valid_range[1] = 255;

      /* Set up QA bitmap information */
      if (allocate_class_metadata (&bmeta[ib], 2) != SUCCESS)
        RETURN_ERROR("allocating 2 classes", "OpenOutput", NULL); 

      bmeta[ib].class_values[0].class = 0;     /* off */
      bmeta[ib].class_values[1].class = 255;   /* on */
      switch (ib - nband_out_extra + 2) {
        case (FILL):
          strcpy (bmeta[ib].class_values[0].description, "not fill");
          strcpy (bmeta[ib].class_values[1].description, "fill");
          break;
        case (DDV):
          strcpy (bmeta[ib].class_values[0].description,
            "not dark dense vegetation");
          strcpy (bmeta[ib].class_values[1].description,
            "dark dense vegetation");
          break;
        case (CLOUD):
          strcpy (bmeta[ib].class_values[0].description, "not cloud");
          strcpy (bmeta[ib].class_values[1].description, "cloud");
          break;
        case (CLOUD_SHADOW):
          strcpy (bmeta[ib].class_values[0].description, "not cloud shadow");
          strcpy (bmeta[ib].class_values[1].description, "cloud shadow");
          break;
        case (SNOW):
          strcpy (bmeta[ib].class_values[0].description, "not snow");
          strcpy (bmeta[ib].class_values[1].description, "snow");
          break;
        case (LAND_WATER):
          strcpy (bmeta[ib].class_values[0].description, "land");
          strcpy (bmeta[ib].class_values[1].description, "water");
          break;
        case (ADJ_CLOUD):
          strcpy (bmeta[ib].class_values[0].description, "not adjacent cloud");
          strcpy (bmeta[ib].class_values[1].description, "adjacent cloud");
          break;
      }
    }

    /* Set up the filename with the scene name and band name and open the
       file for write access */
    sprintf (bmeta[ib].file_name, "%s_%s.img", scene_name,
      bmeta[ib].name);
    this->fp_bin[ib] = open_raw_binary (bmeta[ib].file_name, "w");
    if (this->fp_bin[ib] == NULL)
      RETURN_ERROR("unable to open output band file", "OpenOutput", NULL);
  }  /* for ib */
  this->open = true;

  /* Successful completion */
  return this;
}


bool CloseOutput(Output_t *this)
/* 
!C******************************************************************************

!Description: 'CloseOutput' the output files which are open.
 
!Input Parameters:
 this           'output' data structure

!Output Parameters:
 this           'output' data structure; the following fields are modified:
                   open
 (returns)      status:
                  'true' = okay
                  'false' = error return

!END****************************************************************************
*/
{
  int ib;

  if (!this->open)
    RETURN_ERROR("image files not open", "CloseOutput", false);

  for (ib = 0; ib < this->nband_out; ib++)
    close_raw_binary (this->fp_bin[ib]);

  this->open = false;
  return true;
}


bool FreeOutput(Output_t *this)
/* 
!C******************************************************************************

!Description: 'FreeOutput' frees the 'output' data structure memory.
 
!Input Parameters:
 this           'output' data structure for which the fields are freed

!Output Parameters:
 this           'output' data structure
 (returns)      status:
                  'true' = okay
                  'false' = error occurred

!END****************************************************************************
*/
{
  if (this->open) 
    RETURN_ERROR("file still open", "FreeOutput", false);

  free(this);
  this = NULL;

  return true;
}


bool PutOutputLine(Output_t *this, int iband, int iline, int16 *line)
/* 
!C******************************************************************************

!Description: 'PutOutputLine' writes a line of data to the output file.
 
!Input Parameters:
 this           'output' data structure
 iband          index (within Output_t struct) of output band to be written
 iline          output line number (used for validation only)
 line           buffer of data to be written (int16); if it's QA data then
                it will get converted to uint8 before writing

!Output Parameters:
 (returns)      status:
                  'true' = okay
                  'false' = error return

!END****************************************************************************
*/
{
  int ib;              /* looping variable */
  int nbytes = 0;      /* number of bytes in each pixel */
  uint8 *qabuf = NULL; /* buffer for QA data */
  Espa_band_meta_t *bmeta = this->metadata.band;  /* pointer to band metadata */
  void *void_buf = NULL;

  /* Check the parameters */
  if (this == NULL) 
    RETURN_ERROR("invalid input structure", "PutOutputLine", false);
  if (!this->open)
    RETURN_ERROR("file not open", "PutOutputLine", false);
  if (iband < 0 || iband >= this->nband_out)
    RETURN_ERROR("invalid band number", "PutOutputLine", false);
  if (iline < 0 || iline >= this->size.l)
    RETURN_ERROR("invalid line number", "PutOutputLine", false);

  /* Write the data, only the current line (i.e. one line at a time). If the
     output band is UINT8, then convert the input line to UINT8 before
     writing. */
  if (bmeta[iband].data_type == ESPA_INT16) {
    nbytes = sizeof (int16);
    void_buf = line;
  }
  else {
    nbytes = sizeof (uint8);
    qabuf = (uint8 *)calloc(this->size.s, nbytes);
    if (qabuf == NULL) 
      RETURN_ERROR("allocating output line buffer (uint8)", "PutOutputLine",
        false);
    for (ib = 0; ib < this->size.s; ib++)
      qabuf[ib] = line[ib];
    void_buf = qabuf;
  }

  if (write_raw_binary (this->fp_bin[iband], 1, this->size.s, nbytes, void_buf)
      != SUCCESS)
    RETURN_ERROR("writing output line", "PutOutputLine", false);

  return true;
}

