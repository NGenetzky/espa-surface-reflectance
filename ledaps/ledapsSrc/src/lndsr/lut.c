
#include <stdlib.h>
#include "lndsr.h"
#include "lut.h"
#include "input.h"
#include "mystring.h"
#include "error.h"

#define OUTPUT_FILL (-9999)
#define INPUT_SATU (20000)
#define OUTPUT_SATU (20000)
#define MIN_VALID_SR (-2000)
#define MAX_VALID_SR (16000)
#define AEROSOL_FILL (-9999)
#define AEROSOL_REGION_NLINE (40)
#define AEROSOL_REGION_NSAMP (AEROSOL_REGION_NLINE)
#define LONG_NAME_PREFIX ("band %d reflectance")
#define UNITS            ("reflectance")
#define SCALE_FACTOR     (0.0001)
#define ATMOS_OPACITY_SCALE_FACTOR (0.001)
#define ADD_OFFSET       (0.0)
#define SCALE_FACTOR_ERR (0.0)
#define ADD_OFFSET_ERR   (0.0)
#define CALIBRATED_NT    (DFNT_FLOAT32)

Lut_t *GetLut(int nband, Input_meta_t *meta, Img_coord_int_t *input_size) {
  Lut_t *this;

  /* Create the lookup table data structure */

  this = (Lut_t *)malloc(sizeof(Lut_t));
  if (this == NULL) 
    RETURN_ERROR("allocating Input data structure", "OpenInput", NULL);

  /* Populate the data structure */
  this->nband = nband;
  this->in_fill = meta->fill;
  this->output_fill = OUTPUT_FILL;
  this->in_satu = INPUT_SATU;
  this->output_satu = OUTPUT_SATU;
  this->aerosol_fill = AEROSOL_FILL;
  this->ar_region_size.l = AEROSOL_REGION_NLINE;
  this->ar_region_size.s = AEROSOL_REGION_NSAMP;
  this->ar_size.l = ((input_size->l - 1) / this->ar_region_size.l) + 1;
  this->ar_size.s = ((input_size->s - 1) / this->ar_region_size.s) + 1;
  this->min_valid_sr = MIN_VALID_SR;
  this->max_valid_sr = MAX_VALID_SR;
  this->atmos_opacity_scale_factor= ATMOS_OPACITY_SCALE_FACTOR;
  this->scale_factor=     SCALE_FACTOR;     /* scale factor            */
  this->scale_factor_err= SCALE_FACTOR_ERR; /* scale factor error      */
  this->add_offset=       ADD_OFFSET;       /* add offset              */
  this->add_offset_err=   ADD_OFFSET_ERR;   /* add offset error        */
  this->calibrated_nt=    CALIBRATED_NT;    /* calibrated nt           */

  this->long_name_prefix = DupString(LONG_NAME_PREFIX);
  if (this->long_name_prefix == NULL) {
    free(this);
    RETURN_ERROR("duplicating long name prefix", "GetLut", NULL);
  }

  this->units = DupString(UNITS);
  if (this->units == NULL) {
    free(this);
    RETURN_ERROR("duplicating ref units", "GetLut", NULL);
  }

  if (!InputMetaCopy(meta, nband, &this->meta)) {
    free(this);
    RETURN_ERROR("copying input metadata", "GetLut", NULL);
  }

  return this;
}

bool FreeLut(Lut_t *this) {

  if (this != NULL) {
    free(this);
  }

  return true;
}
