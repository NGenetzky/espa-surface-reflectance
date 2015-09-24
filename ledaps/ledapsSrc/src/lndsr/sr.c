#include "sr.h"
#include "ar.h"
#include "const.h"
#include "sixs_runs.h"

/* !Revision:
 *
 * revision 1.2.1 3/22/2013  Gail Schmidt, USGS
 * - writing UL and LR corners to the output metadata to be able to detect
 *   ascending scenes or scenes where the image is flipped North to South
 * revision 8/11/2015  Gail Schmidt, USGS
 * - input saturated pixels are flagged as such and output as saturated
 */

extern atmos_t atmos_coef;
int SrInterpAtmCoef(Lut_t *lut, Img_coord_int_t *input_loc, atmos_t *atmos_coef,atmos_t *interpol_atmos_coef); 

bool Sr(Lut_t *lut, int nsamp, int il, int16 **line_in, int16 **line_out,
        Sr_stats_t *sr_stats) 
{
  int is;
  bool is_fill;      /* is the pixel fill */
  bool is_satu;      /* is the pixel saturated */
  int ib;
  Img_coord_int_t loc;
  float rho,tmpflt;
  atmos_t interpol_atmos_coef;

  allocate_mem_atmos_coeff(1,&interpol_atmos_coef);
  loc.l = il;

  for (is = 0; is < nsamp; is++) {
    loc.s = is;
    is_fill = false;
    is_satu = false;

/*
NAZMI 6/2/04 : correct even cloudy pixels

*/
	 SrInterpAtmCoef(lut, &loc, &atmos_coef,&interpol_atmos_coef);

    for (ib = 0; ib < lut->nband; ib++) {
      if (line_in[ib][is] == lut->in_fill) {
        /* fill pixel */
        is_fill = true;
        line_out[ib][is] = lut->output_fill;
        sr_stats->nfill[ib]++;
      }
      else if (line_in[ib][is] == lut->in_satu) {
        /* saturated pixel */
        is_satu = true;
        line_out[ib][is] = lut->output_satu;
        sr_stats->nsatu[ib]++;
      }
      else {
        rho=(float)line_in[ib][is]/10000.;
        rho=(rho/interpol_atmos_coef.tgOG[ib][0]-interpol_atmos_coef.rho_ra[ib][0]);
        tmpflt=(interpol_atmos_coef.tgH2O[ib][0]*interpol_atmos_coef.td_ra[ib][0]*interpol_atmos_coef.tu_ra[ib][0]);
        rho /= tmpflt;
        rho /= (1.+interpol_atmos_coef.S_ra[ib][0]*rho);

        line_out[ib][is] = (short)(rho*10000.);  /* scale for output */

        if (line_out[ib][is] < lut->min_valid_sr) {
          sr_stats->nout_range[ib]++;
          line_out[ib][is] = lut->min_valid_sr;
        }
        if (line_out[ib][is] > lut->max_valid_sr) {
          sr_stats->nout_range[ib]++;
          line_out[ib][is] = lut->max_valid_sr;
        }
      }

      if (is_fill || is_satu) continue;

      if (sr_stats->first[ib]) {
        sr_stats->sr_min[ib] = sr_stats->sr_max[ib] = line_out[ib][is];
        sr_stats->first[ib] = false;
      } else {
        if (line_out[ib][is] < sr_stats->sr_min[ib])
          sr_stats->sr_min[ib] = line_out[ib][is];

        if (line_out[ib][is] > sr_stats->sr_max[ib])
          sr_stats->sr_max[ib] = line_out[ib][is];
      } 
    }  /* end for */

  }

  free_mem_atmos_coeff(&interpol_atmos_coef);
  return true;
}


int SrInterpAtmCoef(Lut_t *lut, Img_coord_int_t *input_loc, atmos_t *atmos_coef,atmos_t *interpol_atmos_coef) 
/* 
  Point order:

    0 ---- 1    +--> sample
    |      |    |
    |      |    v
    2 ---- 3   line

 */

{
  Img_coord_int_t p[4];
  int i, n,ipt, ib;
  double dl, ds, w;
  double sum[7][13], sum_w;
  Img_coord_int_t ar_region_half;

  ar_region_half.l = (lut->ar_region_size.l + 1) / 2;
  ar_region_half.s = (lut->ar_region_size.s + 1) / 2;

  p[0].l = (input_loc->l - ar_region_half.l) / lut->ar_region_size.l;

  p[2].l = p[0].l + 1;
  if (p[2].l >= lut->ar_size.l) {
    p[2].l = lut->ar_size.l - 1;
    if (p[0].l > 0) p[0].l--;
  }    
      
  p[1].l = p[0].l;
  p[3].l = p[2].l;

  p[0].s = (input_loc->s - ar_region_half.s) / lut->ar_region_size.s;
  p[1].s = p[0].s + 1;

  if (p[1].s >= lut->ar_size.s) {
    p[1].s = lut->ar_size.s - 1;
    if (p[0].s > 0) p[0].s--;
  }    

  p[2].s = p[0].s;
  p[3].s = p[1].s;

  n = 0;
  sum_w = 0.0;
  for (ipt=0;ipt<9;ipt++)
	 for (ib=0;ib<7;ib++)
		sum[ib][ipt]=0.;

  for (i = 0; i < 4; i++) {
    if (p[i].l != -1  &&  p[i].s != -1) {
		ipt=p[i].l * lut->ar_size.s + p[i].s;
		if (!(atmos_coef->computed[ipt])) continue; 

      dl = (input_loc->l - ar_region_half.l) - (p[i].l * lut->ar_region_size.l);
      dl = fabs(dl) / lut->ar_region_size.l;
      ds = (input_loc->s - ar_region_half.s) - (p[i].s * lut->ar_region_size.s);
      ds = fabs(ds) / lut->ar_region_size.s;
      w = (1.0 - dl) * (1.0 - ds);

      n++;
      sum_w += w;

		for (ib=0;ib<6;ib++) {
      	sum[ib][0] += (atmos_coef->tgOG[ib][ipt] * w);
      	sum[ib][1] += (atmos_coef->tgH2O[ib][ipt] * w);
      	sum[ib][2] += (atmos_coef->td_ra[ib][ipt] * w);
      	sum[ib][3] += (atmos_coef->tu_ra[ib][ipt] * w);
       	sum[ib][4] += (atmos_coef->rho_mol[ib][ipt] * w);
      	sum[ib][5] += (atmos_coef->rho_ra[ib][ipt] * w);
      	sum[ib][6] += (atmos_coef->td_da[ib][ipt] * w);
      	sum[ib][7] += (atmos_coef->tu_da[ib][ipt] * w);
     		sum[ib][8] += (atmos_coef->S_ra[ib][ipt] * w);
     		sum[ib][9] += (atmos_coef->td_r[ib][ipt] * w);
     		sum[ib][10] += (atmos_coef->tu_r[ib][ipt] * w);
     		sum[ib][11] += (atmos_coef->S_r[ib][ipt] * w);
     		sum[ib][12] += (atmos_coef->rho_r[ib][ipt] * w);
		}
    }
  }

  if (n > 0) {
	for (ib=0;ib<6;ib++) {
		interpol_atmos_coef->tgOG[ib][0]=sum[ib][0]/sum_w;
		interpol_atmos_coef->tgH2O[ib][0]=sum[ib][1]/sum_w;
		interpol_atmos_coef->td_ra[ib][0]=sum[ib][2]/sum_w;
		interpol_atmos_coef->tu_ra[ib][0]=sum[ib][3]/sum_w;
		interpol_atmos_coef->rho_mol[ib][0]=sum[ib][4]/sum_w;
		interpol_atmos_coef->rho_ra[ib][0]=sum[ib][5]/sum_w;
		interpol_atmos_coef->td_da[ib][0]=sum[ib][6]/sum_w;
		interpol_atmos_coef->tu_da[ib][0]=sum[ib][7]/sum_w;
		interpol_atmos_coef->S_ra[ib][0]=sum[ib][8]/sum_w;
		interpol_atmos_coef->td_r[ib][0]=sum[ib][9]/sum_w;
		interpol_atmos_coef->tu_r[ib][0]=sum[ib][10]/sum_w;
		interpol_atmos_coef->S_r[ib][0]=sum[ib][11]/sum_w;
		interpol_atmos_coef->rho_r[ib][0]=sum[ib][12]/sum_w;

	}
  }

  return 0;
}

