/************************************************************************
 * $Id: ambralsfor.c 2623 2011-12-23 10:52:38Z robert.buras $
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include "ambralsfor.h"
#include "f77-uscore.h"

#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

#ifndef M_PI_2
#define M_PI_2		1.57079632679489661923	/* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4		0.78539816339744830962	/* pi/4 */
#endif


/* All the following functions are taken from ambmodels.c.reciprocal.
 *The functions have been changed to ANSI C from K&R and trig functions
 * have been replaced by their float equivs
 */

void LiKernel(float, float, float, float, float, float, float *);
void GetPhaang(float, float, float, float, float, float *, float *,float *);
void GetpAngles(float, float, float *, float *, float *);
void GetDistance(float , float , float ,float *);
void GetOverlap(float, float, float, float, float, float, float, float *,
                float *);

/******************************************************************************************/ 
/* AMBRALS - Algorithm for MODIS Bidirectional Reflectance Anisotropy of the Land Surface */
/******************************************************************************************/ 

double ambrals_brdf (double iso, double vol, double geo, 
		     double mu1, double mu2, double phi)
{
  /* this function runs the RossThickLiSparseReciprocal model in the 
   * forward mode and returns the calculated reflectance.
   *
   * Parameters:
   * mu1 mu2 phi                   view geometry
   * iso vol geo                   BRDF kernel weights parameters 
   */

  float cosphi, sinphi, costv, sintv, costi, sinti;
  float cosphaang, phaang, sinphaang, rosskernel, tantv, tanti;
  float likernel, refl;

  /* need to change phi convention, 180 degree = backward */
  phi = 180.0 - phi;
  while (phi<0)
    phi+=360.0;

  while (phi>=360)
    phi-=360.0;
    
  phi = D2R(phi);

  cosphi = cos(phi); /* cosphi = cosf(phi) */
  sinphi = sin(phi); /* sinphi = sinf(phi) */
  costv  = mu1;
  costi  = mu2;
  sintv  = sqrt(1.0-costv*costv);
  sinti  = sqrt(1.0-costi*costi);
        
  GetPhaang(costv, costi, sintv, sinti, cosphi, &cosphaang,
	    &phaang, &sinphaang);
        
  rosskernel = ((M_PI_2 - phaang) * cosphaang + sinphaang)/(costi + costv)
                - M_PI_4;
  
  tantv = sintv/costv;
  tanti = sinti/costi;
  
  LiKernel(HB,BR,tantv,tanti,sinphi,cosphi, &likernel);
  
  refl = iso + vol * rosskernel + geo * likernel;
  
  /* prevent negative reflectivities */
  if (refl < 0)
    refl = 0;

  return (double) refl;
}

void LiKernel(float hbratio, float brratio, float tantv, float tanti,
              float sinphi, float cosphi, float *result)
{

  /* This func calculates the LiSparseModis kernel */

  float sintvp, costvp, tantvp, sintip, costip, tantip;
  float phaangp, cosphaangp, sinphaangp, distancep, overlap, temp;

  GetpAngles(brratio,tantv,&sintvp,&costvp,&tantvp);
  GetpAngles(brratio,tanti,&sintip,&costip,&tantip);
  GetPhaang(costvp,costip,sintvp,sintip,cosphi,&cosphaangp,&phaangp,
             &sinphaangp);
  GetDistance(tantvp,tantip,cosphi,&distancep);
  GetOverlap(hbratio,distancep,costvp,costip,tantvp,tantip,sinphi,
              &overlap,&temp);

  *result = overlap - temp + 0.5 * (1.+cosphaangp)/costvp/costip;
}

void GetPhaang(float cos1, float cos2, float sin1, float sin2, float cos3, 
               float *cosres, float *res,float *sinres)
{
  /* This func calculates the phase angle */

  *cosres = cos1*cos2 + sin1*sin2*cos3;
  if (*cosres > 1.0)
    *cosres=1.0;
  if (*cosres < -1.0)
    *cosres=-1.0;
  *res = acos( *cosres ); /* *res = acosf( MAX(-1., MIN(1.,*cosres)) ) */
  *sinres = sin(*res); /* *sinres = sinf(*res) */
}

void GetpAngles(float brratio, float tan1, float *sinp, float *cosp,
                float *tanp)
{
  /* This func calculates the 'prime' angles */

  float angp;

    *tanp = brratio*tan1;
  if(*tanp < 0) *tanp = 0.;
  angp = atan(*tanp); /* angp = atanf(*tanp) */
  *sinp = sin(angp);  /* *sinp = sinf(angp)  */
  *cosp = cos(angp);  /* *cosp = cosf(angp)  */
}

void GetDistance(float tan1, float tan2, float cos3,float *res)
{
  /* This func calculates the D distance */

  float temp;

  temp = tan1*tan1+tan2*tan2-2.*tan1*tan2*cos3;
  if (temp < 0.0)
    temp=0.0;
  *res = sqrt(temp);  /* *res = sqrtf(MAX(0.,temp)) */
}

void GetOverlap(float hbratio, float distance, float cos1, float cos2,
                float tan1, float tan2, float sin3, float *overlap,
                float *temp)
{
  /* This func calculates the Overlap distance */

  float cost, sint, tvar;

  *temp = 1./cos1 + 1./cos2;
  
   cost = hbratio*sqrt(distance*distance+tan1*tan1*tan2*tan2*sin3*sin3)/
          *temp;
  /*
   cost = hbratio*sqrtf(distance*distance+tan1*tan1*tan2*tan2*sin3*sin3)/
          *temp;
   */
   if (cost > 1.0)
     cost=1.0;
   if (cost < -1.0)
     cost=-1.0;
   tvar = acos(cost); /* tvar = acosf(cost) */
   sint = sin(tvar); /* sint = sinf(tvar) */
   *overlap = 1./M_PI *(tvar-sint*cost)*(*temp);
   if (*overlap < 0.0)
     *overlap = 0.0;
}






void F77_FUNC (ambralsfort, AMBRALSFORT) (float *iso, float *vol, float *geo,
			    float *mu1, float *mu2, float *phi, float *bdref)
{
  *bdref = (float) ambrals_brdf ((double) *iso, (double) *vol, (double) *geo, (double) *mu1, (double) *mu2, (double) *phi);
}


