/*--------------------------------------------------------------------
 * $Id: sofi.c 2623 2011-12-23 10:52:38Z robert.buras $
 * 
 * This file is part of libRadtran.
 * Copyright (c) 1997-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License   
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.        
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
 * GNU General Public License for more details.                    
 * 
 * You should have received a copy of the GNU General Public License          
 * along with this program; if not, write to the Free Software                
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#include "ascii.h"
#include "numeric.h"
#include "sofi.h"

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif


#define PROGRAM "sofi"
#define VERS "0.1"

/*************************************************************/
/* print usage information                                   */
/*************************************************************/


static void print_usage (char *filename)
{
  fprintf (stderr, "%s %s - Calculate irradiance or diffuse flux during a total eclipse. \n", PROGRAM, VERS );
  fprintf(stderr, "Per default it is assumed that the centre of the eclipse is at the \n");
  fprintf(stderr, "observer position. \n\n");
  fprintf (stderr, "written by Claudia Emde,\n");
  fprintf (stderr, "           DLR, e-Mail claudia.emde@dlr.de\n");
  fprintf (stderr, "Version %s finished May 30, 2006\n\n", VERS);
  fprintf (stderr, "Optional arguments:\n");
  fprintf (stderr, "-h            display help message\n");
  fprintf (stderr, "-z <degrees>  solar zenith angle\n");
  fprintf (stderr, "-a <degrees>  solar azimuth angle \n");
  fprintf (stderr, "-l <nm>       wavelength \n");
  fprintf (stderr, "-x <km>       horizontal shift of eclipse shadow \n");
  fprintf (stderr, "-y <km>       vertical shift of eclipse shadow \n");
  fprintf (stderr, "-b <km>       height of atmosphere \n");
  fprintf (stderr, "-c            calculate clearsky \n");
  fprintf (stderr, "-s            calculate standard deviation \n");
  fprintf (stderr, "USAGE: %s [options] <filename>\n\n", filename);
  fprintf (stderr, "%s calculates the irradiance or the diffus flux \n", PROGRAM); 
  fprintf (stderr, "during the solar eclipse at 29-03-2006.\n");
  fprintf (stderr, "The input file is the photon distribution\n");
  fprintf (stderr, "at the top of the atmosphere calculated using  \n");
  fprintf (stderr, "mystic (backward), e.g. mc35.rad. \n");
  fprintf (stderr, "The tool can also be used to calculte clearsky \n");
  fprintf (stderr, "irradiance or flux provided the mystic file above. \n");
  fprintf (stderr, "\n");
}

static int get_options (int argc, char **argv,
                        char *programname,
                        char *infilename,
                        double *sza,
                        double *saa,
                        double *lambda, 
                        double *z,
                        double *dx,
                        double *dy, 
                        int *clear, 
                        int *std)
{
  int c=0;
  char *dummy=NULL;
  
  /* save name of program */
  strncpy (programname, argv[0], FILENAME_MAX);
  
  /* set defaults */
  strcpy (infilename, "");
  *sza = 0.0;
  *saa = 0.0;
  *lambda = 400.0;
  *z=50.0;
  *dx=0.0;
  *dy=0.0;
  /* By default it is assumed that sofi should be calculated*/
  *clear=0;
  *std=0;
  
  while ((c=getopt (argc, argv, "hz:a:l:b:x:y:cs")) != EOF)  {
    switch(c)  {
    case 'h':  /* help */
      print_usage (programname);
      return (-1);
      break;
    case 'z':  /* number of legendre polynomials */
      *sza = strtod(optarg, &dummy);
      break;
    case 'a':
      *saa = strtod(optarg, &dummy);
      break;
    case 'l':
      *lambda = strtod(optarg, &dummy);
      break;
    case 'b':
      *z = strtod(optarg, &dummy);
      break;
    case 'x':
      *dx = strtod(optarg, &dummy);
      break; 
    case 'y':
      *dy = strtod(optarg, &dummy);
      break;
    case 'c':
      *clear = 1;
      break;
    case 's':
      *std = 1;
      break;
    default:
      print_usage (programname);
      return (-1);
    }
  }


  /* check number of remaining command line arguments */
  if (argc - optind != 1)  {
    print_usage (programname);
    return -1;
  }
  
  strncpy (infilename, argv[optind+0], FILENAME_MAX);
  
  return 0;  /* if o.k. */
}


int main(int argc, char **argv)
{
  int i=0, status=0;

  char programname[FILENAME_MAX], infilename[FILENAME_MAX]="";
  
  double sza=0.0;
  double saa=0.0;
  double lambda=0.0;
  double z_toa=50.0;
  double dx_moon=0.0;
  double dy_moon=0.0;
  int clear=0;
  int std=0; 
  double r_earth=6371211.0;
  double irr_tot=0.0;
  double r=0.0, t=0.0, alpha=0.0, x=0.0, x2=0.0, y=0.0, y2=0.0; 
  int i_re=0, dx=0.0, dy=0.0;
  double wvnmlo=0.0, wvnmhi=0.0;
  double *pd, *irr; 
  double cossaa=0.0, sinsaa=0.0, cossza=0.0, tansza=0.0, sinsza=0.0;
  int Nx=10000;
  
  status = get_options (argc, argv, programname, infilename, &sza,
                        &saa, &lambda, &z_toa, &dx_moon, &dy_moon,  &clear, &std);

  if (status!=0)
    return status;
  
  /* Read mystic output */

  double *fi=NULL;
  double *se=NULL;
  double *th=NULL;
  double *fo=NULL;
  double *fv=NULL;
  double *si=NULL;
  double *sv=NULL;
  double *ei=NULL;
  int L=0;
  
  read_8c_file (infilename, &fi, &se, &th, &fo, &fv, &si,&sv, &ei, &L); 
  fprintf (stderr, " ... read %d data points from %s\n", L, infilename);
  
  if (status!=0) {
    fprintf (stderr, "error %d reading %s\n", status, infilename);
    return status;
  }

  if (sza >= 90.0) {
    fprintf (stderr, "Solar eclipse calculations not implemented for \n");
    fprintf (stderr, "solar zenith angles greater or equal 90 degrees. \n");
    return -1;
  }
  
  if (!clear) {
    /* Calculate probability density function. */
   
    wvnmhi=1.0/(lambda*100.0)*1e9;
    wvnmlo=1.0/(lambda*100.0)*1e9;
  
    
    pd = (double*) calloc(Nx,sizeof(double));
    
    sample_photons_sofi(wvnmlo,wvnmhi,pd);

    /*for (i=0; i<5000; i++)
      fprintf(stdout, "%d   %g \n", i, pd[i]);*/ 
    
    /*height of the model atmosphere, convert to meter*/
    z_toa*=1000;
    dx_moon*=1000;
    dy_moon*=1000; 
    
    /*often needed expressions*/ 
    cossaa=cos(saa*PI/180.0);
    sinsaa=sin(saa*PI/180.0);
    tansza=tan(sza*PI/180.0);
    cossza=cos(sza*PI/180.0);
    sinsza=sin(sza*PI/180.0);
    
    /* calculate centre of umbral shadow at TOA */
    r=z_toa*tansza; 
    dx=r*sinsaa; 
    dy=r*cossaa;
    
    /* Sampled irradiance */
    irr = (double*) calloc(L,sizeof(double));
    
    for(i=0; i<L; i++)
      { 
        /* Coordinates on spherical sampling area */
        x=fi[i]-fi[0]-(fi[L-1]+fi[0])/2;
        y=se[i]-se[0]-(se[L-1]+se[0])/2;

        /* If the shift is calculated on the Earth surface, i.e. 
           on the sphere, the shift given as x y inputs must be 
           applied here.
           x+=dx_moon;
           y+=dy_moon;
        */
        
        /* Calculate shift from spherical sampling domain to plane*/
        r=z_toa+r_earth;
        alpha=asin(sqrt(x*x+y*y)/r);
        t=r*(1.0-cos(alpha))/cossza;
        
        /* coordinates on plane */
        x+=t*sinsza*sinsaa;
        y+=t*sinsza*cossaa;
        
        /* Shift shadow to intersection point sun - TOA, 
           needs to be done in plane, since dx, dy are calculated in plane */ 
        x+=(dx+dx_moon);
        y+=(dy+dy_moon);
   
        /* Rotation and dilation in y direction */
        x2=cossaa*x-sinsaa*y;
        y2=(sinsaa*x+cossaa*y)*cossza;
        
        i_re= (int)fabs((sqrt(x2*x2+y2*y2))/1000.0);
        
        /* Out of range of weighting function, set i_re to last value (1) */
        if(i_re > Nx-1)
          i_re = Nx-1;
        
        if(std)
          irr[i]=ei[i]*ei[i]*pd[i_re]*pd[i_re];
        else
          irr[i]=ei[i]*pd[i_re];
        
        irr_tot+=irr[i];
        
        /* Use this if you want to plot the sofi weighted photon distribution
           or the shadow */
        /* fprintf(stdout,  "%f  %f %g \n " ,fi[i], se[i], pd[i_re]);
           fprintf(stdout,  "%f  %f %g \n " ,fi[i], se[i], irr[i]);*/
      }
  }
  else
    for(i=0; i<L; i++)
      irr_tot+=ei[i];
  
  if(std)
    irr_tot=sqrt(irr_tot);
  
  fprintf(stdout,  "Total irradiance/flux: %g \n", irr_tot);
  return 0;
}  


