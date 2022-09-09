/*--------------------------------------------------------------------
 * $Id: tipa.h 2623 2011-12-23 10:52:38Z robert.buras $
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

/*--------------------------------------------------------------------
 * tipa.c contains all functions concerning the tilted independent
 * column (pixel) approximation TICA (TIPA)
 *
 * Author: Ulrike Wissmeier
 *--------------------------------------------------------------------*/

#ifndef __tipa_h
#define __tipa_h

#include "../src/cloud.h"

int tipa_dirdiff ( caoth3d_out_struct *caoth3d,
		   tipa_struct        *tipa,
		   int                 tipaswitch,
		   double              sza,
		   double              phi0,
		   float              *zout,
		   int                 mctipa );

int tipa_xyzintersec ( caoth3d_out_struct *caoth3d,
		       tipa_struct        *tipa,
		       int                 mctipa );

int tipa_dirtilt ( caoth3d_out_struct *caoth3d,
		   atm_out_struct      atm,
		   alt_out_struct      alt,
		   tipa_struct        *tipa,
		   int                 tipaswitch,
		   double              sza,
		   double              phi0,
		   int                 lower_wl_id,
		   int                 upper_wl_id,
		   int                 mctipa );

int tipa_calcdtau ( caoth_inp_struct    input_caoth,
		    caoth3d_out_struct  output_caoth3d,
		    int                 iipa,
		    int                 jipa,
		    int                 iv,
		    input_struct        input,
		    wl_out_struct       wl_out,
		    /* Output */
		    caoth_out_struct   *output_caoth,
		    tipa_struct        *tipa );

int mc_tipa_dir ( caoth_out_struct   *caoth,
		  caoth3d_out_struct *caoth3d,
		  tipa_struct        *tipa,
		  int                 iv,
		  input_struct        input,
		  wl_out_struct      *wl );

#endif



