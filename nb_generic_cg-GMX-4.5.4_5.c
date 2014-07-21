/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "types/simple.h"
#include "vec.h"
#include "typedefs.h"
#include "nb_generic_cg.h"

 void
 gmx_nb_generic_cg_kernel(t_nblist *           nlist,
                          t_forcerec *         fr,
                          t_mdatoms *          mdatoms,
                          real *               x,
                          real *               f,
                          real *               fshift,
                          real *               Vc,
                          real *               Vvdw,
                          real                 tabscale,  
                          real *               VFtab,
                          int *                outeriter,
                          int *                inneriter)
{
     int           nri,ntype,table_nelements,icoul,ivdw;
     real          facel,gbtabscale;
     int           n,is3,i3,k,nj0,nj1,j3,ggid,nnn,n0;
     int           ai0,ai1,ai,aj0,aj1,aj;
     real          shX,shY,shZ;
     real          fscal,tx,ty,tz;
     real          rinvsq;
     real          iq;
     real          qq,vcoul,krsq,vctot;
     int           nti,nvdwparam;
     int           tj;
     real          rt,r,eps,eps2,Y,F,Geps,Heps2,VV,FF,Fp,fijD,fijR;
     real          rinvsix;
     real          Vvdwtot;
     real          Vvdw_rep,Vvdw_disp;
     real          ix,iy,iz,fix,fiy,fiz;
     real          jx,jy,jz;
     real          dx,dy,dz,rsq,rinv;
     real          c6,c12,cexp1,cexp2,br;
     real *        charge;
     real *        shiftvec;
     real *        vdwparam;
     int *         shift;
     int *         type;
     t_excl *      excl;

     int           ljangle_switch;
     int           ljangle_type1;
     int           ljangle_type2;
     real          ipx,ipy,ipz;
     real          jpx,jpy,jpz;
     int           ip3, jp3;
     real          d_iip_x,d_iip_y,d_iip_z;
     real          d_jjp_x,d_jjp_y,d_jjp_z;
     real          r_iip_sq,r_jjp_sq;
     real          r_iip_inv,r_jjp_inv;
     real          scalar_iipij,scalar_jjpji;
     real          cos_theta1,cos_theta1_sq,cos_theta1_pr;
     real          cos_theta2,cos_theta2_sq,cos_theta2_pr;
     real          radial_fac,radial_fac_prime;
     real          f_ljangle_1,f_ljangle_2,f_ljangle_3;
     real          ljangle_cap;
       
     icoul               = nlist->icoul;
     ivdw                = nlist->ivdw;
     
     /* avoid compiler warnings for cases that cannot happen */
     nnn                 = 0;
     vcoul               = 0.0;
     eps                 = 0.0;
     eps2                = 0.0;
     
     /* 3 VdW parameters for buckingham, otherwise 2 */
     nvdwparam           = (nlist->ivdw==2) ? 3 : 2;
     table_nelements     = (icoul==3) ? 4 : 0;
     table_nelements    += (ivdw==3) ? 8 : 0;
     
     charge              = mdatoms->chargeA;
     type                = mdatoms->typeA;
     facel               = fr->epsfac;
     shiftvec            = fr->shift_vec[0];
     vdwparam            = fr->nbfp;
     ntype               = fr->ntype;

     ljangle_switch      = fr->userint1;
     ljangle_type1       = fr->userint2;
     ljangle_type2       = fr->userint3;
     ljangle_cap         = 1000.;
     
     for(n=0; (n<nlist->nri); n++)
     {
         is3              = 3*nlist->shift[n];     
         shX              = shiftvec[is3];  
         shY              = shiftvec[is3+1];
         shZ              = shiftvec[is3+2];
         nj0              = nlist->jindex[n];      
         nj1              = nlist->jindex[n+1];    
         ai0              = nlist->iinr[n];
         ai1              = nlist->iinr_end[n];
         vctot            = 0;              
         Vvdwtot          = 0;              
         fix              = 0;
         fiy              = 0;
         fiz              = 0;
         
         for(k=nj0; (k<nj1); k++)
         {
             aj0              = nlist->jjnr[k];
             aj1              = nlist->jjnr_end[k];
             excl             = &nlist->excl[k*MAX_CGCGSIZE];

             for(ai=ai0; (ai<ai1); ai++)
             {
                 i3               = ai*3;
                 ix               = shX + x[i3+0];
                 iy               = shY + x[i3+1];
                 iz               = shZ + x[i3+2];
                 iq               = facel*charge[ai];
                 nti              = nvdwparam*ntype*type[ai];

                 /* Note that this code currently calculates
                  * all LJ and Coulomb interactions,
                  * even if the LJ parameters or charges are zero.
                  * If required, this can be optimized.
                  */

                 for(aj=aj0; (aj<aj1); aj++)
                 {
                     /* Check if this interaction is excluded */
                     if (excl[aj-aj0] & (1<<(ai-ai0)))
                     {
                         continue;
                     }

                     j3               = aj*3;
                     jx               = x[j3+0];
                     jy               = x[j3+1];
                     jz               = x[j3+2];
                     dx               = ix - jx;      
                     dy               = iy - jy;      
                     dz               = iz - jz;      
                     rsq              = dx*dx+dy*dy+dz*dz;
                     rinv             = gmx_invsqrt(rsq);
                     rinvsq           = rinv*rinv;  
                     fscal            = 0;

                     if (icoul==3 || ivdw==3)
                     {
                         r                = rsq*rinv;
                         rt               = r*tabscale;     
                         n0               = rt;             
                         eps              = rt-n0;          
                         eps2             = eps*eps;        
                         nnn              = table_nelements*n0;           				
                     }
                     
                     /* Coulomb interaction. icoul==0 means no interaction */
                     if (icoul > 0)
                     {
                         qq               = iq*charge[aj]; 
                         
                         switch(icoul)
                         {
                         case 1:
                             /* Vanilla cutoff coulomb */
                             vcoul            = qq*rinv;      
                             fscal            = vcoul*rinvsq; 
                             break;
                             
                         case 2:
                             /* Reaction-field */
                             krsq             = fr->k_rf*rsq;      
                             vcoul            = qq*(rinv+krsq-fr->c_rf);
                             fscal            = qq*(rinv-2.0*krsq)*rinvsq;
                             break;
                             
                         case 3:
                             /* Tabulated coulomb */
                             Y                = VFtab[nnn];     
                             F                = VFtab[nnn+1];   
                             Geps             = eps*VFtab[nnn+2];
                             Heps2            = eps2*VFtab[nnn+3];
                             nnn             += 4;
                             Fp               = F+Geps+Heps2;   
                             VV               = Y+eps*Fp;       
                             FF               = Fp+Geps+2.0*Heps2;
                             vcoul            = qq*VV;          
                             fscal            = -qq*FF*tabscale*rinv;
                             break;
                             
                         case 4:
                             /* GB */
                           gmx_fatal(FARGS,"Death & horror! GB generic interaction not implemented.\n");
                           break;
                           
                         default:
                             gmx_fatal(FARGS,"Death & horror! No generic coulomb interaction for icoul=%d.\n",icoul);
                             break;
                         }
                         vctot            = vctot+vcoul;    
                     } /* End of coulomb interactions */
                     
                     
                     /* VdW interaction. ivdw==0 means no interaction */
                     if (ivdw > 0)
                     {
                         tj               = nti+nvdwparam*type[aj];

                         /* Angular Lennard-Jones interaction */
                         /* Look for ljangle_type1 and ljangle_type2 */
                         if (ljangle_switch == 9999 && ((type[ai] == ljangle_type1 && type[aj] == ljangle_type2) ||
                                                        (type[ai] == ljangle_type2 && type[aj] == ljangle_type1))) 
                         {
                             /* Assume virtual sites have ID ai+1 and aj+1.
                              * Retrieve particle positions. */
                             ip3              = (ai+1)*3;
                             ipx              = x[ip3+0];
                             ipy              = x[ip3+1];
                             ipz              = x[ip3+2];
                             d_iip_x          = ix - ipx;
                             d_iip_y          = iy - ipy;
                             d_iip_z          = iz - ipz;
                             r_iip_sq         = d_iip_x*d_iip_x+d_iip_y*d_iip_y+d_iip_z*d_iip_z;
                             r_iip_inv        = gmx_invsqrt(r_iip_sq);

                             jp3              = (aj+1)*3;
                             jpx              = x[jp3+0];
                             jpy              = x[jp3+1];
                             jpz              = x[jp3+2];
                             d_jjp_x          = jx - jpx;
                             d_jjp_y          = jy - jpy;
                             d_jjp_z          = jz - jpz;
                             r_jjp_sq         = d_jjp_x*d_jjp_x+d_jjp_y*d_jjp_y+d_jjp_z*d_jjp_z;
                             r_jjp_inv        = gmx_invsqrt(r_jjp_sq);

                             /* Calculate two cos angles of theta1 (ip,i,j) and theta2(i,j,jp) */
                             scalar_iipij     = d_iip_x*dx+d_iip_y*dy+d_iip_z*dz;
                             cos_theta1       = scalar_iipij*r_iip_inv*rinv;

                             scalar_jjpji     = -d_jjp_x*dx-d_jjp_y*dy-d_jjp_z*dz;
                             cos_theta2       = scalar_jjpji*r_jjp_inv*rinv;
                                 
                             /* Condition: angles smaller than 90deg. */
                             if (cos_theta1 > 0. && cos_theta2 > 0.)
                             {
                                 cos_theta1_sq    = sqr(cos_theta1);
                                 cos_theta1_pr    = -2.0*cos_theta1;
                                 cos_theta2_sq    = sqr(cos_theta2);
                                 cos_theta2_pr    = -2.0*cos_theta2;
                                     
                                 /* Use vanilla Lennard-Jones coefficients for ljangle interaction. */
                                 c6               = vdwparam[tj];
                                 c12              = vdwparam[tj+1];

                                 /* Ljangle is a 12-10 interaction. */
                                 rinvsix          = rinvsq*rinvsq*rinvsq;
                                 Vvdw_disp        = c6*rinvsix*rinvsq*rinvsq;
                                 Vvdw_rep         = c12*rinvsix*rinvsix;
                                 radial_fac       = Vvdw_rep-Vvdw_disp;
                                 /* Cap radial_fac to ljangle_cap. */
                                 if (radial_fac > ljangle_cap)
                                     radial_fac = ljangle_cap;
                                 Vvdwtot          = Vvdwtot+radial_fac*cos_theta1_sq*cos_theta2_sq;
                                 /* printf("ljangle: %3d-%3d: %6.2f - rad %6.5f - rinv %6.5f %6.5f - coeffs %6.5f %6.5f - rsq %6.5f - c %6.5f %6.5f\n",
                                    ai,aj,Vvdwtot,radial_fac,rinvsq,rinvsix,Vvdw_rep,Vvdw_disp,rsq,c6,c12); */
                                 radial_fac_prime = (10.0*Vvdw_disp-12.0*Vvdw_rep)*rinv;
                                 /* Update forces for all 4 particles.
                                  * f_ljangle_X correspond to different
                                  * terms in the force expression. */
                                 f_ljangle_1      = cos_theta1_sq*cos_theta2_sq*radial_fac_prime*rinv;
                                 f_ljangle_2      = cos_theta2_sq*cos_theta1_pr*radial_fac;
                                 f_ljangle_3      = cos_theta1_sq*cos_theta2_pr*radial_fac;
                                 /* particle i */
                                 f[i3+0]         += 
                                     f_ljangle_1*(-dx) +
                                     f_ljangle_2*( (dx+d_iip_x)*rinv*r_iip_inv -
                                                   cos_theta1*(d_iip_x*r_iip_inv*r_iip_inv + dx*rinvsq) ) -
                                     f_ljangle_3*(d_jjp_x*rinv*r_jjp_inv + cos_theta2*dx*rinvsq);
                                 f[i3+1]         += 
                                     f_ljangle_1*(-dy) +
                                     f_ljangle_2*( (dy+d_iip_y)*rinv*r_iip_inv -
                                                   cos_theta1*(d_iip_y*r_iip_inv*r_iip_inv + dy*rinvsq) ) -
                                     f_ljangle_3*(d_jjp_y*rinv*r_jjp_inv + cos_theta2*dy*rinvsq);
                                 f[i3+2]         += 
                                     f_ljangle_1*(-dz) +
                                     f_ljangle_2*( (dz+d_iip_z)*rinv*r_iip_inv -
                                                   cos_theta1*(d_iip_z*r_iip_inv*r_iip_inv + dz*rinvsq) ) -
                                     f_ljangle_3*(d_jjp_z*rinv*r_jjp_inv + cos_theta2*dz*rinvsq);
                                 /* particle j */
                                 f[j3+0]         += 
                                     f_ljangle_1*dx +
                                     f_ljangle_2*( -d_iip_x*rinv*r_iip_inv + cos_theta1*dx*rinvsq ) +
                                     f_ljangle_3*( (-dx+d_jjp_x)*rinv*r_jjp_inv + 
                                                   cos_theta2*(-d_jjp_x*r_jjp_inv*r_jjp_inv + dx*rinvsq) );
                                 f[j3+1]         += 
                                     f_ljangle_1*dy +
                                     f_ljangle_2*( -d_iip_y*rinv*r_iip_inv + cos_theta1*dy*rinvsq ) +
                                     f_ljangle_3*( (-dy+d_jjp_y)*rinv*r_jjp_inv + 
                                                   cos_theta2*(-d_jjp_y*r_jjp_inv*r_jjp_inv + dy*rinvsq) );
                                 f[j3+2]         += 
                                     f_ljangle_1*dz +
                                     f_ljangle_2*( -d_iip_z*rinv*r_iip_inv + cos_theta1*dz*rinvsq ) +
                                     f_ljangle_3*( (-dz+d_jjp_z)*rinv*r_jjp_inv + 
                                                   cos_theta2*(-d_jjp_z*r_jjp_inv*r_jjp_inv + dz*rinvsq) );
                                 /* particle ip */
                                 f[ip3+0]        +=
                                     f_ljangle_2*(-dx*rinv*r_iip_inv + cos_theta1*d_iip_x*r_iip_inv*r_iip_inv);
                                 f[ip3+1]        +=
                                     f_ljangle_2*(-dy*rinv*r_iip_inv + cos_theta1*d_iip_y*r_iip_inv*r_iip_inv);
                                 f[ip3+2]        +=
                                     f_ljangle_2*(-dz*rinv*r_iip_inv + cos_theta1*d_iip_z*r_iip_inv*r_iip_inv);
                                 /* particle jp */
                                 f[jp3+0]        +=
                                     f_ljangle_3*(dx*rinv*r_jjp_inv + cos_theta2*d_jjp_x*r_jjp_inv*r_jjp_inv);
                                 f[jp3+1]        +=
                                     f_ljangle_3*(dy*rinv*r_jjp_inv + cos_theta2*d_jjp_y*r_jjp_inv*r_jjp_inv);
                                 f[jp3+2]        +=
                                     f_ljangle_3*(dz*rinv*r_jjp_inv + cos_theta2*d_jjp_z*r_jjp_inv*r_jjp_inv);
                                     
                             }
                         } else { 
                         
                         switch(ivdw)
                         {
                         case 1:
                             /* Vanilla Lennard-Jones cutoff */
                             c6               = vdwparam[tj];   
                             c12              = vdwparam[tj+1]; 
                             
                             rinvsix          = rinvsq*rinvsq*rinvsq;
                             Vvdw_disp        = c6*rinvsix;     
                             Vvdw_rep         = c12*rinvsix*rinvsix;
                             fscal           += (12.0*Vvdw_rep-6.0*Vvdw_disp)*rinvsq;
                             Vvdwtot          = Vvdwtot+Vvdw_rep-Vvdw_disp;
                             break;
                             
                         case 2:
                             /* Buckingham */
                             c6               = vdwparam[tj];   
                             cexp1            = vdwparam[tj+1]; 
                             cexp2            = vdwparam[tj+2]; 
                             
                             rinvsix          = rinvsq*rinvsq*rinvsq;
                             Vvdw_disp        = c6*rinvsix;     
                             br               = cexp2*rsq*rinv;
                             Vvdw_rep         = cexp1*exp(-br); 
                             fscal           += (br*Vvdw_rep-6.0*Vvdw_disp)*rinvsq;
                             Vvdwtot          = Vvdwtot+Vvdw_rep-Vvdw_disp;
                             break;
                             
                         case 3:
                             /* Tabulated VdW */
                             c6               = vdwparam[tj];   
                             c12              = vdwparam[tj+1]; 
                             
                             Y                = VFtab[nnn];     
                             F                = VFtab[nnn+1];   
                             Geps             = eps*VFtab[nnn+2];
                             Heps2            = eps2*VFtab[nnn+3];
                             Fp               = F+Geps+Heps2;   
                             VV               = Y+eps*Fp;       
                             FF               = Fp+Geps+2.0*Heps2;
                             Vvdw_disp        = c6*VV;          
                             fijD             = c6*FF;          
                             nnn             += 4;          
                             Y                = VFtab[nnn];     
                             F                = VFtab[nnn+1];   
                             Geps             = eps*VFtab[nnn+2];
                             Heps2            = eps2*VFtab[nnn+3];
                             Fp               = F+Geps+Heps2;   
                             VV               = Y+eps*Fp;       
                             FF               = Fp+Geps+2.0*Heps2;
                             Vvdw_rep         = c12*VV;         
                             fijR             = c12*FF;         
                             fscal           += -(fijD+fijR)*tabscale*rinv;
                             Vvdwtot          = Vvdwtot + Vvdw_disp + Vvdw_rep;						
                             break;
                             
                         default:
                             gmx_fatal(FARGS,"Death & horror! No generic VdW interaction for ivdw=%d.\n",ivdw);
                             break;
                         }
                         }   
                     } /* end VdW interactions */
                     
                     
                     tx               = fscal*dx;     
                     ty               = fscal*dy;     
                     tz               = fscal*dz;     
                     f[i3+0]         += tx;
                     f[i3+1]         += ty;
                     f[i3+2]         += tz;
                     f[j3+0]         -= tx;
                     f[j3+1]         -= ty;
                     f[j3+2]         -= tz;
                     fix             += tx;
                     fiy             += ty;
                     fiz             += tz;
                 }
             }
         }

         fshift[is3]     += fix;
         fshift[is3+1]   += fiy;
         fshift[is3+2]   += fiz;
         ggid             = nlist->gid[n];         
         Vc[ggid]        += vctot;
         Vvdw[ggid]      += Vvdwtot;
     }
     
     *outeriter       = nlist->nri;            
     *inneriter       = nlist->jindex[n];          	
}

