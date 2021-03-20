/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Eduardo Bringa (LLNL)
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   The buck/coul/cut  potential was used as a template and modified into 
   the SHIK potential form with a wolf term for long range interactions. 
   An extra repulsive term has been added to the short range part of the 
   interaction. The long range interaction has been modified to use the 
   wolf method. Both the short range part and long range part have then 
   been shifted to offset the truncation and a smooothening function has 
   been applied to both of them. Length scale of the smoothening function 
   can be changed using the gamma parameters. This potential file has been
   tested upto lammps ver Mar 2018.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Modifying author: Siddharth Sundararaman (RPI) (SHIK/wolf)
 ------------------------------------------------------------------------- */



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "pair_SHIK_wolf.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairSHIKWolf::PairSHIKWolf(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairSHIKWolf::~PairSHIKWolf()
{
  if (!copymode) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(cut_coul);
    memory->destroy(cut_coulsq);
    memory->destroy(a);
    memory->destroy(b);
    memory->destroy(c);
    memory->destroy(d);  //Changed
    memory->destroy(gammalj); //Changed
    memory->destroy(gammacoul); //Changed
    memory->destroy(buck1);
    memory->destroy(buck2);
    memory->destroy(buck3); //Changed

  }
}

/* ---------------------------------------------------------------------- */

void PairSHIKWolf::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,r2inv,r6inv,r24inv,r2invcutlj,r6invcutlj,r24invcutlj,forcecoul,forcebuck,factor_coul,factor_lj; //change
  double r,rexp,rexpcutlj,rcutcoul,rcutlj,ljsmooth,coulsmooth,v_ljcut,v_lj; //change
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r = sqrt(rsq);

        if (rsq < cut_coulsq[itype][jtype]){
          rcutcoul = sqrt(cut_coulsq[itype][jtype]); //change
          coulsmooth = exp(-(gammacoul[itype][jtype]*gammacoul[itype][jtype])/((r-rcutcoul)*(r-rcutcoul))); //Smoothening Function
          forcecoul = qqrd2e * qtmp*q[j]*coulsmooth*((1.0/r-r/(rcutcoul*rcutcoul))-(1.0/r-1.0/rcutcoul+(r-rcutcoul)/(rcutcoul*rcutcoul))*r*2.0*(gammacoul[itype][jtype]*gammacoul[itype][jtype])/((r-rcutcoul)*(r-rcutcoul)*(r-rcutcoul))); //Major Change
        }else forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
        rcutlj = sqrt(cut_ljsq[itype][jtype]); //change
          r2invcutlj = 1.0/cut_ljsq[itype][jtype]; //change
          r6invcutlj = r2invcutlj*r2invcutlj*r2invcutlj; //change
          r24invcutlj = r6invcutlj*r6invcutlj*r6invcutlj*r6invcutlj; //change
          r6inv = r2inv*r2inv*r2inv;
          r24inv = r6inv*r6inv*r6inv*r6inv; //change
          rexp = exp(-r*b[itype][jtype]);
          rexpcutlj = exp(-rcutlj*b[itype][jtype]); //change
          ljsmooth = exp(-(gammalj[itype][jtype]*gammalj[itype][jtype])/((r-rcutlj)*(r-rcutlj))); //change
          v_lj = a[itype][jtype]*rexp - c[itype][jtype]*r6inv + d[itype][jtype]*r24inv; //change
          v_ljcut = a[itype][jtype]*rexpcutlj - c[itype][jtype]*r6invcutlj + d[itype][jtype]*r24invcutlj; //change
          forcebuck = ((buck1[itype][jtype]*r*rexp - buck2[itype][jtype]*r6inv + buck3[itype][jtype]*r24inv) - (v_lj-v_ljcut)*r*2.0*(gammalj[itype][jtype]*gammalj[itype][jtype])/((r-rcutlj)*(r-rcutlj)*(r-rcutlj)))*ljsmooth ; //Major change
	} else forcebuck = 0.0;

        fpair = (factor_coul*forcecoul + factor_lj*forcebuck) * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          if (rsq < cut_coulsq[itype][jtype])
            ecoul = factor_coul * qqrd2e * qtmp*q[j]*(1.0/r-1.0/rcutcoul+(r-rcutcoul)/(rcutcoul*rcutcoul))*coulsmooth; //change
          else ecoul = 0.0;
          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = (v_lj-v_ljcut)*ljsmooth; //change
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSHIKWolf::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");

  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_lj,n+1,n+1,"pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  memory->create(cut_coul,n+1,n+1,"pair:cut_coul");
  memory->create(cut_coulsq,n+1,n+1,"pair:cut_coulsq");
  memory->create(a,n+1,n+1,"pair:a");
  memory->create(b,n+1,n+1,"pair:b");
  memory->create(c,n+1,n+1,"pair:c");
  memory->create(d,n+1,n+1,"pair:d"); //change
  memory->create(gammalj,n+1,n+1,"pair:gammalj"); //change
  memory->create(gammacoul,n+1,n+1,"pair:gammacoul"); //change
  memory->create(buck1,n+1,n+1,"pair:buck1");
  memory->create(buck2,n+1,n+1,"pair:buck2");
  memory->create(buck3,n+1,n+1,"pair:buck3"); //change
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSHIKWolf::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2) error->all(FLERR,"Illegal pair_style command");

  cut_lj_global = force->numeric(FLERR,arg[0]);
  if (narg == 1) cut_coul_global = cut_lj_global;
  else cut_coul_global = force->numeric(FLERR,arg[1]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_lj[i][j] = cut_lj_global;
          cut_coul[i][j] = cut_coul_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSHIKWolf::coeff(int narg, char **arg)
{
  if (narg < 8 || narg > 10)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double a_one = force->numeric(FLERR,arg[2]);
  double b_one = force->numeric(FLERR,arg[3]);
  if (b_one <= 0) error->all(FLERR,"Incorrect args for pair coefficients");
  double c_one = force->numeric(FLERR,arg[4]);
  double d_one = force->numeric(FLERR,arg[5]); //change
  double gammalj_one = force->numeric(FLERR,arg[6]); //change
  double gammacoul_one = force->numeric(FLERR,arg[7]); //chnage
  double cut_lj_one = cut_lj_global;
  double cut_coul_one = cut_coul_global;
  if (narg >= 9) cut_coul_one = cut_lj_one = force->numeric(FLERR,arg[8]);
  if (narg == 10) cut_coul_one = force->numeric(FLERR,arg[9]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a[i][j] = a_one;
      b[i][j] = b_one;
      c[i][j] = c_one;
      d[i][j] = d_one; //change
      gammalj[i][j] = gammalj_one; //change
      gammacoul[i][j] = gammacoul_one; //change
      cut_lj[i][j] = cut_lj_one;
      cut_coul[i][j] = cut_coul_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSHIKWolf::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style shik/wolf requires atom attribute q");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSHIKWolf::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  double cut = MAX(cut_lj[i][j],cut_coul[i][j]);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];
  cut_coulsq[i][j] = cut_coul[i][j] * cut_coul[i][j];

  buck1[i][j] = a[i][j]*b[i][j]; //change
  buck2[i][j] = 6.0*c[i][j];
  buck3[i][j] = 24.0*d[i][j]; //change

  cut_ljsq[j][i] = cut_ljsq[i][j];
  cut_coulsq[j][i] = cut_coulsq[i][j];
  a[j][i] = a[i][j];
  b[j][i] = b[i][j];
  c[j][i] = c[i][j];   //change
  d[j][i] = d[i][j];   //change
  gammalj[j][i] = gammalj[i][j];
  gammacoul[j][i] = gammacoul[i][j];
  buck1[j][i] = buck1[i][j];
  buck2[j][i] = buck2[i][j];
  buck3[j][i] = buck3[i][j];   //change
  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce
 // Tail correction not implemented for SHIK potential implementation yet
  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double b1 = b[i][j];
    double b2 = b1*b1;
    double b3 = b2*b1;
    double rc = cut_lj[i][j];
    double rc2 = rc*rc;
    double rc3 = rc2*rc;
    etail_ij = 2.0*MY_PI*all[0]*all[1]*
      (a[i][j]*exp(-rc/b1)*b1*(rc2 + 2.0*b1*rc + 2.0*b2) -
       c[i][j]/(3.0*rc3));
    ptail_ij = (-1/3.0)*2.0*MY_PI*all[0]*all[1]*
      (-a[i][j]*exp(-rc/b1)*
       (rc3 + 3.0*b1*rc2 + 6.0*b2*rc + 6.0*b3) + 2.0*c[i][j]/rc3);
  }

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSHIKWolf::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&a[i][j],sizeof(double),1,fp);
        fwrite(&b[i][j],sizeof(double),1,fp);
        fwrite(&c[i][j],sizeof(double),1,fp);
	fwrite(&d[i][j],sizeof(double),1,fp);
        fwrite(&gammalj[i][j],sizeof(double),1,fp);
        fwrite(&gammacoul[i][j],sizeof(double),1,fp);
        fwrite(&cut_lj[i][j],sizeof(double),1,fp);
        fwrite(&cut_coul[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSHIKWolf::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&a[i][j],sizeof(double),1,fp);
          fread(&b[i][j],sizeof(double),1,fp);
          fread(&c[i][j],sizeof(double),1,fp);
          fread(&d[i][j],sizeof(double),1,fp);
          fread(&gammalj[i][j],sizeof(double),1,fp);
          fread(&gammacoul[i][j],sizeof(double),1,fp);
          fread(&cut_lj[i][j],sizeof(double),1,fp);
          fread(&cut_coul[i][j],sizeof(double),1,fp);
	}
        MPI_Bcast(&a[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&d[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gammalj[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gammacoul[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_coul[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSHIKWolf::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSHIKWolf::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul_global,sizeof(double),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairSHIKWolf::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
	fprintf(fp,"%d %g %g %g %g %g %g\n",i,a[i][i],b[i][i],c[i][i],d[i][i],gammalj[i][i],gammacoul[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairSHIKWolf::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g %g %g\n",i,j,a[i][j],b[i][j],c[i][j],d[i][j],gammalj[i][j],gammacoul[i][j],cut_lj[i][j],cut_coul[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairSHIKWolf::single(int i, int j, int itype, int jtype,
                               double rsq,
                               double factor_coul, double factor_lj,
                               double &fforce)
{
  double r2inv,r6inv,r24inv,r2invcutlj,r6invcutlj,r24invcutlj,rexpcutlj,ljsmooth,v_ljcut,v_lj,coulsmooth,rcutcoul,rcutlj,r,rexp,forcecoul,forcebuck,phicoul,phibuck;

  r2inv = 1.0/rsq;
  rcutcoul = sqrt(cut_coulsq[itype][jtype]); //change
  rcutlj = sqrt(cut_ljsq[itype][jtype]); //change
  r = sqrt(rsq);
  if (rsq < cut_coulsq[itype][jtype]){
    coulsmooth = exp(-(gammacoul[itype][jtype]*gammacoul[itype][jtype])/((r-rcutcoul)*(r-rcutcoul))); //change
    forcecoul = force->qqrd2e * atom->q[i]*atom->q[j]*coulsmooth*((1.0/r-r/(rcutcoul*rcutcoul))-(1.0/r-1.0/rcutcoul+(r-rcutcoul)/(rcutcoul*rcutcoul))*(2.0*r*(gammacoul[itype][jtype]*gammacoul[itype][jtype])/((r-rcutcoul)*(r-rcutcoul)*(r-rcutcoul))));
 } else forcecoul = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    r24inv = r6inv*r6inv*r6inv*r6inv; //change
    r2invcutlj = 1.0/cut_ljsq[itype][jtype]; //change
    r6invcutlj = r2invcutlj*r2invcutlj*r2invcutlj; //change
    r24invcutlj = r6invcutlj*r6invcutlj*r6invcutlj*r6invcutlj; //change
    rexpcutlj = exp(-rcutlj*b[itype][jtype]); //change
    ljsmooth = exp(-(gammalj[itype][jtype]*gammalj[itype][jtype])/((r-rcutlj)*(r-rcutlj))); //change
    rexp = exp(-r*b[itype][jtype]);
    v_lj = a[itype][jtype]*rexp - c[itype][jtype]*r6inv + d[itype][jtype]*r24inv; //change
    v_ljcut = a[itype][jtype]*rexpcutlj - c[itype][jtype]*r6invcutlj + d[itype][jtype]*r24invcutlj; //change
  forcebuck = ((buck1[itype][jtype]*r*rexp - buck2[itype][jtype]*r6inv + buck3[itype][jtype]*r24inv) - (v_lj-v_ljcut)*r*2.0*(gammalj[itype][jtype]*gammalj[itype][jtype])/((r-rcutlj)*(r-rcutlj)*(r-rcutlj)))*ljsmooth ; //Major change
  } else forcebuck = 0.0;
  fforce = (factor_coul*forcecoul + factor_lj*forcebuck) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq[itype][jtype]) {
    phicoul = force->qqrd2e * atom->q[i]*atom->q[j]*(1.0/r-1.0/rcutcoul+(r-rcutcoul)/(rcutcoul*rcutcoul))*coulsmooth;
    eng += factor_coul*phicoul;
  }
  if (rsq < cut_ljsq[itype][jtype]) {
    	phibuck = (v_lj-v_ljcut)*ljsmooth; //change
	eng += factor_lj*phibuck;
  }
  return eng;
}
