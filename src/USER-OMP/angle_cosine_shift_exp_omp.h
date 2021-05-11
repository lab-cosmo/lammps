/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS

AngleStyle(cosine/shift/exp/omp,AngleCosineShiftExpOMP)

#else

#ifndef LMP_ANGLE_COSINE_SHIFT_EXP_OMP_H
#define LMP_ANGLE_COSINE_SHIFT_EXP_OMP_H

#include "angle_cosine_shift_exp.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class AngleCosineShiftExpOMP : public AngleCosineShiftExp, public ThrOMP {

 public:
  AngleCosineShiftExpOMP(class LAMMPS *lmp);
  virtual void compute(int, int);

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_BOND>
  void eval(int ifrom, int ito, ThrData * const thr);
};

}

#endif
#endif
