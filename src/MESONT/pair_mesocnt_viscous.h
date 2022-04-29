/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
PairStyle(mesocnt/viscous, PairMesoCNTViscous);
#else

#ifndef LMP_PAIR_MESOCNT_VISCOUS_H
#define LMP_PAIR_MESOCNT_VISCOUS_H

#include "pair_mesocnt.h"

namespace LAMMPS_NS {
class PairMesoCNTViscous : public PairMesoCNT {
 public:
  using PairMesoCNT::PairMesoCNT;

  void compute(int, int) override;
  void coeff(int, char **) override;

 protected:
  double a1, a2, b2, vswitch;
};

}    // namespace LAMMPS_NS

#endif
#endif

