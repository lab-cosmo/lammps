# Preset that turns on all existing packages. Using the combination
# of this preset followed by the nolib.cmake preset should configure
# a LAMMPS binary, with as many packages included, that can be compiled
# with just a working C++ compiler and an MPI library.

set(ALL_PACKAGES
  ADIOS
  AMOEBA
  ASPHERE
  ATC
  AWPMD
  BOCS
  BODY
  BPM
  BROWNIAN
  CG-DNA
  CG-SPICA
  CLASS2
  COLLOID
  COLVARS
  COMPRESS
  CORESHELL
  DIELECTRIC
  DIFFRACTION
  DIPOLE
  DPD-BASIC
  DPD-MESO
  DPD-REACT
  DPD-SMOOTH
  DRUDE
  ELECTRODE
  EFF
  EXTRA-COMPUTE
  EXTRA-DUMP
  EXTRA-FIX
  EXTRA-MOLECULE
  EXTRA-PAIR
  FEP
  GPU
  GRANULAR
  H5MD
  INTEL
  INTERLAYER
  KIM
  KOKKOS
  KSPACE
  LATBOLTZ
  LEPTON
  MACHDYN
  MANIFOLD
  MANYBODY
  MC
  MDI
  MEAM
  MESONT
  MGPT
  MISC
  ML-HDNNP
  ML-IAP
  ML-PACE
  ML-POD
  ML-QUIP
  ML-RANN
  ML-SNAP
  MOFFF
  MOLECULE
  MOLFILE
  MPIIO
  NETCDF
  OPENMP
  OPT
  ORIENT
  PERI
  PHONON
  PLUGIN
  PLUMED
  POEMS
  PTM
  PYTHON
  QEQ
  QMMM
  QTB
  REACTION
  REAXFF
  REPLICA
  RIGID
  SCAFACOS
  SHOCK
  SMTBQ
  SPH
  SPIN
  SRD
  TALLY
  UEF
  VORONOI
  VTK
  YAFF)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
