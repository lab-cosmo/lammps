Design taken from pace_kokkos with accessory files in their own directory.
Will probably need some cmake magic to copy them here from somewhere else.

To be compiled as
cmake ../cmake/ -DPKG_ML-METATENSOR=ON -DPKG_KOKKOS=ON -DKokkos_ENABLE_CUDA=ON

Run the example with
../../build/lmp -k on g 1 -pk kokkos newton on -in in.metatensor
