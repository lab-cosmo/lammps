from lammps.mliap.mliap_unified_abc import MLIAPUnified
import numpy as np
import jax
import jax.numpy as jnp
from jax import jit
from functools import partial

#               /- ensure epsilon and sigma are treated as compile-time constants
@partial(jit, static_argnums=(0, 1))
def lj_potential(epsilon: float, sigma: float, rij):
    # a pure function we can differentiate:
    def _tot_e(rij):
        r2inv = 1.0 / jnp.sum(rij ** 2, axis=1)
        r6inv = r2inv * r2inv * r2inv

        lj1 = 4.0 * epsilon * sigma**12
        lj2 = 4.0 * epsilon * sigma**6

        eij = r6inv * (lj1 * r6inv - lj2)
        return jnp.sum(eij)
    #                   /-  construct a function computing _tot_e and its derivative
    tot_e, fij = jax.value_and_grad(_tot_e)(rij)
    return tot_e, fij


class MLIAPUnifiedJAX(MLIAPUnified):
    """Test implementation for MLIAPUnified."""

    epsilon: float
    sigma: float

    def __init__(self, element_types, epsilon=1.0, sigma=1.0, rcutfac=1.25):
        # ARGS: interface, element_types, ndescriptors, nparams, rcutfac
        super().__init__(None, element_types, 1, 3, rcutfac)
        # Mimicking the LJ pair-style:
        # pair_style lj/cut 2.5
        # pair_coeff * * 1 1
        self.epsilon = epsilon
        self.sigma = sigma
        # TODO: Take this from the LAMMPS Cython side.
        self.npair_max = 250000

    def compute_gradients(self, data):
        """Test compute_gradients."""

    def compute_descriptors(self, data):
        """Test compute_descriptors."""

    def compute_forces(self, data):
        """Test compute_forces."""
        rij  = data.rij

        # TODO: Take max npairs from the LAMMPS Cython side.
        if (data.npairs > self.npair_max):
            self.npair_max = data.npairs

        npad = self.npair_max - data.npairs
        # TODO: Take pre-padded rij from the LAMMPS Cython side.
        #       This might account for ~2-3x slowdown compared to original LJ.
        rij = np.pad(rij, ((0,npad), (0,0)), 'constant')

        e_tot, fij = lj_potential(rij)

        data.energy = e_tot.item()
        data.update_pair_forces(np.array(fij, dtype=np.float64))
