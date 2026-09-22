import numpy as np
from constants import *

def built_in_voltage(NA, ND):
    """Built-in voltage (V) of a PN junction.
    NA = acceptor doping on the p-side (cm^-3)
    ND = donor doping on the n-side (cm^-3)
    """
    return Vt * np.log(NA * ND / ni**2)

def depletion_width(NA, ND, V=0):
    """Total depletion width W (cm). V = applied voltage (+ forward, - reverse)."""
    Vbi = built_in_voltage(NA, ND)
    return np.sqrt(2 * eps_s * (Vbi - V) / q * (NA + ND) / (NA * ND))


def depletion_edges(NA, ND, V=0):
    """How far the depletion region reaches into each side: (xp, xn) in cm."""
    W = depletion_width(NA, ND, V)
    xn = W * NA / (NA + ND)
    xp = W * ND / (NA + ND)
    return xp, xn


def max_field(NA, ND, V=0):
    """Peak electric field at the junction (V/cm)."""
    xp, xn = depletion_edges(NA, ND, V)
    return q * ND * xn / eps_s

def saturation_current(NA, ND, A=1e-4, Dn=20, Dp=10, tau=1e-6):
    """Reverse saturation current I0 (A).
    A = area (cm^2), Dn/Dp = diffusion coefficients (cm^2/s), tau = lifetime (s)
    """
    Ln = np.sqrt(Dn * tau)   # electron diffusion length (cm)
    Lp = np.sqrt(Dp * tau)   # hole diffusion length (cm)
    return q * A * ni**2 * (Dn / (Ln * NA) + Dp / (Lp * ND))


def diode_current(V, NA, ND, n=1, **kwargs):
    """Ideal diode (Shockley) current (A) at applied voltage V."""
    I0 = saturation_current(NA, ND, **kwargs)
    return I0 * (np.exp(V / (n * Vt)) - 1)

if __name__ == "__main__":
    NA, ND = 1e17, 1e16
    xp, xn = depletion_edges(NA, ND)
    print(f"Built-in voltage: {built_in_voltage(NA, ND):.3f} V")
    print(f"Depletion width:  {depletion_width(NA, ND) * 1e4:.3f} um")
    print(f"  p-side (xp):    {xp * 1e4:.3f} um")
    print(f"  n-side (xn):    {xn * 1e4:.3f} um")
    print(f"Max field:        {max_field(NA, ND) / 1e3:.1f} kV/cm")