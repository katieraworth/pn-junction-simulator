import numpy as np
import matplotlib.pyplot as plt
from constants import *
from physics import built_in_voltage, depletion_edges, diode_current, saturation_current


def equilibrium_profiles(NA, ND, points=2000):
    """Charge, field, potential and band edges across the junction at zero bias."""
    Vbi = built_in_voltage(NA, ND)
    xp, xn = depletion_edges(NA, ND)
    margin = 0.5 * (xp + xn)
    x = np.linspace(-xp - margin, xn + margin, points)   # cm

    p_dep = (x >= -xp) & (x < 0)     # depleted p-side
    n_dep = (x >= 0) & (x <= xn)     # depleted n-side

    # Charge density divided by q (cm^-3)
    rho = np.zeros_like(x)
    rho[p_dep] = -NA
    rho[n_dep] = ND

    # Electric field (V/cm)
    E = np.zeros_like(x)
    E[p_dep] = -q * NA * (x[p_dep] + xp) / eps_s
    E[n_dep] = -q * ND * (xn - x[n_dep]) / eps_s

    # Potential (V), 0 on the p-side and Vbi on the n-side
    psi = np.zeros_like(x)
    psi[p_dep] = q * NA * (x[p_dep] + xp)**2 / (2 * eps_s)
    psi[n_dep] = Vbi - q * ND * (xn - x[n_dep])**2 / (2 * eps_s)
    psi[x > xn] = Vbi

    # Band edges (eV), with the Fermi level at 0
    Ec = Eg / 2 + Vt * np.log(NA / ni) - psi
    Ev = Ec - Eg
    Ei = Ec - Eg / 2

    return x, rho, E, psi, Ec, Ev, Ei, xp, xn


def plot_equilibrium(NA, ND, filename="equilibrium.png"):
    x, rho, E, psi, Ec, Ev, Ei, xp, xn = equilibrium_profiles(NA, ND)
    x_um = x * 1e4

    fig, ax = plt.subplots(4, 1, sharex=True, figsize=(8, 10))
    fig.suptitle(f"Si PN junction at equilibrium (NA = {NA:.0e}, ND = {ND:.0e} cm$^{{-3}}$)")

    ax[0].plot(x_um, rho)
    ax[0].set_ylabel("Charge / q (cm$^{-3}$)")

    ax[1].plot(x_um, E / 1e3)
    ax[1].set_ylabel("Electric field (kV/cm)")

    ax[2].plot(x_um, psi)
    ax[2].set_ylabel("Potential (V)")

    ax[3].plot(x_um, Ec, label="Ec")
    ax[3].plot(x_um, Ev, label="Ev")
    ax[3].plot(x_um, Ei, "--", label="Ei")
    ax[3].axhline(0, color="k", linestyle=":", label="Ef")
    ax[3].set_ylabel("Energy (eV)")
    ax[3].set_xlabel("Position (µm)")
    ax[3].legend(loc="upper right")

    for a in ax:
        a.axvline(-xp * 1e4, color="gray", linewidth=0.5)
        a.axvline(xn * 1e4, color="gray", linewidth=0.5)
        a.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.show()

def plot_iv(NA, ND, filename="iv_curve.png"):
    V = np.linspace(-1, 0.8, 500)
    I = diode_current(V, NA, ND)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    fig.suptitle("Ideal diode I-V characteristic")

    ax[0].plot(V, I * 1e3)
    ax[0].set_xlabel("Voltage (V)")
    ax[0].set_ylabel("Current (mA)")
    ax[0].set_title("Linear scale")

    ax[1].semilogy(V, np.abs(I))
    ax[1].set_xlabel("Voltage (V)")
    ax[1].set_ylabel("|Current| (A)")
    ax[1].set_title("Log scale")

    for a in ax:
        a.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.show()

if __name__ == "__main__":
    NA, ND = 1e17, 1e16
    print(f"I0 = {saturation_current(NA, ND):.2e} A")
    plot_equilibrium(NA, ND)
    plot_iv(NA, ND)