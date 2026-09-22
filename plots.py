import numpy as np
import matplotlib.pyplot as plt
from constants import *
from physics import *


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
def plot_sweeps(filename="sweeps.png"):
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle("Parameter sweeps")
    ND_range = np.logspace(14, 18, 200)

    # (a) and (b): built-in voltage and depletion width vs. doping
    for NA in [1e15, 1e16, 1e17, 1e18]:
        ax[0, 0].semilogx(ND_range, built_in_voltage(NA, ND_range), label=f"NA = {NA:.0e}")
        ax[0, 1].loglog(ND_range, depletion_width(NA, ND_range) * 1e4, label=f"NA = {NA:.0e}")
    ax[0, 0].set_xlabel("ND (cm$^{-3}$)")
    ax[0, 0].set_ylabel("Built-in voltage (V)")
    ax[0, 1].set_xlabel("ND (cm$^{-3}$)")
    ax[0, 1].set_ylabel("Depletion width (µm)")

    # (c): depletion width vs. applied voltage
    V = np.linspace(-10, 0.5, 300)
    ax[1, 0].plot(V, depletion_width(1e17, 1e16, V) * 1e4)
    ax[1, 0].set_xlabel("Applied voltage (V)")
    ax[1, 0].set_ylabel("Depletion width (µm)")
    ax[1, 0].set_title("NA = 1e17, ND = 1e16")

    # (d): I-V curve at different temperatures
    V = np.linspace(0.05, 0.9, 300)
    for T in [250, 300, 350]:
        ax[1, 1].semilogy(V, diode_current(V, 1e17, 1e16, T=T), label=f"{T} K")
    ax[1, 1].axhline(1e-3, color="gray", linestyle=":", label="1 mA")
    ax[1, 1].set_xlabel("Voltage (V)")
    ax[1, 1].set_ylabel("Current (A)")

    for a in ax.flat:
        a.grid(alpha=0.3)
        if a.get_legend_handles_labels()[0]:
            a.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.show()

def plot_cv(NA, ND, A=1e-4, filename="cv_curve.png"):
    V = np.linspace(-10, 0, 50)
    C = junction_capacitance(NA, ND, V, A)

    # Simulated "measurement": add 1% random noise
    rng = np.random.default_rng(seed=1)
    C_meas = C * (1 + 0.01 * rng.standard_normal(V.size))

    # Fit a straight line to 1/C^2 vs V, then work backward
    y = 1 / C_meas**2
    slope, intercept = np.polyfit(V, y, 1)
    Vbi_fit = -intercept / slope                       # where the line crosses zero
    N_eff = -2 / (q * eps_s * A**2 * slope)            # = NA*ND / (NA + ND)
    ND_fit = 1 / (1 / N_eff - 1 / NA)                  # assumes NA is known

    Vbi_true = built_in_voltage(NA, ND)
    print(f"Vbi: true = {Vbi_true:.3f} V, fitted = {Vbi_fit:.3f} V "
          f"({(Vbi_fit - Vbi_true) / Vbi_true * 100:+.1f}%)")
    print(f"ND:  true = {ND:.2e}, fitted = {ND_fit:.2e} cm^-3 "
          f"({(ND_fit - ND) / ND * 100:+.1f}%)")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    fig.suptitle("C-V characteristic and doping extraction")

    ax[0].plot(V, C * 1e12, label="Model")
    ax[0].plot(V, C_meas * 1e12, "o", markersize=3, label="Simulated measurement")
    ax[0].set_xlabel("Voltage (V)")
    ax[0].set_ylabel("Capacitance (pF)")

    V_line = np.linspace(-10, Vbi_fit, 100)
    ax[1].plot(V, y, "o", markersize=3, label="Simulated measurement")
    ax[1].plot(V_line, slope * V_line + intercept, "--", label=f"Fit (Vbi = {Vbi_fit:.2f} V)")
    ax[1].axhline(0, color="k", linewidth=0.5)
    ax[1].set_xlabel("Voltage (V)")
    ax[1].set_ylabel("1/C$^2$ (1/F$^2$)")

    for a in ax:
        a.grid(alpha=0.3)
        a.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.show()

if __name__ == "__main__":
    NA, ND = 1e17, 1e16
    # plot_equilibrium(NA, ND)
    # plot_iv(NA, ND)
    # plot_sweeps()
    plot_cv(NA, ND)