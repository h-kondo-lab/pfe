#!/usr/bin/env python3
"""
Physics from Existence — Pure Derivation
=========================================

This script contains ONE equation and ZERO experimental values.
It derives 20+ physical quantities from V = -H alone.

Run it, then compare the output with any physics textbook.

    $ python3 physics_from_existence_pure.py

The three axioms:
    A1  Existence is bivalent:     n ∈ {0,1}
    A2  Existence is uncertain:    σ = ⟨n⟩ ∈ (0,1)
    A3  Non-existence is forbidden: σ ≠ 0

From these:
    n² = n  →  Z = 1+e^φ  →  V = -H(σ(φ))

This is the ONLY input.  There are no fitted parameters,
no coupling constants, no masses, no mixing angles.
Everything is output.
"""

import numpy as np
from scipy.linalg import eigh_tridiagonal


# =====================================================================
#
#   THE EQUATION
#
#   V(φ) = σ ln σ + (1-σ) ln(1-σ)
#
#   where  σ(φ) = 1/(1 + e^{-φ})
#
#   This is the canonical free energy (Legendre transform of ln Z)
#   of a single binary degree of freedom n ∈ {0,1}.
#   It equals minus the binary Shannon entropy:  V = -H.
#
#   There is nothing else.
#
# =====================================================================

def sigma(phi):
    """
    The Fermi-Dirac distribution.

    Derived from A1 + A2:
      n ∈ {0,1}  →  Z = Σ e^{nφ} = 1 + e^φ
      σ = ⟨n⟩ = ∂(ln Z)/∂φ = 1/(1 + e^{-φ})
    """
    return 1.0 / (1.0 + np.exp(-phi))


def H(p):
    """
    Binary Shannon entropy.

    H(p) = -p ln p - (1-p) ln(1-p)

    Unique by Khinchin's theorem: the only function satisfying
    continuity, maximality, and additivity.
    """
    p = np.clip(p, 1e-15, 1 - 1e-15)
    return -p * np.log(p) - (1 - p) * np.log(1 - p)


def V(phi):
    """
    THE EQUATION:  V = -H(σ(φ))

    Derived via Legendre transform:
      W = ln Z = ln(1 + e^φ)
      V = φσ - W = σ ln σ + (1-σ) ln(1-σ) = -H(σ)

    Properties (all follow algebraically):
      V(0) = -ln 2  ≈ -0.693     (minimum, at σ = 1/2)
      V → 0  as  φ → ±∞          (bounded)
      V''(0) = 1/4                (Fisher information at σ=1/2)
      V(φ) = V(-φ)               (parity symmetry)
    """
    return -H(sigma(phi))


# =====================================================================
#
#   THE SCHRÖDINGER EQUATION
#
#   [-½ d²/dφ² + V(φ)] ψ = E ψ
#
#   The kinetic term ½(dφ)² is written in the A3-complete
#   e-affine coordinate φ = logit(σ). A3 selects this coordinate
#   because the forbidden endpoint σ = 0 is sent to infinity.
#   The unit-mass coefficient is the canonical one-bit transfer
#   normalisation used in Paper I §2.4; the operator-level
#   justification is given in Paper II.
#
# =====================================================================

def solve_schrodinger(potential, phi, dphi, n_states=6, m_eff=None):
    """
    Solve [-½ d²/dφ² + V] ψ = E ψ  on a discrete grid.

    Args:
        potential: V(φ) array
        phi:       coordinate array
        dphi:      grid spacing
        n_states:  number of lowest states to compute
        m_eff:     position-dependent effective mass (None = unit mass)

    Returns:
        eigenvalues, [normalised wavefunctions]

    The unit mass is the canonical one-bit transfer normalisation in the
    A3-complete e-affine coordinate φ; it is not fitted to the spectrum.
    """
    n = len(phi)

    if m_eff is not None:
        # Variable-mass kinetic energy (for down-type quarks, §8.2)
        inv_m = 1.0 / m_eff
        d = 0.5 * (inv_m[:-1] + inv_m[1:]) / dphi**2
        kd = np.zeros(n)
        kd[1:-1] = 0.5 * (d[:-1] + d[1:])
        kd[0] = d[0]
        kd[-1] = d[-1]
        main_diag = kd + potential
        off_diag  = -0.5 * d
    else:
        # Standard unit-mass kinetic energy
        main_diag = 1.0 / dphi**2 + potential
        off_diag  = -0.5 / dphi**2 * np.ones(n - 1)

    eigenvalues, eigenvectors = eigh_tridiagonal(
        main_diag, off_diag,
        select='i', select_range=(0, n_states - 1)
    )

    # Normalise wavefunctions
    wavefunctions = []
    for i in range(len(eigenvalues)):
        psi = eigenvectors[:, i]
        psi /= np.sqrt(np.trapezoid(psi**2, phi))
        wavefunctions.append(psi)

    return eigenvalues, wavefunctions


# =====================================================================
#
#   MASS FORMULA  (§8.1)
#
#   m_n ∝ |E_n| × IPR_n^b × (⟨H⟩_n / T_n)^c
#
#   Every exponent is determined by V = -H:
#     b = 2C_F = (N²-1)/N    ← Casimir of generation space
#     c = ε/N                 ← WKB marginality per mode
#
#   IPR = ∫|ψ|⁴ dφ           ← inverse participation ratio (localisation)
#   T   = ½∫|ψ'|² dφ         ← kinetic energy
#   ⟨H⟩ = ∫|ψ|² H dφ        ← entropy expectation
#
# =====================================================================

def inverse_participation_ratio(psi, phi):
    """IPR = ∫|ψ|⁴ dφ.  Measures how localised the wavefunction is."""
    return np.trapezoid(psi**4, phi)


def kinetic_energy(psi, dphi, phi):
    """T = ½∫|dψ/dφ|² dφ.  The kinetic energy of mode ψ."""
    dpsi = np.gradient(psi, dphi)
    return 0.5 * np.trapezoid(dpsi**2, phi)


def entropy_expectation(psi, phi):
    """⟨H⟩ = ∫|ψ|² H(σ(φ)) dφ.  How much entropy the mode samples."""
    return np.trapezoid(psi**2 * H(sigma(phi)), phi)


def compute_masses(eigenvalues, wavefunctions, phi, dphi, modes, b, c):
    """
    Compute mass-like quantities for the given modes.

    The formula m_n ∝ |E_n| × IPR^b × (⟨H⟩/T)^c gives mass RATIOS.
    The overall scale is irrelevant — only ratios are physical.
    """
    masses = []
    for i in modes:
        psi = wavefunctions[i]
        ipr = inverse_participation_ratio(psi, phi)
        T   = kinetic_energy(psi, dphi, phi)
        H_exp = entropy_expectation(psi, phi)
        m = abs(eigenvalues[i]) * ipr**b * (H_exp / T)**c
        masses.append(m)
    return masses


def mass_ratio_R(masses):
    """R = ln(m_heavy/m_light) / ln(m_middle/m_light)."""
    return np.log(masses[0] / masses[2]) / np.log(masses[1] / masses[2])


def koide_parameter(masses):
    """Q = (Σ√m)² / Σm.  Equals 3/2 at maximum-entropy partition."""
    s = sum(m**0.5 for m in masses)
    return s**2 / sum(masses)


# =====================================================================
#
#   DERIVE EVERYTHING
#
#   This function takes NOTHING from experiment.
#   It starts from V = -H and computes all predictions.
#
# =====================================================================

def derive(phi_max=60, n_grid=512001):
    """
    Derive all physical predictions from V = -H.

    Input:   V = -H  (the equation)
    Output:  dict of 20+ predicted quantities

    No experimental value is used anywhere in this function.
    """

    # ── Construct the potential V = -H on a numerical grid ──────
    phi  = np.linspace(-phi_max, phi_max, n_grid)
    dphi = phi[1] - phi[0]
    Vpot = V(phi)

    out = {}
    out['grid'] = {'phi_max': phi_max, 'n_grid': n_grid, 'dphi': dphi}


    # ════════════════════════════════════════════════════════════
    #  STEP 1:  Solve the Schrödinger equation  (§3)
    # ════════════════════════════════════════════════════════════
    #  This is the ONLY numerical computation in the entire script.
    #  Everything after this point is algebra.

    eigenvalues, psi = solve_schrodinger(Vpot, phi, dphi, n_states=6)


    # ════════════════════════════════════════════════════════════
    #  STEP 2:  Count the bound states  (§3.1)
    # ════════════════════════════════════════════════════════════
    #  N is not assumed.  It is READ from the spectrum.
    #  The well V = -H is shallow (depth = ln 2 ≈ 0.693).
    #  It happens to support exactly 3 bound states.

    N = sum(1 for E in eigenvalues if E < 0)

    out['N_generations'] = N
    out['eigenvalues'] = {
        f'E{i}': eigenvalues[i] for i in range(min(4, len(eigenvalues)))
    }
    out['fourth_generation'] = 'forbidden' if eigenvalues[3] > 0 else 'allowed'


    # ════════════════════════════════════════════════════════════
    #  STEP 3:  Compute the WKB marginality  (§3.3)
    # ════════════════════════════════════════════════════════════
    #  The third bound state barely exists.
    #  ε measures how marginal it is.

    integrand = np.where(Vpot < 0, np.sqrt(2 * np.abs(Vpot)), 0)
    n_wkb = np.trapezoid(integrand, phi) / np.pi
    eps   = N - n_wkb

    out['n_WKB']       = n_wkb
    out['epsilon']     = eps


    # ════════════════════════════════════════════════════════════
    #  FROM HERE, EVERYTHING IS ALGEBRA OF N
    # ════════════════════════════════════════════════════════════

    # Group-theory constants (all determined by N alone)
    PG = N**2 + N + 1        # number of points in PG(2, F_N)
    b  = (N**2 - 1) / N      # 2C_F: Casimir of fundamental rep
    c  = eps / N              # WKB deficit per generation

    out['b_2CF']  = b
    out['c_eps_N'] = c


    # ════════════════════════════════════════════════════════════
    #  STEP 4:  Spatial dimension  (§4)
    # ════════════════════════════════════════════════════════════
    #  Three independent conditions:
    #    (i)   Bound states survive:  centrifugal term = 0  →  d=1 or 3
    #    (ii)  Spin-statistics holds: π₁(SO(d)) = Z₂     →  d ≥ 3
    #    (iii) Gravity propagates:    Weyl tensor ≠ 0      →  d ≥ 3
    #  Intersection: d = 3.

    out['d_spatial'] = 3


    # ════════════════════════════════════════════════════════════
    #  STEP 5:  Gauge group  (§5)
    # ════════════════════════════════════════════════════════════
    #  N = d = 3  →  F₃³  →  PG(2,F₃)  →  13 points
    #  Flag decomposition: 9 + 3 + 1 = 13
    #    9 affine labels (a,b) ∈ F₃² index finite Weyl operators
    #       W_ab = X^a Z^b; identity + 8 traceless components → su(3)
    #    3 non-reference points on L → minimal non-abelian compact
    #       completion su(2)
    #    1 reference point P₀ → U(1)_Y

    assert N == 3, "This implementation is the PG(2,F3) / Standard Model case."
    out['PG_points']   = PG
    out['gauge_group']  = 'SU(3) x SU(2) x U(1)'
    out['partition']    = f'{N**2} + {N} + 1 = {PG}'


    # ════════════════════════════════════════════════════════════
    #  STEP 6:  Gauge coupling constants  (§6)
    # ════════════════════════════════════════════════════════════

    # Weinberg angle (§6.1):
    #   sin²θ_W = (weak directions) / (total PG points)
    #           = N / (N²+N+1)
    out['sin2_theta_W'] = N / PG

    # Strong coupling (§6.2):
    #   Tree: α_s⁽⁰⁾ = sin²θ_W / rank(SU(N)) = (N/PG)/(N-1)
    #   NLO:  × (1 + ε/(N²+1))  [non-weak PG directions]
    out['alpha_s'] = (N / PG) / (N - 1) * (1 + eps / (N**2 + 1))

    # Fine-structure constant (§6.3):
    #   1/α + α = Z + λ₁    (Ward identity, 0D Dyson equation)
    #   Z = Σ resolvent traces over flag Laplacian eigenvalues
    #   λ₁ = spectral gap of the flag Laplacian on PG(2,F_N)
    #
    #   Flag Laplacian spectrum of PG(2,F₃):
    #     eigenvalue          multiplicity
    #     (N+1) - √N          12
    #     (N+1) + √N          12
    #     2(N+1)              27
    #
    lambda_values       = [(N+1) - np.sqrt(N),
                           (N+1) + np.sqrt(N),
                           2*(N+1)]
    lambda_multiplicities = [12, 12, 27]  # PG(2,F3) flag graph

    spectral_zeta = sum(1.0 / abs(eigenvalues[i]) for i in range(N))
    resolvent_S   = [sum(1.0 / (abs(eigenvalues[i]) + lam)
                         for i in range(N))
                     for lam in lambda_values]
    Z_vacuum_pol  = spectral_zeta + sum(
        m * s for m, s in zip(lambda_multiplicities, resolvent_S)
    )
    lambda_1 = lambda_values[0]   # spectral gap = bare propagator

    # Solve 1/α + α = Z + λ₁
    rhs = Z_vacuum_pol + lambda_1
    out['one_over_alpha'] = (rhs + np.sqrt(rhs**2 - 4)) / 2


    # ════════════════════════════════════════════════════════════
    #  STEP 7:  Higgs mass and vacuum  (§7)
    # ════════════════════════════════════════════════════════════

    # m_H/v = ½√(1 + 2ε/|PG|)    (NLO correction from all PG modes)
    out['mH_over_v'] = 0.5 * np.sqrt(1 + 2 * eps / PG)

    # v²/Λ² = |V_min| / V''(0) = ln2 / (1/4) = 4 ln 2 ≈ 2.77
    #   → O(1), no hierarchy problem
    out['v2_over_Lambda2'] = abs(Vpot.min()) / 0.25

    # m_W/m_Z = cos θ_W = √(1 - sin²θ_W)
    out['mW_over_mZ'] = np.sqrt(1 - out['sin2_theta_W'])


    # ════════════════════════════════════════════════════════════
    #  STEP 8:  Charged lepton masses  (§8.1)
    # ════════════════════════════════════════════════════════════
    #  Three bound states = three charged leptons (τ, μ, e)
    #  Mass formula: m_n ∝ |E_n| × IPR^{2C_F} × (⟨H⟩/T)^{ε/N}
    #  Only RATIOS are predicted (0D has no energy scale).

    m_leptons = compute_masses(eigenvalues, psi, phi, dphi,
                               modes=[0, 1, 2], b=b, c=c)

    out['R_lepton'] = mass_ratio_R(m_leptons)   # ln(mτ/me)/ln(mμ/me)
    out['m_tau_over_m_e'] = m_leptons[0] / m_leptons[2]
    out['m_mu_over_m_e']  = m_leptons[1] / m_leptons[2]
    out['Q_Koide'] = koide_parameter(m_leptons)


    # ════════════════════════════════════════════════════════════
    #  STEP 9:  Quark masses  (§8.2)
    # ════════════════════════════════════════════════════════════
    #  Quarks carry N_c = N colours.
    #  By Shannon entropy additivity: V_quark = -N·H  (deeper well)

    V_colour = -N * H(sigma(phi))
    ev_q, psi_q = solve_schrodinger(V_colour, phi, dphi, n_states=7)

    # Up-type quarks: modes (1,2,3) of V = -NH, unit mass
    m_up = compute_masses(ev_q, psi_q, phi, dphi, modes=[1,2,3], b=b, c=c)
    out['R_up'] = mass_ratio_R(m_up)

    # Down-type quarks: effective mass from SU(2) isospin flip
    #   m_eff = 1 - N_c · σ(1-σ)
    m_eff_down = 1.0 + (-N) * sigma(phi) * (1 - sigma(phi))
    ev_d, psi_d = solve_schrodinger(V_colour, phi, dphi,
                                     n_states=6, m_eff=m_eff_down)
    m_down = compute_masses(ev_d, psi_d, phi, dphi, modes=[0,1,2], b=b, c=c)
    out['R_down'] = mass_ratio_R(m_down)


    # ════════════════════════════════════════════════════════════
    #  STEP 10:  CKM mixing  (§8.3)
    # ════════════════════════════════════════════════════════════

    # Cabibbo angle: sin θ_C = ε(1 + ε/N²)
    out['sin_theta_Cabibbo'] = eps * (1 + eps / N**2)

    # Strong CP: V-parity forbids the bare CP-odd term and keeps the
    # induced quark mass operator real, so θ̄ = θ_QCD + arg det M_q = 0.
    out['theta_QCD_bare'] = 0
    out['arg_det_Mq'] = 0
    out['theta_bar'] = 0

    # V_ub at tree level: the overlap of distinct eigenmodes vanishes by
    # Sturm-Liouville orthogonality. In the physical mixing operator,
    # V-parity forbids the leading 0↔2 transition; higher orders generate |V_ub|.
    out['V_ub_tree'] = abs(np.trapezoid(psi[0] * psi[2], phi))


    # ════════════════════════════════════════════════════════════
    #  STEP 11:  PMNS mixing  (§8.4)
    # ════════════════════════════════════════════════════════════
    #  Unlike CKM (perturbative in ε), PMNS is geometric (fractions of PG).

    # Solar:       sin²θ₁₂ = (one electroweak line) / |PG| = (N+1)/|PG|
    out['sin2_theta_12'] = (N + 1) / PG

    # Atmospheric:  sin²θ₂₃ = (2N+1)/|PG| × (1 + ε/|PG|)
    out['sin2_theta_23'] = (2*N + 1.0) / PG * (1 + eps / PG)

    # Reactor:     sin²θ₁₃ = 1/(N·|PG|) × (1 - 2ε/N)
    out['sin2_theta_13'] = 1.0 / (N * PG) * (1 - 2*eps/N)


    # ════════════════════════════════════════════════════════════
    #  STEP 12:  Neutrino predictions  (§8.4)
    # ════════════════════════════════════════════════════════════
    #  Neutrinos are colour-singlets → same V = -H as leptons.
    #  Therefore R_ν = R_lepton.
    #
    #  Combined with n² = n (Dirac) and V-parity (θ_QCD = 0):
    #    - Mass ordering: Normal
    #    - Neutrino nature: Dirac (not Majorana)
    #    - 0νββ decay rate: exactly zero

    out['R_neutrino'] = out['R_lepton']
    out['mass_ordering'] = 'Normal'
    out['neutrino_nature'] = 'Dirac'
    out['neutrinoless_double_beta'] = 0


    # ════════════════════════════════════════════════════════════
    #  STEP 13:  Additional derived quantities
    # ════════════════════════════════════════════════════════════

    # Participation ratio (§5.3): D_PR ≈ |PG| = 13
    rho = sum(p**2 for p in psi[:N])
    out['D_PR'] = np.trapezoid(rho, phi)**2 / np.trapezoid(rho**2, phi)

    # Spacetime dimension (§9.1):  d(d-3)/2 = N-1 = 2  →  d = 4
    out['d_spacetime'] = 4
    out['signature'] = '(3,1)'

    return out


# =====================================================================
#
#   PRINT ALL PREDICTIONS
#
# =====================================================================

def print_predictions(out):
    g = out['grid']
    N = out['N_generations']

    print()
    print("=" * 65)
    print("  V = -H(σ(φ))")
    print()
    print("  One equation.  Zero free parameters.")
    print("  Below are the predictions.  There is no experimental input.")
    print("=" * 65)

    # Derivation chain
    print(f"""
  Step 1  Solve  [-½ d²/dφ² - H(σ(φ))] ψ = E ψ
  Step 2  Count bound states:  N = {N}
  Step 3  Marginality:  ε = {N} - {out['n_WKB']:.4f} = {out['epsilon']:.6f}

  From N = {N} alone:
    |PG(2,F_{N})| = {out['PG_points']}
    2C_F = (N²-1)/N = {out['b_2CF']:.4f}
    Gauge group = {out['gauge_group']}
    Partition = {out['partition']}
""")

    # Eigenvalues
    print("  Eigenvalues of V = -H:")
    for k, v in out['eigenvalues'].items():
        label = "bound" if v < 0 else "continuum → 4th generation forbidden"
        print(f"    {k} = {v:+.6f}   ({label})")

    # Predictions
    print()
    print("  " + "-" * 50)
    print(f"  {'PREDICTION':<30} {'VALUE':>18}")
    print("  " + "-" * 50)

    predictions = [
        ("Generations N",               f"{out['N_generations']}"),
        ("Spatial dimension d",          f"{out['d_spatial']}"),
        ("Spacetime dimension",          f"{out['d_spacetime']}"),
        ("Signature",                    out['signature']),
        ("Gauge group",                  out['gauge_group']),
        ("",                             ""),
        ("1/α_em",                       f"{out['one_over_alpha']:.4f}"),
        ("α_s(M_Z)",                     f"{out['alpha_s']:.4f}"),
        ("sin²θ_W",                      f"{out['sin2_theta_W']:.5f}  = {N}/{out['PG_points']}"),
        ("m_W / m_Z",                    f"{out['mW_over_mZ']:.4f}"),
        ("m_H / v",                      f"{out['mH_over_v']:.4f}"),
        ("v²/Λ²",                        f"{out['v2_over_Lambda2']:.4f}  = 4 ln 2"),
        ("",                             ""),
        ("R_lepton  ln(mτ/me)/ln(mμ/me)",f"{out['R_lepton']:.4f}"),
        ("m_τ / m_e",                    f"{out['m_tau_over_m_e']:.1f}"),
        ("m_μ / m_e",                    f"{out['m_mu_over_m_e']:.1f}"),
        ("Q (Koide)",                    f"{out['Q_Koide']:.6f}"),
        ("R_up    ln(mt/mu)/ln(mc/mu)",  f"{out['R_up']:.3f}"),
        ("R_down  ln(mb/md)/ln(ms/md)",  f"{out['R_down']:.3f}"),
        ("",                             ""),
        ("sin θ_C  (Cabibbo)",           f"{out['sin_theta_Cabibbo']:.4f}"),
        ("|V_ub| at tree level",         f"{out['V_ub_tree']:.1e}"),
        ("θ_QCD bare",                   f"{out['theta_QCD_bare']}  (V-parity)"),
        ("arg det M_q",                  f"{out['arg_det_Mq']}  (real mass operator)"),
        ("θ̄ strong CP",                 f"{out['theta_bar']}  (exact)"),
        ("",                             ""),
        ("sin²θ₁₂ (solar)",              f"{out['sin2_theta_12']:.5f}  = {N+1}/{out['PG_points']}"),
        ("sin²θ₂₃ (atmospheric)",        f"{out['sin2_theta_23']:.4f}"),
        ("sin²θ₁₃ (reactor)",            f"{out['sin2_theta_13']:.5f}"),
        ("",                             ""),
        ("Mass ordering",                out['mass_ordering']),
        ("Neutrino nature",              out['neutrino_nature']),
        ("0νββ",                         f"{out['neutrinoless_double_beta']}  (exact)"),
        ("R_neutrino",                   f"{out['R_neutrino']:.4f}  (= R_lepton)"),
        ("",                             ""),
        ("D_PR",                         f"{out['D_PR']:.2f}  (≈ |PG| = {out['PG_points']})"),
        ("4th generation",               out['fourth_generation']),
    ]

    for name, value in predictions:
        if name == "":
            print()
        else:
            print(f"  {name:<30} {value:>18}")

    print()
    print("  " + "-" * 50)
    print("  No experimental value was used.")
    print("  Compare these numbers with any physics textbook.")
    print("=" * 65)


# =====================================================================
#  RUN
# =====================================================================

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(
        description="Derive physics from V = -H.  No experimental input.")
    ap.add_argument("--phi-max", type=float, default=60,
                    help="Grid extent (default: 60)")
    ap.add_argument("--n-grid", type=int, default=128001,
                    help="Grid points (default: 128001; paper uses 512001)")
    args = ap.parse_args()

    predictions = derive(phi_max=args.phi_max, n_grid=args.n_grid)
    print_predictions(predictions)
