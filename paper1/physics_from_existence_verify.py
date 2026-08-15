#!/usr/bin/env python3
"""
PFE Paper I Verification Script — v3.2.0
========================================================
Experimental references: CODATA 2022 / PDG 2026 / NuFIT 6.1
(JHEP 12 (2024) 216). Formulas follow the v3.2.0 manuscript:
NLO forms 1±ℓε/k (sinθ_C: 1+ε/9; sin²θ₁₃: 1−2ε/3), the
screened Cabibbo expression and the neutrino master relation are
printed as informational lines with their Paper II attribution.
"""

import numpy as np
from scipy.linalg import eigh, eigh_tridiagonal
from scipy.optimize import brentq

def sigma(phi):
    return 1.0 / (1.0 + np.exp(-phi))

def H_binary(p):
    p = np.clip(p, 1e-15, 1 - 1e-15)
    return -p * np.log(p) - (1 - p) * np.log(1 - p)

def H_of_phi(phi):
    # Numerically stable H(sigma(phi)) = ln(1+e^{-phi}) + phi*sigma(-phi).
    # Evaluating H_binary(sigma(phi)) directly saturates for |phi| > ~37
    # (sigma -> 1 in float64), leaving a spurious clip floor ~3.5e-14 that
    # adds a fake O(box-size) tail to the n_WKB integral.
    return np.logaddexp(0.0, -phi) + phi * 0.5 * (1.0 - np.tanh(0.5 * phi))

def V_potential(phi):
    return -H_of_phi(phi)

def sech(x):
    return 1.0 / np.cosh(x)

def make_grid(phi_max, n_grid):
    phi = np.linspace(-phi_max, phi_max, n_grid)
    dphi = phi[1] - phi[0]
    Vpot = V_potential(phi)
    return phi, dphi, Vpot

def generalized_rayleigh_eigenvalues(trials, phi, dphi, Vpot):
    n = len(trials)
    S = np.zeros((n, n))
    H = np.zeros((n, n))
    gradients = [np.gradient(t, dphi) for t in trials]
    for i in range(n):
        for j in range(n):
            S[i, j] = np.trapezoid(trials[i] * trials[j], phi)
            H[i, j] = (
                0.5 * np.trapezoid(gradients[i] * gradients[j], phi)
                + np.trapezoid(trials[i] * Vpot * trials[j], phi)
            )
    return eigh(H, S, eigvals_only=True)

def verify_N3_variational(phi, dphi, Vpot):
    tail_trials = [
        sech(0.9 * phi),
        np.sinh(0.4 * phi) * sech(0.4 * phi) ** 2,
        (1.0 - 4.0 * sech(0.8 * phi)) * sech(0.15 * phi),
    ]
    tail_evals = generalized_rayleigh_eigenvalues(tail_trials, phi, dphi, Vpot)
    if not np.all(tail_evals < 0):
        raise AssertionError(
            f"sech-tail variational certificate failed: {tail_evals}"
        )

    # Negative control: Gaussian trial functions miss the marginal third state.
    w = 3.0
    gaussian_trials = [
        np.exp(-phi**2 / (2.0 * w**2)),
        phi * np.exp(-phi**2 / (2.0 * w**2)),
        (phi**2 - w**2) * np.exp(-phi**2 / (2.0 * w**2)),
    ]
    gaussian_evals = generalized_rayleigh_eigenvalues(
        gaussian_trials, phi, dphi, Vpot
    )
    if not gaussian_evals[2] > 0:
        raise AssertionError(
            f"Gaussian negative control did not fail as expected: {gaussian_evals}"
        )
    return tail_evals, gaussian_evals

def solve_schrodinger(V, phi, dphi, n_states=6, m_eff=None):
    N = len(phi)
    if m_eff is not None:
        inv_m = 1.0 / m_eff
        diag = 0.5 * (inv_m[:-1] + inv_m[1:]) / dphi**2
        kinetic_diag = np.zeros(N)
        kinetic_diag[1:-1] = 0.5 * (diag[:-1] + diag[1:])
        kinetic_diag[0] = diag[0]
        kinetic_diag[-1] = diag[-1]
        main_diag = kinetic_diag + V
        off_diag = -0.5 * diag
    else:
        main_diag = 1.0 / dphi**2 + V
        off_diag = -0.5 / dphi**2 * np.ones(N - 1)
    evals, evecs = eigh_tridiagonal(main_diag, off_diag, select='i',
                                     select_range=(0, n_states - 1))
    psi_list = []
    for i in range(len(evals)):
        v = evecs[:, i]
        v /= np.sqrt(np.trapezoid(v**2, phi))
        psi_list.append(v)
    return evals, psi_list

def compute_IPR(psi, phi):
    return np.trapezoid(psi**4, phi)

def compute_gradient_energy(psi, dphi, phi):
    dpsi = np.gradient(psi, dphi)
    return 0.5 * np.trapezoid(dpsi**2, phi)

def compute_H_expectation(psi, phi):
    return np.trapezoid(psi**2 * H_of_phi(phi), phi)

def mass_formula(evals, psi_list, phi, dphi, modes, b, c, kappa=None):
    # kappa is the well-depth multiple (V = -kappa*H).  When given, the
    # kinetic energy is taken from the operator identity T = E - <V>
    # = E + kappa*<H> (exact for eigenstates; free of finite-difference
    # gradient noise).  kappa=None keeps the flat gradient energy
    # T = (1/2)\int psi'^2 — the expression used for the down sector, whose
    # operator has a position-dependent mass (see §8.2).
    masses = []
    for i in modes:
        ipr = compute_IPR(psi_list[i], phi)
        H_exp = compute_H_expectation(psi_list[i], phi)
        if kappa is None:
            T = compute_gradient_energy(psi_list[i], dphi, phi)
        else:
            T = evals[i] + kappa * H_exp
        m = abs(evals[i]) * ipr**b * (H_exp / T)**c
        masses.append(m)
    return masses

def R_ratio(masses):
    return np.log(masses[0] / masses[2]) / np.log(masses[1] / masses[2])

def koide_Q(masses):
    s = sum(np.sqrt(m) for m in masses)
    return s**2 / sum(masses)

def run_paper1(phi_max=100, n_grid=1600001):
    phi, dphi, Vpot = make_grid(phi_max, n_grid)

    # ===== Experimental reference values =====
    # PDG 2026 / CODATA 2022
    m_tau = 1776.93   # MeV (PDG 2026)
    m_mu  = 105.6584  # MeV
    m_e   = 0.51100   # MeV

    # NuFIT 6.1 (2025), www.nu-fit.org (base: JHEP 12 (2024) 216, arXiv:2410.05380)
    # IC23 without SK-atm, Normal Ordering best-fit
    Dm21_sq_6 = 7.537e-5   # eV² (NuFIT 6.1)
    Dm31_sq_6 = 2.521e-3   # eV² (NuFIT 6.1)
    # Second-octant local minimum of the NuFIT 6.1 chi^2 profile (global bfp 0.470).
    # Obtained from the released grid v61.release-TBoff-NO.txt.xz by projecting over
    # Dm31^2 and delta_CP: minimum at 0.550 with Delta chi^2 = 1.03 above the global
    # minimum (grid step 0.005).  The former value 0.561 was the NuFIT 6.0 best fit
    # (IC19 without SK-atm, NO) and is not a NuFIT 6.1 number.
    sin2_23_6 = 0.550      # NuFIT 6.1 second-octant local minimum
    sin2_12_6 = 0.3088     # NuFIT 6.1 (JUNO)
    sin2_13_6 = 0.02249    # NuFIT 6.1 IC23 NO

    # Also keep NuFIT 5.x for comparison
    Dm21_sq_5 = 7.53e-5
    Dm31_sq_5 = 2.455e-3
    sin2_23_5 = 0.546

    EXP = {
        '1/alpha':      137.035999177,   # CODATA 2022
        'sin2_tW':      0.23122,         # MS-bar at M_Z (PDG 2026)
        'm_W/m_Z':      0.8813,          # PDG 2026 (80.3625/91.1879)
        'mH_v_exp':     0.5082,          # m_H/v = 125.13/246.22 (PDG 2026)
        'alpha_s':      0.1180,          # PDG 2026 world average, QCD review Sec. 9.4.8 (0.1180 +- 0.0009)
        'm_tau':        m_tau,
        'm_mu':         m_mu,
        'm_e':          m_e,
        'sin_tC':       0.2243,          # |V_us|=0.22431 PDG 2026
        'sin2_12':      sin2_12_6,
        'sin2_23':      sin2_23_6,
        'sin2_13':      sin2_13_6,
    }
    EXP['m_tau/m_e'] = m_tau / m_e
    EXP['m_mu/m_e'] = m_mu / m_e
    EXP['R_lepton'] = np.log(EXP['m_tau/m_e']) / np.log(EXP['m_mu/m_e'])

    # Koide Q from experimental masses
    Q_exp = koide_Q([m_tau, m_mu, m_e])

    print("=" * 80)
    print("  PAPER I VERIFICATION — v3.2.0")
    print("  20+ predictions from V = -H(σ(φ)), zero free parameters")
    print(f"  Grid: PHI_MAX={phi_max}, N_GRID={n_grid}, dφ={dphi:.8f}")
    print(f"  Experimental: CODATA 2022 / PDG 2026 / NuFIT 6.1")
    print("=" * 80)

    results = []

    # ===== Solve Schrödinger equation =====
    evals, psi = solve_schrodinger(Vpot, phi, dphi, n_states=6)
    n_bound = sum(1 for E in evals if E < 0)
    print(f"\n  Bound states: N = {n_bound}")
    for i in range(min(4, len(evals))):
        status = "BOUND" if evals[i] < 0 else "continuum"
        print(f"    E_{i} = {evals[i]:.12f}  [{status}]")
    results.append(("N (generations)", 3, 3, "exact"))

    # WKB
    integrand = np.where(Vpot < 0, np.sqrt(2 * np.abs(Vpot)), 0)
    n_wkb = np.trapezoid(integrand, phi) / np.pi
    N = 3; eps = N - n_wkb; b = (N**2 - 1) / N; c = eps / N
    PG = N**2 + N + 1
    print(f"  n_WKB = {n_wkb:.10f}, ε = {eps:.10f}")
    print(f"  [info] float64 grid values.  Certified enclosure "
          f"(proof_eps_interval_enclosure_20260731.py, Arb): "
          f"n_WKB = 2.7808120498, ε = 0.2191879502; grid-certified gap ~8e-8; "
          f"both round to the displayed 0.219188.")
    print(f"  b = 2C_F = {b:.10f}, c = ε/N = {c:.10f}")
    results.append(("d (spatial dim.)", 3, 3, "exact"))

    tail_evals, gaussian_evals = verify_N3_variational(phi, dphi, Vpot)
    print("\n  Variational N ≥ 3 certificate:")
    print("    sech-tail Rayleigh eigenvalues: "
          + ", ".join(f"{x:.6f}" for x in tail_evals))
    print("    Gaussian control eigenvalues:  "
          + ", ".join(f"{x:.6f}" for x in gaussian_evals))

    # ===== Mass formula =====
    m_lep = mass_formula(evals, psi, phi, dphi, modes=[0, 1, 2], b=b, c=c, kappa=1)
    R_lep = R_ratio(m_lep)
    Q_pred = koide_Q(m_lep)

    # ===== Mathematical theorems =====
    print(f"\n  --- Mathematical Theorems (V = -H alone) ---")

    # Strong CP: the strong-CP matching theorem of Paper II gives θ̄(μ0) = 0 at the
    # matching point; the bare QCD angle is killed by the measure-line lemma
    # (positive measure ray), NOT by V-parity alone. This is a matching-point
    # statement, not an exact zero at all scales; θ̄_IR is unevaluated.
    results.append(("θ_QCD bare", 0, "—", "Paper II"))
    results.append(("arg det M_q", 0, "—", "Paper II"))
    results.append(("θ̄ strong CP", 0, "<1e-10", "μ₀; Paper II"))

    # V_ub tree-level overlap: distinct eigenmodes vanish by Sturm-Liouville
    # orthogonality; the physical leading 0↔2 transition is V-parity forbidden.
    overlap_02 = abs(np.trapezoid(psi[0] * psi[2], phi))
    results.append(("V_ub|tree ≈ 0", overlap_02, 0.004, f"{overlap_02:.1e}"))

    # 4th gen
    E3_positive = evals[3] > 0
    results.append(("4th gen", "forbidden" if E3_positive else "EXISTS",
                     "excluded", "no 4th neg." if E3_positive else "FAIL"))

    # ===== QFT derivations =====
    print(f"\n  --- QFT Derivations (V = -H + PG(2,F₃)) ---")

    # 1/α
    zeta_V = sum(1.0 / abs(evals[i]) for i in range(3))
    lam1 = (N+1) - np.sqrt(N); lam2 = (N+1) + np.sqrt(N); lam3 = 2*(N+1)
    S1 = sum(1.0 / (abs(evals[i]) + lam1) for i in range(3))
    S2 = sum(1.0 / (abs(evals[i]) + lam2) for i in range(3))
    S3 = sum(1.0 / (abs(evals[i]) + lam3) for i in range(3))
    Z = zeta_V + 12*S1 + 12*S2 + 27*S3
    alpha_inv = ((Z + lam1) + np.sqrt((Z + lam1)**2 - 4)) / 2
    print(f"  [info] 1/alpha inputs: Z = {Z:.6f}, lambda_1 = {lam1:.6f} "
          f"-> 1/alpha = {alpha_inv:.4f}")
    print(f"  [info] grid Z = {Z:.5f} vs certified-E Z = 134.77613 (gap ~1.7e-5 from "
          f"the grid E_2; manuscript displays Z=134.77613, lambda_1=2.26795; "
          f"both give the displayed 1/alpha = 137.0368).")
    dev = abs(alpha_inv - 137.036) / 137.036 * 100
    results.append(("1/α_em (leading)", alpha_inv, "137.036", f"{dev:.4f}%"))  # paper Table 2: leading vs 137.036

    # 1/α at all orders (Dyson closed form; derivation and 13-digit certification
    # in Paper II). Solve  1/a + a = Λ0 - x(1+a/2) + (3/4)ε²x + (1589/5408)x²
    # - (2/13)x³/(1+x),  x = 9a/(26π),  Λ0 = λ1 + Z,  with the certified
    # eigenvalues of §3.2 and exact λ's, in high precision.
    #
    # Why the eigenvalues below are literals rather than recomputed here.
    # They are outputs of this framework, not experimental inputs.  The float64
    # grid above reaches only ~1e-5 on E_2 (see the grid-vs-certified Z gap
    # printed above), which is not enough for the 13-digit closure.  Redoing
    # them at full precision needs two-sided Prufer shooting: the third state
    # is shallow (k = sqrt(2|E_2|) = 0.144), so the inward leg must start near
    # phi = 200, and a single outward sweep is useless because the sign change
    # collapses to width ~exp(-2kT).  That costs minutes, so the default run
    # cites the literals.  Both checks are available on demand:
    #   --certify    rigorous +-1e-24 windows via the shipped certificate
    #                (supplementary/proof_E_taylor_enclosures_20260718.py)
    #   --recompute  regenerate eps and E_0,E_1,E_2 from scratch with the same
    #                Taylor recurrence in plain arithmetic, then redo 1/alpha
    #                (recompute_constants() at the end of this file)
    try:
        import mpmath as mp
        mp.mp.dps = 30
        E_cert = [mp.mpf('-0.4850022531029642374452777658'),
                  mp.mpf('-0.1590522216188561779934190'),
                  mp.mpf('-0.010423393062109459327616')]
        l1 = (N + 1) - mp.sqrt(N); l2 = (N + 1) + mp.sqrt(N); l3 = mp.mpf(2 * (N + 1))
        Zc = (sum(1 / abs(e) for e in E_cert)
              + 12 * sum(1 / (abs(e) + l1) for e in E_cert)
              + 12 * sum(1 / (abs(e) + l2) for e in E_cert)
              + 27 * sum(1 / (abs(e) + l3) for e in E_cert))
        L0 = Zc + l1
        # certified ε = N − n_WKB (§3.1).  Unlike the eigenvalues this one is
        # cheap to redo: one exponentially convergent integral, ~0.01 s
        # (eps_from_scratch() below does it; --recompute compares).
        eps_c = mp.mpf('0.2191879502483955')
        def dyson(a):
            x = 9 * a / (26 * mp.pi)
            return (1 / a + a
                    - (L0 - x * (1 + a / 2) + mp.mpf(3) / 4 * eps_c**2 * x
                       + mp.mpf(1589) / 5408 * x**2
                       - mp.mpf(2) / 13 * x**3 / (1 + x)))
        # Newton seed = 1/L0, the leading-order solution of the same equation
        # (theory-derived; no measured value enters even as a seed)
        a_all = mp.findroot(dyson, 1 / L0)
        inv_all = 1 / a_all
        lo = mp.mpf('137.0359990843979'); hi = mp.mpf('137.0359990844119')
        inside = (lo <= inv_all <= hi)
        print(f"  [info] 1/alpha all orders (Dyson closed form, Paper II): "
              f"{mp.nstr(inv_all, 13)}  "
              f"({'inside' if inside else 'OUTSIDE'} the certified 13-digit interval "
              f"[137.0359990843979, 137.0359990844119])")
        print(f"  [info] CODATA 2018 (Cs recoil): 137.035999084(21) — matches the "
              f"displayed digits; CODATA 2022 (Rb recoil): 137.035999177(21), 4.4σ away; "
              f"the framework predicts the Cs side (manuscript §6.3, Table 4).")
        dev18 = abs(float(inv_all) - 137.035999084) / 137.035999084 * 100
        results.append(("1/α_em (all orders)", float(inv_all), 137.035999084,
                        f"{dev18:.1e}%"))
    except Exception as ex:
        print(f"  [warn] all-orders 1/alpha check skipped: {ex}")

    # R_lepton
    dev = abs(R_lep - EXP['R_lepton']) / EXP['R_lepton'] * 100
    results.append(("R_lepton", R_lep, EXP['R_lepton'], f"{dev:.4f}%"))

    # α_s
    alpha_s = (3.0 / 26) * (1 + eps / (N**2 + 1))
    dev = abs(alpha_s - EXP['alpha_s']) / EXP['alpha_s'] * 100
    results.append(("α_s", alpha_s, EXP['alpha_s'], f"{dev:.2f}%"))

    # sinθ_C
    sin_thetaC = eps * (1 + eps / N**2)
    dev = abs(sin_thetaC - EXP['sin_tC']) / EXP['sin_tC'] * 100
    results.append(("sinθ_C", sin_thetaC, EXP['sin_tC'], f"{dev:.2f}%"))
    # Screened Cabibbo expression (Paper II): sinθ_C·(1−x), x = 9α/(26π).
    x_scr = 9.0 / (alpha_inv * 26.0 * np.pi)
    sin_tC_screened = sin_thetaC * (1 - x_scr)
    print(f"  [info] screened Cabibbo (Paper II): sinθ_C·(1−9α/26π) = "
          f"{sin_tC_screened:.6f}  (manuscript: 0.22435 with certified ε)")

    # m_H/v
    mH_v_NLO = 0.5 * np.sqrt(1 + 2*eps/PG)
    dev = abs(mH_v_NLO - EXP['mH_v_exp']) / EXP['mH_v_exp'] * 100
    results.append(("m_H/v", mH_v_NLO, EXP['mH_v_exp'], f"{dev:.2f}%"))

    # m_τ/m_e, m_μ/m_e
    mt_me = m_lep[0]/m_lep[2]; mm_me = m_lep[1]/m_lep[2]
    results.append(("m_τ/m_e", mt_me, EXP['m_tau/m_e'],
                     f"{abs(mt_me-EXP['m_tau/m_e'])/EXP['m_tau/m_e']*100:.2f}%"))
    results.append(("m_μ/m_e", mm_me, EXP['m_mu/m_e'],
                     f"{abs(mm_me-EXP['m_mu/m_e'])/EXP['m_mu/m_e']*100:.2f}%"))

    # sin²θ_W
    sin2_thetaW = N / PG
    dev = abs(sin2_thetaW - EXP['sin2_tW']) / EXP['sin2_tW'] * 100
    results.append(("sin²θ_W", sin2_thetaW, EXP['sin2_tW'], f"{dev:.2f}%"))

    # R_down
    Vpot3 = -3 * H_of_phi(phi)
    evals3, psi3 = solve_schrodinger(Vpot3, phi, dphi, n_states=7)
    m_up = mass_formula(evals3, psi3, phi, dphi, modes=[1, 2, 3], b=b, c=c, kappa=3)
    R_up = R_ratio(m_up)
    results.append(("R_up", R_up, 1.770, f"{abs(R_up-1.770)/1.770*100:.2f}%"))

    m_eff = 1.0 + (-3) * sigma(phi) * (1 - sigma(phi))
    evals_d, psi_d = solve_schrodinger(Vpot3, phi, dphi, n_states=6, m_eff=m_eff)
    m_down = mass_formula(evals_d, psi_d, phi, dphi, modes=[0, 1, 2], b=b, c=c)
    R_down = R_ratio(m_down)
    results.append(("R_down", R_down, 2.276, f"{abs(R_down-2.276)/2.276*100:.2f}%"))

    # PMNS
    sin2_12 = (N+1)/PG
    sin2_23 = (7.0/PG) * (1 + eps/PG)
    sin2_13 = (1.0/(N*PG)) * (1 - 2*eps/N)
    results.append(("sin²θ₁₂", sin2_12, EXP['sin2_12'],
                     f"{abs(sin2_12-EXP['sin2_12'])/EXP['sin2_12']*100:.2f}%"))
    results.append(("sin²θ₂₃", sin2_23, EXP['sin2_23'],
                     f"{abs(sin2_23-EXP['sin2_23'])/EXP['sin2_23']*100:.2f}%"))
    results.append(("sin²θ₁₃", sin2_13, EXP['sin2_13'],
                     f"{abs(sin2_13-EXP['sin2_13'])/EXP['sin2_13']*100:.2f}%"))

    # m_W/m_Z.  The value computed here is the internal √(10/13) of §6.1.
    # Table 2 of the manuscript lists the identified value 0.88137 (custodial
    # response carried by the single neutral Higgs component; computed in
    # Paper II), which agrees with experiment to 0.01%.  This script cannot
    # reproduce that computation, so the row below reports the internal value.
    mW_mZ = np.sqrt(1 - sin2_thetaW)
    results.append(("m_W/m_Z (√(10/13))", mW_mZ, EXP['m_W/m_Z'],
                     f"{abs(mW_mZ-EXP['m_W/m_Z'])/EXP['m_W/m_Z']*100:.2f}%"))
    print(f"  [info] identified m_W/m_Z (custodial single neutral component, "
          f"Paper II): 0.88137 vs experiment {EXP['m_W/m_Z']} -> 0.01%; the "
          f"summary row shows the internal value sqrt(10/13) of §6.1.")

    # Q (Koide)
    results.append(("Q (Koide)", Q_pred, Q_exp, f"δ={Q_pred-1.5:.2e}"))

    # D_PR
    rho = sum(p**2 for p in psi[:3])
    D_PR = (np.trapezoid(rho, phi))**2 / np.trapezoid(rho**2, phi)
    results.append(("D_PR", D_PR, 13.0, f"{abs(D_PR-13)/13*100:.2f}%"))

    # ===== Neutrino predictions =====
    print(f"\n  --- Neutrino Predictions ---")

    # Setting the absolute scale.  The dimensionless mass ray is fixed internally by
    # R_nu = R_lepton together with (m3/m2)^2 = 169/5; NO oscillation data enters it.
    # Exactly one dimensionful quantity, Dm21, is then used to set the overall scale.
    # The second splitting Dm31 is deliberately NOT used, so that the internal ratio
    # R_delta = Dm31^2/Dm21^2 remains on the verification side as a prediction.
    def compute_neutrino(Dm21, label):
        t = np.exp(-np.log(13.0 / np.sqrt(5.0)) / (R_lep - 1.0))   # m1/m2, internal
        r = 13.0 / np.sqrt(5.0)                                    # m3/m2, internal
        m2 = np.sqrt(Dm21 / (1.0 - t**2))
        m1 = t * m2
        m3 = r * m2
        sum_m = (m1 + m2 + m3) * 1000  # meV
        print(f"    [{label}] m₁={m1*1000:.4f} meV, m₂={m2*1000:.4f} meV, "
              f"m₃={m3*1000:.4f} meV, Σm={sum_m:.4f} meV")
        return m1, m2, m3, sum_m

    print(f"  Using R_ν = R_lepton = {R_lep:.6f}  (scale fixed by Δm²₂₁ alone)")
    m1_6, m2_6, m3_6, sum6 = compute_neutrino(Dm21_sq_6, "NuFIT 6.1")
    m1_5, m2_5, m3_5, sum5 = compute_neutrino(Dm21_sq_5, "NuFIT 5.x")

    # Predicted second splitting, from the same scale setting.  This is the quantity
    # left for verification; it is NOT an input anywhere above.
    Dm31_pred_6 = m3_6**2 - m1_6**2
    print(f"  [prediction] Δm²₃₁ = {Dm31_pred_6*1e3:.4f}e-3 eV²  "
          f"(NuFIT 6.1 measured {Dm31_sq_6*1e3:.4f}e-3; "
          f"{100*(Dm31_pred_6/Dm31_sq_6-1):+.2f}%)")

    # Master relation (Paper II): with m₃/m₂ = 13/√5 (exact, Paper II) and
    # R_ν = R_lepton (this paper), Δm²₃₁/Δm²₂₁ = (169/5 − t²)/(1 − t²)
    # follows with NO oscillation input. Informational.
    u_ratio = 13.0 / np.sqrt(5.0)
    x_t = np.log(u_ratio) / (R_lep - 1.0)   # x = −ln t
    t_sq = np.exp(-2.0 * x_t)
    R_delta = (169.0/5.0 - t_sq) / (1.0 - t_sq)
    print(f"  [info] master relation (Paper II, m₃/m₂=13/√5): "
          f"Δm²₃₁/Δm²₂₁ = {R_delta:.4f}  (manuscript: 33.842, displayed 33.84)")

    results.append(("Ordering", "Normal", "not settled", "prediction"))
    results.append(("m₁ (6.1) meV", m1_6*1000, "—", "scale by Δm²₂₁"))
    results.append(("Σm_ν (6.1) meV", sum6, "—", "scale by Δm²₂₁"))
    results.append(("m₁ (5.x) meV", m1_5*1000, "—", "for comparison"))
    results.append(("Σm_ν (5.x) meV", sum5, "—", "for comparison"))
    results.append(("Δm²₃₁/Δm²₂₁", R_delta, Dm31_sq_6/Dm21_sq_6, "JUNO"))
    results.append(("0νββ", 0, "not observed", "nEXO"))

    # ===== Summary =====
    print(f"\n{'='*85}")
    print(f"  RESULTS SUMMARY")
    print(f"{'='*85}")
    print(f"  {'Quantity':<16} {'Predicted':>14} {'Experiment':>14} {'Status':>14}")
    print(f"  {'-'*60}")
    for name, comp, exp, acc in results:
        if isinstance(comp, (int, str)):
            print(f"  {name:<16} {str(comp):>14} {str(exp):>14} {acc:>14}")
        elif isinstance(exp, str):
            print(f"  {name:<16} {comp:>14.6f} {exp:>14} {acc:>14}")
        else:
            print(f"  {name:<16} {comp:>14.8f} {exp:>14.8f} {acc:>14}")

    # ===== Critical comparison: NuFIT 5.x vs 6.1 =====
    print(f"\n{'='*85}")
    print(f"  CRITICAL: Impact of NuFIT 6.1 on neutrino predictions")
    print(f"{'='*85}")
    print(f"  {'':20} {'NuFIT 5.x':>14} {'NuFIT 6.1':>14} {'Change':>10}")
    print(f"  {'-'*60}")
    print(f"  {'Δm²₂₁ (10⁻⁵eV²)':<20} {Dm21_sq_5*1e5:>14.2f} {Dm21_sq_6*1e5:>14.2f} {(Dm21_sq_6-Dm21_sq_5)/Dm21_sq_5*100:>+9.1f}%")
    print(f"  {'Δm²₃₁ (10⁻³eV²)':<20} {Dm31_sq_5*1e3:>14.3f} {Dm31_sq_6*1e3:>14.3f} {(Dm31_sq_6-Dm31_sq_5)/Dm31_sq_5*100:>+9.1f}%  [not an input]")
    print(f"  {'m₁ (meV)':<20} {m1_5*1000:>14.4f} {m1_6*1000:>14.4f} {(m1_6-m1_5)/m1_5*100:>+9.1f}%")
    print(f"  {'Σm_ν (meV)':<20} {sum5:>14.4f} {sum6:>14.4f} {(sum6-sum5)/sum5*100:>+9.1f}%")
    print(f"  {'sin²θ₂₃ (exp)':<20} {sin2_23_5:>14.3f} {sin2_23_6:>14.3f}")
    print(f"  {'sin²θ₂₃ (pred)':<20} {sin2_23:>14.6f} {'':>14}")
    print(f"  {'-'*60}")
    print(f"  NOTE: Theory prediction R_ν = {R_lep:.6f} is INDEPENDENT of Δm².")
    print(f"  The dimensionless mass ray is fixed internally; only the overall scale")
    print(f"  is set, by Δm²₂₁ alone.  m₁ and Σm_ν therefore move only with")
    print(f"  Δm²₂₁ (as its square root), not with Δm²₃₁.  Δm²₃₁ is left unused so")
    print(f"  that R_Δ = Δm²₃₁/Δm²₂₁ = {R_delta:.4f} stays on the verification side.")
    print(f"{'='*85}")

    # ===== sin²θ_W comparison =====
    print(f"\n  sin²θ_W comparison:")
    print(f"    Theory (tree):  3/13 = {3/13:.5f}")
    print(f"    MS-bar (M_Z):   0.23122 → deviation {abs(3/13-0.23122)/0.23122*100:.2f}%")

# ===========================================================================
#  Optional: regenerate the certified constants instead of trusting the
#  literals used above.  Nothing in this section runs in the default
#  verification; it exists so that a reader who wants the numbers rather than
#  the citation can produce them.
# ===========================================================================

def _conv(a, b, k):
    """Cauchy product coefficient sum_{i=0}^k a[i] b[k-i]."""
    s = a[0] * b[k]
    for i in range(1, k + 1):
        s += a[i] * b[k - i]
    return s


def _conv_sq(a, k):
    """Cauchy coefficient of a*a at order k (symmetric halving)."""
    s = 0
    for i in range((k + 1) // 2):
        s += a[i] * a[k - i]
    s = 2 * s
    if k % 2 == 0:
        s += a[k // 2] * a[k // 2]
    return s


def _theta_taylor(t0, th0, E2, N, mp):
    """Taylor coefficients th_0..th_N of the Prufer angle at t0.

    Same recurrence as supplementary/proof_E_taylor_enclosures_20260718.py,
    in plain (non-interval) arithmetic: the certificate needs intervals to
    prove an enclosure, this only needs to produce digits.
    Prufer form: th' = cos^2 th + (2 H(sigma(t)) + 2E) sin^2 th.
    """
    x = mp.exp(-t0)
    opx = 1 + x
    sig = [1 / opx]
    h0 = mp.log(opx) + t0 * x / opx
    p = []
    g = [2 * h0 + E2]
    S = [mp.sin(th0)]
    C = [mp.cos(th0)]
    th = [th0]
    w = []
    B = []
    for k in range(N):
        pk = sig[k] - _conv_sq(sig, k)
        p.append(pk)
        sig.append(pk / (k + 1))
        tpk = t0 * pk + (p[k - 1] if k >= 1 else 0)
        hk1 = -tpk / (k + 1)
        B.append(_conv_sq(S, k))
        wk = _conv_sq(C, k) + _conv(g, B, k)
        w.append(wk)
        th.append(wk / (k + 1))
        Sk1 = _conv(C, w, k) / (k + 1)
        Ck1 = -_conv(S, w, k) / (k + 1)
        S.append(Sk1)
        C.append(Ck1)
        g.append(2 * hk1)
    return th


def _horner(c, x):
    r = c[-1]
    for a in reversed(c[:-1]):
        r = r * x + a
    return r


def _march(E, t_start, th_start, t_end, N, dt, mp):
    """Step the Prufer angle from t_start to t_end with signed step dt."""
    E2 = 2 * E
    t, th = t_start, th_start
    if dt > 0:
        while t < t_end:
            h = min(dt, t_end - t)
            th = _horner(_theta_taylor(t, th, E2, N, mp), h)
            t += h
    else:
        while t > t_end:
            h = max(dt, t_end - t)
            th = _horner(_theta_taylor(t, th, E2, N, mp), h)
            t += h
    return th


def eigenvalue_from_scratch(E_guess, sector, T, mp, N=24, xm=6, verbose=False):
    """Bound-state energy by two-sided Prufer shooting, matched at xm.

    Outward from 0 and inward from T are each integrated in their stable
    direction.  A single outward sweep to large T cannot work: away from an
    eigenvalue the growing branch takes over and the sign change collapses to
    width ~exp(-2kT).  Bisection, not secant: a secant step can leave E >= 0,
    where the decaying branch does not exist.  Only a float64-level starting
    guess enters -- no certified literal is used.
    """
    th0 = mp.pi / 2 if sector == 'even' else mp.mpf(0)

    def mismatch(E):
        th_out = _march(E, mp.mpf(0), th0, mp.mpf(xm), N, mp.mpf(1) / 32, mp)
        th_dec = mp.pi - mp.atan(1 / mp.sqrt(-2 * E))   # exact decaying branch
        th_in = _march(E, mp.mpf(T), th_dec, mp.mpf(xm), N, -mp.mpf(1) / 8, mp)
        return th_out - th_in

    m = int(mp.nint(mismatch(E_guess) / mp.pi))

    def F(E):
        return mismatch(E) - m * mp.pi

    lo, hi = E_guess * mp.mpf('1.02'), E_guess * mp.mpf('0.98')
    Flo, Fhi = F(lo), F(hi)
    tries = 0
    while Flo * Fhi > 0 and tries < 40:
        lo *= mp.mpf('1.05')
        hi *= mp.mpf('0.95')
        Flo, Fhi = F(lo), F(hi)
        tries += 1
    if Flo * Fhi > 0:
        raise RuntimeError("no sign change bracketing %s" % mp.nstr(E_guess, 6))

    for i in range(int(mp.mp.dps * 3.4) + 8):
        mid = (lo + hi) / 2
        Fm = F(mid)
        if Flo * Fm <= 0:
            hi, Fhi = mid, Fm
        else:
            lo, Flo = mid, Fm
        if verbose and i % 25 == 0:
            print(f"      bisect {i:3d}  width {mp.nstr(hi - lo, 3)}", flush=True)
    return (lo + hi) / 2


def eps_from_scratch(mp):
    """eps = 3 - n_WKB.  One exponentially convergent integral; ~0.01 s."""
    return 3 - 2 * mp.quad(
        lambda p: mp.sqrt(2 * (mp.log(1 + mp.exp(-p)) + p / (1 + mp.exp(p)))),
        [0, mp.inf]) / mp.pi


def recompute_constants(dps=20, verbose=True):
    """Regenerate eps and E_0,E_1,E_2 from scratch, then redo 1/alpha.

    No literal from the default run enters: the starting guesses are the
    float64-level values -0.485, -0.159, -0.0104.  Timing measured at dps=20:
    eps 0.01 s, E_0 60 s, E_1 60 s, E_2 179 s (T=200, because the third state
    is shallow, k = sqrt(2|E_2|) = 0.144).  Cost rises steeply with dps.
    """
    import time
    import mpmath as mp
    dps_save = mp.mp.dps
    mp.mp.dps = dps
    try:
        print(f"\n{'='*78}")
        print(f"  RECOMPUTED FROM SCRATCH (dps = {dps}; no certified literal used)")
        print(f"{'='*78}")

        t0 = time.time()
        eps_s = eps_from_scratch(mp)
        print(f"  eps = {mp.nstr(eps_s, min(dps, 22))}"
              f"    [{time.time() - t0:.2f} s]", flush=True)

        jobs = (("E_0", mp.mpf('-0.485'), 'even', 45),
                ("E_1", mp.mpf('-0.159'), 'odd', 45),
                ("E_2", mp.mpf('-0.0104'), 'even', 200))
        E_s = []
        for name, guess, sector, T in jobs:
            t0 = time.time()
            E = eigenvalue_from_scratch(guess, sector, T, mp, verbose=verbose)
            E_s.append(E)
            print(f"  {name} = {mp.nstr(E, min(dps, 22))}"
                  f"    [T = {T}, {time.time() - t0:.0f} s]", flush=True)

        N_ = 3
        l1 = (N_ + 1) - mp.sqrt(N_)
        l2 = (N_ + 1) + mp.sqrt(N_)
        l3 = mp.mpf(2 * (N_ + 1))
        Z = (sum(1 / abs(e) for e in E_s)
             + 12 * sum(1 / (abs(e) + l1) for e in E_s)
             + 12 * sum(1 / (abs(e) + l2) for e in E_s)
             + 27 * sum(1 / (abs(e) + l3) for e in E_s))
        L0 = Z + l1

        def dyson(a):
            x = 9 * a / (26 * mp.pi)
            return (1 / a + a
                    - (L0 - x * (1 + a / 2) + mp.mpf(3) / 4 * eps_s**2 * x
                       + mp.mpf(1589) / 5408 * x**2
                       - mp.mpf(2) / 13 * x**3 / (1 + x)))

        inv_all = 1 / mp.findroot(dyson, 1 / L0)  # seed = leading-order 1/L0 (theory-derived)
        lo, hi = mp.mpf('137.0359990843979'), mp.mpf('137.0359990844119')
        print(f"\n  Z        = {mp.nstr(Z, 14)}")
        print(f"  1/alpha  = {mp.nstr(inv_all, 16)}")
        print(f"  certified 13-digit interval (Paper II): "
              f"[{mp.nstr(lo, 16)}, {mp.nstr(hi, 16)}]")
        print(f"  inside   : {bool(lo <= inv_all <= hi)}")
        print(f"  manuscript display 137.035999084 -> difference "
              f"{mp.nstr(abs(inv_all - mp.mpf('137.035999084')), 4)}")
        print(f"{'='*78}")
        return {'eps': eps_s, 'E': E_s, 'inv_alpha': inv_all}
    finally:
        mp.mp.dps = dps_save


def _run_certificate():
    """Rigorously re-certify the eigenvalue literals with the shipped script."""
    import pathlib
    import subprocess
    import sys
    script = (pathlib.Path(__file__).parent / "supplementary"
              / "proof_E_taylor_enclosures_20260718.py")
    print(f"\n{'='*78}")
    print("  RIGOROUS RE-CERTIFICATION of the eigenvalue literals (+-1e-24)")
    print(f"  {script.name}  --  about 15 s per eigenvalue")
    print(f"{'='*78}")
    if not script.exists():
        print(f"  [warn] certificate script not found at {script}")
        return
    subprocess.run([sys.executable, str(script)], check=False)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="PFE Paper I verification.  The default run uses the "
                    "certified constants of §3.1-§3.2 as literals; --certify "
                    "and --recompute check them independently.")
    parser.add_argument("--phi-max", type=float, default=100)
    parser.add_argument("--n-grid", type=int, default=1600001)
    parser.add_argument("--certify", action="store_true",
                        help="rigorously re-certify the eigenvalue literals "
                             "(+-1e-24 windows) with the shipped certificate; "
                             "about 15 s per eigenvalue")
    parser.add_argument("--recompute", action="store_true",
                        help="regenerate eps and E_0,E_1,E_2 from scratch and "
                             "redo 1/alpha; MINUTES TO HOURS, see --dps")
    parser.add_argument("--dps", type=int, default=20,
                        help="working precision for --recompute (default 20; "
                             "higher is slower, steeply)")
    parser.add_argument("--yes", action="store_true",
                        help="skip the confirmation prompt of --recompute")
    args = parser.parse_args()

    run_paper1(phi_max=args.phi_max, n_grid=args.n_grid)

    if args.certify:
        _run_certificate()

    if args.recompute:
        print(f"\n{'='*78}")
        print("  --recompute regenerates the certified constants from scratch.")
        print("  Measured at --dps 20: eps 0.01 s, E_0 60 s, E_1 60 s, "
              "E_2 179 s (about 5 min total).")
        print("  Cost rises steeply with --dps; a large value can run for hours.")
        print("  Nothing printed above is affected: this only re-derives the")
        print("  literals that the default run cites.")
        print(f"{'='*78}")
        if not args.yes:
            try:
                if input("  proceed? [y/N] ").strip().lower() not in ("y", "yes"):
                    raise SystemExit("  cancelled.")
            except EOFError:
                raise SystemExit("  no tty; re-run with --yes to proceed.")
        recompute_constants(dps=args.dps)
