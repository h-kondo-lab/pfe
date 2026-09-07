# Physics from Existence I: The Equation — The Binary-Entropy Equation $V = -H$ and the Structure of Standard-Model Parameters

**Hidekazu Kondo** — Kondo Research Institute, Tokyo, Japan (April 2026)

---

## Abstract

No theory has simultaneously derived all 26 free parameters of the Standard Model. In this paper, extending the quantum indefiniteness of existence to the existence of the universe itself, we posit three axioms about existence — existence is bivalent (n ∈ {0,1}, i.e. n² = n), existence is uncertain (σ ≠ 1), and non-existence is forbidden (σ ≠ 0) — and show that they uniquely determine the single equation V = −H as an algebraic identity, and that the mathematical structure derived from its exact algebraic expansion corresponds to the structures of the Standard Model. As a result, the 26 physical constants of the Standard Model are reproduced with no adjustable parameters. Of these, the nineteen derived in this paper agree with the established measured values to 0.0006%–2.6%.

Derived from the axioms of existence, V = −H is the canonical free energy of a single binary degree of freedom n ∈ {0,1}, determined as the Legendre transform of the Fermi-Dirac distribution. From this potential and the kinetic term of the canonical normalisation (g = 1, §3.1), the Schrödinger equation is derived as a transfer matrix, and we show that it has exactly three bound states (discrete solutions, N = 3), corresponding to the three generations of fermions. The spatial dimension d = 3 follows from two conditions — preservation of the well and propagating gravity — and the Lie-algebra structure unfolding on the finite geometry constructed from N = 3 corresponds to the structure of the gauge group SU(3)×SU(2)×U(1) (the coincidence of the two numbers, N = d = 3, is an independent consistency).

Twenty-six physical constants are derived from this structure, including 1/α = 137.035999084 (all orders, matching CODATA 2018 to all displayed digits; the leading term is 137.0368 at 0.0006%), α_s = 0.1179 (0.07%), sin²θ_W = 0.23122 (<0.01%), the charged-lepton mass ratio R = 1.5293 (0.006%), and sinθ_C = 0.2245 (0.10%). The results further predict the absolute neutrino masses (Σm_ν = 59.37 meV, with no oscillation input) and strictly zero neutrinoless double-beta decay, among other testable predictions. These quantities are not independent input parameters — they are all computed from the same mathematical structure fixed by V = −H.

All of these numerical values are reproducible with the accompanying verification code.

# 1. Introduction

The Standard Model is the most successful theory of particle physics. It describes all known particle interactions and has withstood extensive experimental tests, including the discovery of the Higgs boson [1, 2] (2012). At present, however, the theory contains 26 free parameters[^sm26] — fermion masses (Yukawa couplings), mixing angles, gauge couplings, the Higgs self-coupling, and CP phases. Their values can only be put in from measurement; why these physical constants take the values they do is not explained within the Standard Model. The generation number N = 3 is likewise an unexplained input.

Despite many attempts — grand unified theories (GUTs), supersymmetry (SUSY), string theory — no theory has simultaneously derived all 26 parameters with zero free parameters [3]. Indeed, even for the single constant α, the fine-structure constant, attempts at a first-principles derivation [35, 36, 37] have a long history, yet no theory has derived it from an equation with a unifying mechanism.

Approaches that derive physics from information theory have also been attempted. Wheeler's "it from bit" [4] proposed the vision that physical reality arises from information, but had no mathematical framework. Jaynes [5] derived statistical mechanics from entropy maximisation, but did not reach quantum field theory. Verlinde [7] derived gravity from entropic reasoning, but gravity alone. There are also Frieden's Fisher-information approach [6], Caticha's information geometry [8], and the information-theoretic axiomatisations of Hardy [9] and Chiribella et al. [10]. The prior work closest in spirit is the Connes–Chamseddine spectral Standard Model [39, 40, 41], which derives the structure of the gauge group and the Higgs sector from noncommutative geometry (a related algebraic reconstruction is Furey's division-algebra approach [42]). All these approaches share a common limitation: they can derive the *form* of physical equations, but not the *values* of the constants appearing in them.

Now, "Why do I exist? Why does the universe exist?" — who has not, at some point, let their thoughts dwell on such questions? Traditionally, they have been regarded as purely philosophical questions, beyond the scope of physics. In fact, however, quantum theory has already established a decisive fact about the very notion of existence. As is well known, the existence of a quantum object — whether a mode is occupied — is in general not definite, and is given only as a probability. In other words, at the microscopic level, existence is not a settled given. If so, might the "existence" of the universe itself perhaps behave in the same way? What happens if that question is written down as an equation? That is the starting point of this paper.

This paper first formalises this question about the existence of the universe as three axioms — Existence (n ∈ {0,1}), Uncertainty (σ ≠ 1), and Non-existence is forbidden (σ ≠ 0). Together A2 and A3 determine σ ∈ (0,1), a partition function Z = 1 + e^φ, and, via the Legendre transform of the Fermi–Dirac distribution, the canonical free energy V = −H(σ(φ)) as an algebraic identity, where H is the binary Shannon entropy, σ(φ) is the Fermi–Dirac distribution (sigmoid function), and φ = logit(σ) is the natural parameter [14, 15].

Combining this potential V = −H(σ(φ)) with the kinetic term of the canonical normalisation (g = 1, §3.1), the Schrödinger equation is derived as a transfer matrix. Carrying the computation further, this Schrödinger equation has exactly three bound states (N = 3), and a structure corresponding to the three generations of fermions appears. At the same time, the spatial dimension d = 3 follows from two conditions — preservation of the well and propagating gravity. The spacetime dimension is D = d + 1 = 4, supported as consistency checks under stated premises by the gravitational degree-of-freedom count, by ℂP², and by the connection of the signature (3,1) to the M₂(ℂ)_sa structure (§3.4, §9). From N = 3 the finite geometry PG(2, F₃) is fixed (the coincidence of the two numbers, N = d = 3, is an independent consistency), and the Lie-algebra structure unfolding from its 9+3+1 partition corresponds uniquely to the structure of the gauge group SU(3)×SU(2)×U(1).

From the representations of this gauge group, a structure yielding all the particles — 48 fermions, 12 gauge bosons, and 1 Higgs — is derived (§4). From this structure, 19 parameters are derived — including 1/α_em = 137.0368 (0.0006%; 137.035999084 at all orders — proved in Paper II), sin²θ_W = 3/13 (0.2%), the charged-lepton mass ratio R = 1.5293 (0.006%), Koide's Q ≈ 3/2 (δ = −1.5 × 10⁻⁵), and CKM/PMNS mixing angles — all at 0.0006%–2.6% accuracy with zero adjustable parameters. The seven derived in Paper II (|V_cb|, |V_ub|, δ_CKM, δ_PMNS, m₁/m_e, m₃/m₂, Δm²₃₁/Δm²₂₁) are also displayed in the tables (§10). The associated Paper II reconstruction additionally gives the matching-point condition θ̄(μ₀) = 0.

Thus the algebraic expansion of "existence" alone yields numbers that agree with the Standard-Model constants to high accuracy; beyond this numerical agreement, the mathematical structure derived from V = −H also corresponds to the structural problems of the Standard Model themselves. For example, the role played in the Standard Model by μ² > 0 — moving the vacuum from the symmetric point to a finite field value — is played here by the shape of the well of V (§6); the bounded internal potential (V ∈ [−ln 2, 0]) supplies no instability channel, so neither the metastability problem nor the hierarchy problem arises within this construction (the correspondence to 4D vacuum stability and the hierarchy problem requires the reconstruction matching). V = −H is uniquely determined as the canonical free energy of n ∈ {0,1} and fixes the shape of the well. That the Higgs is unique (weak-doublet scalar multiplicity 1) follows not from the uniqueness of V but from the classification of the minimal operator algebra, whose proof is given in Paper II (§6.2).

What this paper shows is that, when the mathematical structure generated by V = −H is placed in correspondence with the structure of the Standard Model, the derived numbers agree with experiment — this paper itself does not carry out a rigorous first-principles derivation of the Standard Model; that derivation is treated in Paper II as a theorem under the stated realisation conditions (§9.2). The method of computing the physical constants from the derived equations alone is given in this paper, and can be reproduced independently with the accompanying verification code.

Comparing the derived numbers with the experimental observables requires making two things explicit. First, dimensionful quantities require a correspondence to humanly chosen units of measurement — setting the absolute scale (in this paper the electron mass m_e; §7.3, §10). Second, dimensionless quantities require a correspondence that identifies which part of the structure is which physical structure or observable — for example, reading and identifying sin²θ_W as the fraction of the 13 directions occupied by the 3 weak ones. In other words, one must identify which known physical structure each derived mathematical structure corresponds to.

Each of these identifications is provisional at the stage of its declaration. The text declares each identification where it enters (§3 generations, §4.2 gauge sectors, §6 Higgs, §7 masses, §8 mixings) and argues under that hypothesis; the agreement of the full set of derived values with experiment is what corroborates the identifications. This paper states all of these correspondences (the mapping rules onto the observables) explicitly and computes the values. That each mapping rule is not a free parameter — the proof of its uniqueness — together with the operator-level derivations of the individual constructions, would require an enormous body of argument and is therefore deferred to Paper II [20]. Such claims are marked with a superscript ^II^ in the text (indicating rigorisation in Paper II) and are not restated individually.

The paper is organized as follows. §2 derives V = −H from three axioms. §3 establishes the Schrödinger equation, proves the bound-state count N = 3 as a mathematical theorem, and derives the spatial dimension d = 3 as the unique common solution of two independent realisation-consistency conditions. §4 constructs the gauge structure from PG(2,3). §5 derives the gauge coupling constants. §6 derives the Higgs mass and shows that the vacuum is pinned at an interior point with the internal potential bounded (no instability channel). §7 identifies the bound states with the three generations of charged leptons and derives the fermion mass ratios. §8 derives the inter-generation mixing and the neutrino properties. §9 confirms the spacetime dimension D = 4 and signature (3,1) by consistency checks with stated premises, and establishes the tensor-product structure ℋ_int(φ) ⊗ ℋ_space(x) from the independence of the construction. §10 summarises all results and presents falsifiable predictions. §11 gives the conclusion.

[^sm26]: 19 in the minimal SM with massless neutrinos; 26 with three neutrino masses and four PMNS parameters. Some references count 25, absorbing v_EW into the energy scale.

# 2. From Existence to the Equation

This chapter formalises three axioms about existence mathematically and shows that the fundamental equation V = −H is derived as an algebraic identity. From this equation, the geometric representation corresponding to the internal-symmetry structure of the Standard Model is uniquely determined^II^.

## 2.1. The Axioms

Paring the question of the introduction — might the existence of the universe itself, like quantum existence, be undetermined? — down to a form that can become an equation leaves three short sentences. What the three sentences speak of is a single proposition, prior to any property: that something exists. Which particles, what space, what laws — nothing of that is said yet. These three sentences are the entire input of the theory; no physical constant, no spacetime, and no postulate of quantum mechanics enters here.

**A1 (Existence):** Something exists. "Exists" is a proposition; it is either true or false — by the law of excluded middle, there is no third value. Writing the state of existence as n, the state space is n ∈ {0,1}. A1 defines the state space but does not determine the **value** of n.

What is bivalent is the *proposition* of existence, not a claim that the world is deterministic. To ask about the existence of something is, before asking about any of its properties, the single question: is it, or is it not? A1 fixes only the minimal structure this question carries — a one-bit state space. Everything that follows is built on this one bit.

**A2 (Uncertainty):** Existence is uncertain. Whatever exists does not exist with certainty: the occupation probability σ = ⟨n⟩ satisfies σ ≠ 1. A2 excludes the upper endpoint (certain existence) but says nothing about the lower endpoint.

If existence were certain (σ = 1), it would be an immovable background: nothing left to ask, nothing left to happen. As quantum theory teaches, microscopic existence is given only as a probability — A2 lifts this fact to the existence of the universe itself. Existence is not a given but something that *occurs*, carrying uncertainty. This one sentence is what later generates the statistical ensemble and the dynamics.

**A3 (Prohibition of non-existence):** Non-existence is not a state (σ ≠ 0). A way of being that consists in not-being is a contradiction.

"Nothing" is not a way for something to be. If σ = 0 were a realisation state, non-existence would be *there*, carrying data — the absence of being would be. A3 is the hardest-working of the three axioms: here it only removes the lower endpoint, but later it selects the coordinate by forbidding model data at the point of non-existence (§3.1) and selects the spatial dimension as the prohibition of boundary conditions (§3.4).

The three axioms are mutually independent — A1 removes the third truth value, A2 the upper endpoint, A3 the lower, each removing what the other two do not. And the three are exhaustive: together they leave only the open interval 0 < σ < 1 (§2.2), and no further axiom enters anywhere in this paper. The axioms are these three sentences; how they are realised mathematically is constructed in minimal form in the following sections, and the proof that no freedom of choice remains in that realisation is given^II^.

## 2.2. Bivalence and the Partition Function

Starting from the state space n ∈ {0,1} established by A1, we construct the partition function.

Since n takes only the values 0 and 1,

$$n^2 = n$$

follows — because the solutions of n(n − 1) = 0 are exactly {0,1}. However, A1 alone still admits the endpoint pure states n = 0 and n = 1, and no dynamics arises. A2 (σ ≠ 1) excludes certain existence; A3 (σ ≠ 0) excludes non-existence. Together they determine σ ∈ (0,1): the system is described as a canonical ensemble of a binary variable — an exponential family with mean parameter σ — neither certainly present nor certainly absent. The tilt weight e^{nφ} has its tilt coefficient fixed to 1 (the scale of φ; §3.1)^II^. n² = n truncates the sum to two terms, and the partition function is

$$Z = \sum_{n=0}^{1} e^{n\varphi} = 1 + e^{\varphi}$$

where φ = ln[σ/(1−σ)] (the logit transform) is the canonical coordinate of the exponential family.

Were n able to take values 0, 1, 2, … (Bose statistics), the sum would be an infinite series. As a consequence of the law of excluded middle, n² = n forbids this and restricts the occupation number of each binary mode to {0,1} — the exclusion of occupations 2, 3, … of that same mode, the Fermi–Dirac occupation statistics. This one-mode relation alone does not determine the exchange law between distinct modes — exchange antisymmetry (spin statistics) is fixed uniquely by the many-body realisation (§4.3)^II^.

The expectation σ = ⟨n⟩ = ∂lnZ/∂φ coincides with the Fermi–Dirac distribution, and its variance is the source of all subsequent structure:

$$\sigma = \frac{1}{1+e^{-\varphi}}, \qquad g_F = \sigma(1-\sigma) = \mathrm{Var}(n) > 0$$

These are not assumed but derived, from the combination of A1 (existence: n² = n), A2 (uncertainty: σ ≠ 1), and A3 (non-existence forbidden: σ ≠ 0). g_F is the gap between the binary observable (n² = n) and its uncertain expectation (σ² ≠ σ), and vanishes only at the pure values σ ∈ {0,1}.

## 2.3. The Effective Potential V = −H

§2.2 established the partition function Z = 1+e^φ. The next step is to construct an effective potential from it. The naive choice −ln Z diverges as φ → ∞, so a Legendre transform to the canonical free energy is required.

The logarithm W = ln Z of the partition function generates the classical field σ = dW/dφ, yielding the effective potential:

$$V(\varphi) = \varphi\sigma - W(\varphi) = \sigma\ln\sigma + (1-\sigma)\ln(1-\sigma) = -H(\sigma)$$

where H(p) = −p ln p − (1−p) ln(1−p) is the binary Shannon entropy [12].

That is, written out with its arguments, V(φ) = −H(σ(φ)): at each point of the canonical coordinate φ, the value of the potential equals minus the binary entropy of the occupation probability σ(φ). This is an algebraic identity — the canonical free energy of a single binary degree of freedom n ∈ {0,1} equals the negative of the Shannon entropy — and there is no freedom of choice. The bare form V = −H, used throughout, abbreviates this identity of functions.

The potential

$$V = -H$$

is thereby uniquely determined by axioms A1–A3 together with the exponential-family realisation above.

From this single equation, a mathematical structure corresponding to the structure of the Standard Model emerges, and values agreeing with physical constants to high precision are derived (§3–§10).

The derivation chain of this chapter is summarised as:

$$A1\text{–}A3 \;\Rightarrow\; n^2 = n,\; 0<\sigma<1 \;\Rightarrow\; Z = 1+e^{\varphi} \;\xrightarrow{\text{Legendre}}\; V = -H$$

# 3. Derivation of the Schrödinger Equation and N = d = 3

V = −H is a potential; by itself it does not yet say how many discrete states exist. This chapter first derives the eigenvalue equation — the Schrödinger equation — from the same exponential-family structure as the potential (§3.1). The spectral question — how many bound states the well holds — then becomes well-posed, and its answer is N = 3 (§3.2–3.3). At the end of §3.3 the three bound states are identified with the three generations of fermions — the one physical identification this chapter makes. Only under that identification does the second question arise: the particles the bound states represent exist in real space, so what determines the dimension of that space? Realisation consistency answers it, and d = 3 follows (§3.4). The derivation of d = 3 makes no use of the value N = 3 — the two derivations share the operator and nothing else. The coincidence N = d = 3 becomes the input of §4.
## 3.1. Derivation of the Schrödinger Equation

To write the eigenvalue equation, a kinetic term is needed in addition to the potential of §2.3, and it follows from the same structure as the potential. Dynamics enters by reading the canonical ensemble of §2 as sequential transfer: a family of transfer operators that, at each step, exponentiate the quadratic variation of the displacement (the kinetic term) together with V = −H (the action weight); its continuum limit yields the eigenvalue equation at the end of this section^II^. Two data remain: the coordinate carrying the motion, and one dimensionless normalisation. Both are fixed uniquely; the uniqueness is proved^II^.

First, the coordinate. The exponential family of a binary variable has two natural coordinates:

$$u = 2\arcsin\sqrt{\sigma}\quad(\text{Fisher metric flat}),\qquad \phi = \mathrm{logit}(\sigma)\quad(\text{canonical parameter in the exponent})$$

A3 decides between them. The coordinate u places the forbidden point σ = 0 at the finite endpoint u = 0, so dynamics on u would need boundary data exactly at the point of non-existence — non-existence would become a repository of model-defining data. The coordinate φ sends σ = 0 to φ → −∞: no boundary data attaches to non-existence, and under A3 none is admissible. Among the admissible coordinates, φ is the one on which the tilting action of A1 (the tilting weight e^{nφ}, §2.2) is linear — unique up to affine maps — and the binary indicator n ∈ {0,1} fixes its scale. Free motion on φ therefore has a constant-coefficient quadratic Lagrangian; position-dependent coefficients are excluded by displacement homogeneity^II^.

Second, the normalisation. Writing the kinetic coefficient as a²/2 and the action weight as b, the apparent two-parameter family reduces to a single dimensionless ratio:

$$\frac{a^2}{2}\Bigl(\frac{d}{d\phi}\Bigr)^{\!2},\quad b\,V \qquad\longrightarrow\qquad g = \frac{b}{a^2}$$

The reason is that the binary indicator n ∈ {0,1} of A1 fixes the unit of φ, the Legendre identity fixes the unit of V = −H, and rescaling the transfer step moves a² and b in the same proportion. Measuring the quadratic variation per step and the action weight as quantities of one and the same transfer process sets g = 1; the Landauer bit-normalisation (which takes as its unit the minimal work of erasing one bit) gives the same value independently. The uniqueness of both fixings is proved^II^.

The bound-state count established below (§3.3) does not depend on this last choice: it is the same across the whole family −½ d²/dφ² − g H(σ(φ)) for g ∈ [0.80, 1.50]. The number of negative eigenvalues is non-decreasing in g (since H(σ) ≥ 0, raising g lowers the quadratic form uniformly), and the count of three at both endpoints g = 4/5 and g = 3/2 is certified by an accompanying interval-arithmetic certificate — so N = 3 holds throughout the interval. An accompanying script further checks numerically that the count is stable across box sizes L = 60/100/160.

The ε-dependent corrections of §5–§8, by contrast, do rely on g = 1, since ε = N − √g·n_WKB(1) changes sign inside the same window (+0.51 at g = 0.80, −0.41 at g = 1.50); at the lower edge the action count 2.49 dips below the Maslov threshold 2.5 while the exact count is unchanged — the threshold rule of §3.2 is a heuristic only.

Combining V = −H as the potential with the kinetic term of the canonical normalisation (g = 1):

$$\left[-\frac{1}{2}\frac{d^2}{d\phi^2} - H(\sigma(\phi))\right]\Psi = E\,\Psi$$

This is nothing other than a Schrödinger equation: with only the three axioms of existence and the canonical normalisation g = 1 above as input, the eigenvalue equation of quantum mechanics emerges as the continuum limit of the transfer matrix of the canonical action.

The question posed at the head of this chapter — how many discrete states the well holds — is now well-posed, and §3.2–3.3 answer it. Whatever that count turns out to be, it is a consequence of this construction — the coordinate φ fixed uniquely above together with a flat kinetic term. Writing the flat kinetic term in a different coordinate yields a different operator and a different spectrum, so the uniqueness of φ is essential.

The derivation chain so far is summarised as:

$$A1\text{--}A3 \;\Rightarrow\; n^2 = n,\; 0<\sigma<1 \;\Rightarrow\; Z = 1+e^{\varphi} \;\xrightarrow{\text{Legendre}}\; V = -H \;\xrightarrow{\text{transfer matrix}}\; \left[-\tfrac{1}{2}\partial_\phi^2 + V\right]\Psi = E\Psi$$

## 3.2. Three Bound States

Why do fermions come in exactly three generations — why τ, μ, e, and why not a fourth? This and the following subsection prove rigorously that the potential V = −H has exactly three bound states and no fourth (hereafter N = 3), and identify them with the three generations of fermions. The identification is corroborated by the agreement of the mass ratios and mixing angles in §7.

The potential V = −H is a well. It is deepest at the centre (φ = 0) where V = −ln 2, returning to zero at both ends φ → ±∞. In quantum mechanics, such a well confines discrete energy levels — by the same principle that the hydrogen atom confines electrons to discrete orbitals. In general, how many levels a well holds is set by two free parameters, its depth and its width — a deep, wide well holds many levels, a shallow, narrow well holds few.

The well of V = −H, however, has no such free parameters. Its depth is fixed as the number ln 2 ≈ 0.693 directly by the functional form of the Shannon entropy itself (§2.3), and the coefficient of the kinetic term that plays the role of width is likewise fixed to the unique value g = 1 (the uniqueness theorem of §3.1). In an ordinary well problem, depth and width are free to choose; here neither is — both are completely fixed by the form of the equations alone.

How many levels, then, does this uniquely fixed well hold? That can only be found by actually solving it and counting.

![The potential V = −H(σ(φ)) and the probability densities |ψ_n|² of its three bound states (left)](img/fig_potential.png)

*Figure 1. The potential V = −H(σ(φ)) and the probability densities |ψ_n|² of its three bound states (left). More localised states correspond to heavier particles (right). No fourth negative eigenvalue: the prediction that no fourth sequential chiral generation exists.*

We denote the number of levels confined in this well by N. A semiclassical (WKB) approximation estimates how many half-wavelengths fit inside the well. Each integer half-wavelength corresponds to one bound state:

$$n_\mathrm{WKB} = \frac{1}{\pi}\int_{-\infty}^{\infty}\sqrt{2H(\sigma(\phi))}\,d\phi = 2.781$$

In the standard convention including the Maslov correction, the threshold for the k-th state is n_WKB > k − 1/2. The value n_WKB = 2.781 exceeds the third-state threshold 2.5 by 0.281, while falling short of the fourth-state threshold 3.5 by 0.719. The WKB estimate therefore suggests N = 3. We define ε ≡ N − n_WKB = 3 − 2.781 = 0.219: the gap between the exact count and the semiclassical action count.[^epsenc] That ε > 0 says the third state is barely bound.[^epsconv] This ε is the source of every correction exponent in this paper: it controls the NLO corrections to the couplings and mixing angles of §5 and §8.
[^epsenc]: The rigorous interval enclosure 0.219188 is certified by the accompanying certificate script (rigorous Arb-ball integration with an analytic tail bound); the certified interval lies inside the certified interval of Paper II.

[^epsconv]: Measured against the Maslov-corrected threshold, the same marginality reads 2.781 − 2.5 = 0.281. We use ε because it compares two unambiguously defined quantities — the exact count and the action integral — whereas the threshold rule is an asymptotic heuristic; the g = 0.80 example of §3.1 shows it can fail at finite depth.

Nothing in this count is adjustable. The three ingredients of the construction are each unique: the coordinate and the normalisation are the unique ones of §3.1, and Shannon entropy is the unique entropy satisfying the Khinchin axioms, which include strong additivity (Khinchin's theorem [16, 17]). Under this construction, N = 3 is rigorously proved by interval arithmetic (§3.3).

## 3.3. Rigorous Proof

**Theorem (N = 3).** The Schrödinger equation of §3.1 has exactly three bound states, and their eigenvalues are rigorously enclosed by interval arithmetic:

$$\begin{aligned}
E_0 &= -0.4850022531029642374452777658 \pm 10^{-24},\\
E_1 &= -0.1590522216188561779934190 \pm 10^{-24},\\
E_2 &= -0.010423393062109459327616 \pm 10^{-22}
\end{aligned}$$

No fourth negative eigenvalue exists; the essential spectrum begins at 0 (the small positive values seen on finite grids are discretisation artifacts of the box).

*Proof (computer-assisted, interval arithmetic).* The proof has three independent layers.

*Counting.* By parity the problem splits into a Neumann problem (even) and a Dirichlet problem (odd) on the half-line [0,∞). Propagating the Prüfer angle by interval enclosures (using the monotonicity of the potential), combined with analytic tail estimates, gives rigorous two-sided bounds on the angle at the trial energies. By Sturm oscillation theory the rotation of the angle counts the negative eigenvalues exactly: two in the even sector, one in the odd — N = 3.

*Enclosures.* The eigenvalue intervals above are certified by two-sided shooting with a verified high-order Taylor integrator over interval coefficients.

*Independent checks.* A separate implementation in Arb ball arithmetic re-certifies N = 3, and variational upper bounds from Rayleigh quotients of three sech-type (Pöschl–Teller-type [18]) trial functions confirm the existence of three negative eigenvalues. Three certificate scripts are included in the ancillary files. □

**Identification with generations.** Here we hypothesise that these three bound states ψ₀, ψ₁, ψ₂ are the mathematical structure of the three generations of fermions. The ordering then fixes which state is which generation: binding energy orders localisation (|E₀| > |E₁| > |E₂|, from most to least localised), and the generations carry the mass hierarchy — the most localised ψ₀ is the heaviest generation. For charged leptons, ψ₀, ψ₁, ψ₂ ↔ τ, μ, e.

At this stage the identification is provisional. In §7–§8 the mass ratios and mixing angles are derived from it and agree with experiment; that agreement is what makes it credible. If the identification is correct, the absence of a fourth negative eigenvalue is the prediction that no fourth sequential chiral generation exists.

## 3.4. Three Spatial Dimensions

The preceding subsections proved N = 3. That proof lives in the internal space (the φ coordinate), and φ is not a spatial coordinate — as §9.2 shows, the internal and spacetime factors are independent. Under the provisional identification of §3.3, however, the fermions represented by the three bound states are particles that exist in real space. What, then, determines the dimension of that space? The answer is the consistency of realisation: the internal well must be realised, unbroken, on a particle in space. This subsection shows that the two conditions of this consistency — preservation of the well (the operator must not be deformed under the radial transfer, the A3-based criterion: d = 1, 3 only) and propagating gravity (d ≥ 3) — have d = 3 as their unique common solution.

Consider the same well shape V = −H realised as a radial profile in d-dimensional space. Here σ(r) is the same well shape as σ(φ) on φ ∈ ℝ, transferred to the radial variable r ≥ 0; the difference in domain is addressed below in the comparison with full-line quantisation. Reducing the d-dimensional Schrödinger equation to a one-dimensional radial equation by standard partial-wave decomposition, an additional effective force (d−1)(d−3)/(8r²) appears in the zero-angular-momentum (spherically symmetric) radial channel:

$$V_\mathrm{eff}(r) = -H(\sigma(r)) + \frac{(d-1)(d-3)}{8r^2}.$$

The numerator (d−1)(d−3) vanishes for d = 1 and d = 3, leaving V = −H unchanged. In other dimensions, this geometric term deforms the operator:

| d | (d-1)(d-3)/8 | Effect on the operator |
|---|---|---|
| 1 | 0 | No deformation |
| 2 | -1/8 | An attractive term is added (tending to create extra states) |
| 3 | **0** | **No deformation: the well shape is preserved** |
| 4 | +3/8 | A repulsive term is added (tending to remove the marginal third state) |
| ≥5 | ≥1 | Strong repulsion |

Full-line quantisation (φ ∈ ℝ) is essentially self-adjoint at the endpoints and requires no boundary condition. Under A3 (non-existence is not a state — no additional data is admissible at the boundary; the same conclusion as the exclusion of u-type coordinates in §3.1)^II^, this canonical full-line quantisation is the reference, and a spatial realisation is required not to add an extra inverse-square term to the canonical operator. This requirement is met only when the coefficient (d−1)(d−3)/8 vanishes, i.e. only for d = 1, 3. Note that a vanishing inverse-square term does not make the full-line and radial problems identical — a half-line radial model adds endpoint structure at the origin (a choice of self-adjoint extension) and is a different model, outside the comparison fixed by this correspondence^II^.

A realisation-consistency condition excludes d = 1. The realisation of a spin-2 field (gravity) — taken as a premise, as in §9.1, with its derivation left to a subsequent paper — propagating dynamically requires d ≥ 3 (for d ≤ 2 the Weyl tensor vanishes identically, leaving no local degrees of freedom). The fermionic structure (CAR) derived from A1–A3^II^ has, in its spatial realisation, finite-dimensional spinor representations for d ≥ 3 (π₁(SO(d)) = Z₂, the double cover Spin(d)) and is consistent with this choice; but since fermionic statistics can be realised in lower dimensions as well, this is not a condition that excludes d = 1 on its own, and it is not used for the selection.

We summarise this as a theorem.

**Theorem (d = 3).** Under the above choice of full-line quantisation, d = 3 is the unique spatial dimension simultaneously satisfying: *(i)* the well operator is preserved with no extra geometric term: d = 1, 3; *(ii)* gravity propagates (Weyl tensor non-trivial): d ≥ 3. Intersection: {3}.

*Proof.* (i) The centrifugal term (d−1)(d−3)/(8r²) = 0 if and only if d = 1 or 3; (ii) the Weyl tensor vanishes identically in spacetime dimension ≤ 3, i.e. spatial d ≤ 2; propagating gravity requires d ≥ 3. □

The three spatial dimensions are also obtained independently from the fact that the traceless part of the self-adjoint sector of the complex carrier M₂(ℂ) is three-dimensional^II^; the well-preservation argument of this section is also a consistency check of that result.




# 4. Gauge Structure from the Projective Plane

The previous chapter derived two numbers: the bound-state count of the Schrödinger equation, N = 3 (§3.3), and the spatial dimension fixed by its realisation consistency, d = 3 (§3.4). The two were derived independently — yet they take the same value, N = d = 3. In this chapter we construct, from N = 3, the three-dimensional space F₃³ over the three-element field F₃ — why "3" becomes a field, and not a mere count, and why the rank 3 comes from the three bound states (d = 3 is not an input), is shown in §4.1. This chapter shows that the geometric structure arising from it has the same mathematical construction as the Standard Model gauge structure.

## 4.1. Origin of the Gauge Structure

What gauge structure arises from N = 3? The key is degeneracy. Within this construction, a continuous internal rotation symmetry is carried by degeneracy — directions with the same internal energy.

On the Standard-Model side, quark colour (red, green, blue) is the archetype: every colour has the same energy, so the continuous rotation symmetry SU(3) stands. On the derived side, the generations — the three bound states of V = −H — all have distinct energies (E₀ ≠ E₁ ≠ E₂, by the Sturm-Liouville theorem), so a preferred eigenbasis exists. The symmetry preserving this basis is discrete and cannot support a continuous gauge group.[^gencas] The gauge directions therefore cannot arise from the interchange of generations — they must come from the geometry of the space F₃³, which we now construct.

[^gencas]: The contact index b = 8/3 in the mass formula of §7.1 is derived from the finite incidence geometry and is compatible with this (its numerical agreement with the SU(3) Casimir is an N = 3-only consistency check, §7.1).

$$N = 3 \xrightarrow{\text{generation cycling}} C_3 \xrightarrow{\;\mathrm{End}\;} \mathbb{F}_3 \xrightarrow{\text{free module}} \mathbb{F}_3^3 \xrightarrow{\text{projectivisation}} \mathrm{PG}(2,\mathbb{F}_3)\;(13\text{ points}) \xrightarrow{\;\text{flag}\;} 9{+}3{+}1$$


From N = 3, the three-dimensional space F₃³ over the field F₃ = {0, 1, 2} is determined. The origin of the field structure lies in the exchange of generations. The cyclic shifts of the three bound states (1→2→3→1 and its iterates) form C₃, the unique non-trivial proper normal subgroup of S₃, and carry any state to any other by exactly one shift. The maps preserving the composition of shifts are only three — multiplication by 0, by 1, and by 2 (the endomorphism ring) — and they form the field F₃ = ℤ/3ℤ with its addition and multiplication: here 3 becomes a field, not merely a count.

$$C_3 = \{e,\ c,\ c^2\},\qquad \mathrm{End}(C_3) = \{\,c \mapsto c^k : k = 0, 1, 2\,\} \cong \mathbb{Z}/3\mathbb{Z} = \mathbb{F}_3$$

Assigning one basis direction to each of the three bound states and completing freely, admitting only F₃-linear combinations, gives F₃³ (the free F₃-module on the three states; |F₃³| = 3³ = 27 elements; the theorem-level construction is given in Paper II). The spatial dimension d = 3 is an independent theorem carrying the same number, not an input to the rank of F₃³; N = d = 3 remains a coincidence of consistency. The projective construction that follows requires a field (not merely a group); since N = 3 is prime, F₃ = ℤ/3ℤ is the unique field with N elements — indeed, the additive group of a three-element field has order 3, so it coincides with ℤ/3ℤ generated by 1, and multiplication is then determined from addition by distributivity, making the operation tables unique.

Removing the origin and taking a two-dimensional cross-section gives F₃² = {0,1,2} × {0,1,2} — a 3×3 grid. These 9 points are the starting point for the construction.

Drawing lines on the 3×3 grid, lines with the same slope are parallel and do not intersect. In projective geometry, parallel lines meet at a single "point at infinity" — the same principle by which parallel railway tracks converge to a point on the horizon in perspective drawing. The lines on F₃² have |F₃| + 1 = 4 slopes (0, 1, 2, ∞ — horizontal, right-diagonal, left-diagonal, vertical), each generating one point at infinity. These four points form the "line at infinity" L — the horizon.

Adding the 4 points at infinity to the original 9 gives 9 + 4 = 13 points. This is the projective plane PG(2, F₃) — the space obtained by defining "points" as one-dimensional subspaces and "lines" as two-dimensional subspaces through the origin of F₃³, with a symmetric structure of 13 points, 13 lines, 4 points per line, and 4 lines through each point (Figure 2).

$$|\mathrm{PG}(2,\mathbb{F}_3)| = \frac{3^3 - 1}{3 - 1} = 13$$

![The structure of the projective plane PG(2, F₃)](img/fig_pg23_v2.png)

*Figure 2. The structure of the projective plane PG(2, F₃). Parallel lines on the 3×3 grid (9 points, filled circles) meet at directions on the line at infinity L (open circles). Choosing P₀ (filled diamond) produces the partition 9+3+1 = 13, corresponding to the sector structure of SU(3)×SU(2)×U(1).*

Under the full symmetry PGL(3, F₃) of PG(2, F₃), all 13 points are equivalent and no partition arises. We now fix one of the 4 points on the line at infinity L as a reference point P₀ — the pair (P₀, L) is called a flag. Fixing the flag is a choice of reference within the construction: PGL(3, F₃) acts transitively on flags, so every choice yields an isomorphic partition. The fixing splits the 13 points into three orbits that do not mix under the transformations preserving the flag: 9 points off L (the 3×3 grid), 3 points on L excluding P₀, and P₀ itself. Choosing P₀ alone gives only two groups, 1+12; choosing a flag produces 9+3+1. This is the minimal structure distinguishing three forces:

$$N^2_{(9\;\mathrm{colour})} + N_{(3\;\mathrm{weak})} + 1_{(1\;\mathrm{hyp})} = N^2 + N + 1.$$

N = d is satisfied only by N = 3 (N = 2, 4, 5 all have d = 3 ≠ N) — a coincidence of consistency, not an input. All numbers in the partition 9+3+1 are determined by N alone: N² = 9, N + 1 = 4 (excluding the reference point gives N = 3), reference point 1. Total: N² + N + 1 = 13. The 52 flags (point-line incidence pairs) provide the basis for the Laplacian (§5.3).

## 4.2. From 13 Points to the Standard Model Gauge Group

The partition 9+3+1 has the same form as the sector structure of the Standard Model gauge group SU(3)×SU(2)×U(1)_Y. Here we hypothesise that the three blocks of this partition are the gauge sectors of the Standard Model. At this stage this is a provisional identification; but the coupling constants (§5) and the mixing angles (§8) are derived from it and agree with experiment, which corroborates the identification.

**9 points → SU(3).** The 9 off-line points constitute the affine plane AG(2, F₃). These 9 points are labelled by (a,b) ∈ F₃² and naturally index the finite Weyl operators W_{a,b} = XᵃZᵇ in the three-dimensional internal representation. The identity W_{0,0} is the central component; the remaining 8 non-trivial Weyl components form the traceless sector. By the Killing-Cartan classification, the unique rank-2, dimension-8 compact simple Lie completion is su(3). The 9th, central component corresponds to the global baryon number U(1)_B (an ungauged global charge; the exact anomaly-free combination is B−L, §8.3).

**3 points → SU(2).** The 3 points on L excluding P₀ give three non-trivial weak directions relative to the reference point. The minimal non-abelian compact Lie completion of these three directions is the 3-generator su(2), whose doublet representation gives the weak SU(2) action.

**1 point → U(1)_Y.** The reference point P₀ gives the abelian central module phase, and physically corresponds to hypercharge U(1)_Y.

**Direct product guarantee.** Translations of the affine plane fix the line at infinity L — shifting the grid does not change the "directions." The colour sector (9 points) transformations do not mix with the electroweak sector (3+1 points). This is a geometric property following from the incidence structure of PG(2, F₃): it yields a direct-product structure of three mutually non-mixing sectors — the same shape as the Standard-Model gauge group. The global form of the group (in the Standard Model, including the Z₆ quotient [38]) is not fixed by the geometry alone; it is determined by the faithful action on the particle carrier^II^.

$$n^2{=}n \xrightarrow{V=-H} N{=}3 \xrightarrow{\text{free module}} \mathbb{F}_3^3 \xrightarrow{\text{count directions}} 13\text{ points} \xrightarrow{\text{flag}} 9{+}3{+}1 \xrightarrow{\text{Killing-Cartan}} \mathfrak{su}(3)\oplus\mathfrak{su}(2)\oplus\mathfrak{u}(1).$$

This mathematical result agrees with the experimentally established gauge structure SU(3)×SU(2)×U(1).

The 4 on-line points correspond to the four real degrees of freedom of the SU(2) doublet — the three Goldstone bosons required for SU(2)×U(1) → U(1)_em breaking and the physical Higgs. V = −H takes its minimum −ln 2 at φ = 0 (σ = 1/2) and returns to zero as φ → ±∞ — a stable well (V''(0) = 1/4 > 0). The vacuum sits at the interior point σ = 1/2 rather than at the vanishing-field points σ ∈ {0,1} — the role played in the Standard Model by the negative quadratic term (μ² > 0), moving the vacuum from the symmetric point to a finite field value, is here played by the shape of the well.


## 4.3. Origin of Fermion Diversity

§4.2 established that the 13 directions of PG(2, F₃) give a structure with the same shape as the Standard-Model gauge group. But the 13 directions do not correspond to fermions themselves. They define the gauge symmetry — the structure of forces — and the fermion structure emerges as **representations** of this gauge group. Here we show why the single existence variable n gives rise to a diverse matter-particle structure.

Before the gauge structure emerges, n has no quantum numbers — no generation label, no colour, no isospin. Without distinguishing properties, multiple instances of n are identical (Leibniz's principle of the identity of indiscernibles), and n remains unique. The gauge structure of PG(2, F₃) changes this: the representations of the gauge group create diverse "slots" — configurations — in the internal space, with each configuration carrying an occupation number n_i ∈ {0,1}. The configuration decomposition is:

$$\sum_{i} \hat{n}_i = \hat{n}, \qquad \hat{n}_i^2 = \hat{n}_i, \qquad \hat{n}_i \hat{n}_j = \delta_{ij}\,\hat{n}_i$$

The three equations are one-particle statements: the first states that existence decomposes into configurations; the second that bivalence (n² = n) is inherited by each configuration; the third that the unique unit of existence occupies exactly one configuration. Taken together, the n̂ᵢ are mutually orthogonal projections summing to n̂ — the orthogonal resolution of the one-particle state space.

The Pauli exclusion principle itself — no double occupancy of a single mode, together with exchange antisymmetry — is carried by the many-body realisation (the CAR construction; CAR = canonical anticommutation relations)^II^, with each configuration as a mode. In the many-body realisation distinct configurations can be simultaneously occupied, but the occupation of each mode is restricted to {0,1} — this is the occupation-number representation of the Pauli exclusion principle [11]. Note that n_i = 0 does not mean non-existence (the σ = 0 forbidden by A3); it only means existence occupies a different configuration.

Each configuration i is characterised by two labels — generation k (which bound state of V = −H) and gauge representation α (which sector of PG(2, F₃)):

$$n_i = n_{k,\alpha},\qquad k = 1, 2, 3$$

The representations of SU(3)×SU(2)×U(1) give 16 states per generation — 48 configurations across the three generations (table below).

The chiral asymmetry (SU(2) acting only on left-handed fields) is part of this particle-species specification and is not derived in this paper. In Paper II, under the stated enumeration conditions, mixed configurations are excluded and the left/right assignment is unique up to a global orientation convention. The hypercharge values in the table below are the unique anomaly-free assignments compatible with the 9+3+1 sector structure, Yukawa closure, the inclusion of ν_R, and the specification of a single local abelian gauge factor (B−L remains an ungauged global charge, §4.2); the operator-level derivation is likewise given in Paper II:

| Representation | Content | States |
|---|---|---|
| (3, 2, 1/6) | Left-handed quarks Q_L | 6 |
| (3, 1, 2/3) | Right-handed up-type u_R | 3 |
| (3, 1, −1/3) | Right-handed down-type d_R | 3 |
| (1, 2, −1/2) | Left-handed leptons L_L | 2 |
| (1, 1, −1) | Right-handed charged lepton e_R | 1 |
| (1, 1, 0) | Right-handed neutrino ν_R | 1 |

Total: 16 × 3 generations = 48 configurations (counted as Weyl components, antiparticles not included, ν_R included). The inclusion of ν_R is required for Yukawa closure and for the B−L-protected Dirac-neutrino structure (see §8.3 of this paper). The generators of SU(3)×SU(2)×U(1) correspond to 8+3+1 = 12 gauge bosons (8 gluons, 3 weak bosons W⁺W⁻Z, 1 photon). Together with the Higgs boson from §4.2, the 48 fermions, 12 gauge bosons, and 1 Higgs constitute the complete particle content of the Standard Model.

Thus, diverse matter (n_i) arises from unique existence (n) — the diversity of matter is a consequence of gauge structure.

The internal consistency of this construction is independently confirmed by the eigenfunction participation ratio

$$D_{PR} = \frac{(\int \rho\, d\varphi)^2}{\int \rho^2\, d\varphi} = 13.06 \approx |\mathrm{PG}| = 13,\qquad \rho = \sum_n \psi_n^2$$

with ρ the total density of the three bound states (reproduced by the accompanying verification code), and by the spectral decomposition of the 52×52 flag Laplacian.[^flagspec] No step in the derivation contains a free parameter. The validity of this construction is supported by the fact that all independent physical quantities derived in §5–§10 agree with experiment to 0.0006%–2.6% accuracy.

[^flagspec]: The finite-kernel theorem and its spectral gauge-decomposition corollary in Paper II.


# 5. Gauge Coupling Constants

§4.2 derived the structure corresponding to the gauge group (which forces exist). This chapter derives their strengths (coupling constants).

## 5.1. Weinberg Angle

The partition 9+3+1 of §4.2 gives values corresponding to the coupling constants. The weight of each point is fixed by the symmetry prior to choosing a flag: PGL(3, F₃) acts transitively on PG(2,3), so the 13 directions carry equal weight (invariant measure). Fixing the flag does not alter these weights; it only classifies the 13 points into 9+3+1. On this basis, the Weinberg angle [21] (weak mixing angle) sin²θ_W is read and identified as the fraction that the 3 weak directions occupy among all 13 — taking all internal directions as the denominator^II^:

$$\sin²θ_W = 3/13 = 0.23077 \quad (\text{expt: } 0.23122\;\overline{\mathrm{MS}}\text{ at }M_Z,\; 0.2\%)$$

Alternatives such as taking only the four electroweak directions in the denominator correspond to a different observable type; that the type-preserving correspondence is unique to this fraction is proved^II^.

This angle also fixes the tree-level W/Z mass ratio (§5.4).

**Radiative corrections.** The predictions above are algebraic values, and their comparison with M_Z-scale experiments requires the 4D RG evolution and radiative corrections of the running quantities (§9.3). The next correction to the weak angle — the radiative correction, with expansion parameter α, that carries the tree value 3/13 to the MS-bar value at M_Z — is fixed in closed form by the framework itself:

$$\Delta_\mathrm{rad} = \alpha\cdot\frac{5N^3}{|\mathrm{PG}|^3} = \alpha\cdot\frac{135}{2197} = 0.00045$$

(N and |PG| are quantities of this paper; the uniqueness of this product form, including the coefficient 5, is proved as a theorem^II^.) This gives sin²θ_W = 3/13 + Δ_rad = 0.23077 + 0.00045 = 0.23122, in agreement with the MS-bar experimental value 0.23122 ± 0.00003 [22] within errors.

## 5.2. Strong Coupling Constant

The strong coupling constant follows by distributing the weak-sector fraction over the rank of the colour sector. Of the 9 colour-sector directions, the independent directions number rank(SU(3)) = 2 (Cartan generators); equidistributing sin²θ_W over these 2 directions: α_s^(0) = sin²θ_W / (N−1) = 3/26 = 0.1154 (the uniqueness of this distribution rule is proven^II^). The spectral marginality ε = 0.219 from §3 provides the NLO correction. The correction is distributed over N²+1 = 10 non-weak PG directions (colour 9 + hypercharge 1):

$$α_s = (3/26)(1 + ε/(N²+1)) = 0.1179 \quad (\text{expt: } 0.1180\;\overline{\mathrm{MS}}\text{ at }M_Z,\; 0.07\%)$$

The value k = N²+1 = 10 follows from the number of non-weak directions in PG(2, F₃). Throughout this paper the NLO corrections take the form 1 ± ℓε/k, where the denominator k is the size of the finite support over which the correction is distributed — k = N²+1 = 10 for the non-weak sector (α_s; affine + hypercharge), k = N² = 9 for generation mixing (sinθ_C; affine plane), k = |PG| = 13 for the full set of PG modes (m_H/v, sin²θ₂₃; full projective space), and k = N = 3 for the three weak directions (sin²θ₁₃). The assignments of (ℓ, k) to each observable, including the sign and coefficient ℓ, are each proven unique^II^.

## 5.3. Fine-Structure Constant

The Weinberg angle and strong coupling constant were determined by simple point ratios on PG(2,3). The fine-structure constant requires more — the full spectral structure of the internal space.

The inverse electromagnetic coupling 1/α measures the weakness of the electromagnetic interaction. The observed charge is the charge screened by the matter fields — the stronger the screening, the weaker the coupling and the larger 1/α. In this framework 1/α is given (up to the self-consistent correction −α) by the sum of the bare photon propagation cost λ₁ and the matter screening quantity Z, and the Thomson-point value is dominated by the screening term.

In 4D QED this structure is formalised as the Dyson equation. In 0D the same structure holds but becomes algebraically exact. The propagator is a constant, so the vertex correction and self-energy are determined by the same resolvent, and the 4D perturbative cancellation Z₁ = Z₂ holds as an exact algebraic identity^II^. The self-consistency equation is:

$$\frac{1}{\alpha} + \alpha = \lambda_1 + Z$$

The right-hand side has two terms: λ₁ is the bare photon inverse propagator (propagation cost) and Z corresponds to the scalar-channel spectral quantity (screening). This correspondence is a spectral one — Z is the resolvent trace whose energy-denominator structure matches the matter vacuum polarisation.

The form of the equation (the left-hand side 1/α + α with unit coefficient) is derived^II^ — what remains external is the comparison convention identifying the incidence current with the Thomson-point electromagnetic current. The full justification of this correspondence as a self-consistent prediction, and the proof that the physical root is unique under this identification, are both given in Paper II. The dominant contribution to Z is 1/|E₂| = 95.94; since E₂ is rigorously certified by the interval arithmetic of §3.3, the certified enclosure of E₂ does not change the displayed digits of this spectral contribution.

The minimum cost λ₁ for gauge-field propagation through internal space is the smallest non-zero eigenvalue of the Laplacian on internal space. Since gauge interactions involve both positions (points) and directions (lines), the natural object for diffusion is a flag (point-line incidence pair). The spectrum of the graph Laplacian ℒ on the flag graph of PG(2, F₃), a 52×52 matrix, is:

| Eigenvalue | Value at N = 3 | Multiplicity |
|---|---|---|
| 0 | 0 | 1 |
| (N+1) − √N | 2.268 | 12 |
| (N+1) + √N | 5.732 | 12 |
| 2(N+1) | 8 | 27 |

This spectrum follows elementarily from the incidence counts of the projective plane alone. Each point lies on 4 lines and each line carries 4 points, so flag adjacency (sharing the point or sharing the line) forms a 6-regular graph, and the Laplacian can be written

$$\mathcal{L} = 8\,\mathbf{1} - P - \Lambda$$

where P (Λ) is the operator summing over the 4 flags with the same point (the same line). On constant functions P = Λ = 4, so ℒ = 0 (multiplicity 1). On the space of functions whose sum vanishes within every point group and within every line group — of dimension 52 − 13 − 13 + 1 = 27 — one has P = Λ = 0, so ℒ = 8 (multiplicity 27).

The remaining 24 dimensions are spanned by pairs of a zero-sum point function u and a zero-sum line function v, on which ℒ acts as

$$(u, v) \mapsto (4u - Mv,\ 4v - M^{\mathsf T}u)$$

with M the 13×13 point-line incidence matrix. Since exactly one line passes through two distinct points, M Mᵀ = 3I + J (J the all-ones matrix), which on zero-sum functions reduces to M Mᵀ = 3I. Hence

$$(4-\lambda)^2 = 3,\qquad \lambda = 4 \mp \sqrt{3}\quad(\text{multiplicity 12 each})$$ The √3 comes from the uniqueness of the line joining two points; the 4 = N + 1 from the number of lines through a point.

The smallest non-zero eigenvalue λ₁ = 4 − √3 = 2.268 is the quantity that corresponds to the bare photon inverse propagator in the self-consistent equation 1/α + α = λ₁ + Z above.

The screening term Z corresponds to the scalar-channel spectral quantity from three generations of matter fields, determined by the resolvent of the internal kinetic operator. The internal kinetic operator is a tensor sum of the generation sector (K_gen, whose eigenvalues are the binding energies |E_n|) and the gauge sector (ℒ, eigenvalues λ_k):

$$K_\mathrm{int} = K_\mathrm{gen} \otimes \mathbf{1} + \mathbf{1} \otimes \mathcal{L}$$

This form — the tensor sum, with both factors measured in the same unit at relative scale 1 — is established as a uniqueness theorem^II^ under the single-clock unit clause (product couplings, twisted products, and rescalings of the relative unit are excluded). Decomposing the trace of the resolvent by multiplicity:

$$Z = \mathrm{Tr}[1/K_\mathrm{int}] = \zeta_V(1) + 12\,S_1 + 12\,S_2 + 27\,S_3 = 134.77613$$

where ζ_V(1) = Σ_n |E_n|⁻¹ and S_i = Σ_n (|E_n| + λ_i)⁻¹. The intermediate values are

$$\zeta_V(1) = 104.28714,\qquad 12S_1 = 14.57025,\qquad 12S_2 = 6.05684,\qquad 27S_3 = 9.86190$$

with ζ_V(1) dominated by the E₂⁻¹ contribution 95.94. The 27-dimensional sector (multiplicity 27 = 3³) contributes 9.86 ≈ N², providing the spectral structure corresponding to the QCD vacuum-polarisation correction (all intermediate values are reproducible with the accompanying verification code).

Substituting λ₁ = 4 − √3 = 2.26795 and Z = 134.77613 into the self-consistency equation:

$$\frac{1}{\alpha} = \frac{(Z + \lambda_1) + \sqrt{(Z + \lambda_1)^2 - 4}}{2} = 137.0368 \quad (\text{expt: } 137.036,\; 0.0006\%)$$

Adding the correction terms^II^:

$$\frac{1}{\alpha} + \alpha = \Lambda_0 - x\Bigl(1+\frac{\alpha}{2}\Bigr) + \frac{3}{4}\varepsilon^2 x + \frac{1589}{5408}\,x^2 - \frac{2}{13}\,\frac{x^3}{1+x},\qquad x = \frac{9\alpha}{26\pi}$$

where Λ₀ = λ₁ + Z is the right-hand side of the leading term and x is the same expansion parameter as the colour screening of §8.1. The origin of each term — the self-consistent screening correction, the generation branching for ε²x, the second-order screening count for x², and the first return of the 52-flag orbit for x³ — and the proof that the series closes in this form are given in Paper II. The physical root of this equation is

$$\frac{1}{\alpha}\Big|_{\text{all orders}} = 137.035999084$$

Here "all orders" refers to all orders of the typed finite Dyson closure defined in Paper II. The value is certified to 13 digits by the interval [137.0359990843979, 137.0359990844119]^II^. We note that the experimental recommendations split by measurement system: CODATA 2018, whose principal input is the Cs photon-recoil measurement, gives 137.035999084(21) — agreement to all displayed digits — while CODATA 2022, based on the Rb recoil, gives 137.035999177(21), 4.4σ away. The Cs and Rb measurements themselves disagree by more than 5σ; on this open experimental dispute the framework gives the prediction that the fine-structure constant lies on the Cs side (Table 4).

## 5.4. The W/Z Mass Ratio

From the weak angle the W/Z mass ratio follows: m_W/m_Z = cosθ_W = √(10/13) = 0.8771 — 0.5% below experiment; this residual has approximately the size and sign of the standard SM top-quark Δρ correction. That the custodial response is carried by exactly one neutral component of the one Higgs doublet is a theorem; what remains as identification is the naming of the experimental definition points (the on-shell convention) and the readout word. Under that identification, m_W/m_Z = 0.88137 follows uniquely and agrees with experiment at +1.0σ^II^.

# 6. The Higgs Sector

§5 derived the gauge coupling constants. The Standard Model has one more free parameter in the scalar sector — the Higgs self-coupling λ. This chapter derives the mathematical structure corresponding to it and shows the pinning of the internal vacuum and the boundedness of the internal potential (no instability channel).

## 6.1. Symmetry Breaking and Higgs Mass

Electroweak symmetry breaking is the statement that, while the dynamical laws remain symmetric, the vacuum (ground state) departs from the symmetric point. Since the W and Z are massive while the photon is massless, the vacuum cannot sit at the field-vanishing point (the symmetric point) — breaking is required by experiment. The Standard Model realises this through the assumption μ² > 0 in the Higgs potential V = −μ²|Φ|² + λ|Φ|⁴: it destabilises the origin, rolling the field to a finite value. But "why μ² > 0" is not explained within the Standard Model.

The V = −H(σ(φ)) of this paper answers this question through structure. The binary Shannon entropy H vanishes at both endpoints σ ∈ {0,1} (the trivial configurations corresponding to the field-vanishing point) and is maximal at the interior point σ = 1/2, so V = −H is a stable well that is 0 at the endpoints and reaches its minimum −ln 2 at the interior point (V''(0) = 1/4 > 0). Moreover A2 (σ ≠ 1) and A3 (σ ≠ 0) forbid the endpoints themselves. The unique minimum of the classical internal potential therefore necessarily lies at the interior point σ = 1/2 — the condition that the vacuum cannot sit at the symmetric point, i.e. the condition for breaking, follows from the shape of the entropy and the axioms rather than from destabilising the origin. What the Standard Model assumes as μ² > 0 becomes here a consequence of the shape of V.

Under the identification of §4.2, the line structure of PG(2, F₃) corresponds to the Higgs SU(2) doublet. The curvature at the centre — the Fisher metric at σ = 1/2 — gives the Higgs mass:

$$V''(0) = g_F(0) = \sigma(1-\sigma)\big|_{\sigma=1/2} = \tfrac{1}{4}$$

$$m_H/v = \sqrt{V''(0)} = 1/2 \quad (LO)$$

The NLO correction enters from all |PG| = 13 modes of the Laplacian matrix ℒ:

$$m_H/v = \frac{1}{2}\sqrt{1 + 2\varepsilon/|PG|} = 0.5084 \quad (\text{expt: } 0.5082,\; 0.03\%)$$

This corresponds to

$$\lambda = \tfrac{1}{2}(m_H/v)^2 = 0.1292\quad(\text{expt: } 0.1291,\ 0.06\%)$$

## 6.2. Vacuum Stability and the Hierarchy Problem

The Standard Model Higgs potential V_SM = −μ²h² + λh⁴ diverges as h → ∞, and the UV correction δμ² ∝ Λ² requires ~10³⁴ fine-tuning between the electroweak and Planck scales (hierarchy problem). Furthermore, the Higgs quartic coupling λ(μ) runs negative at μ ~ 10¹⁰ GeV, making the electroweak vacuum cosmologically metastable.

In the internal sector of V = −H, the instabilities corresponding to these problems do not arise. V = −H is bounded (V ∈ [−ln 2, 0]), with no high-energy divergence. The well depth |V_min| = ln 2 and curvature V''(0) = 1/4 are both O(1) algebraic constants; their ratio |V_min|/V''(0) = 4 ln 2 ≈ 2.77 is a ratio of internal quantities, an order-unity number requiring no fine-tuning:

$$V \in [-\ln 2,\ 0],\qquad |V_{\min}|/V''(0) = 4\ln 2 \approx 2.77$$

Relating this O(1) ratio to the 4D v²/Λ², and interpreting the boundedness as the reason the running coupling λ(μ) in 4D does not turn negative, both require the reconstruction matching (likewise the full 4D vacuum-stability interpretation and the correspondence to thermal properties such as symmetry restoration at high temperature)^II^. The shape of V is fixed by the axioms, and scale dependence does not alter the algebraic structure of V.

Furthermore, V = −H is uniquely determined as the canonical free energy of n ∈ {0,1} (§2). What the uniqueness of V fixes is the shape of the internal potential. Within the minimal operator algebra of the theory, no second fundamental scalar potential is generated, and that the weak-doublet scalar has multiplicity 1 — that the minimal reconstruction predicts a single Higgs doublet — is proved as a classification of one-forms on the finite internal carrier^II^.

**Summary of §5–§6.** sin²θ_W(M_Z) = 0.23122 (<0.01%), α_s(M_Z) = 0.1179 (0.07%), 1/α_em = 137.0368 (0.0006%; 137.035999084 at all orders^II^), m_H/v = 0.5084 (0.03%). All four coupling constants are derived from V = −H and PG(2,3) and agree with experiment. The vacuum is pinned at an interior point, the internal potential is bounded (no instability channel), and there is only one Higgs.

# 7. Fermion Masses

§4–6 established the structures corresponding to gauge structure, coupling constants, and the Higgs mass. Based on the provisional identification of the generations made in §3, this chapter derives formulas for the mass ratios from the localisation of the bound states; the inter-generation mixing is derived in the next chapter. The derived mass ratios agree to high precision with experiment for charged leptons, quarks, and neutrinos — this agreement corroborates the provisional identification of §3.

## 7.1. Mass Formula

The three bound states of V = −H differ in localisation. The deepest ψ₀ is localised at the well bottom and couples most strongly to the gauge field, corresponding to the heaviest τ:

| Bound state | Energy | Localisation | Particle | Mass (MeV) |
|---|---|---|---|---|
| ψ₀ | E₀ = −0.485 | High (well bottom) | τ | 1777 |
| ψ₁ | E₁ = −0.159 | Medium | μ | 106 |
| ψ₂ | E₂ = −0.010 | Low (well edge) | e | 0.511 |

The mass ratios are determined by the localisation of the wave functions, not as independent parameters. The hierarchy that the Standard Model inputs by hand as three Yukawa couplings is given here by the difference in localisation of the three modes of one well. The key quantity is R = ln(m_τ/m_e)/ln(m_μ/m_e), which compresses the three-body mass hierarchy into a single number.

$$V = -H \xrightarrow{\text{§3}} \psi_n \xrightarrow{\text{IPR}} \mathrm{localisation} \xrightarrow{b} m_n \xrightarrow{} R = 1.5293$$

Each eigenstate has a binding energy |E_n| and a localisation, measured by the inverse participation ratio IPR_n = ∫|ψ_n|⁴dφ. Since the coupling to the gauge field is proportional to the wave-function amplitude at the interaction point, a more localised state acquires a larger mass: m_n ∝ IPR_n^b.

The exponent b is not a free parameter. In 0D there are no loop integrals, hence no running that would generate an anomalous dimension dynamically. By the Peter-Weyl theorem for finite groups, the heat kernel on the group is the finite sum K(t) = Σ_ρ dim(ρ) χ_ρ(g) e^{−λ_ρ t}, admitting only exponential decay with no power-law terms — the exponent is a structurally fixed quantity. The exponent b is the contact index derived from the normalisation of the finite incidence geometry PG(2, F₃): with the normalised single-leg index w_leg = (N+1)/N, b = 2w_leg = 8/3[^casimir]. The uniqueness of this mapping rule is a theorem^II^[^massuniq].

To this, a small correction from the anharmonicity of V = −H is added. The virial ratio |⟨V⟩_n|/T_n (potential to kinetic energy) measures how differently each mode samples the well, and its exponent is read as the per-mode marginality ε/N = 0.219/3 = 0.073 (ε from §3.2). The reading takes the total deficit to be carried by the three generation terminals (Tr D = ε) and reads the generation-terminal weight as the normalised evaluation on the three-terminal support, (1/3)Tr D = ε/3, independently of how the insertion is distributed. The value itself agrees with the Koide inversion below to 0.07%.

The complete mass formula is:

$$m_n \propto |E_n| \times \mathrm{IPR}_n^{\,b} \times (|\langle V \rangle_n| / T_n)^{\varepsilon/N}$$

The three factors have different origins:

| Factor | Physical meaning | Determined by |
|---|---|---|
| \|E_n\| | Energy scale | Binding energy |
| IPR^{b} | Localisation → coupling strength | Finite incidence geometry (contact index) |
| (|⟨V⟩|/T)^{ε/N} | Anharmonic correction | Spectral marginality |

All exponents are determined by V = −H. The result:

$$R = \frac{\ln(m_\tau/m_e)}{\ln(m_\mu/m_e)} = 1.5293 \quad (\text{expt: } 1.5294,\; 0.006\%)$$

(All spectral quantities computed with |φ|_max = 100, Δφ = 1.25×10⁻⁴ (N_grid = 1600001), verified stable under grid-resolution doubling and box enlargement.)

Koide's relation Q ≈ 3/2 [23] (this paper's convention Q = (Σ√m)²/Σm is the reciprocal of the standard 2/3) is not an input that determines b, but a cross-check — a verification on the output side. Parametrising the mass formula in Koide form yields the pair (b, c), where c is the exponent of the virial correction. Back-solving, for each b, the c required to force Q = 3/2, only b = 2w_leg = 8/3 returns a value, 0.0730, matching the spectral value ε/N = 0.073 (to 0.07%):

| b | c (Q=3/2) | Deviation |
|---|---|---|
| 2.0 | 0.437 | 500% |
| 2.5 | 0.163 | 120% |
| **8/3 = 2w_leg** | **0.0730** | **0.07%** |
| 2.78 | 0.012 | 84% |
| 3.0 | -0.107 | 246% |

| | Prediction | Experiment | Accuracy |
|---|---|---|---|
| m_τ/m_e | 3481 | 3477 | 0.1% |
| m_μ/m_e | 207.0 | 206.8 | 0.1% |
| R | 1.5293 | 1.5294 | 0.006% |
| Q (Koide) | 1.499985 | 1.500005(11) | 2σ (from PDG masses) |

The 0.006% agreement of R = 1.5293 quantitatively supports the identification of the bound states with the lepton generations. Because this quantity is measured to 6 × 10⁻⁶, the 0.006% difference amounts to about 10σ of the experimental uncertainty; it is a consequence of the truncation order used here, and the evaluation including the next order gives R = 1.529365 (0.0008%, −1.3σ), and the evaluation including the feedback term gives R = 1.529371 (−0.7σ)^II^. The stability of the Koide relation under radiative corrections is discussed by Sumino [43].

[^massuniq]: The functional form of the mass response and the contact exponent b = 2w_leg = 8/3 are given as uniqueness theorems in Paper II; the latter is a theorem once the correspondence carrying the finite-side eigenvalue into a power of the contact coordinate is fixed.

[^casimir]: This is not the physical colour charge of leptons. The exponent b is the contact index of the incidence geometry; its numerical agreement with the SU(3) fundamental Casimir, 2C_F = 8/3, holds only at N = 3 — (N²−1)/2N = (N+1)/N only for N = 3 — and is a consistency check at N = 3. See Paper II [20].

## 7.2. Quark Masses

The mass formula of §7.1 reproduced R = 1.5293 for leptons. Can the same formula be applied to quarks? The formula itself does not change — only the well depth changes.

Leptons are colour singlets and feel V = −H as a single copy. Quarks carry three colours. The colour-triplet well composes the three colour contributions — isomorphic copies of the same Bernoulli degree of freedom — as independent factors pulled back to the common existence coordinate φ: Z₃ = Z₁³, additivity of the logarithm gives log Z₃ = 3 log Z₁, and additivity of the Legendre dual gives V = −3H^II^; weak isospin does not deepen the well, because the doublet is an orientation structure within a single copy — it fixes which component is called I₃ = +1/2 (the g_c correction below) — and adds no second copy. The deeper well supports 5 bound states (n_WKB = 4.82). In the tables below, κ denotes the well-depth multiplier (the κ in V = −κH):

| Sector | Operator | κ | Bound states | Modes used |
|---|---|---|---|---|
| Lepton | −H | 1 | 3 | (0,1,2) = τ,μ,e |
| Up-type | −3H | 3 | 5 | (1,2,3) = t,c,u |
| Down-type | −3H (correction in kinetic term) | 3 | 4 | (0,1,2) = b,s,d |

**Up-type.** Mode 0 is deeply confined and belongs to the I_3 = −1/2 sector. Up-type quarks (I_3 = +1/2) use modes (1,2,3), giving R_up = 1.777 (expt: 1.7697(25), +0.41% = +2.9σ). The evaluation including the colour feedback term gives R_up = 1.767 (−0.17%, −1.2σ), and the same evaluation gives the light-quark mass ratio m_u/m_d = 0.477 (expt: 0.460(16))^II^.

**Down-type.** Up-type and down-type belong to the same SU(2) doublet and feel the same V = −3H, but differ in the sign of weak isospin (I_3 = +1/2 vs −1/2). This sign supplies a correction to the mass formula.

The total gauge charge of the SU(2) doublet (u, d), 2C_F + 1/N_c = N_c, is an exact algebraic identity. The physical mechanism originates from PG(2, F₃) geometry: among the N+1 = 4 points on line L, the reference point P₀ is the reference that orients the doublet (which component is called I_3 = +1/2), and the remaining 3 points generate SU(2). Reversing this orientation flips the sign of I_3 and reverses the SU(2) action, coupling the total charge N_c to g_F(φ) = σ(1−σ) (the Bernoulli Fisher metric g_F = Var(n) of §2.2; its uniqueness is Čencov's theorem [13]) with a negative sign. The effective mass is:

$$m_\mathrm{eff}(φ) = 1 − N_c σ(1−σ)$$

This is the static Fisher-limit form of the down-sector mass formula. The operator actually solved is the symmetric position-dependent-mass form −½ ∂_φ (1/m_eff) ∂_φ − 3H (identical to the discretisation in the accompanying verification code).

The correction coefficient g_c = −N_c = −3 (distinct from the normalisation g of §3.1) is fixed by the algebraic structure of the SU(2) doublet. The corrected operator has 4 negative eigenvalues (−1.5709, −0.5455, −0.1952, −0.00733; stable for box size L_box = 40/80/120), of which 3 modes (0,1,2) correspond to the physical window b, s, d — the physical window is uniquely determined by type exclusivity and order preservation^II^ (the window differs sector by sector: all 3 modes for leptons, (1,2,3) for up-type, (0,1,2) for down-type). With modes (0,1,2): R_down = 2.273 (expt: 2.276, 0.1%).

The mass ratios of all three charged-fermion sectors are reproduced:

| Sector | κ | g_c | R (pred.) | R (expt.) | Accuracy |
|---|---|---|---|---|---|
| Lepton | 1 | 0 | 1.529 | 1.529 | 0.006% |
| Up | N_c = 3 | 0 | 1.777 | 1.770 | 0.4% |
| Down | N_c = 3 | −N_c | 2.273 | 2.276 | 0.1% |

Here R captures the shape of the log-hierarchy---it is invariant under a uniform stretching of the mass ratios (raising every m_i/m_j to a common power)---not the individual masses. The accuracy column gives the deviation of the prediction from the experimental central value. The experimental R inherits the light-quark mass uncertainties (PDG 2026: m_u = 2.16 ± 0.07 MeV, etc.), which propagate to a width of about ±0.15–0.25%, and the comparison further depends on the scheme and scale conventions for each mass (light quarks MS-bar at 2 GeV, c and b at m(m), t from direct measurement).[^absquark]

[^absquark]: The individual quark magnitudes are not derived in this paper. Paper II derives the within-sector ratios and then gives a closed form for the six absolute masses written in derived quantities alone — the fixings of meaning are explicit physical identifications, with the conversion between schemes and definition points remaining as continuing computation. Agreement with the six measured masses used in assessing those identifications is not counted as independent validation.

## 7.3. Neutrino Masses and the Absolute Scale

§7.1–7.2 derived the charged-fermion (lepton and quark) mass ratios. What remains are the neutrinos.

**Mass predictions.** Neutrinos are colour singlets with κ = 1, coupling to V = −H — the same potential as charged leptons. The mass ratio R depends only on the well shape (IPR and virial ratios), not on particle species. The same well gives the same R:

$$R_\nu = R_\mathrm{lepton} = 1.529$$

This is the input-free prediction of this section. Together with the spectral ratio (m₃/m₂)² = 169/5 established below, it fixes the *dimensionless* mass ray completely, with no oscillation input at all.

Absolute masses in meV require one quantity that carries the dimensionless internal quantities into human units. This paper takes it to be the electron mass m_e = 0.51099895069(16) MeV (CODATA 2022 [34]) — the most precisely measured particle mass (relative 3 × 10⁻¹⁰), and a member of the linked mass system of charged leptons and neutrinos. This one point fixes the internal mass unit, and the electroweak scale v becomes a quantity predicted in that unit: v = 246.2360 GeV (+6.7 × 10⁻⁵ against the 246.2196 GeV fixed by the Fermi constant G_F = 1.1663788(6) × 10⁻⁵ GeV⁻²)^II^. No second scale is needed for the neutrino sector, and m_e itself is not counted as an independent prediction. The mass ratio between the lightest neutrino and the electron is derived from the specified operator word and read as a physical mass ratio in the same unit^II^:

$$\frac{m_1}{m_e} = x^3\,\frac{7}{6}\Bigl(1 + \frac{\varepsilon^2}{9}\Bigr) = 6.097 \times 10^{-10},\qquad x = \frac{9\alpha}{26\pi}$$

Here x is the expansion quantity of §5.3, which reappears in §8.1 as the colour screening of the Cabibbo angle; the ten orders of magnitude by which the neutrino masses fall below the charged-lepton masses follow from that quantity entering at the third power. The unit that displays the charged leptons therefore carries over to the neutrinos. No oscillation data enter the choice of any formula, coefficient or branch. Conditional on the normal-ordering branch, this gives[^numass]:

| Quantity | Prediction | Status |
|---|---|---|
| m_1 | 0.3116 meV | Unmeasured |
| m_2 | 8.667 meV | Consistent |
| m_3 | 50.39 meV | Consistent |
| Σm_i | 59.37 meV | Below Planck (<120) |

[^numass]: The displayed digits depend on where the higher-order charged-lepton evaluation used for the absolute normalisation is truncated; an identification concerning corrections beyond the displayed order is not yet fixed, and the displayed digits carry that condition^II^. Since the mass ray is fixed internally, all three masses scale with the common unit. The theory has no free parameters. Normal ordering is an input condition in this section; its independent derivation is given as a theorem^II^ (a consequence of R_ν > 1 and (m₃/m₂)² = 169/5 > 1). These values are consequences of setting the absolute scale m_e, and are not counted as independent predictions.

In that same unit, both oscillation splittings become quantities to be compared only after the computation is complete: Δm²₂₁ = 7.5016 × 10⁻⁵ eV² (against the measured 7.537 × 10⁻⁵ eV², −0.47%) and Δm²₃₁ = 2.5387 × 10⁻³ eV² (against the measured 2.521 × 10⁻³ eV², +0.70%). One degree of freedom tests this common unit, however, so the two must not be counted as two independent tests; their ratio tests the mass ray discussed next.

What the flag-Laplacian spectrum fixes exactly is

$$\Bigl(\frac{m_3}{m_2}\Bigr)^{\!2} = \frac{169}{5},\qquad \frac{m_3}{m_2} = \frac{|\mathrm{PG}|}{\sqrt{N+2}} = \frac{13}{\sqrt{5}} = 5.814$$

The numerator 169 is fixed by two elementary facts that follow from the spectrum of §5.3 alone. First, the non-zero eigenvalues 4∓√3 form a conjugate pair whose product reproduces the size of the projective plane exactly,

$$(4-\sqrt3)(4+\sqrt3) = 16 - 3 = 13 = |\mathrm{PG}(2,\mathbb F_3)|$$

so no separate normalisation enters. Second, write ς for the involution exchanging the eigenvectors v₁, v₂ of that pair and consider the doublet ψ = cos θ_H·v₁ + sin θ_H·v₂ (the coefficient angle θ_H is distinct from the atmospheric angle θ₂₃). Exchange invariance ςψ = ±ψ forces equal weights |cos θ_H| = |sin θ_H|, so the coherent amplitude has magnitude |13 sin 2θ_H| = 13 regardless of the choice of sign. The squaring of this amplitude to the intensity 13² = 169 is not a second identification; it is a theorem under the identification^II^.

The denominator 5 = N+2 = 3+2 is the rank of the solar-channel inventory, and its derivation — including the exhaustiveness of that inventory — is given as an operator-trace theorem^II^. The mass-squared-difference ratio, obtained by eliminating t = m₁/m₂ using R_ν = R_lepton, is given by the full relation^II^

$$\frac{\Delta m^2_{31}}{\Delta m^2_{21}} = \frac{169/5 - t^2}{1 - t^2} = 33.842\ldots\quad(\text{displayed as } 33.84)$$

With non-zero m₁, 169/5 is not the exact value of the mass-squared-difference ratio; the full relation above gives 33.84, with no oscillation input. The measured ratio is 2.521/0.07537 = 33.45, so the theoretical and measured values differ by 1.2%. Because no oscillation data were used to set the absolute scale, this difference is not absorbed anywhere: it stands as the quantity under test, and the unrounded masses behind the table above satisfy the internal ratio m₃/m₂ = 13/√5 = 5.814 exactly (recomputing from the displayed digits leaves rounding residue).

**Scope of the mass formula.** The mass formula of §7.1 determines mass ratios within a sector, and its overall proportionality constant is sector-dependent. Only the ratio joining the lightest neutrino to the electron is derived from the specified operator word, and with it the lepton mass system closes on the single absolute scale m_e^II^ — the dimensionless inter-sector ratio m(ν₃)/m(τ) ~ 3×10⁻¹¹ is fixed along with it. The absolute quark masses are given by the closed form of Paper II (footnote of §7.2), with the conversion between schemes and definition points remaining as continuing computation.

# 8. Mixing Angles — CKM and PMNS

§7 derived the fermion mass ratios. This chapter derives the inter-generation mixing — the quark CKM matrix and the neutrino PMNS matrix — and, with it, the Dirac nature of the neutrinos.

## 8.1. CKM Mixing

§7.1–7.2 established the mass ratios. With masses determined, we now derive the inter-generation mixing [24]. Why is V_ub small, and why does the Cabibbo angle [19] take its particular value? The Standard Model does not explain these. In V = −H, the potential symmetry corresponds to the former and the marginality of the third generation to the latter.

The potential V = −H is symmetric: V(φ) = V(−φ). This fixes the parity of each eigenfunction: ψ_0 even, ψ_1 odd, ψ_2 even. The direct overlap between distinct modes vanishes by Sturm-Liouville orthogonality: ⟨ψ_0|ψ_2⟩ = 0. Since ψ_0 and ψ_2 are both even, ⟨ψ_0|M|ψ_2⟩ = 0 for a V-parity-odd mixing operator M — the leading 0↔2 transition is forbidden:

$$V_{ub}\big|_\mathrm{tree} = 0$$

The measured value |V_{ub}| ≈ 0.004 arises at higher orders of the Wolfenstein hierarchy [25]; the first non-zero term is given by the closed form below.

The Cabibbo angle follows from spectral marginality. The third bound state ψ₂ is marginally bound — its binding energy |E₂| = 0.010 is only 1.5% of the well depth. The parameter measuring this marginality is the gap between the exact bound-state count and the semiclassical action count:

$$\varepsilon = N - n_\mathrm{WKB} = 3 - 2.781 = 0.219$$

Mixing angles reflect the mismatch between mass and weak eigenstates; the marginality of the third bound state (ε > 0) is the source of mixing, and in the formal limit ε → 0 each transition formula tends to zero, leaving the three generations intact (that this formal limit does not define a family of unitary CKM matrices is noted in §8.2). The expression giving the angle as a series in ε (which operator carries the transition, and what fixes its order and denominator) is proven unique per observable^II^ (see the list in §5.2). To second order:

$$\sin\theta_C = \varepsilon(1 + \varepsilon/N^2) = 0.2245 \quad (\text{expt: } 0.2243,\; 0.10\%)$$

This value is the unscreened internal value. The value composed with colour screening is sinθ_C·(1−x) = 0.22435 (x = N²α/(2π|PG|) = 9α/(26π) is the colour-screening expansion parameter; derived in Paper II).

The Wolfenstein hierarchy follows from the Cabibbo angle: |V_{us}| ∼ ε, |V_{cb}| ∼ ε², and |V_{ub}| is suppressed to at least O(ε³) with first non-zero term 2ε⁴ (experiment: 0.225, 0.041, 0.004). NLO corrections take the form 1 ± ℓε/k. The assignment of the denominator k and of the sign and coefficient ℓ is as listed in §5.2 (1+ε/9 for sinθ_C, 1−2ε/3 for sin²θ₁₃). The full closed forms, proved as a theorem^II^, are

$$|V_{cb}| = \varepsilon^2\left(1-\frac{2\varepsilon}{3}\right) = 0.0410,\qquad |V_{ub}| = 2\varepsilon^4\left(1-\frac{2\varepsilon}{3}\right)\left(1-\frac{3\varepsilon}{26}\right) = 0.00384$$

The origin of these exponents and of the coefficient 2 is discrete. Generation transitions are carried only by parity-odd permutations of the three generations (a consequence of V-parity). The adjacent transitions 1↔2 and 2↔3 are realised by the transpositions t₁₂ and t₂₃. With respect to the flag V₁ ⊂ V₂ fixed by the generation order (the line spanned by the first generation and the plane spanned by the first two), count the cost of a wall as the shallowest flag level it moves: t₁₂ moves V₁, cost 1; t₂₃ fixes V₁ and moves only V₂, cost 2 — these are the exponents of |V_us| ∼ ε and |V_cb| ∼ ε².

Parity makes the direct one-step 1↔3 element vanish, but it does not by itself forbid the two-step composite 1→2→3: two odd insertions compose to an even one, and the product of the two matrix elements is nonzero. What separates them is the generative type—the composite carries a different composition tree, so it is a different observable rather than a neglected term of the same one^II^.

The odd permutation carrying 1↔3 directly is the longest permutation w₀, which reverses the three generations. Its gallery crosses all three walls t₁₂, t₂₃, t₁₃, with total cost 1 + 2 + 1 = 4 — the reason the first non-zero term is of order ε⁴. The coefficient 2 arises because w₀ has exactly two shortest decompositions (galleries), exchanged by point-line duality. The assignment of ε per wall and this count are proved as theorems under the specified mapping rule^II^.

The parity symmetry V(φ) = V(−φ) suppresses the direct 0↔2 flavour transition, but it does not by itself fix the QCD vacuum angle: parity alone leaves θ ∈ {0, π}.

The observable strong-CP angle is

$$\bar\theta = \theta_\mathrm{QCD} + \arg\det M_q$$

That both terms vanish is shown by the strong-CP matching theorem of Paper II: applying A1–A3 so as to preserve the action, the fermion determinant and the realisation of large gauge transformations gives arg det(M_u M_d) = 0 and a vanishing bare θ_QCD, hence θ̄(μ₀) = 0 at the matching point. What removes the bare θ_QCD is that the only unitary scalar preserving the positive measure ray is 1 (the measure-line triviality^II^) — not V-parity on its own.

However, this is a statement at the matching point and not an exact zero at all scales. For the infrared θ̄_IR, which includes weak CKM loops and thresholds, Paper II gives the conditional bound |θ̄_IR| < 6.0×10⁻¹⁷ at the declared orders — weak O(G_F²), chiral O(p²), in the large-N_c chiral-perturbation framework; the all-orders evaluation is a continuing computation. Consequently, beyond a PQ-type axion — excluded by the absence of a PQ-type global U(1) — no claim is made about the non-existence of axion-like particles in general.

Yet CP violation is observed in the CKM sector (J = 3.16 × 10⁻⁵ [22]). This is consistent with the bare θ_QCD = 0: the CKM phase arises at higher order in ε:

$$J = \frac{\pi\varepsilon^5}{N_\mathrm{flags}} = \frac{\pi\varepsilon^5}{52} = 3.06\times10^{-5}$$

(N_flags = 52 is the total number of flags, §4.1; a theorem under the specified mapping rule^II^.) The matching-point condition remains exactly zero.

The strong-CP attributions are summarised:

| Quantity | Prediction | Experiment |
|---|---|---|
| Bare θ_QCD | **Not a consequence of V-parity** (matching-point value θ̄(μ₀)=0 from the strong-CP matching theorem^II^) | — (the bound \|θ̄_IR\| ≲ 10⁻¹⁰ applies to the infrared quantity; not directly comparable) |
| CKM CP violation | J = πε⁵/N_flags (theorem under the mapping rule^II^) | 3.06 × 10⁻⁵ (expt 3.16 × 10⁻⁵) |

## 8.2. PMNS Mixing

In §8.1, CKM mixing angles were obtained from ε = 0.219 — small values (sinθ_C ≈ 0.22). Yet PMNS mixing angles are large (sin²θ₁₂ ≈ 1/3). Why does the same V = −H produce both small and large mixing?

The answer lies in colour charge. Quarks carry colour and reside in the deep well V = −3H, which places their mixing in the perturbative regime of small mode overlap — a localisation picture consistent with the smallness of quark mixing, whose expansion parameter is the generation-sector marginality ε (§3.2); the formulae themselves rest on the generating flag orders and screening words of Paper II.

Neutrinos are colour-neutral; instead, the geometric allocation of gauge directions on PG(2,3) directly gives values corresponding to the mixing angles. The solar angle is the fraction of one complete electroweak line (N+1 = 4 points) in the 13 directions of PG:

$$\sin²θ_{12} = 4/13 = 0.3077 \quad (\text{expt: } 0.3088,\; 0.4\%)$$

At LO the atmospheric angle is the fraction that the disjoint union of the N pure-weak directions and another line through the reference point (N+1 points) occupies among the 13 directions of PG, namely 7/13. With the PG-scale NLO correction (k = |PG| = 13):

$$\sin²θ_{23} = (7/13)(1 + ε/13) = 0.5475 \quad (\text{second octant})$$

The LO fraction 7/13 exceeds 1/2, predicting the second octant. In NuFIT 6.1 [28] the global best fit has moved to the first octant (sin²θ₂₃ = 0.470), but the octant is unresolved and the prediction 0.5475 lies within the 3σ range (0.432–0.587). The second-octant local minimum of the NuFIT 6.1 χ² profile lies at 0.550, and the prediction deviates from it by 0.45%.[^nufitgrid]

[^nufitgrid]: The local minimum is obtained from the released grid by projecting over Δm²₃₁ and δ_CP; it lies at Δχ² ≈ 1.03 above the global minimum, with a grid step of 0.005 in sin²θ₂₃. Unless stated otherwise, all NuFIT 6.1 values quoted in this paper are the without-SK-atm, normal-ordering set, consistent with the Δm² inputs used above.

The reactor angle is sin²θ_{13} = 1/(N·|PG|) = 1/39 with NLO correction (k = N = 3):

$$\sin²θ_{13} = (1/39)(1 - 2ε/3) = 0.0219 \quad (\text{expt: } 0.02249,\; 2.6\%)$$

Sum rule: 7/13 = 3/13 + 4/13 (at LO, sin²θ₂₃ = sin²θ_W + sin²θ₁₂) is a finite-counting identity with denominator 13. This is a formal sum, not a partition of point sets: the 3 pure-weak points are contained in the 4 points of the complete electroweak line (the two sets overlap; their union is 4 points). The geometric origin of the atmospheric 7 directions is the disjoint union above. NLO corrections break this sum rule; Paper II gives its form.

**The PMNS and CKM phases.** The atmospheric holonomy ratio q₂₃ = (7/13)(1+ε/13) = sin²θ₂₃ also fixes the PMNS phase directly, proved as a theorem^II^:

$$\delta_{\rm PMNS} = 2\pi q_{23} = 197.1^\circ$$

The CKM phase does not reduce to a single closed-form expression of ε alone: it is fixed by embedding the three mixing quantities of §8.1 — the screened sinθ_C, |V_cb|, and |V_ub| — together with the CP-odd invariant J = πε⁵/N_flags in one unitary CKM matrix and solving sinδ = J/(c₁₂s₁₂c₂₃s₂₃c₁₃²s₁₃) together with cos δ > 0 (the sign fixed by the initial path lifting of A1). This is proved as a theorem^II^ and gives

$$\delta_{\rm CKM} = 62.595^\circ,\qquad \beta = 23.972^\circ,\qquad \gamma = 62.560^\circ$$

with the experimental values δ_CKM = 66.1° and γ = 66.4° [22]. Because β and γ follow rigidly from the same unitary matrix, their deviations from experiment are one and the same phase difference and are not counted as two independent tests^II^. The three angles and J are derived independently, so |sinδ| ≤ 1 is a non-trivial constraint; it holds for ε ≳ 0.204, and the derived value 0.2192 lies inside. Each transition formula formally tends to zero as ε → 0, but the four fit into a single unitary CKM matrix only in that range; they do not define a common ε → 0 CKM family.

## 8.3. Neutrino Properties

Whether neutrinos are Dirac or Majorana is an unresolved experimental question. This framework gives a clear prediction. Combining the phase covariance of the number-operator algebra (a → e^{iθ}a) with the derived particle census, the single Higgs, the four non-vanishing Dirac edges, the anomaly-free integral character lattice, and the naturality of the generation structure fixes the residual global phase direction uniquely as B−L^II^. The Majorana bilinear aa + a†a† (Δ(B−L) = 2) carries non-zero weight under this B−L and therefore has no invariant component. V-parity then pairs the distinct states ν and ν^c (lepton number L = ±1) into a Dirac fermion (operator-level proof: the Dirac-selection theorem of Paper II [20]). Therefore neutrinoless double beta decay is exactly zero:

$$0\nu\beta\beta = 0$$

With this, fermion mass ratios (§7), mixing angles, and neutrino properties have been obtained from the structure of V = −H and PG(2, F₃) (the absolute masses with normal ordering assumed and m_e as the single absolute scale). Where the Standard Model treats these as independent free parameters, this framework treats them all as different aspects of the algebraic expansion of a single equation V = -H.

# 9. Connection to Spacetime

All calculations in §2–8 were performed in internal space (the φ coordinate). Neither space nor time appeared. How, then, does this zero-dimensional theory connect to four-dimensional spacetime? This chapter answers in three steps: the dimension and signature (3,1) of spacetime (§9.1), the independence of internal space and spacetime (§9.2), and why the zero-dimensional algebraic values correspond to the four-dimensional experimental values (§9.3).

## 9.1. Spacetime Dimension and Gravity

§3.4 derived the spatial dimension d_space = 3. But "why four-dimensional spacetime?" is a separate question — it is not obvious that time is one-dimensional. We show that two independent paths give the same D = 4 — both are consistency checks with stated premises (hereafter we distinguish the spatial dimension d from the spacetime dimension D = d + 1).

The first path comes from gravitational degrees of freedom. A massless spin-2 particle (graviton) in D-dimensional spacetime has D(D−3)/2 physical degrees of freedom (following from Poincaré group representation theory, without assuming any particular D or Einstein's equations). Meanwhile, the N = 3 bound states of V = −H provide N−1 = 2 independent internal degrees of freedom (since three probabilities sum to one). As a consistency check that takes Poincaré symmetry, a spin-2 field, and gauge reduction as premises, matching the spacetime dynamical degrees of freedom (graviton polarisations) one-to-one with the internal degrees of freedom gives:

$$D(D−3)/2 = N − 1 = 2 \quad \Longrightarrow \quad D = 4$$

The second path comes from projective geometry. The projective plane PG(2, ℂ) = ℂP² is a real 4-dimensional manifold, and for N = 3:

$$D = \dim_\mathbb{R}(\mathbb{C}P^{N-1}) = 2(N-1) = 4$$

Both are consistent only for N = 3.

These two paths are independent — one from field-theoretic representation theory (spin-2 polarisations), the other from pure algebraic geometry (complex carrier ℂP²). That both give D = 4 is a non-trivial consistency check of N = 3.

Even with D = 4 fixed, the signature could be (4,0), (3,1), or (2,2). The finite field F₃ (from the primality of N = 3) and the complex carrier ℂ³ (from the unitary realisation) are independent constructions moored to the same N = 3 — since the characteristics differ, no field homomorphism F₃ → ℂ exists, and no lift via F₉ succeeds either. What links the two is the correspondence between the exponent labels of the Heisenberg–Weyl system (F₃ side) and the representation carrier (ℂ³ side); what is preserved is the Weyl commutation relations and the MUB structure. ℂ³ carries a complex structure from the outset. Furthermore, the transfer action S = ∫[½φ̇² + V]dτ of the internal coordinate φ opened by A2 (§2.2) has a single evolution parameter τ, so one of the four real directions is distinguished as the evolution direction.

The complex structure and the single evolution direction alone do not uniquely fix the signature. (3,1) is fixed as the connection to the determinant form on M₂(ℂ)_sa and the structure SL₂(ℂ) → SO⁺(1,3).

Indeed, an element of M₂(ℂ)_sa is parametrised as

$$X = \begin{pmatrix} x^0 + x^3 & x^1 - i x^2 \\ x^1 + i x^2 & x^0 - x^3 \end{pmatrix},\qquad \det X = (x^0)^2 - |\boldsymbol{x}|^2$$

and the determinant is precisely a quadratic form of signature (+,−,−,−). The SL₂(ℂ) action X ↦ AXA† preserves this determinant, giving SO⁺(1,3).

For a four-dimensional metric theory satisfying the locality, naturality, diffeomorphism-covariance, divergence-free and second-order assumptions of Lovelock's theorem [30], the field equation is restricted to the Einstein form with a cosmological term

$$G_{\mu\nu} + \Lambda g_{\mu\nu} = \kappa T_{\mu\nu}.$$

Jacobson [31] showed that, assuming local Rindler horizons, the Unruh temperature, an entropy density proportional to area, and local equilibrium, the same equation follows from δQ = T_H dS on a causal horizon. The gravitational action and its coefficient are derived from the four-dimensional local manifestation in a subsequent paper. Jacobson’s relation is an independent consistency check of that result: when the assumptions above are realised on the same local manifestation, δQ = T_H dS reproduces the same Einstein equation. What stands here is therefore not a second derivation from the identity S = H(σ) alone, but agreement between the derived gravitational structure and its thermodynamic description.

## 9.2. Why Internal Space and Spacetime Are Independent

The coordinate φ = logit(σ) lives in internal space; the coordinate x lives in spacetime. Why are they independent?

The answer lies in the construction of V = −H. The potential V = −H(σ(φ)) is built as the internal zero-mode of the Bernoulli degree of freedom and contains no explicit spacetime coordinate x. However, the coordinate-independence ∂V/∂x = 0 alone does not exclude a mixed kinetic term. What makes the decomposition hold is the construction in which φ (internal) and x (spacetime) are independent degrees of freedom, with generators acting on separate factors and no cross term. Under this, the full Hamiltonian takes the form Ĥ = Ĥ_φ + Ĥ_x, and the Hilbert space decomposes as a tensor product ℋ_int(φ) ⊗ ℋ_space(x).

$$\hat{H} = \hat{H}_\phi \otimes \mathbf{1} + \mathbf{1} \otimes \hat{H}_x,\qquad \mathcal{H} = \mathcal{H}_\mathrm{int}(\phi) \otimes \mathcal{H}_\mathrm{space}(x)$$

Since this composite space is separable with dim ≥ 3 (ℋ_int(φ) is 3-dimensional from the generation factor alone, ℋ_space(x) is infinite-dimensional), Gleason's theorem [32] restricts every non-contextual probability assignment on all projections that is countably additive over orthogonal families to the trace form p(P) = Tr(ρP) — under those assumptions the form of the Born rule is a theorem rather than an axiom. What is unique is the form of the assignment; the single detection probability σ = ⟨n⟩ of §2 does not by itself determine the state ρ.

Two observational facts confirm this. First, all gauge couplings are independent of generation — direct evidence that φ and x are not coupled. Second, the 0D theory computes only dimensionless quantities (mass ratios, mixing angles), while dimensionful quantities (absolute masses) require spacetime — precisely the prediction of the tensor-product structure.

The analogy with spin illuminates this structure. Spin algebra is finite-dimensional and independent of the Lagrangian, yet determines 4D physical quantities (magnetic moment, fine structure). Generations have the same structure — the finite-dimensional structure V = −H corresponds to 4D mass ratios and mixing angles.

Once both factors are in place, the form of the 4D field theory is determined as well. The spacetime factor carries the local Lorentz carrier of §9.1 (spinor fields); the internal factor carries the finite data of §4 — the gauge algebra with its representations, the single Higgs doublet, and the four Yukawa edges. The first-order local operator coupling the two is narrowed to a single form compatible with the tensor-product structure and the internal data:

$$D_{A,\Phi} = i\gamma^\mu(\nabla_\mu + A_\mu) + \gamma_5 \otimes \Phi$$

where A_μ is the gauge connection and Φ contains only the single Higgs and the four Yukawa edges (the classification of §4).

The interactions then decompose uniquely as the curvature components of this single operator — the field strength F_A (gauge self-interaction), the covariant derivative D_AΦ (the gauge–Higgs coupling), and the zero-form component Φ² (whose relativisation to the vacuum point, Φ² − Φ₀², gives the Higgs potential) — whose magnitude gives the bosonic action, while ⟨Ψ, D_{A,Φ}Ψ⟩ gives the fermionic action. Each term of the Standard-Model Lagrangian is thus not a separately postulated assumption but a curvature component of a single first-order coupling.

The rigorous uniqueness of this form (the exclusion of other coupling types) and the values of the dimensionless coefficients are given, under the stated realisation conditions, by the manifestation theorem of Paper II. The electromagnetic-clock correspondence is unique, as a continuous map preserving semigroup composition, up to a single scale-setting quantity^II^. The β functions of individual renormalisation schemes and the numerical conversion between schemes are the continuing computation that displays the same physical quantity in the coordinates of each scheme, and the absolute scale (m_e) is what carries the derived quantities into physical units — the conventional layer (both in Paper II).

## 9.3. Correspondence between 0D Values and 4D Experiment

The tensor-product structure explains why 0D algebraic values coincide with tree-level values of the 4D renormalised theory. The 0D values are algebraic constants. In 4D, tree level is the limit of zero radiative corrections — i.e. the limit in which spacetime dynamics does not affect internal structure — which is precisely the independence guaranteed by the tensor product ℋ_int(φ) ⊗ ℋ_space(x). Radiative corrections arise from loop momenta in the spacetime factor ℋ_space and enter as sub-percent corrections to the internal algebraic values (§5.1: sin²θ_W = 3/13 + Δ_rad = 0.23122). The 0D values are neither RG fixed points nor UV boundary conditions, but the unique algebraically determined point serving as the reference for comparison with 4D theory — running quantities acquire the 4D RG evolution and radiative corrections, while quantities such as discrete state counts do not run.

Concretely, each observable carries a stated definition point for the comparison (α at the Thomson point; sin²θ_W and α_s at M_Z; §5), and the 0D value is compared with the renormalised 4D value at that point. Changing the scale is described on the 4D side as RG evolution; on the internal side it corresponds to reading the same finite resolvent that gave Z in §5.3 at a different evaluation point, and the closed form of this response kernel is derived^II^. The closed form of the β function in any renormalisation scheme is obtained by running between definition points.

For the Thomson-normalised physical electromagnetic effective charge, by contrast, the identification of the electromagnetic response together with the unique clock correspondence yields, as a theorem, uniform boundedness on the whole Euclidean-momentum half-line and the absence of a Landau pole^II^ — α can be followed along it from its Thomson-point value to a finite ultraviolet value. Poles appearing in perturbative displays are coordinate singularities created by a singular reparametrisation of the coupling, not poles of the physical current response.

With this, the relationship between internal structure (φ) and spacetime (x) is established. V = −H gives internal values corresponding to mass ratios, mixing angles, and coupling constants. The spacetime factor supplies, in addition to the overall energy scale, the 4D RG evolution and radiative corrections for running quantities.

# 10. Results

This chapter collects the results. Nineteen parameters are derived in this paper. A further seven (|V_cb|, |V_ub|, δ_CKM, δ_PMNS, m₁/m_e, m₃/m₂, Δm²₃₁/Δm²₂₁) require the specified operator words and the spectral structure of the flag Laplacian, and are derived in Paper II — their values are displayed in the tables of this paper, with the attribution marked in the § column.

Of the 26 Standard-Model parameters, one — in this paper the electron mass m_e — serves as the unit in which the dimensionless internal theory is displayed in human units; its specification is not a residue of the derivation but the setting of the absolute scale, the conventional layer (§9.2, Paper II). All predictions of this paper — the dimensionless quantities and the absolute masses of §7.3 alike — are displayed under this single unit m_e, and the electroweak scale v is predicted in that unit^II^.

The full list of results is given below (experimental values from PDG 2026 [22], NuFIT 6.1 [28], and CODATA 2022 [34]; the comparison value for 1/α_em alone is the Cs-side CODATA 2018 value, per the Cs/Rb dispute — §5.3, Table 4). The largest deviation among the quantities derived in this paper (sin²θ₁₃ at 2.6%) lies within about 1σ of the corresponding experimental uncertainty (§8.2). The 0.006% for the charged-lepton mass ratio R amounts to about 10σ because that quantity is measured so precisely, but the difference is one of truncation order ; the evaluation including the feedback term gives −0.7σ (§7.1). The sin²θ₂₃ comparison is against the second-octant local minimum 0.550 of the NuFIT 6.1 χ² profile; the NuFIT 6.1 global best fit lies in the first octant (0.470), the octant is unresolved, and the prediction lies within the 3σ range (§8.2). The 5.3% for δ_CKM is a genuine residual of this sector.

**Table 1.** The derived equations and their values.[^t1provenance]

| Quantity | Equation | Value |
|---|---|---|
| 1/α (leading) | physical root of 1/α + α = λ₁ + Z (λ₁ = 4 − √3, Z = 134.77613) | 137.0368 |
| 1/α (all orders) | physical root of the all-orders Dyson series (§5.3) | 137.035999084 |
| sin²θ_W | 3/13 (next correction +α·135/2197) | 0.23077 → 0.23122 |
| α_s | (3/26)(1 + ε/10) | 0.1179 |
| m_W/m_Z | √(10/13) | 0.8771 |
| m_H/v | ½√(1 + 2ε/13) | 0.5084 |
| mass ratios R (three sectors) | m_n ∝ \|E_n\|·IPR_n^{b}·(\|⟨V⟩_n\|/T_n)^{ε/N} (κ = 1, 3) | 1.529, 1.777, 2.273 |
| down-type correction | m_eff = 1 − N_c σ(1−σ) | (enters R_down) |
| sinθ_C | ε(1 + ε/9), after screening ×(1 − 9α/(26π)) | 0.2245 → 0.22435 |
| \|V_cb\| | ε²(1 − 2ε/3) | 0.0410 |
| \|V_ub\| | 2ε⁴(1 − 2ε/3)(1 − 3ε/26) | 0.00384 |
| J | πε⁵/52 | 3.06 × 10⁻⁵ |
| δ_CKM | sinδ = J/(c₁₂s₁₂c₂₃s₂₃c₁₃²s₁₃), cosδ > 0 | 62.6° |
| sin²θ₁₂ | 4/13 | 0.3077 |
| sin²θ₂₃ | (7/13)(1 + ε/13) | 0.5475 |
| sin²θ₁₃ | (1/39)(1 − 2ε/3) | 0.0219 |
| δ_PMNS | 2π·(7/13)(1 + ε/13) | 197.1° |
| m₁/m_e | x³·(7/6)(1 + ε²/9) (x = 9α/(26π)) | 6.097 × 10⁻¹⁰ |
| m₃/m₂ | 13/√5 | 5.814 |
| Δm²₃₁/Δm²₂₁ | (169/5 − t²)/(1 − t²) (t = m₁/m₂) | 33.84 |
| 0νββ | exactly zero | 0 |
| θ̄(μ₀) | zero at the matching point | 0 |

[^t1provenance]: All theoretical values in Table 1 are computed solely from the spectral and incidence data derived from V = −H and PG(2, F₃); no experimental value enters their calculation. The correspondence between the resulting mathematical structures and physical structures, and the proof of its uniqueness under the stated realisation conditions, are given in Paper II.

**Table 2.** Comparison with Standard Model constants.[^t2count]

| Parameter | Prediction | Experiment | Rel. difference | § |
|---|---|---|---|---|
| **Gauge couplings** | | | | |
| 1/α_em (→ g₁) (leading; Thomson point) | 137.0368 | 137.036 | 0.0006% | 5.3 |
| 1/α_em (all orders) | 137.035999084 | 137.035999084(21) | all digits (Cs) | II |
| α_s (→ g₃) (MS-bar, M_Z) | 0.1179 | 0.1180 | 0.07% | §5.2 |
| sin²θ_W (→ g₁/g₂) (MS-bar, M_Z) | 0.23122 | 0.23122 | <0.01% | §5.1 |
| m_W/m_Z (on-shell identification) | 0.88137 | 0.8813 | +1.0σ | II |
| **Higgs** | | | | |
| m_H/v (→ λ) | 0.5084 | 0.5082 | 0.03% | §6 |
| v (electroweak scale, GeV) | 246.2360 | 246.2196 | +0.0067% | II |
| **Fermion masses** | | | | |
| R_lepton (→ y_e, y_μ, y_τ) | 1.5293 | 1.5294 | 0.006% | §7.1 |
| Q (Koide) | 1.499985 | 1.500005(11) | 2σ | §7.1 |
| R_up (up-type hierarchy shape) | 1.777 | 1.7697(25) | +0.41% (+2.9σ) | §7.2 |
| R_up (with feedback term) | 1.767 | 1.7697(25) | −0.17% (−1.2σ) | II |
| R_down (down-type hierarchy shape) | 2.273 | 2.276 | 0.1% | §7.2 |
| **CKM / PMNS mixing** | | | | |
| sinθ_C (→ CKM θ₁₂) | 0.2245 | 0.2243 | 0.10% | §8.1 |
| \|V_cb\| | 0.0410 | 0.0407(13) | 0.8% | II |
| \|V_ub\| | 0.00384 | 0.00389(16) | 1.2% | II |
| δ_CKM | 62.6° | 66.1° | 5.3% | II |
| sin²θ₁₂ | 4/13 | 0.3088 | 0.4% | §8.2 |
| sin²θ₂₃ | 0.5475 | 0.550 (local min.; global 0.470) | 0.45% | §8.2 |
| sin²θ₁₃ | 0.0219 | 0.02249 | 2.6% | §8.2 |

[^t2count]: II in the § column denotes Paper II. The 26: three gauge couplings, λ, nine charged-fermion masses, four CKM parameters, θ_QCD, v, three neutrino masses, four PMNS parameters (of the nine charged-fermion masses, the electron mass m_e is the unit). The mass formula fixes within-sector ratios; the absolute scales of the two quark sectors are treated in Paper II. Relative differences are taken with respect to the comparison central value; rows quoted in σ compare against the experimental uncertainty.

**Table 3.** Structural results (derivations in the sections; D and the signature in §9 are consistency checks under stated premises).

| Structure | Result | § |
|---|---|---|
| Number of generations N | 3 (a fourth sequential chiral generation is excluded) | §3 |
| Spatial dimension d | 3 | §3.4 |
| Spacetime dimension D and signature | 4, (3,1) | §3.4, §9 |
| Gauge group | SU(3)×SU(2)×U(1) | §4, Paper II |
| Neutrino properties | fixed by U(1)_{B−L} and V-parity | §8.3 |
| θ̄(μ₀) | 0 | §8.1, Paper II |

**Table 4.** Testable predictions.[^frunit]

| Prediction | Value | Verification experiment | § |
|---|---|---|---|
| Δm²₃₁/Δm²₂₁ | 33.84, from the full relation (leading m₃/m₂ = 13/√5, non-zero m₁ included) | JUNO [26] (2026–27) | II |
| 1/α (all orders) | 137.035999084 — the Cs/Rb recoil discrepancy resolves on the Cs side | next-generation photon-recoil / g−2 comparisons | II |
| θ₂₃ octant | Second (>45°) | DUNE / HK (~2030) | §8.2 |
| δ_PMNS | 2π(7/13)(1+ε/13) = 197.1° | DUNE / HK (~2030) | II |
| R_ν (mass ratio; input-free) | 1.529 (= R_lepton) | osc. + absolute mass | §7.3 |
| Σm_ν | 59.37 meV (m_e as the single absolute scale, no oscillation input; ≈ NO min. 58.9) | CMB-S4 [27] (~2030) | §7.3, Paper II |
| m₁ | 0.3116 meV (m_e as the single absolute scale) | CMB-S4 (via Σm_ν) | §7.3, Paper II |
| 0νββ | Exactly zero | nEXO [29] (~2030) | §8.3 |
| Koide deviation δ | −1.5 × 10⁻⁵ | Belle II [33] (2027–30) | §7.1 |

[^frunit]: One further prediction on an open experimental dispute, of the same kind as the 1/α entry: first-row CKM unitarity is exact^II^, and the |V_ud| dispute — global fit 0.97431(16) versus the direct superallowed 0⁺→0⁺ value 0.97367(32) — resolves on the unitarity side.

The Standard Model treats these physical quantities as independent parameters to be measured and input, but the derivation in this paper has no input parameters. All numerical results above can be independently reproduced with the accompanying verification-code suite, which prints for every entry of Tables 1–4 whether it is evaluated, interval-certified, comparison-only, or not implemented here — some quantities derived in Paper II are not evaluated by this paper’s code.

# 11. Conclusion

This paper has shown that structures and numerical values normally treated as independent inputs in the Standard Model are reproduced, with no adjustable parameters, by a single equation fixed by the three axioms of existence. This agreement supports the identifications between mathematical and physical structure made throughout this paper. If so, this suggests that the fundamental laws of matter and the Universe are different aspects of the single equation

$$V = -H(\sigma(\phi))$$

and can be described as a necessity from this equation derived from the axioms of existence.

Paper II [20] treats the derivation of the seven quantities listed in §10, the spectral proof of the gauge decomposition 9+3+1, and a more detailed account of the uniqueness of the construction used in this paper. Subsequent papers extend this framework to gravity, the Universe, and the unification of forces.

*To be, or not to be — that is the equation.*

# Acknowledgements

To the information theorists, mathematicians, and physicists whose profound insights have been revealed anew through the equations of this work, I express my deepest respect and gratitude. Their legacy lives in every line of the derivation. I also thank the friends, colleagues, and family who provided constant encouragement. Above all, I thank my wife Kana, who supported this work with dedication throughout its entire duration. Thank you always.

# References

[1] ATLAS Collaboration (G. Aad et al.), Phys. Lett. B **716**, 1 (2012).

[2] CMS Collaboration (S. Chatrchyan et al.), Phys. Lett. B **716**, 30 (2012).

[3] Xing, Z.-Z. "Flavor structures of charged fermions and massive neutrinos," Phys. Rept. **854**, 1 (2020) [arXiv:1909.09610].

[4] J. A. Wheeler, "Information, physics, quantum: The search for links," in *Complexity, Entropy, and the Physics of Information* (Addison-Wesley, 1990).

[5] E. T. Jaynes, Phys. Rev. **106**, 620 (1957).

[6] B. R. Frieden, *Physics from Fisher Information* (Cambridge Univ. Press, 1998).

[7] E. Verlinde, JHEP **04**, 029 (2011) [arXiv:1001.0785].

[8] A. Caticha, Entropy **17**, 6110 (2015) [arXiv:1509.03222].

[9] L. Hardy, arXiv:quant-ph/0101012 (2001).

[10] G. Chiribella, G. M. D'Ariano, and P. Perinotti, Phys. Rev. A **84**, 012311 (2011) [arXiv:1011.6451].

[11] Pauli, W. Z. Phys. **31**, 765 (1925).

[12] Shannon, C. E. Bell Syst. Tech. J. **27**, 379 (1948).

[13] Čencov, N. N. *Statistical Decision Rules and Optimal Inference* (AMS, 1982).

[14] Amari, S. *Information Geometry and Its Applications* (Springer, 2016).

[15] Ay, N., Jost, J., Lê, H. V. & Schwachhöfer, L. *Information Geometry* (Springer, 2017).

[16] Aczél, J. *Lectures on Functional Equations and Their Applications* (Academic Press, 1966).

[17] Khinchin, A. I. *Mathematical Foundations of Information Theory* (Dover, 1957).

[18] Pöschl, G. & Teller, E. Z. Phys. **83**, 143 (1933).

[19] Cabibbo, N. Phys. Rev. Lett. **10**, 531 (1963).

[20] H. Kondo, "Physics from Existence II: The Proofs", in preparation.

[21] Weinberg, S. Phys. Rev. Lett. **19**, 1264 (1967).

[22] Particle Data Group (F. Takahashi et al.), to be published in Int. J. Mod. Phys. A **41**, 2630011 (2026); 2026 web edition, https://pdg.lbl.gov/2026/

[23] Koide, Y. Lett. Nuovo Cim. **34**, 201 (1982).

[24] Kobayashi, M. & Maskawa, T. Prog. Theor. Phys. **49**, 652 (1973).

[25] Wolfenstein, L. Phys. Rev. Lett. **51**, 1945 (1983).

[26] JUNO Collaboration (A. Abusleme et al.), Prog. Part. Nucl. Phys. **123**, 103927 (2022) [arXiv:2104.02565].

[27] CMB-S4 Collaboration (K. N. Abazajian et al.), arXiv:1610.02743 (2016).

[28] Esteban, I. et al. "NuFit-6.0: updated global analysis of three-flavor neutrino oscillations," JHEP **12** (2024) 216 [arXiv:2410.05380]; NuFIT 6.1 (2025), http://www.nu-fit.org/

[29] nEXO Collaboration (J. B. Albert et al.), Phys. Rev. C **97**, 065503 (2018).

[30] Lovelock, D. J. Math. Phys. **12**, 498 (1971).

[31] Jacobson, T. Phys. Rev. Lett. **75**, 1260 (1995).

[32] Gleason, A. M. J. Math. Mech. **6**, 885 (1957).

[33] Belle II Collaboration (E. Kou, P. Urquijo et al.), Prog. Theor. Exp. Phys. **2019**, 123C01 (2019).

[34] P. J. Mohr et al., J. Phys. Chem. Ref. Data **54**, 033105 (2025); CODATA 2022 recommended values.

[35] A. S. Eddington, *Fundamental Theory* (Cambridge Univ. Press, 1946).

[36] A. Wyler, C. R. Acad. Sci. Paris A **269**, 743 (1969).

[37] U. D. Jentschura and I. Nándori, "Attempts at a determination of the fine-structure constant from first principles: A brief historical overview," Eur. Phys. J. H **39**, 591 (2014) [arXiv:1411.4673].

[38] Baez, J. and Huerta, J. "The algebra of grand unified theories," Bull. Amer. Math. Soc. **47**, 483 (2010).

[39] Connes, A. "Gravity coupled with matter and the foundation of non-commutative geometry," Commun. Math. Phys. **182**, 155 (1996).

[40] Chamseddine, A. H. and Connes, A. "Universal Formula for Noncommutative Geometry Actions: Unification of Gravity and the Standard Model," Phys. Rev. Lett. **77**, 4868 (1996).

[41] Chamseddine, A. H., Connes, A. and Marcolli, M. "Gravity and the standard model with neutrino mixing," Adv. Theor. Math. Phys. **11**, 991 (2007).

[42] Furey, C. Phys. Lett. B **785**, 84 (2018) [arXiv:1910.08395].

[43] Sumino, Y. Phys. Lett. B **671**, 477 (2009).
