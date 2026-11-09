## ============================================================================
## lib/FP_finder_v6.R
##
## Fixed-point finder for lib/ODE_SYSTEM_v6.R, Tregs OFF or ON (theta_treg = 0/1),
## deterministic (C_rat_perceived = C_rat), full stability (phi_instability = 0).
## Both the sterile and pathogenic finders now solve the Treg subsystem when
## theta_treg = 1 (sterile via sterile_reduction, pathogenic via the Tr,S slaving
## in pathogenic_reduction); theta_treg = 0 recovers the Tregs-off system exactly.
##
## Provides two finders:
##   find_fixed_points_sterile(...)     PLP = PL = 0.  Reduces to ONE scalar
##                                      equation G(L) = 0 scanned on a bounded
##                                      interval -> the count is complete.
##   find_fixed_points_pathogenic(...)  PLP, PL live.  Reduces to a bounded
##                                      2-D root-find in (M1, PLP), solved by
##                                      grid-seeded Newton.
##
## Both label each fixed point "stable" or "saddle/unstable" using a REDUCED
## Jacobian.  The full Jacobian carries two structural zero eigenvalues from the
## conservation laws  M0 + M1 = M_tot  and  T0 + Tr = T_tot; those zeros mask
## the sign of the real dynamical eigenvalues.  We remove them by slaving
## M0 = M_tot - M1 and T0 = T_tot - Tr and differentiating only the independent
## variables.
##
## The method is described in docs/methods.pdf, section "Fixed points and their
## stability".
##
## This file assumes make_ode_system() (from lib/ODE_SYSTEM_v6.R) is already loaded;
## it sources it defensively if not.
## ============================================================================

if (!exists("make_ode_system")) source("./lib/ODE_SYSTEM_v6.R")

## STATE_ORDER of the RHS output (CumC is a quadrature, not a real state).
FP_STATE_ORDER = c("PL","PLP","CLP","Eh","Ea","En","M0","M1","R","T0","Tr","D","S")

## ----------------------------------------------------------------------------
## Fixed structural constants.  These MUST match params_base in
## pipeline/00_parameter_grids.R.  M_tot and T_tot are NOT here: they are the
## conserved pools set by build_state() (M0(0) = K_M, M1(0) = 0  =>  M_tot = K_M;
## T0(0) = K_T, Tr(0) = 0  =>  T_tot = K_T).  They are passed in by the driver.
## ----------------------------------------------------------------------------
FP_FIXED = list(
  K_E = 100, K_M = 100, K_T = 100, K_R = 1,
  K_L_P = 1000, K_LP_P = 100, K_LP_C = 100, CL = 1000,
  a0 = 0.0045, a_floor = 0.0002, n_g = 4
)

## ----------------------------------------------------------------------------
## Central-difference Jacobian (avoids a numDeriv dependency on the cluster).
## f: a function R^n -> R^m.  Returns the m-by-n Jacobian at x.
## ----------------------------------------------------------------------------
fp_jacobian = function(f, x, eps = 1e-6) {
  n = length(x)
  m = length(f(x))
  J = matrix(0, m, n)
  for (j in seq_len(n)) {
    h = eps * max(abs(x[j]), 1e-4)
    x_plus = x; x_plus[j] = x[j] + h
    x_minus = x; x_minus[j] = x[j] - h
    J[, j] = (f(x_plus) - f(x_minus)) / (2 * h)
  }
  J
}

## ----------------------------------------------------------------------------
## The efferocytosis switch g(L): a DECREASING Hill function of corpse load L.
## Falls from a0 (clearance intact) to a_floor (clearance collapsed).
## ----------------------------------------------------------------------------
fp_g = function(L, Kg) {
  L = pmax(L, 0)
  FP_FIXED$a_floor +
    (FP_FIXED$a0 - FP_FIXED$a_floor) * Kg^FP_FIXED$n_g / (Kg^FP_FIXED$n_g + L^FP_FIXED$n_g)
}

## ----------------------------------------------------------------------------
## The DAMP saturation constant K_En (needs K_Ea first).  Verbatim from the RHS.
## ----------------------------------------------------------------------------
fp_K_En = function(par) {
  K_Ea = (FP_FIXED$K_E * par$tau_r * FP_FIXED$K_R) / (par$mu_Ea + FP_FIXED$a0 * FP_FIXED$K_M)
  (FP_FIXED$K_E * par$tau_P_B + par$mu_Ea * K_Ea) / (FP_FIXED$a0 * FP_FIXED$K_M)
}

## ----------------------------------------------------------------------------
## Build the parameter list the RHS expects, for the deterministic Jacobian.
## ----------------------------------------------------------------------------
fp_rhs_params = function(par, Kg, beta_ros, r_L_P, theta_treg = 0) {
  c(par, FP_FIXED,
    list(Kg = Kg, beta_ros = beta_ros, theta_treg = theta_treg,
         r_L_P = r_L_P, r_LP_P = r_L_P, phi_instability = 0))
}

## ----------------------------------------------------------------------------
## A deterministic RHS (noise switched off: C_rat_perceived = C_rat).
## ----------------------------------------------------------------------------
fp_make_rhs = function() {
  make_ode_system(list(use_crat_fun = function(t) 1, u_fun = function(t) 0))
}

## ============================================================================
## ============================================================================
## SECTION 1 -- STERILE (PLP = PL = 0), Tregs off
## ============================================================================
## ============================================================================
##
## At a sterile Tregs-off equilibrium everything can be written as a function of
## the corpse load L = Ea + En:
##
##   g(L)   = fp_g(L)
##   A_req  = M_tot * L * g(L)                         (the damage drive)
##   M1     = A_req * kappa_E / ( aa * (kappa_E*K_E - A_req) ),  aa = tau_r*beta_ros/mu_R
##   R,Eh,Ea,En,CLP,D  follow in closed form.
##   G(L)   = (dM1/dt at that state)                   -> zero at a fixed point.
##
## Fixed points are the roots of the scalar G(L).  We scan L on a bounded
## interval and bracket every sign change: nothing can be missed.
## ============================================================================

## Build the closures (state-at-L, residual G, full-state reconstruction) for
## one sterile parameter set.  Kept together so the algebra reads top to bottom.
sterile_reduction = function(par, Kg, beta_ros, M_tot, T_tot = 0, theta_treg = 0) {

  aa   = par$tau_r * beta_ros / par$mu_R      # so that tau_r * R = aa * M1
  Amax = par$kappa_E * FP_FIXED$K_E           # A(M1) saturates at this value
  K_En = fp_K_En(par)

  ## Reconstruct the full reduced state at a given corpse load L.
  ## Returns NULL if this L cannot be realised by any admissible M1.
  state_at = function(L) {
    g_L   = fp_g(L, Kg)
    A_req = M_tot * L * g_L
    if (A_req >= Amax) return(NULL)                      # would need M1 = Inf
    M1 = A_req * par$kappa_E / (aa * (Amax - A_req))
    if (M1 < 0 || M1 > M_tot) return(NULL)               # M0 = M_tot - M1 >= 0

    M0 = M_tot - M1
    R  = beta_ros * M1 / par$mu_R
    Eh = par$kappa_E * FP_FIXED$K_E / (par$kappa_E + aa * M1)
    Ea = A_req / (par$mu_Ea + g_L * M_tot)
    En = par$mu_Ea * Ea / (g_L * M_tot)                  # Ea + En = L by construction

    ## commensals: dCLP = 0 is linear in CLP
    clearance = par$alpha_0 * M0 + par$alpha_1_e * M1 + par$alpha_1_r * R
    source_C  = FP_FIXED$CL * par$gamma_B_C * (1 - Eh / FP_FIXED$K_E)
    CLP = source_C / (clearance + source_C / FP_FIXED$K_LP_C)

    ## danger signal D at steady state
    f_En = En / (K_En + En)
    f_PC = CLP / (FP_FIXED$K_LP_P + FP_FIXED$K_LP_C)      # sterile: PLP = 0
    D = (par$sigma_D / par$tau_D) * (par$w_E_over_M * f_En + (1 - par$w_E_over_M) * f_PC)

    ## Tregs (theta_treg = 1): sterile => C_rat = 1 (no pathogen), so dTr = 0 and
    ## dS = 0 give closed forms slaved to M1.  Both are 0 when Tregs are off.
    Tr = if (theta_treg == 1)
           par$tau_t_a * M1 * T_tot / (par$tau_t_d + par$tau_t_a * M1) else 0
    S  = if (theta_treg == 1) par$sigma_S * Tr / par$tau_S else 0

    list(L = L, g_L = g_L, M1 = M1, M0 = M0, R = R, Eh = Eh,
         Ea = Ea, En = En, CLP = CLP, D = D, Tr = Tr, S = S)
  }

  ## The macrophage-balance residual: this is dM1/dt, the one condition NOT
  ## already enforced by state_at.  Its roots are the FPs.  With Tregs on it
  ## gains the -tau_10*S*M1 suppression term (S = 0 when Tregs are off).
  G = function(L) {
    s = state_at(L)
    if (is.null(s)) return(NA_real_)
    par$tau_01 * s$D * (M_tot - s$M1) -
      par$tau_10 * s$S * s$M1 -
      par$tau_1_d * s$M1 -
      par$tau_eff * s$g_L * s$Ea * s$M1
  }

  ## Assemble the full 14-long state vector (STATE_ORDER + CumC) for the Jacobian.
  full_state = function(s) {
    c(PL = 0, PLP = 0, CLP = s$CLP, Eh = s$Eh, Ea = s$Ea, En = s$En,
      M0 = s$M0, M1 = s$M1, R = s$R, T0 = T_tot - s$Tr, Tr = s$Tr, D = s$D, S = s$S, CumC = 0)
  }

  list(state_at = state_at, G = G, full_state = full_state,
       Amax = Amax, M_tot = M_tot, T_tot = T_tot)
}

## Independent variables for the sterile stability Jacobian.  M0, T0 are slaved
## by conservation; PL = PLP = 0 holds the sterile subspace; so they are excluded.
STERILE_INDEP    = c("CLP","Eh","Ea","En","M1","R","D","Tr","S")
PATHOGENIC_INDEP = c("CLP","PLP","Eh","Ea","En","M1","R","D","Tr","S")

## Classify one sterile fixed point as stable / saddle via the reduced Jacobian.
classify_sterile = function(s, red, rhs_params) {
  rhs = fp_make_rhs()
  M_tot = red$M_tot; T_tot = red$T_tot
  y_full = red$full_state(s)                 # named length-14
  y_indep = y_full[STERILE_INDEP]

  ## reduced vector field on the independent variables only
  reduced_rhs = function(y) {
    y = setNames(y, STERILE_INDEP)
    full = y_full
    full[STERILE_INDEP] = y
    full["M0"] = M_tot - y["M1"]             # conservation law 1
    full["T0"] = T_tot - y["Tr"]             # conservation law 2
    full["PL"] = 0; full["PLP"] = 0          # sterile subspace
    full["CumC"] = 0
    dy = rhs(0, full, rhs_params)[[1]]
    names(dy) = c(FP_STATE_ORDER, "CumC")
    dy[STERILE_INDEP]
  }

  eig = eigen(fp_jacobian(reduced_rhs, y_indep), only.values = TRUE)$values
  max_re = max(Re(eig))
  type = if (max_re < -1e-8) "stable" else if (max_re > 1e-8) "saddle/unstable" else "marginal"
  list(max_re = max_re, type = type)
}

## ----------------------------------------------------------------------------
## PUBLIC: find every sterile Tregs-off fixed point for one parameter set.
## Returns a data.frame, one row per fixed point.
## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
## PUBLIC: find every sterile Tregs-off fixed point for one parameter set.
##
## TWO CORRECTIONS relative to the previous version this replaces, both of which cost
## real results before they were found:
##
##  1. THE SCAN STARTS AT L = 0, not 1e-4. The sterile HEALTHY equilibrium sits at
##     L = Ea + En = 0 EXACTLY -- on the boundary -- so the old domain excluded it
##     by construction and the finder returned NOTHING for every system that heals
##     completely. That was 1,353 of 1,434 excluded sterile combos in one batch,
##     and it hit both seeding paths because grid_info() and ray_one() call this
##     same routine. L = 0 is an exact root: at Eh = K_E the (1 - Eh/K_E) factor
##     kills the CLP source, so CLP = D = M1 = R = Ea = En = 0 and every term of
##     G is identically zero -- in exact arithmetic and in floating point.
##     A second sub-scan covers (0, dL] because G(0) = 0 makes the first cell's
##     sign-change test 0 * G(dL), which can never fire.
##
##  2. IT NOW ACCEPTS extra_seeds, the same endpoint matrix ray_one() passes to
##     the pathogenic finder. The scan is complete only to grid resolution; a
##     root the integration demonstrably sits on must not be missed because two
##     crossings shared a cell. Seeds can only ADD roots, never move or remove
##     one. n_seed_recovered is returned as an attribute -- if it is ever
##     non-zero, the scan resolution is the thing to look at.
##
## Returns a data.frame, one row per fixed point.
## ----------------------------------------------------------------------------
find_fixed_points_sterile = function(par, Kg, beta_ros, M_tot,
                                        T_tot = 0, theta_treg = 0,
                                        extra_seeds = NULL, ...) {

  red = sterile_reduction(par, Kg, beta_ros, M_tot, T_tot, theta_treg)

  L_hi   = (red$Amax / (M_tot * FP_FIXED$a_floor)) * 1.2
  n_scan = max(6000, round(L_hi * 12))
  ## ---- THE FIX: the grid starts at 0, not at 1e-4 -------------------------
  L_grid = seq(0, L_hi, length.out = n_scan)
  ## ------------------------------------------------------------------------
  G_vals = vapply(L_grid, red$G, numeric(1))
  finite = is.finite(G_vals)

  G_ACCEPT = 1e-6
  roots = numeric(0)
  accept = function(a, b) {
    u = tryCatch(uniroot(red$G, c(a, b), tol = .Machine$double.eps), error = function(e) NULL)
    if (!is.null(u) && is.finite(u$f.root) && abs(u$f.root) <= G_ACCEPT)
      roots[[length(roots) + 1L]] <<- u$root
  }

  ## Boundary root, checked explicitly rather than left to the bracketing loop:
  ## a root sitting exactly ON the left endpoint is not a sign change, so the
  ## interior logic below can only catch it via the G_vals[i] == 0 branch, which
  ## is an exact-equality test.  Belt and braces -- accept it on the same
  ## |G| <= G_ACCEPT criterion every other root has to pass.
  if (is.finite(G_vals[1]) && abs(G_vals[1]) <= G_ACCEPT &&
      !is.null(red$state_at(L_grid[1]))) roots = c(roots, L_grid[1])

  ## ...and a SECOND root can sit between 0 and the first interior grid point.
  ## That interval is invisible to the loop below for two compounding reasons:
  ##   * G(0) = 0 exactly, so the sign-change test is 0 * G(dL), never < 0;
  ##   * the G_vals[i] == 0 branch `next`s past the interval entirely.
  ## Caught in regression: 2 of 300 random combos have their only non-zero root
  ## at Eh ~ 99.8-99.9, i.e. L well inside the first cell (dL = L_hi/(n_scan-1)
  ## ~ 0.1), and the old 1e-4-start grid DID bracket them. Sub-scan the cell
  ## log-spaced so a root arbitrarily close to 0 is still bracketed. ~400 extra
  ## G evaluations, against ~6000-8000 for the main grid.
  if (length(L_grid) > 1L && L_grid[1] == 0 && L_grid[2] > 0) {
    sub = 10^seq(-10, log10(L_grid[2]), length.out = 400)
    gs  = vapply(sub, red$G, numeric(1))
    for (k in seq_len(length(gs) - 1L))
      if (is.finite(gs[k]) && is.finite(gs[k + 1L]) && gs[k] * gs[k + 1L] < 0)
        accept(sub[k], sub[k + 1L])
  }

  for (i in seq_len(length(G_vals) - 1L)) {
    fa = finite[i]; fb = finite[i + 1L]
    if (fa && fb) {
      if (G_vals[i] == 0) { roots = c(roots, L_grid[i]); next }
      if (G_vals[i] * G_vals[i + 1L] < 0) accept(L_grid[i], L_grid[i + 1L])
    } else if (fa != fb) {
      finL = if (fa) L_grid[i] else L_grid[i + 1L]
      naL  = if (fa) L_grid[i + 1L] else L_grid[i]
      lo = finL; hi = naL
      for (b in 1:60) { mid = (lo + hi) / 2; if (is.finite(red$G(mid))) lo = mid else hi = mid }
      sub = seq(min(finL, lo), max(finL, lo), length.out = 200)
      gs  = vapply(sub, red$G, numeric(1))
      for (k in seq_len(length(gs) - 1L))
        if (is.finite(gs[k]) && is.finite(gs[k + 1L]) && gs[k] * gs[k + 1L] < 0)
          accept(sub[k], sub[k + 1L])
    }
  }
  ## --------------------------------------------------------------------------
  ## ENDPOINT CROSS-CHECK.  The scan above is complete IN PRINCIPLE -- it brackets
  ## every sign change of G on [0, L_hi] -- but only to the resolution of the
  ## grid.  Two roots inside one cell (an even number of crossings) are invisible
  ## to it, and so is any root the bracketing logic mishandles.  The L = 0 bug was
  ## precisely that failure mode: the integration sat on a root the scan could not
  ## see, and nothing in the finder noticed the contradiction.
  ##
  ## extra_seeds is the same matrix ray_one() passes to the pathogenic finder: one
  ## row per integrated endpoint, with L = Ea + En.  Each seed is a place the
  ## dynamics demonstrably went, so a root had better be there.  For each seed we
  ##   (a) skip it if the scan already found a root at that L;
  ##   (b) accept the seed itself if |G| is already at root level;
  ##   (c) otherwise sub-scan a fine local window around it.
  ## This can only ADD roots, never move or remove one, so it cannot change any
  ## result the scan alone would have produced.  n_seed_recovered is attached to
  ## the output as an attribute: if it is ever non-zero the scan resolution needs
  ## looking at, not just the seed.
  ## --------------------------------------------------------------------------
  n_seed_rec = 0L
  if (!is.null(extra_seeds)) {
    sm = as.matrix(extra_seeds)
    Ls = if ("L" %in% colnames(sm)) sm[, "L"] else numeric(0)
    Ls = unique(Ls[is.finite(Ls) & Ls >= 0])
    for (Lseed in Ls) {
      if (length(roots) && min(abs(roots - Lseed)) <= 1e-3 * (1 + Lseed)) next
      gseed = tryCatch(red$G(Lseed), error = function(e) NA_real_)
      if (is.finite(gseed) && abs(gseed) <= G_ACCEPT &&
          !is.null(red$state_at(Lseed))) {
        roots = c(roots, Lseed); n_seed_rec = n_seed_rec + 1L; next }
      w   = max(1e-6, 0.05 * max(Lseed, 1))
      sub = seq(max(0, Lseed - w), Lseed + w, length.out = 400)
      gv  = vapply(sub, red$G, numeric(1))
      before = length(roots)
      for (k in seq_len(length(gv) - 1L))
        if (is.finite(gv[k]) && is.finite(gv[k + 1L]) && gv[k] * gv[k + 1L] < 0)
          accept(sub[k], sub[k + 1L])
      if (length(roots) > before) n_seed_rec = n_seed_rec + 1L
    }
  }

  roots = sort(unique(round(roots, 8)))
  if (length(roots) == 0) return(structure(data.frame(), n_seed_recovered = n_seed_rec))

  rhs_params = fp_rhs_params(par, Kg, beta_ros, r_L_P = 0, theta_treg = theta_treg)
  rows = lapply(roots, function(L) {
    s = red$state_at(L)
    if (is.null(s)) return(NULL)
    cl = classify_sterile(s, red, rhs_params)
    data.frame(M1 = s$M1, PLP = 0, PL = 0, L = s$L,
               Eh = s$Eh, Ei = FP_FIXED$K_E - s$Eh,
               Ea = s$Ea, En = s$En, R = s$R, D = s$D, CLP = s$CLP,
               Tr = s$Tr, S = s$S,
               max_Re_eig = cl$max_re, type = cl$type,
               theta_treg = theta_treg, branch = "sterile",
               stringsAsFactors = FALSE)
  })
  structure(do.call(rbind, rows), n_seed_recovered = n_seed_rec)
}

## ============================================================================
## ============================================================================
## SECTION 2 -- PATHOGENIC (PLP, PL live), Tregs off OR on (theta_treg = 0/1)
## ============================================================================
## ============================================================================
##
## theta_treg couples in only through the macrophage balance: Tr and S are
## slaved to (M1, C_rat) at steady state (dTr=dS=0), and S adds a -tau_10*S*M1
## suppression to dM1.  theta_treg=0 => Tr=S=0 and the reduction is identical to
## the original Tregs-off system.  The stability Jacobian carries the Tr,S
## directions (PATHO_INDEP) and is evaluated with the matching theta_treg.
##
## Now Eh depends on PLP, and En has an extra necrotic source from PLP.  The
## inner corpse solve keeps its form  L*g(L) = kappa_E*(K_E - Eh)/M_tot, but Eh
## is now a function of (M1, PLP).  The problem is a bounded 2-D root-find in
## (M1, PLP); we carry L as a third unknown so the efferocytosis fold is handled
## naturally, and slave PL by its own branch (dPL = 0 gives PL = 0 or a positive
## logistic-balance value).  Solved by grid-seeded Newton over both PL branches.
## The sterile fixed points reappear here as the PLP = 0 solutions.
## ============================================================================

## Build the reduced-quantity and residual closures for one pathogenic set.
pathogenic_reduction = function(par, Kg, beta_ros, r_L_P, M_tot, T_tot = 0, theta_treg = 0) {

  aa   = par$tau_r * beta_ros / par$mu_R
  K_En = fp_K_En(par)
  r_LP_P = if (!is.null(par$r_LP_P)) par$r_LP_P else r_L_P   # equal per build_params()

  ## Reduced state at (M1, PLP, L) for the chosen PL branch ("zero" or "pos").
  quant = function(M1, PLP, L, pl_branch) {
    if (M1 < 0 || M1 > M_tot || PLP < 0 || L <= 0) return(NULL)
    M0 = M_tot - M1
    R  = beta_ros * M1 / par$mu_R
    ## Eh now feels the pathogen load PLP directly (extra term in dEh)
    Eh = par$kappa_E * FP_FIXED$K_E /
         (par$kappa_E + aa * M1 + par$tau_P_B * PLP / FP_FIXED$K_LP_P)
    g_L = fp_g(L, Kg)
    A1  = Eh * aa * M1                        # apoptotic drive
    Ea  = A1 / (par$mu_Ea + g_L * M_tot)
    En  = L - Ea                              # Ea + En = L (holds at the solution)

    clearance = par$alpha_0 * M0 + par$alpha_1_e * M1 + par$alpha_1_r * R
    source_C  = FP_FIXED$CL * par$gamma_B_C * (1 - Eh / FP_FIXED$K_E)
    CLP = source_C / (clearance + source_C / FP_FIXED$K_LP_C)

    ## Tregs slaved to (M1, C_rat) at steady state.  From the RHS:
    ##   dTr = theta_treg*tau_t_a*(M1*C_rat)*T0 - tau_t_d*Tr,   T0 = T_tot - Tr
    ##   dS  = sigma_S*Tr - tau_S*S
    ## Setting both to 0 gives the closed forms below.  theta_treg = 0 (or
    ## T_tot = 0) => Tr = S = 0, which recovers the Tregs-off reduction EXACTLY,
    ## so the pathogenic-off results are unchanged bit-for-bit.  C_rat is the
    ## geometric commensal fraction f_C/(f_P+f_C) the deterministic RHS perceives.
    f_P_treg = PLP / FP_FIXED$K_LP_P
    f_C_treg = CLP / FP_FIXED$K_LP_C
    C_rat = if ((f_P_treg + f_C_treg) > 0) f_C_treg / (f_P_treg + f_C_treg) else 0
    a_treg = theta_treg * par$tau_t_a * M1 * C_rat
    Tr = a_treg * T_tot / (par$tau_t_d + a_treg)
    S  = par$sigma_S * Tr / par$tau_S

    f_En = En / (K_En + En)
    f_PC = (PLP + CLP) / (FP_FIXED$K_LP_P + FP_FIXED$K_LP_C)
    D = (par$sigma_D / par$tau_D) * (par$w_E_over_M * f_En + (1 - par$w_E_over_M) * f_PC)

    ## dPL = 0 has two branches:
    PL = if (pl_branch == "zero") 0
         else FP_FIXED$K_L_P * (1 - par$alpha_AMP * D / r_L_P)

    list(M0 = M0, R = R, Eh = Eh, g_L = g_L, Ea = Ea, En = En,
         clearance = clearance, CLP = CLP, D = D, PL = PL,
         Tr = Tr, S = S, C_rat = C_rat)
  }

  ## The three residuals (corpse balance, macrophage balance, LP-pathogen balance).
  ## Their simultaneous zero is a fixed point.
  resid = function(x, pl_branch) {
    M1 = x[1]; PLP = x[2]; L = x[3]
    s = quant(M1, PLP, L, pl_branch)
    if (is.null(s) || s$En < 0 || s$PL < 0) return(rep(1e6, 3))

    r_inner = L * s$g_L * M_tot - par$kappa_E * (FP_FIXED$K_E - s$Eh)
    r_M1 = par$tau_01 * s$D * (M_tot - M1) -
           par$tau_10 * s$S * M1 -
           par$tau_1_d * M1 -
           par$tau_eff * s$g_L * s$Ea * M1
    r_PLP = s$PL * par$gamma_B_P * (1 - s$Eh / FP_FIXED$K_E) * (1 - PLP / FP_FIXED$K_LP_P) +
            PLP * r_LP_P * (1 - PLP / FP_FIXED$K_LP_P) -
            PLP * s$clearance
    c(r_inner, r_M1, r_PLP)
  }

  list(quant = quant, resid = resid, M_tot = M_tot, T_tot = T_tot, r_LP_P = r_LP_P)
}

## Damped Newton for the 3-equation residual system.
##
## ACCEPT_TOL: the residual at which a point counts as a fixed point. The final
## line of this function has always used 1e-6 for that; the constant just names
## it so the line search can apply the SAME standard (see below).
##
## BUG FIXED (endpoint seeding): the line search demands a STRICTLY smaller
## residual. When x0 is already essentially ON a root -- exactly the case when
## Newton is seeded from a converged trajectory endpoint (grid_seed = FALSE, as
## the whole cluster sweep runs) -- max(abs(fx)) sits between tol (1e-9) and
## ACCEPT_TOL, the finite-difference Jacobian at eps = 1e-6 is dominated by
## round-off, and no lambda in 25 halvings improves anything. The old code then
## did `return(NULL)`, discarding a point that was already converged, and
## jumping over the 1e-6 acceptance test on the last line. Symptom: the finder
## returns 0 fixed points for a combo whose endpoints sit numerically on top of
## roots that grid seeding finds without trouble (row_ind 527: endpoint
## M1 = 70.5393 vs root M1 = 70.539255). Now a non-improving step returns x when
## x already satisfies ACCEPT_TOL, and NULL only when it genuinely does not.
ACCEPT_TOL = 1e-6
fp_newton = function(f, x0, tol = 1e-9, max_iter = 80) {
  x = x0
  for (iter in seq_len(max_iter)) {
    fx = f(x)
    if (any(!is.finite(fx))) return(NULL)
    if (max(abs(fx)) < tol) return(x)
    ## Every "cannot proceed" exit below applies the SAME rule: if x already
    ## satisfies ACCEPT_TOL it IS a fixed point, so return it; only give up when
    ## the residual is genuinely too large. The pathogenic reduction is singular
    ## as PLP -> 0 (cond ~1e17 at a pathogen-free attractor), so a converged
    ## endpoint seed reaches `solve()` with residual ~1e-9..1e-8 -- above tol,
    ## below ACCEPT_TOL -- and the old code discarded it.
    J = fp_jacobian(f, x)
    if (any(!is.finite(J))) return(if (max(abs(fx)) < ACCEPT_TOL) x else NULL)
    step = tryCatch(solve(J, fx), error = function(e) NULL)
    if (is.null(step)) return(if (max(abs(fx)) < ACCEPT_TOL) x else NULL)
    ## backtracking line search: shrink the step until the residual decreases
    lambda = 1; improved = FALSE
    for (bt in 1:25) {
      x_try = x - lambda * step
      f_try = tryCatch(f(x_try), error = function(e) rep(Inf, length(x)))
      if (all(is.finite(f_try)) && max(abs(f_try)) < max(abs(fx))) { improved = TRUE; break }
      lambda = lambda / 2
    }
    ## cannot improve: accept x if it already IS a root, otherwise give up
    if (!improved) return(if (max(abs(fx)) < ACCEPT_TOL) x else NULL)
    x = x - lambda * step
  }
  if (max(abs(f(x))) < ACCEPT_TOL) x else NULL
}

## Independent variables for the pathogenic stability Jacobian.  PL and PLP are
## LIVE now (kept in); only M0, T0 are slaved out by conservation.
PATHO_INDEP = c("PL","PLP","CLP","Eh","Ea","En","M1","R","D","Tr","S")

## Near-marginal stability confirmation. When the linear leading eigenvalue is within
## MARGINAL_TOL of zero, the finite-difference Jacobian cannot resolve its sign (round-off at
## ~1e-6), so a high-En / PLP~0 attractor that a trajectory settles onto can read as a
## marginally-positive saddle. We confirm by a short forward SETTLE TEST: perturb the full
## state slightly and integrate; if it stays in a small neighbourhood of the fixed point it
## is (marginally) stable, if it departs it is a saddle. Needs deSolve; NA => keep "marginal".
## TOL is tight (1e-4) so the test fires ONLY on the genuinely-near-zero attractors (~1e-6);
## clear saddles like the pathogen-free twins (~1e-3) stay saddle by the linear test, no
## integration -- this keeps the per-combo cost low now that the gate classifies every combo.
MARGINAL_TOL = 1e-4

## hmax IS NOT OPTIONAL HERE EITHER, AND ITS ABSENCE WAS A DEFECT. This call ran
## with no hmax on a two-point output grid to T_end = 20000, so deSolve set
## hmax = max(diff(times)) = 20000, i.e. no cap at all. That is precisely the
## configuration pipeline/01_sweep.R documents at its own ode() call as making lsoda
## stride the stiff tail until the local error blows up: PLP crosses zero where
## the logistic term is quadratic and negative and escapes to -Inf (measured
## PLP = -1.5e156), and the failure is not always loud -- it can return finite
## numbers with Eh above K_E, or silently settle in a different basin. The verdict
## this function returns decides stable against saddle for every near-marginal
## fixed point, and that propagates into n_field_stable_*, `stab` in ray_one(),
## n_ray, and hence every class label the sweep emits. The settings come from
## lib/extinction_floor.R when the caller has loaded it, from the environment
## otherwise, so the finder still stands alone. Resolved at CALL time, not at
## source time: this file is sourced before extinction_floor.R in both callers.
fp_settle_stable = function(y_full, rhs, rhs_params, T_end = 20000, pert = 1e-3) {
  if (!requireNamespace("deSolve", quietly = TRUE)) return(NA)
  .hmax = if (exists("HMAX")) HMAX else as.numeric(Sys.getenv("HMAX", "5"))
  .atol = if (exists("ATOL")) ATOL else as.numeric(Sys.getenv("ATOL", "1e-14"))
  .rtol = if (exists("RTOL")) RTOL else as.numeric(Sys.getenv("RTOL", "1e-8"))
  keys  = FP_STATE_ORDER                                  # 13 real states (exclude CumC)
  scale = abs(y_full) + 1
  y0 = y_full; y0[keys] = y_full[keys] + pert * scale[keys]
  odef = function(t, y, p) list(rhs(t, setNames(y, names(y_full)), p)[[1]])
  out = tryCatch(as.data.frame(deSolve::ode(y0, c(0, T_end), odef, rhs_params,
                                            method = "lsoda", maxsteps = 1e5,
                                            hmax = .hmax, atol = .atol, rtol = .rtol)),
                 error = function(e) NULL)
  if (is.null(out) || nrow(out) < 2) return(NA)
  yT = unlist(out[nrow(out), keys])
  d0 = sqrt(sum(((y0[keys] - y_full[keys]) / scale[keys])^2))
  dT = sqrt(sum(((yT       - y_full[keys]) / scale[keys])^2))
  dT <= 5 * d0                                            # stayed near -> stable; departed -> saddle
}

## Classify one pathogenic fixed point as stable / saddle.
classify_pathogenic = function(row, red, rhs_params) {
  rhs = fp_make_rhs()
  M_tot = red$M_tot; T_tot = red$T_tot
  s = red$quant(row$M1, row$PLP, row$L, row$branch)
  y_full = c(PL = s$PL, PLP = row$PLP, CLP = s$CLP, Eh = s$Eh, Ea = s$Ea, En = s$En,
             M0 = s$M0, M1 = row$M1, R = s$R, T0 = T_tot - s$Tr, Tr = s$Tr,
             D = s$D, S = s$S, CumC = 0)
  y_indep = y_full[PATHO_INDEP]

  reduced_rhs = function(y) {
    y = setNames(y, PATHO_INDEP)
    full = y_full
    full[PATHO_INDEP] = y
    full["M0"] = M_tot - y["M1"]
    full["T0"] = T_tot - y["Tr"]
    full["CumC"] = 0
    dy = rhs(0, full, rhs_params)[[1]]
    names(dy) = c(FP_STATE_ORDER, "CumC")
    dy[PATHO_INDEP]
  }

  eig = eigen(fp_jacobian(reduced_rhs, y_indep), only.values = TRUE)$values
  max_re = max(Re(eig))
  type =
    if (max_re < -1e-8)          "stable"
    else if (max_re > MARGINAL_TOL) "saddle/unstable"
    else {                        # near-marginal: linear sign unreliable -> forward settle test
      st = fp_settle_stable(y_full, fp_make_rhs(), rhs_params)
      if (isTRUE(st)) "stable" else if (isFALSE(st)) "saddle/unstable" else "marginal"
    }
  list(max_re = max_re, type = type)
}

## ----------------------------------------------------------------------------
## PUBLIC: find every pathogenic Tregs-off fixed point for one parameter set.
## r_L_P must be > 0 (the "pos" PL branch divides by it; use the sterile finder
## for r_L_P = 0).  Returns a data.frame, one row per fixed point.
##
## Seed resolution: n_M1 x n_PLP grid of (M1, PLP) starts, crossed with L_seeds
## and both PL branches.  Defaults (4 x 4 x 4) were benchmarked to recover the
## same fixed points as a much denser grid; raise them if a count looks off.
## ----------------------------------------------------------------------------
find_fixed_points_pathogenic = function(par, Kg, beta_ros, r_L_P, M_tot, T_tot = 0,
                                        theta_treg = 0, n_M1 = 4, n_PLP = 4,
                                        L_seeds = c(1, 5, 30, 100, 400),
                                        extra_seeds = NULL, grid_seed = TRUE) {  # L=1: low-L attractors
  ## extra_seeds: optional matrix/data.frame with columns M1, PLP, L -- extra Newton starts,
  ## e.g. the observed simulation ENDPOINTS, which sit ON the real attractors and recover
  ## attractors the coarse grid misses.
  ## grid_seed: TRUE = seed Newton from the (M1,PLP,L) grid (finds saddles + all field FPs;
  ## the default, used by the analysis scripts). FALSE = skip the grid and seed ONLY from
  ## extra_seeds -- for the sweep, where we care only about the fan-REACHED attractors (each
  ## endpoint converges back to its own attractor), which is far cheaper (~10 vs ~250 solves).

  stopifnot(r_L_P > 0)
  red = pathogenic_reduction(par, Kg, beta_ros, r_L_P, M_tot, T_tot, theta_treg)

  M1_seeds  = seq(1, M_tot * 0.95, length.out = n_M1)
  ## small PLP seeds: a pathogen-carrying attractor can sit at PLP << 1 (e.g. ~0.26).
  PLP_seeds = c(0, 0.2, 0.5, seq(1, FP_FIXED$K_LP_P * 0.95, length.out = n_PLP - 1))

  push = function(x, pl_branch) {
    if (is.null(x)) return(invisible())
    s = red$quant(x[1], x[2], x[3], pl_branch)
    if (is.null(s) || s$En < -1e-6 || s$PL < -1e-6) return(invisible())
    if (x[1] < -1e-6 || x[1] > M_tot + 1e-6 || x[2] < -1e-6) return(invisible())
    solutions[[length(solutions) + 1]] <<-
      list(M1 = x[1], PLP = x[2], L = x[3], Eh = s$Eh, PL = s$PL,
           Ea = s$Ea, En = s$En, R = s$R, D = s$D, CLP = s$CLP,
           Tr = s$Tr, S = s$S, branch = pl_branch)
  }

  solutions = list()
  if (grid_seed)
    for (pl_branch in c("zero", "pos"))
      for (m in M1_seeds) for (p in PLP_seeds) for (l in L_seeds)
        push(fp_newton(function(x) red$resid(x, pl_branch), c(m, p, l)), pl_branch)

  ## Re-seed the POS branch from every fixed point found so far, nudging PLP off
  ## ~0, to catch the pathogen-carrying attractor that bifurcates off a
  ## pathogen-free saddle (unstable when lambda_PL = r_L_P - alpha_AMP*D > 0)
  ## and can sit in a low-(M1,PLP,L) corner the grid misses.
  for (sf in solutions)
    for (nudge in c(0.2, 0.5, 2))
      push(fp_newton(function(x) red$resid(x, "pos"),
                     c(sf$M1, max(sf$PLP, nudge), sf$L)), "pos")

  ## endpoint-seeded starts: Newton from observed simulation endpoints (M1, PLP, L=Ea+En),
  ## which lie ON the real attractors -> recovers stable fixed points the grid seeds miss.
  if (!is.null(extra_seeds) && NROW(extra_seeds) > 0) {
    es = as.matrix(extra_seeds)
    for (r in seq_len(nrow(es))) {
      s0 = c(es[r, "M1"], es[r, "PLP"], es[r, "L"])
      if (any(!is.finite(s0))) next
      ## NEGATIVE PLP. At a pathogen-free attractor lsoda undershoots PLP below
      ## zero -- observed at -8e-45, -1.7e-16, -4.6e-08 in the same fan. quant()
      ## rejects PLP < 0 and resid falls back to its 1e6 penalty, whose Jacobian
      ## is constant and therefore singular, so Newton returns NULL and the
      ## attractor is silently lost. The endpoint is physically at PLP = 0; the
      ## sign is integrator noise. Try the clamped value AND a small positive
      ## offset (the reduction evaluates cleanly at PLP ~ 1e-40, so 1e-12 is
      ## safely inside the well-behaved region and far below K_LP_P = 100).
      plp_try = unique(c(max(s0[2], 0), 1e-12))
      for (pv in plp_try) {
        sv = c(s0[1], pv, s0[3])
        for (pl_branch in c("zero", "pos"))
          push(fp_newton(function(x) red$resid(x, pl_branch), sv), pl_branch)
      }
    }
  }

  if (length(solutions) == 0) return(data.frame())

  df = do.call(rbind, lapply(solutions, function(v)
    data.frame(v, stringsAsFactors = FALSE)))

  ## de-duplicate fixed points found from different seeds
  dup_key = paste(round(df$Eh, 3), round(df$PLP, 3), round(df$M1, 3))
  df = df[!duplicated(dup_key), , drop = FALSE]
  df = df[order(df$Eh), , drop = FALSE]
  rownames(df) = NULL

  ## classify each
  rhs_params = fp_rhs_params(par, Kg, beta_ros, r_L_P = r_L_P, theta_treg = theta_treg)
  cl = lapply(seq_len(nrow(df)), function(i) classify_pathogenic(df[i, ], red, rhs_params))
  df$Ei = FP_FIXED$K_E - df$Eh
  df$theta_treg = theta_treg
  df$max_Re_eig = vapply(cl, function(z) z$max_re, numeric(1))
  df$type       = vapply(cl, function(z) z$type, character(1))

  df[, c("M1","PLP","PL","L","Eh","Ei","Ea","En","R","D","CLP","Tr","S",
         "max_Re_eig","type","theta_treg","branch")]
}

## ============================================================================
## PUBLIC DISPATCHER: one entry point, routes by r_L_P (0 => sterile, > 0 =>
## pathogenic) and forwards ... (e.g. extra_seeds) to the pathogenic finder.
## Output schema is identical across regimes, so results row-bind directly.
## ============================================================================
find_fixed_points = function(par, Kg, beta_ros, r_L_P = 0,
                             M_tot = FP_FIXED$K_M, T_tot = FP_FIXED$K_T,
                             theta_treg = 0, ...) {
  if (r_L_P == 0) {
    ## ... carries extra_seeds to sterile too. Both branches are now seeded
    ## from the ray endpoints; they differ only in the algebra of the solve.
    find_fixed_points_sterile(par, Kg, beta_ros,
                              M_tot = M_tot, T_tot = T_tot,
                              theta_treg = theta_treg, ...)
  } else {
    find_fixed_points_pathogenic(par, Kg, beta_ros, r_L_P = r_L_P,
                                 M_tot = M_tot, T_tot = T_tot,
                                 theta_treg = theta_treg, ...)
  }
}
