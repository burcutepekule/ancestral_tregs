## ============================================================================
## lib/ODE_SYSTEM_v6.R -- the ODE system of the paper: make_ode_system() returns
## the right-hand side of the 13-state model. The equations are written out in
## docs/methods.pdf, section "Model".
##
## Development history, kept for the record. The vector field is unchanged from
## the version this replaces apart from the non-negativity clamp noted below:
## efferocytosis
## is a load-dependent, collapsing Hill function g(L) instead of a flat
## per-capita clearance rate, and epithelial damage is mediated by a ROS
## compartment R (driven by M1, decaying at mu_R) instead of M1 acting
## directly on Eh/Ea via tau_tnfa.
##
## Relative to the ORIGINAL (pre-Hill) Della version:
##  - REMOVED: tau_tnfa-driven epithelial damage (Eh, Ea no longer see M1
##    directly)
##  - REMOVED: alpha_1_p as a flat per-capita efferocytosis/clearance rate
##    on Ea and En
##  - ADDED:   R (ROS) state variable, with dR = beta_ros*M1 - mu_R*R
##  - ADDED:   g(L) = a_floor + (a0-a_floor)*Kg^n_g/(Kg^n_g + L^n_g), where
##             L = Ea + En (total corpse load). g(L) replaces alpha_1_p as
##             the per-capita clearance rate on Ea and En, and also drives
##             the M1<->M0 efferocytosis feedback term (tau_eff * g * Ea * M1)
##  - ADDED:   alpha_1_r (ROS-augmented clearance of PLP/CLP by R), which
##             REPLACES alpha_1_p in the PLP/CLP clearance terms too
##  - K_Ea / K_En are now defined off tau_r/R and a0, instead of off
##    tau_tnfa/alpha_1_p
##
## Parameters required in `params`:
##    kappa_E, tau_r, tau_P_B, mu_Ea, mu_R, tau_01, tau_10, tau_1_d, tau_eff,
##    tau_t_a, tau_t_d, alpha_0, alpha_1_e, alpha_1_r, alpha_AMP, gamma_B_P,
##    gamma_B_C, sigma_D, sigma_S, tau_D, tau_S, w_E_over_M, beta_ros,
##    theta_treg, r_L_P, r_LP_P, a0, a_floor, Kg, n_g, phi_instability,
##    K_E, K_M, K_T, K_R, K_L_P, K_LP_P, K_LP_C, CL
## ============================================================================

make_noise_fun = function(t_max, treg_eff) {
  n   = as.integer(t_max) + 2L
  p   = runif(n)                          # bernoulli: < treg_eff → use C_rat
  u   = runif(n)                          # uniform draw for noise case
  use_crat = as.numeric(p < treg_eff)    # 1 = use C_rat, 0 = use u
  list(
    use_crat_fun = stepfun(x = seq_len(n - 1L), y = use_crat),
    u_fun        = stepfun(x = seq_len(n - 1L), y = u)
  )
}

make_ode_system = function(noise_fun) {
  function(t, state, params) {
    ## ------------------------------------------------------------------------
    ## NON-NEGATIVITY CLAMP. Every state here is a concentration or a cell count;
    ## none is meaningful below zero. lsoda does not know that, and PLP = 0 is an
    ## INVARIANT MANIFOLD -- the pathogen terms all carry a factor of PLP -- so a
    ## trajectory that should sit on it instead straddles it by roundoff. Once
    ## PLP < 0 the logistic term
    ##      PLP * r_LP_P * (1 - PLP / K_LP_P)   ~   -r * x^2 / K_LP_P  at PLP = -x
    ## is negative and QUADRATIC: not decay toward a carrying capacity but
    ## finite-time escape to -Inf. Measured endpoints reached -1.9e11, and the
    ## resulting relative deviations reached 1e154.
    ##
    ## Evaluating the vector field at the clamped state fixes it at the source:
    ## at PLP <= 0 the field is evaluated at PLP = 0, where dPLP = PL*gamma_B_P*
    ## (1 - Eh/K_E) >= 0, so the manifold is repelling from below instead of a
    ## launch pad. This changes NOTHING for any trajectory that stays in the
    ## positive orthant, which is every trajectory that was ever meaningful.
    ## ------------------------------------------------------------------------
    state = pmax(state, 0)
    with(as.list(c(state, params)), {
      # === State variables
      # Eh  = healthy epithelial cells
      # Ea  = apoptotic cells awaiting clearance (convert to secondary necrosis at rate mu_Ea)
      # En  = necrotic debris (DAMPs source, cleared by efferocytosis Hill switch g(L))
      # M0  = resting macrophages
      # M1  = M1 inflammatory macrophages
      # R   = ROS, produced by M1, decays at mu_R, drives epithelial damage
      # T0  = inactive Tregs
      # Tr  = active Tregs
      # PLP = pathogens in the lamina propria
      # PL  = pathogens in the lumen
      # CLP = commensals in the lamina propria
      # CL  = commensals in the lumen (held constant upstream or as param)
      # D   = DAMPs + MAMPs = danger signals
      # S   = suppression signals (SAMPs from Tregs)
      # CumC = cumulative integral of C_rat_perceived (diff post-hoc)

      f_P  = PLP / K_LP_P
      f_C  = CLP / K_LP_C
      f_PC = (PLP + CLP) / (K_LP_P + K_LP_C)
      C_rat = ifelse((f_P + f_C) > 0, f_C / (f_P + f_C), 0)

      use_crat        = noise_fun$use_crat_fun(floor(t))
      u               = noise_fun$u_fun(floor(t))
      C_rat_perceived = use_crat * C_rat + (1 - use_crat) * u

      ## efferocytosis (Hill, the switch): g(L) decreases as corpse load
      ## L = Ea + En grows, replacing the old flat alpha_1_p clearance rate.
      K_Ea = (K_E * tau_r * K_R) / (mu_Ea + a0 * K_M)
      K_En = (K_E * tau_P_B + mu_Ea * K_Ea) / (a0 * K_M)

      L    = max(Ea + En, 0)
      g    = a_floor + (a0 - a_floor) * Kg^n_g / (Kg^n_g + L^n_g)
      effA = g * (M0 + M1) * Ea
      effN = g * (M0 + M1) * En
      fb   = tau_eff * g * Ea * M1

      # Epithelium -- damage now mediated by ROS (R), not M1 directly
      dEh = kappa_E * (K_E - Eh) - Eh * (tau_r * R + tau_P_B * (PLP / K_LP_P))
      dEa = Eh * tau_r * R              - mu_Ea * Ea - effA
      dEn = Eh * tau_P_B * (PLP/K_LP_P) + mu_Ea * Ea - effN

      # Macrophages — efferocytosis feedback from apoptotic clearance only
      dM0 = - tau_01 * D * M0 + tau_10 * S * M1 + tau_1_d * M1 + fb
      dM1 =   tau_01 * D * M0 - tau_10 * S * M1 - tau_1_d * M1 - fb

      # Danger / suppression signals
      # DAMPs proportional to necrotic debris (bounded Michaelis-Menten form)
      f_En = En / (K_En + En)
      dD   = sigma_D * (w_E_over_M * f_En + (1 - w_E_over_M) * f_PC) - tau_D * D
      dS   = sigma_S * Tr  - tau_S * S

      # Tregs
      tau_t_d_effective = tau_t_d + phi_instability * (D / (sigma_D/tau_D))
      dT0 = - theta_treg * tau_t_a * (M1 * C_rat_perceived) * T0 + tau_t_d_effective * Tr
      dTr =   theta_treg * tau_t_a * (M1 * C_rat_perceived) * T0 - tau_t_d_effective * Tr

      # ROS -- produced by M1, decays at mu_R
      dR = beta_ros * M1 - mu_R * R

      # Luminal pathogens
      dPL = PL * r_L_P * (1 - PL / K_L_P) - alpha_AMP * D * PL

      # LP pathogens -- alpha_1_r (ROS-augmented clearance) replaces alpha_1_p here
      dPLP = PL * gamma_B_P * (1 - Eh/K_E) * (1 - PLP / K_LP_P) +
        PLP * r_LP_P * (1 - PLP / K_LP_P) -
        PLP * (alpha_0 * M0 + alpha_1_e * M1 + alpha_1_r * R)

      # LP commensals -- alpha_1_r (ROS-augmented clearance) replaces alpha_1_p here
      dCLP = CL * gamma_B_C * (1 - Eh/K_E) * (1 - CLP / K_LP_C) -
        CLP * (alpha_0 * M0 + alpha_1_e * M1 + alpha_1_r * R)

      dCumC = C_rat_perceived

      list(c(dPL, dPLP, dCLP, dEh, dEa, dEn, dM0, dM1, dR, dT0, dTr, dD, dS, dCumC))
    })
  }
}
