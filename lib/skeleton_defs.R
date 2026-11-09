## ============================================================================
## lib/skeleton_defs.R -- the shared model definitions the skeleton figures need.
## ----------------------------------------------------------------------------
## WHY THIS FILE EXISTS
## 06b_skeleton_fig_*.R used to obtain rhs_app(), seed(), cmpf(), tg and the
## STRUCTURAL/EFFERO/KIN constants by reading 10_corridor_circularity.R as TEXT,
## truncating it at the "computing frac_rescued" line and source-ing the prefix.
## That prefix also loads ./corridor_bulge_<branch>.rds and a "tau_grid_batches"
## file -- a retired stage output and a filename that never existed -- so the
## figure died on a missing file it never used.
##
## Everything below is copied verbatim from 05b so the two cannot drift, minus the
## two data reads. tg now comes from TAU_FILE with the real filename.
##
## Requires BRANCH to be set before sourcing (PL0/PLP0 depend on it).
## ============================================================================
suppressMessages({library(data.table); library(deSolve)})
source("./lib/ODE_SYSTEM_v6.R"); source("./lib/FP_finder_v6.R")
if (!exists("BRANCH")) BRANCH = Sys.getenv("BRANCH", "sterile")

STRUCTURAL=list(K_E=100,K_M=100,K_T=100,K_R=1,K_L_P=1000,K_LP_P=100,K_LP_C=100,CL=1000)
EFFERO=list(a0=0.0045,a_floor=0.0002,n_g=4)
KIN=c("kappa_E","tau_r","tau_P_B","mu_Ea","mu_R","tau_01","tau_10","tau_1_d","tau_eff",
      "tau_t_a","tau_t_d","alpha_0","alpha_1_e","alpha_1_r","alpha_AMP","gamma_B_P",
      "gamma_B_C","sigma_D","sigma_S","tau_D","tau_S","w_E_over_M")
STATE_ORDER=c("PL","PLP","CLP","Eh","Ea","En","M0","M1","R","T0","Tr","D","S")
M_TOT=STRUCTURAL$K_M; T_TOT=STRUCTURAL$K_T; K_E=STRUCTURAL$K_E
PL0=if(BRANCH=="sterile")0 else 500; PLP0=if(BRANCH=="sterile")0 else 10
## ---- finder(), AND WHY IT TAKES SEEDS ---------------------------------------
## This wrapper used to call the branch finders with the grid seeds only, while
## pipeline/01_sweep.R called the same functions with the integrated final states as
## extra_seeds. Same function, two call sites, different seeding -- and the
## analysis scripts got the weak one. Measured on batch 2 row_ind 1756 Kg 40
## beta_ros 0.005 r_L_P 0.08, Tregs OFF: the grid-only solve returns one stable
## fixed point at composite 80.415, while the Eh0 line settles for Eh0 >= 15 on a
## state at composite 58.550 that is not in its list. 80.415 - 58.550 = 21.865,
## which is hi_lo_gap_off in rescue_pathogenic_help.rds to five decimals, so the
## sweep is right and the missing attractor is the finder's. It sits at
## M1 = 9.7, PLP = 13.5, L ~ 5, a corner no grid seed is near. Every one of the 6
## combinations 14_eh0_star_fine.R called monostable in pathogenic_help is this
## same failure.
##
## extra_seeds is a matrix with columns M1, PLP, L (= Ea + En), one row per
## integrated final state, exactly what ray_one() passes. Seeds start Newton ON
## the attractor, so they can only ADD fixed points, never move or remove one.
## The grid is kept (grid_seed = TRUE): it is what finds the saddles, and
## n_stable is read off the same set.
finder=function(par,KG,BR,RLP,th,extra_seeds=NULL) if(BRANCH=="sterile")
  find_fixed_points_sterile(par,KG,BR,M_tot=M_TOT,T_tot=T_TOT,theta_treg=th,
                            extra_seeds=extra_seeds) else
  find_fixed_points_pathogenic(par,KG,BR,r_L_P=RLP,M_tot=M_TOT,T_tot=T_TOT,theta_treg=th,
                               extra_seeds=extra_seeds,grid_seed=TRUE)


## NON-NEGATIVITY CLAMP -- must match pipeline/01_sweep.R's rhs_app, or a phase
## portrait drawn here would show trajectories the sweep never integrated.
## See lib/ODE_SYSTEM_v6.R for why (PLP = 0 is an invariant manifold).
## eff_ko MULTIPLIES g, AND ONLY g. Absent, it is 1 and this is the sweep's RHS
## unchanged -- every existing caller is unaffected. Passing 0 removes clearance.
##
## Setting a0 = a_floor = 0 would look like the same thing and is not: K_En is
## defined as (K_E*tau_P_B + mu_Ea*K_Ea)/(a0*K_M), so a0 = 0 sends K_En to Inf,
## f_En to 0, and silently switches off DAMP production from necrotic debris too.
## K_Ea and K_En are normalising constants defined off the INTACT a0 and must stay
## that way under knockout.
rhs_app=function(t,y,p){ y=pmax(y,0); if(is.null(p$eff_ko)) p$eff_ko=1
  with(as.list(c(setNames(as.list(y),STATE_ORDER),p)),{
  f_P=PLP/K_LP_P;f_C=CLP/K_LP_C;f_PC=(PLP+CLP)/(K_LP_P+K_LP_C);C_rat=ifelse((f_P+f_C)>0,f_C/(f_P+f_C),0)
  K_Ea=(K_E*tau_r*K_R)/(mu_Ea+a0*K_M);K_En=(K_E*tau_P_B+mu_Ea*K_Ea)/(a0*K_M)
  L=max(Ea+En,0);g=eff_ko*(a_floor+(a0-a_floor)*Kg^n_g/(Kg^n_g+L^n_g));effA=g*(M0+M1)*Ea;effN=g*(M0+M1)*En;fb=tau_eff*g*Ea*M1
  dEh=kappa_E*(K_E-Eh)-Eh*(tau_r*R+tau_P_B*(PLP/K_LP_P));dEa=Eh*tau_r*R-mu_Ea*Ea-effA;dEn=Eh*tau_P_B*(PLP/K_LP_P)+mu_Ea*Ea-effN
  dM0=-tau_01*D*M0+tau_10*S*M1+tau_1_d*M1+fb;dM1=tau_01*D*M0-tau_10*S*M1-tau_1_d*M1-fb;f_En=En/(K_En+En)
  dD=sigma_D*(w_E_over_M*f_En+(1-w_E_over_M)*f_PC)-tau_D*D;dS=sigma_S*Tr-tau_S*S
  dT0=-theta_treg*tau_t_a*(M1*C_rat)*T0+tau_t_d*Tr;dTr=theta_treg*tau_t_a*(M1*C_rat)*T0-tau_t_d*Tr;dR=beta_ros*M1-mu_R*R
  dPL=PL*r_L_P*(1-PL/K_L_P)-alpha_AMP*D*PL
  dPLP=PL*gamma_B_P*(1-Eh/K_E)*(1-PLP/K_LP_P)+PLP*r_LP_P*(1-PLP/K_LP_P)-PLP*(alpha_0*M0+alpha_1_e*M1+alpha_1_r*R)
  dCLP=CL*gamma_B_C*(1-Eh/K_E)*(1-CLP/K_LP_C)-CLP*(alpha_0*M0+alpha_1_e*M1+alpha_1_r*R)
  list(c(dPL,dPLP,dCLP,dEh,dEa,dEn,dM0,dM1,dR,dT0,dTr,dD,dS))}) }
seed=function(e0){y=setNames(numeric(13),STATE_ORDER);y["Eh"]=e0;y["En"]=K_E-e0
  y["M0"]=M_TOT;y["T0"]=T_TOT; if(PL0>0)y["PL"]=PL0; if(PLP0>0)y["PLP"]=PLP0; y}
cmpf=function(Eh,PLP) sqrt((K_E - Eh)^2 + PLP^2)   # COMPOSITE ruler: resolved=LOW burden, chronic=HIGH (matches 03_benefit). sterile PLP==0 -> reduces to the old Eh-only ruler.

## ---- THE EXTINCTION FLOOR AND THE SOLVER SETTINGS ---------------------------
## Defined once in lib/extinction_floor.R and sourced by both this file and
## pipeline/01_sweep.R, so the figures cannot integrate under different dynamics or
## different tolerances than the sweep that produced their labels. It provides
## EXT_TOL, rootfun_ext, EVENTS_EXT, ATOL, RTOL, HMAX and MAXSTEPS, and needs
## STATE_ORDER, which is defined above.
source("./lib/extinction_floor.R")

## ---- TAU GRIDS, ALL OF THEM, KEYED BY (batch, row_ind) ----------------------
## This used to load ONE grid, ./tau_grid_v6_1.rds by default, and the figure
## scripts looked up tg[row_ind == cc$row_ind]. Each grid numbers its own rows
## from 1 (tau_grid[, row_ind := .I]), so once a second batch was ingested that
## lookup silently returned batch 1's parameters for a batch-2 combination: the
## panel drew one system's fixed points with another system's dynamics, and the
## trajectories no longer landed on the plotted attractors.
##
## Every ./tau_grid_v6_<n>.rds present is loaded and tagged with its own n, and
## the table is keyed on (batch, row_ind) so a lookup that omits the batch is a
## visible error rather than a wrong row. TAU_FILES overrides the set.
TAU_FILES = Sys.getenv("TAU_FILES", "")
.tf = if (nzchar(TAU_FILES)) strsplit(TAU_FILES, ",")[[1]] else sort(Sys.glob("./data/tau_grid_v6_*.rds"))
.tf = .tf[file.exists(.tf)]
if (!length(.tf)) stop("no ./tau_grid_v6_<n>.rds found -- set TAU_FILES to the grids you swept")
## THE BATCH NUMBER NOW COMES FROM THE DATA, NOT THE FILENAME. pipeline/00_parameter_grids.R
## stamps `batch`, `tau_id` and `tau_stride` into every grid it writes. Reading
## them here instead of parsing the path means a renamed or copied grid still
## carries the right key, and a grid written before that change is caught rather
## than silently mis-keyed. `tau_id` is the only globally unique combination key
## and nothing downstream should ever reconstruct it by arithmetic.
.tb = as.integer(sub(".*tau_grid_v6_([0-9]+)\\.rds$", "\\1", .tf))
tg = rbindlist(Map(function(p, b) {
  d = as.data.table(readRDS(p))
  if (!"batch" %in% names(d))
    stop(sprintf("%s predates the stamped key -- rerun pipeline/00_parameter_grids.R", p))
  if (any(d$batch != b))
    stop(sprintf("%s carries batch %s but is named batch %d",
                 p, paste(unique(d$batch), collapse = ","), b))
  d }, .tf, .tb), fill = TRUE)
stopifnot(!anyNA(tg$tau_id), !anyDuplicated(tg$tau_id))
setkey(tg, batch, row_ind)
cat(sprintf("[skeleton_defs] BRANCH=%s | PL0=%g PLP0=%g | %d tau grid(s) %s | %d rows keyed on (batch, row_ind)\n",
            BRANCH, PL0, PLP0, length(.tf), paste(.tb, collapse = ","), nrow(tg)))
