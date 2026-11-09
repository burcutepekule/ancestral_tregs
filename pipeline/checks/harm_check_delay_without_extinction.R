## ============================================================================
## pipeline/checks/harm_check_delay_without_extinction.R (diagnostic for Supplementary Fig. 3)
## ----------------------------------------------------------------------------
## THE QUESTION. The delay panel of Supplementary Fig. 3 (panel B) draws only combinations whose Tregs OFF arm
## drives the pathogen to the extinction floor, and every point on it is a
## harm-to-help switch. That raises an obvious objection: is extinction in the
## OFF arm NECESSARY for a delay to improve the outcome, or is it merely the
## case that the combinations we happened to look at were the ones that cleared?
##
## EXTINCTION HERE MEANS BOTH COMPARTMENTS AT THE FLOOR, not P_LP alone. The
## reason is in the block that computes it below; the short version is that the
## absorbing event roots on max(P_L, P_LP) and P_L feeds P_LP, so P_LP under the
## floor with P_L above it is not an eliminated pathogen.
##
## A combination whose OFF arm dips to a small but nonzero load might still be
## improved by delaying Treg induction -- to neutral if not to help. If such
## cases exist, the Supplementary Fig. 3 delay panel is drawing a subset and the claim in the text has to be
## narrowed. If none exist, that is a real finding: the delay treatment is
## all-or-nothing, and the reason is mechanical -- suppression after clearance is
## pure ROS saving with no pathogen cost, whereas suppression after a nonzero dip
## always leaves a population that regrows into the suppressed state.
##
## WHAT IT DOES. Reads ./data/harm_rescue_data.rds, which pipeline/50_harm_rescue_compute.R
## writes. Integrates nothing. Takes the combinations that are HARMED at delay 0
## -- the untreated Tregs ON arm -- and for each asks two independent questions:
## did its Tregs OFF arm reach the floor, and what is the best outcome any delay
## achieved? The cross-tabulation of those two answers is the whole result.
##
## It then prints, individually, every combination that improved WITHOUT its OFF
## arm having cleared. Those are the counterexamples, and if there are any they
## are worth looking at one by one rather than counting: the question is whether
## the mechanism is possible, not how often it happens.
##
## in : ./data/harm_rescue_data.rds; ./data/selected_big.rds if present
## out: prints; writes ./data/harm_check_delay_without_extinction.csv when
##      counterexamples exist and ./data/harm_check_help_harm_off_state.csv when
##      selected_big.rds is present
## RUN: Rscript pipeline/checks/harm_check_delay_without_extinction.R
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
suppressMessages(library(data.table))
RDS = Sys.getenv("IN_RDS", "./data/harm_rescue_data.rds")
TOL = as.numeric(Sys.getenv("EFFECT_TOL", "1.0"))
if (!file.exists(RDS)) stop(RDS, " not found -- run pipeline/50_harm_rescue_compute.R first")

Z  = readRDS(RDS)
SM = as.data.table(Z$sm)
## Only the scalar entries are printed. The stamp also records the treatment
## ladders, which are vectors, and unlisting the whole thing mixes lengths.
if (!is.null(Z$settings)) {
  .st = Z$settings; .sc = .st[vapply(.st, length, 1L) == 1L]
  cat("[chk] settings stamp:",
      paste(sprintf("%s=%s", names(.sc), unlist(.sc)), collapse = "  "), "\n")
}

## The key. cid exists only in files written after it was added; fall back to the
## five-part key so this runs against a data file written before cid existed.
KEY = if ("cid" %in% names(SM)) "cid" else
        c("batch", "row_ind", "Kg", "beta_ros", "r_L_P")
cat(sprintf("[chk] keyed on %s | %d combinations | %d rows\n",
            paste(KEY, collapse = "+"), uniqueN(SM[, ..KEY]), nrow(SM)))

## effect is derived here rather than read: it is assigned in the figure script
## AFTER the file is written, so it is not in the stored object.
SM[, effect := fifelse(dcomp < -TOL, "help", fifelse(dcomp > TOL, "harm", "neutral"))]

## ---- DID THE TREGS OFF ARM DRIVE THE PATHOGEN EXTINCT? ---------------------
## THE DEFINITION MATTERS AND THE OBVIOUS ONE IS TOO WEAK. t_zero is the first
## stored time the FLOORED P_LP series reads zero, so it is defined as soon as
## P_LP alone passes under EXT_TOL. Extinction is a different condition: the
## absorbing event has the root max(P_L, P_LP) - EXT_TOL and zeroes BOTH
## compartments, because P_L feeds P_LP through
##     dPLP = PL*gamma_B_P*(1 - Eh/K_E)*(1 - PLP/K_LP_P) + ...
## a term with no factor of P_LP. While P_L is above the floor that influx is
## still running, so a P_LP that has gone under the floor is not an eliminated
## pathogen and can come back.
##
## The claim being tested is about extinction, so extinction is what is tested:
## both compartments at the floor at the end of the run. Nothing has to be
## re-integrated for it -- run() records PL_end and PLP_end on every arm, and
## once the event has fired both are exactly zero.
##
## t_zero IS KEPT AS A SECOND OPINION rather than dropped. Where the two
## disagree, P_LP went under the floor and the pathogen was not eliminated, and
## that combination would be counted as cleared under the weaker definition and
## is not here. The count is printed, because a large disagreement would mean
## the figures that select on t_zero are selecting on the wrong thing.
EXTTOL = as.numeric(Sys.getenv("EXT_TOL",
           if (!is.null(Z$settings$EXT_TOL)) as.character(Z$settings$EXT_TOL) else "1e-11"))
OFF = SM[treat == "off"]
if (!all(c("PL_end", "PLP_end") %in% names(OFF)))
  stop("the summary rows carry no PL_end/PLP_end, so extinction cannot be tested; re-run pipeline/50_harm_rescue_compute.R")
OFF[, extinct := PL_end < EXTTOL & PLP_end < EXTTOL]
cat(sprintf("[chk] extinction on the Tregs OFF arm (both P_L and P_LP under %g at t = TMAX): %d of %d\n",
            EXTTOL, sum(OFF$extinct), nrow(OFF)))
if ("t_zero" %in% names(OFF)) {
  .w = OFF[, sum((!is.na(t_zero)) != extinct)]
  cat(sprintf("[chk] t_zero (P_LP alone under the floor) and extinction disagree on %d combinations%s\n",
              .w, if (.w) " -- listed below" else ""))
  if (.w) print(OFF[(!is.na(t_zero)) != extinct,
                    c(KEY, "t_zero", "PL_end", "PLP_end"), with = FALSE][seq_len(min(10L, .N))])
}
OFF[, cleared := extinct]
OFF = OFF[, c(KEY, "cleared", "PLP_dip", "t_dip"), with = FALSE]
setnames(OFF, c("PLP_dip", "t_dip"), c("dip_off", "t_dip_off"))

## ---- HARMED UNTREATED, AND THE BEST ANY DELAY ACHIEVED ---------------------
DEL   = SM[treat == "delay"]
HARM0 = DEL[level == min(level) & effect == "harm", ..KEY]
cat(sprintf("[chk] %d of %d combinations are harmed by untreated Tregs (delay 0)\n",
            nrow(HARM0), uniqueN(DEL[, ..KEY])))
if (!nrow(HARM0)) stop("nothing is harmed at delay 0 -- nothing to test")

D = merge(DEL[HARM0, on = KEY], OFF, by = KEY)
RANK = c(harm = 1L, neutral = 2L, help = 3L)
BEST = D[, .(best = names(RANK)[max(RANK[effect])],
             best_level = level[which.max(RANK[effect])],
             n_levels = .N),
         by = c(KEY, "cleared", "dip_off", "t_dip_off")]
BEST[, best := factor(best, levels = c("harm", "neutral", "help"))]

cat("\n[chk] BEST OUTCOME REACHED BY ANY DELAY, against whether the Tregs OFF arm cleared:\n")
print(dcast(BEST[, .N, by = .(cleared, best)], cleared ~ best,
            value.var = "N", fill = 0L, drop = FALSE))

## ---- THE COUNTEREXAMPLES ---------------------------------------------------
## Improvement without clearance. These are the cases that would narrow the claim
## the Supplementary Fig. 3 delay panel makes, so they are listed rather than counted.
CE = BEST[cleared == FALSE & best != "harm"][order(-dip_off)]
cat(sprintf("\n[chk] combinations that improved WITHOUT their OFF arm clearing: %d\n", nrow(CE)))
if (nrow(CE)) {
  print(CE[, c(KEY, "dip_off", "t_dip_off", "best", "best_level"), with = FALSE])
  fwrite(CE, "./data/harm_check_delay_without_extinction.csv")
  cat("[chk] wrote ./data/harm_check_delay_without_extinction.csv\n")
  cat("\n   READ THIS AS: a delay improved the outcome even though the pathogen was\n")
  cat("   never eliminated in the Tregs OFF arm. dip_off is the load it reached at\n")
  cat("   its lowest. The Supplementary Fig. 3 delay panel excludes these, and if the smallest dip_off\n")
  cat("   values dominate the list, the boundary is a matter of degree rather than\n")
  cat("   a clean requirement for extinction.\n")
} else {
  cat("\n   NONE. Every combination a delay improved had already driven the\n")
  cat("   pathogen extinct in its Tregs OFF arm, and no combination that did not\n")
  cat("   reached even neutral at any delay. That is the all-or-nothing result,\n")
  cat("   and it is what justifies the Supplementary Fig. 3 delay panel drawing only the extinct set.\n")
}
## ---- CLEARED, BUT NEVER HELPED ---------------------------------------------
## The other off-diagonal cell. These combinations drove the pathogen extinct
## without Tregs and still could not be moved into help by any delay -- and one
## of them stops at neutral, which is a different statement from stopping at
## harm and is worth separating.
##
## WHY NEUTRAL IS THE FLOOR OF WHAT A DELAY CAN DO. help is dcomp < -EFFECT_TOL,
## and dcomp is measured against the Tregs OFF arm, so help means the delayed
## arm ends BETTER than having no Tregs at all. A long enough delay makes the
## treated arm approach the untreated one, which sends dcomp to zero: neutral.
## Getting past neutral needs something the delay cannot supply on its own, and
## that something is the saving from suppressing an inflammatory response after
## the pathogen is already gone -- no pathogen cost, only ROS damage avoided.
##
## SO A CLEARED COMBINATION THAT STOPS AT NEUTRAL IS ONE WHOSE POST-CLEARANCE
## SAVING IS SMALLER THAN EFFECT_TOL, not one where the mechanism is absent.
## best_dcomp says which: a value just above -EFFECT_TOL is a combination the
## tolerance excluded, and one near zero is a combination the delay genuinely
## only restored to the untreated outcome. Both are printed rather than
## described, because the difference decides how the case is reported.
NH = BEST[cleared == TRUE & best != "help"]
if (nrow(NH)) {
  BD = D[, .(best_dcomp = min(dcomp), at_level = level[which.min(dcomp)],
             dcomp_at_0 = dcomp[which.min(level)]), by = KEY]
  NH = merge(NH, BD, by = KEY)[order(best_dcomp)]
  cat(sprintf("\n[chk] cleared but never helped: %d combinations (%d stop at harm, %d at neutral)\n",
              nrow(NH), NH[best == "harm", .N], NH[best == "neutral", .N]))
  cat(sprintf("      best_dcomp is the most improvement any delay reached; help needs < %g\n", -TOL))
  ## best_level is the level at which the best CLASS was reached; at_level is
  ## the level with the lowest burden difference. They differ whenever the class
  ## never changes, and both are shown so the two are not read as one number.
  print(NH[best != "harm" | seq_len(.N) <= 10L,
           c(KEY, "best", "best_level", "at_level", "dcomp_at_0", "best_dcomp"),
           with = FALSE])
  .near = NH[best_dcomp > -TOL & best_dcomp < -0.5 * TOL, .N]
  cat(sprintf("      %d of them sit between %g and %g, i.e. excluded by the tolerance rather\n",
              .near, -TOL, -0.5 * TOL))
  cat("      than by the absence of the effect. Report that count with the cell.\n")
  ## DELAY CAN ALSO MAKE THINGS WORSE, and this is where that shows. at_level 0
  ## means the least harmful delay was no delay at all: every longer one raised
  ## the burden difference. That is not a failure to rescue, it is the treatment
  ## running the wrong way, and it has a reading -- inducing Tregs early spares
  ## the epithelium at the cost of leaving a pathogen, while inducing them late
  ## clears the pathogen but lets ROS run, and the composite can prefer the
  ## first. Counted rather than described, because it bears on whether a delay
  ## is safe to propose where clearance is not certain.
  .worse = NH[at_level == min(D$level), .N]
  cat(sprintf("      %d of them are LEAST harmful at delay 0, so every delay tested made\n", .worse))
  cat("      them worse rather than better.\n")
}

## WHAT "NO DELAY HELPS" CAN AND CANNOT MEAN. The delays tested are a ladder,
## so the statement is that no delay ON THE LADDER helps, and a class reached
## only between two rungs of it would be missed. Two things bound that. The
## ladder is dense where the thresholds are -- the located ones sit between 10
## and 2400 and the spacing there is a factor of two -- and the outcome is
## monotone in the delay on every combination that has a bracket, so a class
## reached between two levels would have to be reached at the higher one as
## well. The ceiling is the real limit: a combination still harmed at the
## longest delay might switch at a longer one that was never integrated. The
## count is printed so the claim is bounded rather than assumed.
.top = D[level == max(level)]
.topharm = merge(.top[effect == "harm", ..KEY], BEST[, c(KEY, "cleared"), with = FALSE], by = KEY)
cat(sprintf("\n[chk] still harmed at the longest delay tested (%g): %d combinations, %d of them\n",
            max(D$level), nrow(.topharm), .topharm[cleared == TRUE, .N]))
cat("      with an extinct OFF arm. Those are the ones a longer ladder could still\n")
cat("      convert; the rest are the all-or-nothing result and no ceiling affects them.\n")

## ---- HOW CLOSE DO THE NON-CLEARING ONES GET? -------------------------------
## If the answer above is NONE, the next question is whether it is a sharp
## boundary or an artefact of the population sampled. A non-clearing combination
## whose dip is 1e-3 is a near miss; one whose dip is 20 is not. This prints the
## distribution of dip_off among the combinations no delay could improve, so the
## gap between the two groups can be seen rather than assumed.
NC = BEST[cleared == FALSE]
if (nrow(NC)) {
  cat(sprintf("\n[chk] dip_off among the %d non-clearing combinations (min / q25 / med / q75 / max):\n   %s\n",
              nrow(NC), paste(sprintf("%.4g", quantile(NC$dip_off, c(0, .25, .5, .75, 1))),
                              collapse = "  ")))
  cat(sprintf("   of these, %d have dip_off below 1 (i.e. under 1e11 organisms) and\n",
              NC[dip_off < 1, .N]))
  cat(sprintf("   %d of THOSE reached neutral or help at some delay.\n",
              NC[dip_off < 1 & best != "harm", .N]))
}

## ---- AND THE OTHER THREE TREATMENTS, FOR CONTRAST --------------------------
## The delay arm is the one the Supplementary Fig. 3 delay panel is about, but the same question applies to
## the rate treatments, and their answer is the comparison that makes the delay
## result interesting. Reported per treatment, not pooled.
cat("\n[chk] the same question for every treatment -- best outcome by clearance:\n")
ALL = merge(SM[treat != "off"][HARM0, on = KEY], OFF, by = KEY)
AB = ALL[, .(best = names(RANK)[max(RANK[effect])]), by = c(KEY, "treat", "cleared")]
print(dcast(AB[, .N, by = .(treat, cleared, best)], treat + cleared ~ best,
            value.var = "N", fill = 0L))

## ============================================================================
## IS THE DELAY-PROOF GROUP ALREADY IN BAD SHAPE WITHOUT TREGS?
## ----------------------------------------------------------------------------
## The claim to be made in the text is that the combinations a delay cannot help
## are the ones whose Tregs OFF arm is poor to begin with -- so the delay is not
## failing to rescue a recoverable system, it is failing on systems that were
## never going to resolve. That has to be MEASURED, not assumed, and one obvious
## way of measuring it is circular: "cleared" means the OFF arm ended at P_LP = 0,
## so a comparison of final P_LP between the two groups can only restate the
## grouping. It is reported anyway, labelled as the tautology it is, so that the
## non-circular comparison beside it can be read against something.
##
## E_h IS THE NON-CIRCULAR ONE. Nothing about the grouping constrains where the
## epithelium ends up: a combination could clear the pathogen and still finish
## with heavy damage, or hold a chronic load at a nearly intact epithelium. If
## the delay-proof group sits at low E_h_off, the claim holds on its own terms.
## comp_off combines the two and is reported third.
cat("\n[chk] THE TREGS OFF ARM ITSELF, split by whether a delay could improve the\n")
cat("      combination. Eh_off is the informative one -- PLP_off restates the\n")
cat("      grouping, because clearance IS final P_LP = 0.\n")
.q = function(x) sprintf("%7.3g %7.3g %7.3g %7.3g %7.3g",
                         quantile(x, 0, na.rm=TRUE), quantile(x, .25, na.rm=TRUE),
                         quantile(x, .5, na.rm=TRUE), quantile(x, .75, na.rm=TRUE),
                         quantile(x, 1, na.rm=TRUE))
OFFQ = SM[treat == "off"][BEST[, c(KEY, "best"), with = FALSE], on = KEY]
OFFQ[, rescuable := best != "harm"]
for (v in c("Eh_off", "PLP_off", "comp_off")) {
  if (!v %in% names(OFFQ)) next
  cat(sprintf("\n   %s   (min / q25 / med / q75 / max)%s\n", v,
              if (v == "PLP_off") "   <- circular, see above" else ""))
  for (r in c(FALSE, TRUE)) {
    z = OFFQ[rescuable == r][[v]]
    if (length(z)) cat(sprintf("     delay %-13s n=%3d  %s\n",
                               if (r) "CAN help" else "cannot help", length(z), .q(z)))
  }
  z0 = OFFQ[rescuable == FALSE][[v]]; z1 = OFFQ[rescuable == TRUE][[v]]
  if (length(z0) > 2 && length(z1) > 2) {
    w = suppressWarnings(wilcox.test(z0, z1))
    cat(sprintf("     Wilcoxon rank sum p = %.3g (two-sided; a difference in location,\n", w$p.value))
    cat("       not a claim about any individual combination)\n")
  }
}

## ============================================================================
## WHERE HARM LIVES ACROSS THE WHOLE SWEEP
## ----------------------------------------------------------------------------
## The comparison above is inside the harmed set. The wider question is whether
## the combinations Tregs HARM sit at a different Tregs OFF pathogen load from
## the ones Tregs HELP -- a statement about the gated population rather than
## about this experiment. It needs selected_big.rds, which 01_select.R writes, so
## it is skipped rather than required when that file is absent.
##
## lo_off is the composite burden of the LOW-burden Tregs OFF attractor and
## lo_off_PLP is its pathogen load; treg_effect is 01_select.R's own label. This
## is the quantity an overlaid distribution would draw.
SEL_F = Sys.getenv("SELECTED_RDS", "./data/selected_big.rds")
if (!file.exists(SEL_F)) {
  cat(sprintf("\n[chk] %s absent -- skipping the sweep-wide help/harm comparison\n", SEL_F))
} else {
  SEL = as.data.table(readRDS(SEL_F))
  SEL = SEL[patho == TRUE & treg_effect %in% c("help", "harm")]
  cat(sprintf("\n[chk] SWEEP-WIDE, pathogenic branch: %d help, %d harm\n",
              SEL[treg_effect == "help", .N], SEL[treg_effect == "harm", .N]))
  for (v in c("lo_off_PLP", "lo_off_Eh", "lo_off")) {
    if (!v %in% names(SEL)) next
    cat(sprintf("\n   %s at the Tregs OFF low-burden attractor (min / q25 / med / q75 / max)\n", v))
    for (e in c("help", "harm")) {
      z = SEL[treg_effect == e][[v]]
      if (length(z)) cat(sprintf("     %-5s n=%5d  %s\n", e, length(z), .q(z)))
    }
  }
  fwrite(SEL[, c("treg_effect",
                 intersect(c("lo_off_PLP", "lo_off_Eh", "lo_off"), names(SEL))),
             with = FALSE], "./data/harm_check_help_harm_off_state.csv")
  cat("\n[chk] wrote ./data/harm_check_help_harm_off_state.csv -- two columns per quantity, ready\n")
  cat("      for an overlaid density or ECDF without re-reading selected_big.rds\n")
}
