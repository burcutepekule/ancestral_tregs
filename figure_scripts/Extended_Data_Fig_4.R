## ============================================================================
## Extended_Data_Fig_4.R -- Extended Data Fig. 4: without Tregs, what a stronger ROS
## response buys and costs, across the sweep's grid of beta_r and pathogen growth rate.
## ----------------------------------------------------------------------------
## One tile per (beta_r, r_P) cell of the sweep, split along its diagonal into two
## triangles, both on the same colour scale:
##   upper triangle   the chance that the pathogen is held below P_CUT at rest (P_CUT = 1
##                    by default, 1% of the lamina propria carrying capacity; set
##                    P_CUT = 1e-3, the extinction tolerance, for clearance instead)
##   lower triangle   the chance that the epithelium ends above half (E_h > EH_CUT)
## Probabilities are reported as percentages by default; PROB01=1 reports them on 0-1
## (tile labels, legend title and scale all follow).
## These are the two marginal chances, not a joint probability: a steady state can be both.
## BOTH ARE PROPERTIES OF A STEADY STATE, so they are read off the stored fixed points and
## nothing is integrated. The Tregs-OFF system is bistable in every cell drawn; each of its
## stable steady states counts once, and the chance is the share of them with the property,
## averaged over the parameter sets of the cell. It is therefore a statement about the
## steady states a cell can settle on, not about how often each is reached from an injury.
## Data: pipeline/45_beta_effero_tradeoff.R. Tregs are OFF throughout.
##
## in : ./data/fig_beta_effero.rds
## out: ./figures/extended_data/Extended_Data_Fig_4.png
## RUN: Rscript pipeline/45_beta_effero_tradeoff.R   (once, writes the table)
##      Rscript figure_scripts/Extended_Data_Fig_4.R
## ============================================================================
suppressMessages({library(data.table); library(ggplot2)})
SRC = Sys.getenv("SRC", "."); setwd(SRC)
BASE_SZ = as.numeric(Sys.getenv("BASE_SIZE", "30"))
EH_CUT  = as.numeric(Sys.getenv("EH_CUT", "50"))
P_CUT   = as.numeric(Sys.getenv("P_CUT", "1"))      # the pathogen criterion; 1e-3 = cleared
PROB01  = as.integer(Sys.getenv("PROB01", "1"))     # 1 = report probabilities on 0-1, not %
SCALE   = if (PROB01) 1 else 100
FMT     = if (PROB01) "%.2f" else "%.0f"
TXT_CUT = if (PROB01) 0.45 else 45                  # light text below, dark above
LEG     = if (PROB01) "Probability" else "Probability (%)"
F_POP   = Sys.getenv("F_POP", "./data/fig_beta_effero.rds")
OUT_FIG = Sys.getenv("OUT_FIG", "./figures/extended_data/Extended_Data_Fig_4.png")
if (!file.exists(F_POP)) stop(F_POP, " not found -- run: Rscript pipeline/45_beta_effero_tradeoff.R")

B = as.data.table(readRDS(F_POP))
## each stable steady state of a parameter set counts once
C = B[, .(sets = .N,
          clear = SCALE * mean(c(P_lo < P_CUT, P_hi < P_CUT)),
          eh    = SCALE * mean(c(Eh_lo > EH_CUT, Eh_hi > EH_CUT))), keyby = .(r_P = r_L_P, beta)]
BX = sort(unique(C$beta)); RY = sort(unique(C$r_P))
C[, `:=`(x = match(beta, BX), y = match(r_P, RY))]
P_LAB = if (P_CUT <= 1e-2) "pathogen cleared" else sprintf("pathogen below %g", P_CUT)
## the panel title names both quantities as probabilities of a steady-state (SS) condition
TITLE = bquote("Upper triangle: " * P(P[LP]^"SS" < .(P_CUT)) *
               "          Lower triangle: " * P(E[h]^"SS" > .(EH_CUT)))
cat(sprintf("[ED Fig. 4] %d cells, %d parameter sets; chance %s %.1f-%.1f%%, chance E_h > %g: %.1f-%.1f%%\n",
            nrow(C), uniqueN(B[, paste(tau_id, Kg, r_L_P)]), P_LAB, min(C$clear), max(C$clear), EH_CUT, min(C$eh), max(C$eh)))
print(C[, .(r_P, beta, sets, cleared = round(clear, 1 + 2*PROB01), epithelium = round(eh, 1 + 2*PROB01))])

## each tile split along its main diagonal: the upper-left triangle carries the
## clearance chance, the lower-right one the epithelium chance
tri = function(d, which) {
  o = if (which == "upper") list(dx = c(-.5, -.5,  .5), dy = c(-.5, .5, .5))
      else                  list(dx = c(-.5,  .5,  .5), dy = c(-.5, -.5, .5))
  d[, .(id = paste(x, y, which), x = x + o$dx, y = y + o$dy,
        v = if (which == "upper") clear else eh), by = .(x0 = x, y0 = y)]
}
P = rbind(tri(C, "upper"), tri(C, "lower"))
LB = rbind(C[, .(x = x - 1/6, y = y + 1/6, v = clear)], C[, .(x = x + 1/6, y = y - 1/6, v = eh)])

g = ggplot() +
  geom_polygon(data = P, aes(x, y, group = id, fill = v), colour = "white", linewidth = 0.5) +
  geom_segment(data = C, aes(x = x - .5, xend = x + .5, y = y - .5, yend = y + .5), colour = "white", linewidth = 0.5) +
  geom_text(data = LB, aes(x, y, label = sprintf(FMT, v), colour = v > TXT_CUT), size = BASE_SZ * 0.21, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "grey15", `FALSE` = "grey95")) +
  scale_fill_viridis_c(option = "viridis", name = LEG, limits = c(0, NA)) +
  scale_x_continuous(breaks = seq_along(BX), labels = function(z) format(BX[z]), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(RY), labels = function(z) format(RY[z]), expand = c(0, 0)) +
  coord_equal() +
  labs(x = expression(beta[r] * ", ROS production rate"), y = expression(r[P] * ", pathogen growth rate"),
       title = TITLE) +
  theme_minimal(base_size = BASE_SZ) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = BASE_SZ * 0.8),
        axis.text.y = element_text(size = BASE_SZ * 0.8),
        plot.title = element_text(size = BASE_SZ * 0.85, hjust = 0, margin = margin(b = BASE_SZ * 0.5)),
        legend.key.height = unit(BASE_SZ * 0.15, "cm"), legend.title = element_text(size = BASE_SZ * 0.8))
dir.create(dirname(OUT_FIG), showWarnings = FALSE, recursive = TRUE)
png(OUT_FIG, width = 3800, height = 2800, res = 200); print(g); invisible(dev.off())
cat("[ED Fig. 4] wrote", OUT_FIG, "\n")
