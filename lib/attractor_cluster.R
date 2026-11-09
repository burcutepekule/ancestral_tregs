## ============================================================================
## lib/attractor_cluster.R -- when two fixed points are one attractor.
## ----------------------------------------------------------------------------
## ONE DEFINITION, SOURCED BY EVERYTHING THAT CLUSTERS. It lived in
## pipeline/01_sweep.R while the sweep was the only consumer, and pipeline/21_corridor_tzero.R
## carried a different rule of its own: a gap in the SORTED COMPOSITE, breaking
## wherever diff(composite) exceeded the tolerance. That rule is wrong and the
## sweep abandoned it, but the corridor kept it, so the two stages could
## disagree about how many attractors a combination has -- and the corridor is
## the stage that decides which basin a state is in.
##
## WHY THE COMPOSITE GAP FAILS. The composite sqrt((K_E - Eh)^2 + P_LP^2) is one
## number standing for two coordinates, and two states far apart in P_LP can sit
## at nearly the same composite. Measured on the sweep: Eh 38.7 / P_LP 0 and
## Eh 43.7 / P_LP 26.9 have composites 61.30 and 62.40, a gap of 1.1, so any
## tolerance of 2 merges a cleared attractor with one carrying 27 units of
## pathogen and calls the combination monostable.
##
## THE RULE. Two stable fixed points are the same attractor when they are within
## RAY_CLUST_TOL in BOTH coordinates -- a Chebyshev ball, not a distance in the
## composite. The tolerance is 1, not the composite rule's 2.
##
## WHY SINGLE-LINKAGE UNION-FIND. A Chebyshev ball is not transitive: A can be
## within tol of B and B of C while A and C are not. In one dimension the old
## code could sort and break on gaps; in two there is no sort order that does
## the same job, so the grouping is explicit. It is deterministic and
## independent of input order.
##
## RAY_CLUST_TOL is expected in the calling environment, so a script that sets
## it from its own environment variable keeps control of it.
## RUN: source("./lib/attractor_cluster.R")
## ============================================================================
fp_same = function(Eh1, PLP1, Eh2, PLP2, tol = RAY_CLUST_TOL)
  abs(Eh1 - Eh2) <= tol & abs(PLP1 - PLP2) <= tol

## Cluster ids are renumbered by increasing composite, so id 1 is always the
## lowest-burden attractor and callers can select it by identity rather than by
## position in an arbitrary ordering.
fp_cluster = function(Eh, PLP, cmp, tol = RAY_CLUST_TOL) {
  n = length(Eh); if (!n) return(integer(0))
  parent = seq_len(n)
  root = function(i) { while (parent[i] != i) i = parent[i]; i }
  if (n > 1L) for (i in seq_len(n - 1L)) for (j in (i + 1L):n)
    if (fp_same(Eh[i], PLP[i], Eh[j], PLP[j], tol)) {
      ri = root(i); rj = root(j); if (ri != rj) parent[rj] = ri
    }
  lab = vapply(seq_len(n), root, integer(1))
  u   = unique(lab)
  u   = u[order(vapply(u, function(k) min(cmp[lab == k]), numeric(1)))]
  match(lab, u)
}
