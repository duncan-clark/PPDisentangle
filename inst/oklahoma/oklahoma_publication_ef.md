## Main-text LaTeX (E/F only)

We focus on the all-free KDE pair: Fit E (naive labels) and Fit F (SEM
relabelled). The bootstrap distribution of the all-or-nothing total
saved,
*Δ̂*<sub>AoN</sub> = *E*\[*N*(control everywhere)\] − *E*\[*N*(actual regime)\],
shows a qualitatively different conclusion across the two fits. For Fit
E, bootstrap replicate means are centered on positive savings (mean
22.2, median -10.3; 95% quantiles \[-203.2,270.1\]), whereas Fit F is
centered near/under zero (mean 29.0, median 8.3; 95% quantiles
\[-62.5,266.7\]). The cumulative diagnostics make the same point: the
SEM-based F fit does not support robust policy savings, while the naive
E fit can suggest savings that are not stable once relabelling
uncertainty is propagated.

## Appendix LaTeX (full E/F exposition)

The post-treatment analysis window is 100 days. The target estimand is
the all-or-nothing total saved over the window under control everywhere
versus the observed treatment regime.

Fits E/F use the bivariate ETAS with KDE background and all
structural/productivity terms free. The model uses process-specific
background rates and a 2 triggering matrix with shared structural terms
(*c*, *p*, *D*, *γ*, *q*) and process-pair productivities
(*A*<sub>00</sub>, *A*<sub>11</sub>, *A*<sub>01</sub>, *A*<sub>10</sub>)
with magnitude scaling
(*α*<sub>*m*, 00</sub>, *α*<sub>*m*, 11</sub>, *α*<sub>*m*, 01</sub>, *α*<sub>*m*, 10</sub>).

Fit E uses naive (location-based) labels in the post-treatment fit; Fit
F uses adaptive SEM relabelling, then refits ETAS on the SEM-implied
labels. Therefore E/F differ primarily in latent process assignment, not
in the ETAS likelihood family.

Temporal truncation is set by *t*<sub>trunc</sub> = 12.6474266462322
days (source: auto\_from\_pre50(rel=0.050); relative rule=0.05). SEM
settings in this run: outer iterations=1, inner iterations=2000, inner
proposals=20, labellings=20. Bootstrap summaries below are based on the
stored replicate means from the same run.
