# PPDisentangle 0.2.0

This release archives the major-revision analysis code and current Oklahoma
results workflow, superseding the July 29 v0.1.0 release.

- Publish the existing bivariate C/D fitting and bootstrap implementation,
  including the 250-day finite temporal kernel and matching simulation/refit law.
- Make observed AOI versus no treatment the explicit Oklahoma paper contrast;
  reject unavailable contrasts instead of silently substituting all-or-nothing.
- Regenerate partition sensitivity and full-precision count tables from saved fits.
- Freeze Census county geometry for offline paper maps.
- Resolve archived robustness files locally and support output paths with spaces.
- Generate and validate publication assets inside archive staging; include source
  provenance, checksums and a reproduction-session snapshot.
- Document archived settings, recentered bootstrap quantiles, and the different
  simulation/application intervention conventions. No saved fits are changed
  or rerun by this release preparation.
