# JWST Extragalactic Contamination Tool

Public beta release. The hosted application is intended for exploratory
scientific use; please report issues and cite the accompanying Forney et al.
work when using its results.

This tool estimates
the probability that at least one extragalactic background source falls within a
JWST direct-imaging search radius. It supersedes the surface-brightness method in
`../jit` with the aperture-photometry framework in the Forney et al. (2026) draft.

## Method

For an uploaded two-column contrast curve, the tool calculates a limiting
magnitude as a function of angular separation. The field is split at a user-set
region cutoff:

- Inside the cutoff (contrast-limited), all nonstellar COSMOS-3D sources whose
  core-aperture magnitude is brighter than the encircled-energy-corrected limit
  are counted.
- Outside the cutoff (background-limited), only nonstellar COSMOS-3D sources
  consistent with the empirical point-source compactness bounds and brighter
  than the total-magnitude limit are counted.

The radial surface density is integrated in circular annuli to obtain the
expected contaminant count, lambda. The reported probability is
`P(X >= 1) = 1 - exp(-lambda)`.

Contrast-less mode is intended for observations of very faint infrared hosts,
such as white dwarfs. It ignores the contrast curve and host magnitude, treats
the full field as background limited, and integrates the selected filter's full
COSMOS-3D point-source surface density from zero separation.

The second output is the Vega-system COSMOS-3D F200W-F444W versus F444W
color-magnitude diagram. The user can overlay Sonora Bobcat, ATMO 2020, or
Sonora Flame Skimmer planetary tracks shifted to the system distance. The
packaged Flame Skimmer tracks use equilibrium chemistry, [M/H] = 0.0, and
C/O = 0.458 (solar) from the
[official model release](https://doi.org/10.5281/zenodo.20030439).

The interface provides CSV output plus 300 DPI PNG downloads for both plots.
The CMD uses fixed 0.25-mag color by 0.5-mag F444W bins, normalized to the
1282 arcmin² COSMOS-3D area, with empty bins shown in white and a logarithmic
objects-per-arcmin²-per-bin colorbar.

## Run

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
streamlit run app.py
```

Select **Use fiducial Carter F444W curve** to use the packaged Carter et al.
(2023) curve. Its flux-space host-dependent residual is rescaled using the
entered host-star apparent magnitude and system distance; the user selects the
background limiting magnitude added back during recomposition. Other filters
require an uploaded curve, which is used unchanged. Host-star magnitudes must
be AB magnitudes in the selected filter.

## Inputs

- Contrast curve: separation in arcsec plus flux ratio or delta-mag
- Filter corresponding to that curve
- Optional fiducial Carter F444W curve, with a selectable F444W AB background limiting magnitude
- Optional contrast-less mode; the filter remains required but contrast-specific inputs are omitted
- Host-star AB magnitude in that filter
- Optional CMD color range (F200W-F444W, Vega) and magnitude range in the selected contrast-curve filter (AB); either bound may be left blank
- Region cutoff in arcsec
- Core-aperture diameter in arcsec (0.05, 0.1, 0.2, or 0.3)
- System distance in parsecs
- Planetary model for the CMD overlay (Sonora Bobcat, ATMO 2020, or Sonora Flame Skimmer)
- Optional minimum and maximum plotted radii (0 and 10 arcsec by default)

## Reference pack

The app reads `data/cosmos3d_reference.npz`, a compact deployment artifact made
from `catalogs/cosmos_v1.2_merged.fits`, `catalogs/type_in_3d.npy`, and
`EE_dict.pkl`. Rebuild it from the repository root with:

```powershell
python jittles/build_reference_pack.py
```

The pack contains sorted magnitude arrays, per-filter total magnitudes,
encircled energies, survey areas, a binned CMD, and compact source-index
mappings for applying CMD constraints. It does not contain coordinates or
source identifiers.

## Scientific conventions inherited from the notebooks

- COSMOS-3D source detection requires total-flux S/N > 5 for point-source
  classification.
- Compactness is the 0.1 arcsec diameter aperture flux divided by the 0.3
  arcsec diameter aperture flux; a source is point-like when its value or its
  1-sigma interval overlaps the filter-specific empirical bounds.
- Matched stars (`type_in_3d == 1`) are removed from contaminant counts.
- Contrast-less mode uses the complete nonstellar point-source population and the selected filter's survey area, without applying a contrast-derived magnitude limit.
- The fiducial Carter F444W curve is decomposed into a host-dependent flux term
  and its original background floor. The host-dependent term is rescaled for
  the entered host magnitude and distance, then recomposed with the selected
  background limiting magnitude; uploaded curves are not rescaled.
- The final contrast value is held constant beyond the measured curve.
- When enabled, the color-magnitude constraint curve applies the entered Vega F200W-F444W color range and selected-filter AB magnitude range to catalog objects in both counting regimes.
- Survey areas are 1932, 1930, 1282, 1930, 1241, and 1936 arcmin² for F115W
  through F444W, respectively.
