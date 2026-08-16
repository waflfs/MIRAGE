# MIRAGE

**Modeling Interloper Risk from Alignments with Galaxies in Exoplanet searches**

MIRAGE estimates the probability that one or more extragalactic sources could
be mistaken for a directly imaged planetary companion with JWST. It implements
the aperture-photometry method for calculation contamination described by Forney et al. (2026).

This is a public beta for exploratory scientific use. Results should be
interpreted with the assumptions and limitations documented below.

## What the app calculates

For a supplied contrast curve, MIRAGE calculates the cumulative probability of
at least one contaminant as a function of angular separation. It also reports
the expected contaminant count at the selected separation of interest and
within the plotted search radius. An optional second probability curve applies
user-selected color and magnitude constraints to the catalog sources.

The calculation divides the field at the user-selected region cutoff:

- **Contrast-limited region:** all nonstellar sources with finite flux in the
  selected core aperture are counted when their core-aperture magnitude is
  brighter than the corrected contrast limit. No point-source compactness cut
  is applied in this region, because PSF subtraction can make resolved sources
  appear point-like.
- **Background-limited region:** only nonstellar sources classified as
  point-like by the COSMOS-3D compactness selection are counted when their
  total-magnitude measurement is brighter than the contrast limit.

The transition is applied as `radius < region cutoff` for the inner regime and
`radius >= region cutoff` for the outer regime.

## Inputs and units

- A two-column contrast curve: separation in arcseconds followed by either
  flux ratio or delta magnitude. The parser accepts comma- or whitespace-
  separated files, comments beginning with `#`, and one header row.
- The filter associated with the curve: F115W, F150W, F200W, F277W, F356W,
  or F444W.
- Host-star magnitude in **AB magnitudes** in the selected contrast-curve
  filter.
- Region cutoff in arcseconds.
- Core-aperture diameter: 0.05, 0.1, 0.2, or 0.3 arcsec.
- System distance in parsecs. Distance is used for the planetary CMD tracks and
  for rescaling the fiducial Carter curve; an uploaded observation-specific
  curve is not distance-rescaled.
- Minimum and maximum plotted separation. The default is 0 to 10 arcsec,
  subject to the first separation in the supplied curve.
- Optional color range in Vega `F200W-F444W` and optional magnitude range in
  AB magnitudes in the selected contrast-curve filter. Each lower or upper
  bound may be left blank; ranges are inclusive.
- A planetary model for the CMD overlay: Sonora Bobcat, ATMO 2020, or Sonora
  Flame Skimmer.

## Contrast-curve processing

In automatic unit-detection mode, a contrast-curve second column with median
value below 0.1 is interpreted as a positive flux ratio and converted with
`delta_mag = -2.5 log10(flux ratio)`. Otherwise the second column is treated as
delta magnitude. Rows are sorted by separation, duplicate separations are
averaged, and the curve is linearly interpolated onto the calculation grid.
Beyond the curve's outermost measured separation, its final delta magnitude is
held constant. The maximum requested radius must be larger than the curve's
first separation.

For an uploaded curve, the values are assumed to describe the entered host star,
distance, and observing setup already. MIRAGE therefore uses the uploaded curve
unchanged apart from unit conversion, interpolation, and outer-radius
extrapolation.

## COSMOS-3D photometry and source selection

The packaged reference data are derived from the COSMOS-3D v1.2 merged catalog.
The app does not load the original catalog at runtime (that would be too large); it 
reads the `data/cosmos3d_reference.npz` reference file included in this repository.

The reference-pack builder uses the following catalog quantities in each of the
six filters:

- `f_auto_<filter>` and `e_auto_<filter>`: total flux and total-flux error in
  microJy. Total AB magnitude is calculated as
  `m_AB = -2.5 log10[(f_microJy * 10^-6) / 3631]`.
- `f_aper_nat_<filter>` and `e_aper_nat_<filter>`: aperture fluxes and errors
  for native aperture diameters 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, and 1.0
  arcsec.

The `type_in_3d` is the type parameter from COSMOS-Web 2025 matched by source
position to sources in COSMOS-3D. Objects with `type_in_3d == 1` are treated as 
stars and removed from the contamination counts. The remaining nonstellar 
objects (`type_in_3d != 1`) are retained for the inner-aperture arrays and for the 
color-magnitude constraint arrays.

### Point-source selection for the background-limited regime

For each filter, the point-source list is built using total-flux S/N and
compactness:

`S/N = f_auto / e_auto`

`compactness = f_aper(0.1 arcsec) / f_aper(0.3 arcsec)`

The compactness uncertainty is propagated as

`sigma_C = |C| sqrt[(sigma_0.3 / f_0.3)^2 + (sigma_0.1 / f_0.1)^2]`.

An object is retained as point-like when total-flux S/N is greater than 5 and
either its measured compactness lies inside the filter-specific interval or its
1-sigma compactness interval overlaps that interval. The intervals are:

| Filter | Compactness interval |
| F115W | 0.50--0.75 |
| F150W | 0.50--0.75 |
| F200W | 0.50--0.70 |
| F277W | 0.30--0.45 |
| F356W | 0.27--0.37 |
| F444W | 0.23--0.30 |

The resulting point-source magnitudes are total `f_auto` AB magnitudes. 
The inner-regime arrays are different, they use the selected `f_aper_nat` core 
aperture for every nonstellar source without applying the S/N or compactness selection.

### Encircled-energy corrections

The selected core-aperture magnitude limit is corrected using the encircled
energy value for the chosen filter and aperture. The values stored in the
packaged reference data are:

| Filter | EE(0.05) | EE(0.1) | EE(0.2) | EE(0.3) |
|---|---:|---:|---:|---:|
| F115W | 0.4330 | 0.6258 | 0.7752 | 0.8135 |
| F150W | 0.3386 | 0.6260 | 0.7686 | 0.8129 |
| F200W | 0.2339 | 0.5614 | 0.7019 | 0.8118 |
| F277W | 0.1701 | 0.3706 | 0.6554 | 0.7358 |
| F356W | 0.1155 | 0.2702 | 0.5944 | 0.6872 |
| F444W | 0.0827 | 0.2015 | 0.5171 | 0.6739 |

The survey areas used to convert counts to surface densities are 1932, 1930,
1282, 1930, 1241, and 1936 arcmin^2 for F115W, F150W, F200W, F277W, F356W,
and F444W, respectively. This was calculated by binning a position map by RA and Dec
and tallying the area of the boxes that contain sources.

## Probability calculation

For each radius, the code converts the contrast delta magnitude into two
possible limiting magnitudes. In the contrast-limited region, the core limit is

`m_core,lim = m_host,AB + delta_mag - 2.5 log10(EE)`.

In the background-limited region, the total-flux limit is

`m_total,lim = m_host,AB + delta_mag`.

The code counts catalog magnitudes strictly brighter than each limit. 
The cumulative count is divided by the filter-specific COSMOS-3D survey area 
to obtain a surface density in arcmin^-2. The expected number inside radius `r` 
is then calculated by integrating circular annuli:

`lambda(r) = integral[0 to r] 2 pi r' Sigma(r') dr'`,

where radius is converted from arcseconds to arcminutes before integration.
The probability of at least one source is computed assuming a Poisson spatial
process:

`P(N >= 1) = 1 - exp[-lambda(r)]`.

The expected number and probability are evaluated numerically with trapezoidal
integration on a 900-sample radial grid augmented with the supplied curve
separations, the region cutoff, and the requested minimum separation. In normal
contrast-curve mode, the grid begins at the curve's first measured separation;
the cumulative result is therefore zero interior to that starting separation.
Contrast-less mode instead includes radius zero explicitly.

## Fiducial Carter F444W curve

Selecting **Use fiducial Carter F444W curve** fixes the filter to F444W and
provides the Carter et al. (2023) curve. The code uses the Carter fiducial pair
`m_F444W = 9.76` AB at 109 pc, corresponding to the quoted absolute magnitude
4.57 AB. It identifies the faintest point of the magnitude limit curve as the
fiducial background limit, converts the total limit and background limit to
flux, and subtracts the background flux to obtain the host-dependent residual.
That residual is scaled by the target host's apparent-flux offset, the selected
background flux floor is added back, and the resulting flux limit is converted
to the new magnitude limit curve. The user must supply the background-limiting 
F444W magnitude for the new observation. This is an approximate transformation 
for planning or speculation; a measured observation-specific curve
should be used whenever available.

## Contrast-less mode

Contrast-less mode is intended for observations of very faint infrared hosts,
such as white dwarfs, where the field can be treated as background limited. It
ignores the contrast curve, host-star magnitude, region cutoff, and core-aperture
choice. The selected-filter point-source count is divided by that
filter's survey area to make a constant surface density, which is integrated
from zero separation. This mode is not equivalent to setting the region cutoff
to zero in a normal calculation, because a normal calculation still applies its
contrast-derived detection limit to the point-source population.

## Color-magnitude constraints and CMD

When enabled, the optional constraint mask is applied to the nonstellar
source-level arrays using Vega `F200W-F444W` color and total AB magnitude in the
selected filter. The mask is combined with the applicable core-aperture or
point-source detection list in both radial regimes, producing a second expected
number and probability curve. Blank lower or upper bounds are treated as
unbounded; if both bounds for a quantity are blank, that quantity imposes no
restriction.

The CMD itself uses only sources with `type_in_3d == 0` (the catalog galaxy
class). AB magnitudes are converted to Vega using offsets of 1.71 mag for
F200W and 3.25 mag for F444W. The displayed CMD uses 0.25-mag color bins from
-1 to 5 and 0.5-mag F444W bins from 13 to 27, divided by the 1282 arcmin^2
COSMOS-3D F200W area. Empty bins are white and populated bins are shown on a
logarithmic objects-per-arcmin^2 bin scale.

Planetary tracks are shown at 10, 100, and 1000 Myr and shifted from absolute
F444W magnitude to the entered distance using the distance modulus
`5 log10(d / 10 pc)`. The available overlays are Sonora Bobcat, ATMO 2020, and
Sonora Flame Skimmer.

## Packaged reference data

The repository already includes the compact reference data needed by the app in
`data/cosmos3d_reference.npz`. No additional catalogs, calibration files, or
source-data downloads are required to run MIRAGE.

## Run locally (For offline use)

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
streamlit run app.py
```

## Scope and limitations

The calculation assumes that the COSMOS-3D field represents the background
population around the target and that sources follow a spatial Poisson process.
It does not model galaxy clustering, target-field cosmic variance, extinction,
source variability, or candidate-specific astrometry. The transition radius,
compactness intervals, encircled-energy values, catalog completeness, and
contrast-curve calibration are all potential sources of systematic error. The
tool is intended for observation planning and candidate vetting; the probability
is not a classification of an individual candidate as planetary or extragalactic.
