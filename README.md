# MIRAGE

**Modeling Interloper Risk from Alignments with Galaxies in Exoplanet searches**

MIRAGE estimates the probability that one or more extragalactic sources could
be mistaken for a directly imaged planetary companion with JWST. It implements
the aperture-photometry method for calculation contamination described by Forney et al. (2026).

This is a public beta for exploratory scientific use. Results should be
interpreted with the assumptions and limitations documented below.

## Using MIRAGE

The hosted Streamlit application is available at
[mirage.streamlit.app](https://mirage.streamlit.app/):

[![Open MIRAGE](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://mirage.streamlit.app/)

In the hosted app, choose the observing inputs in the sidebar, select
**Calculate**, and download the probability table or figures from the results
page. The repository also contains everything needed for an offline run:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
streamlit run app.py
```

The app does not require the original COSMOS-3D catalog or the raw planetary
model archives at runtime; the compact reference and track files in `data/`
are already packaged with the repository. See [Inputs and units](#inputs-and-units)
before interpreting a result.

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
- In contrast-less mode, a required background limiting magnitude in **AB
  magnitudes** in the selected filter. This is the total-magnitude detection
  threshold applied to the point-source population.
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
from zero separation. Only point-like sources with total AB magnitudes strictly
brighter than the user-entered background limiting magnitude are included. The
optional color-magnitude constraint is then applied as an additional mask to
that detection-limited population. This mode is not equivalent to setting the
region cutoff to zero in a normal calculation, because a normal calculation
still applies its contrast-derived detection limit to the point-source
population.

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

## Planetary-model grids and assumptions

The planetary tracks are diagnostic overlays on the CMD; they do not enter the
extragalactic-contamination probability. MIRAGE displays the 10, 100, and
1000 Myr isochrones and marks the 2000, 1200, and 800 K locations when those 
temperatures exist in a track. The model files provide synthetic Vega magnitudes in 
F200W and F444W. The color is therefore unchanged with distance, while the F444W 
magnitude is shifted from absolute magnitude to the entered system distance with
`5 log10(d / 10 pc)`. No claim is made that an observed companion follows one
of these tracks exactly; age, mass, composition, clouds, chemistry, and model
physics can all move a real object away from the tracks.

### Sonora Bobcat

The Bobcat overlay is derived from the public Sonora Bobcat evolution and
photometry release ([Marley et al. 2021](https://doi.org/10.3847/1538-4357/ac141d);
[Zenodo data release](https://doi.org/10.5281/zenodo.5063476)). The selected
branch is the cloudless, non-irradiated, rainout chemical-equilibrium grid at
solar metallicity (`[M/H] = 0.0`) and solar C/O composition. The published
model family spans approximately 200--2400 K and `3.25 <= log(g) <= 5.5` in
cgs units; the CMD file is a compact subset prepared from its JWST synthetic
photometry and evolution tables.

The packaged NPZ contains precomputed arrays at several reference distances.
At runtime, MIRAGE uses the smallest-distance reference array for each age
(currently the 10 pc array) and applies the distance modulus above. It does not
recompute the atmosphere, interpolate between metallicities, add clouds, or
apply a non-equilibrium chemistry correction.

### ATMO 2020

The ATMO overlay uses the public ATMO 2020 evolutionary/atmosphere grids
([Phillips et al. 2020](https://doi.org/10.1051/0004-6361/201937381)). ATMO
2020 couples a one-dimensional radiative--convective atmosphere to an interior
evolution calculation for roughly `0.001--0.075 M_sun` objects. The source
release provides equilibrium-chemistry (CEQ), weak-mixing, and strong-mixing
non-equilibrium branches; MIRAGE deliberately uses **CEQ only**, from files in
`ATMO_CEQ/JWST_photometry/JWST_phot_NIRCAM_modA`.

`build_atmo_tracks.py` reads the source mass tracks, rejects missing or zero
F200W/F444W magnitudes, interpolates each valid sequence to 10, 100, and
1000 Myr, and samples 300 uniformly spaced masses across the valid mass range.
The F200W and F444W values are the source table's Vega synthetic photometry.
The published atmosphere grid nominally spans `T_eff = 200--3000 K` and
`log(g) = 2.5--5.5`, but the authors caution that the atmosphere is most
reliable below about 2000 K because some high-temperature opacity sources are
not included. MIRAGE preserves the supplied track values (subject to the CMD
plot temperature cutoff) and only applies the distance-modulus shift.

### Sonora Flame Skimmer

The Flame Skimmer overlay is built from the public release documented at
[`doi:10.5281/zenodo.20030439`](https://doi.org/10.5281/zenodo.20030439). MIRAGE
selects the cloud-free, equilibrium-chemistry branch with `[M/H] = 0.0` and
C/O scaling `= 1.0` (the solar C/O branch, approximately C/O = 0.458 in the
release). `build_flameskimmer_tracks.py` interpolates the evolution tables in
log age at 10, 100, and 1000 Myr, then interpolates the corresponding JWST
F200W/F444W magnitudes over the release's `(T_eff, log(g))` photometry grid.
Other Flame Skimmer metallicities, C/O ratios, and non-equilibrium branches
are not used by the current MIRAGE overlay.

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

## Credits

- **Tomas Forney, Caitlin Casey, Brendan Bowler** — scientific design, methodology, and project direction.
- **OpenAI Codex** — software and deployment assistance.
