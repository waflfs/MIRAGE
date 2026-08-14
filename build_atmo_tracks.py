"""Build the ATMO 2020 JWST CMD table used by the JITTLES application.

The original packaged table was sampled only from 1--20 M_J.  That is
insufficient at 1 Gyr, where objects above 20 M_J are needed to cover the
800, 1200, and 2000 K portions of an ATMO isochrone.  This utility reads the
locally stored ATMO CEQ Vega NIRCam grid and samples its entire valid mass
range at 10, 100, and 1000 Myr.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GRID_DIRECTORY = (
    PROJECT_ROOT
    / "data"
    / "atmo_evolutionary_tracks"
    / "ATMO_CEQ"
    / "JWST_photometry"
    / "JWST_phot_NIRCAM_modA"
)
OUTPUT_PATHS = (
    PROJECT_ROOT / "atmo_theoretical_JWST_CMD.csv",
    Path(__file__).resolve().parent / "data" / "atmo_theoretical_JWST_CMD.csv",
)

# The raw ATMO evolution grid gives mass in solar masses and radius in solar
# radii.  These IAU nominal conversion factors place the output in M_J and
# R_J, matching the pre-existing table schema.
SOLAR_TO_JUPITER_MASS = 1047.348644
SOLAR_TO_JUPITER_RADIUS = 9.958
AGES_GYR = (0.01, 0.1, 1.0)
POINTS_PER_TRACK = 300


def _read_mass_tracks() -> tuple[np.ndarray, list[np.ndarray]]:
    """Return valid ATMO source rows, sorted by mass in Jupiter masses."""
    masses: list[float] = []
    tracks: list[np.ndarray] = []

    for path in GRID_DIRECTORY.glob("*_ATMO_CEQ_vega_NIRCAM.txt"):
        try:
            mass_mjup = float(path.name.split("_", maxsplit=1)[0]) * SOLAR_TO_JUPITER_MASS
        except ValueError:  # Ignore macOS resource-fork files, if present.
            continue

        table = np.loadtxt(path, comments="#")
        # Columns 18 and 27 are NIRCam F200W and F444W respectively.  A zero
        # records a missing atmosphere-grid magnitude rather than a magnitude
        # of zero, so it must not participate in interpolation.
        valid = np.isfinite(table).all(axis=1) & (table[:, 18] > 0.0) & (table[:, 27] > 0.0)
        table = table[valid]
        if table.size:
            masses.append(mass_mjup)
            tracks.append(table[np.argsort(table[:, 1])])

    order = np.argsort(masses)
    return np.asarray(masses)[order], [tracks[index] for index in order]


def _track_at_age(
    masses_mjup: np.ndarray, source_tracks: list[np.ndarray], age_gyr: float
) -> pd.DataFrame:
    """Interpolate the evolutionary sequences to one requested age."""
    source_values = np.full((len(source_tracks), 5), np.nan)

    for index, table in enumerate(source_tracks):
        if table[0, 1] <= age_gyr <= table[-1, 1]:
            # Interpolate in age using the same linear age convention used by
            # the original Species-generated table.
            source_values[index] = [
                np.interp(age_gyr, table[:, 1], table[:, column])
                for column in (2, 5, 4, 18, 27)
            ]

    valid = np.isfinite(source_values).all(axis=1)
    valid_masses = masses_mjup[valid]
    values = source_values[valid]
    target_masses = np.linspace(valid_masses.min(), valid_masses.max(), POINTS_PER_TRACK)
    interpolated = np.vstack(
        [np.interp(target_masses, valid_masses, values[:, column]) for column in range(values.shape[1])]
    ).T

    teff, logg, radius_rsun, f200w, f444w = interpolated.T
    return pd.DataFrame(
        {
            "age_myr": age_gyr * 1000.0,
            "mass_mjup": target_masses,
            "teff_K": teff,
            "logg_cgs": logg,
            "radius_rjup": radius_rsun * SOLAR_TO_JUPITER_RADIUS,
            "F200W_minus_F444W": f200w - f444w,
            "absolute_F444W": f444w,
        }
    )


def build_tracks() -> pd.DataFrame:
    """Build all three ATMO CMD tracks."""
    masses_mjup, source_tracks = _read_mass_tracks()
    return pd.concat(
        [_track_at_age(masses_mjup, source_tracks, age_gyr) for age_gyr in AGES_GYR],
        ignore_index=True,
    )


def main() -> None:
    table = build_tracks()
    for output_path in OUTPUT_PATHS:
        table.to_csv(output_path, index=False)
        print(f"Wrote {output_path} ({len(table)} rows)")


if __name__ == "__main__":
    main()
