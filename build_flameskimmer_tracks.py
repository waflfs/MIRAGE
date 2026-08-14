"""Build compact JWST CMD isochrones from the public Flame Skimmer release.

The source archives are available at https://doi.org/10.5281/zenodo.20030439.
This script uses the equilibrium-chemistry, [M/H] = 0.0, C/O-scaling = 1.0
evolution and JWST magnitude tables.
"""

from __future__ import annotations

import argparse
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.interpolate import RegularGridInterpolator


EVOLUTION_MEMBER = "flameskimmer_evolution_tables/eq/mh+0.0_co1.0.txt"
PHOTOMETRY_MEMBER = "flame_skimmer_publication/eq/JWST/mh+0.0_co1.0_mag.txt"
TARGET_AGES_GYR = (0.01, 0.1, 1.0)


def _read_table(archive_path: Path, member: str) -> np.ndarray:
    with ZipFile(archive_path) as archive:
        raw_table = archive.read(member)
    return np.genfromtxt(BytesIO(raw_table), names=True, delimiter="\t", encoding="utf-8")


def _photometry_interpolator(
    photometry: np.ndarray, field: str
) -> RegularGridInterpolator:
    temperatures = np.unique(photometry["Teff"])
    gravities = np.unique(photometry["logg"])
    values = np.full((temperatures.size, gravities.size), np.nan)
    temperature_index = {value: index for index, value in enumerate(temperatures)}
    gravity_index = {value: index for index, value in enumerate(gravities)}
    for row in photometry:
        values[temperature_index[row["Teff"]], gravity_index[row["logg"]]] = row[field]
    if not np.isfinite(values).all():
        raise ValueError(f"The Flame Skimmer {field} photometry grid is incomplete.")
    return RegularGridInterpolator(
        (temperatures, gravities),
        values,
        bounds_error=False,
        fill_value=np.nan,
    )


def build_tracks(evolution_archive: Path, photometry_archive: Path) -> np.ndarray:
    evolution = _read_table(evolution_archive, EVOLUTION_MEMBER)
    photometry = _read_table(photometry_archive, PHOTOMETRY_MEMBER)
    f200w = _photometry_interpolator(photometry, "F200W")
    f444w = _photometry_interpolator(photometry, "F444W")

    output_rows: list[tuple[float, ...]] = []
    for target_age_gyr in TARGET_AGES_GYR:
        for mass in np.unique(evolution["MassMjup"]):
            mass_track = evolution[evolution["MassMjup"] == mass]
            order = np.argsort(mass_track["AgeGyr"])
            mass_track = mass_track[order]
            ages = mass_track["AgeGyr"]
            if not ages.min() <= target_age_gyr <= ages.max():
                continue

            log_ages = np.log10(ages)
            target_log_age = np.log10(target_age_gyr)
            temperature = float(np.interp(target_log_age, log_ages, mass_track["Teff"]))
            gravity = float(np.interp(target_log_age, log_ages, mass_track["logg"]))
            radius = float(np.interp(target_log_age, log_ages, mass_track["RadiusRjup"]))
            interpolation_point = np.array([[temperature, gravity]])
            f200w_magnitude = float(f200w(interpolation_point)[0])
            f444w_magnitude = float(f444w(interpolation_point)[0])
            if not np.isfinite((f200w_magnitude, f444w_magnitude)).all():
                continue
            output_rows.append(
                (
                    target_age_gyr * 1000.0,
                    float(mass),
                    temperature,
                    gravity,
                    radius,
                    f200w_magnitude - f444w_magnitude,
                    f444w_magnitude,
                )
            )

    dtype = [
        ("age_myr", "f8"),
        ("mass_mjup", "f8"),
        ("teff_K", "f8"),
        ("logg_cgs", "f8"),
        ("radius_rjup", "f8"),
        ("F200W_minus_F444W", "f8"),
        ("absolute_F444W", "f8"),
    ]
    return np.array(output_rows, dtype=dtype)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evolution_archive", type=Path)
    parser.add_argument("photometry_archive", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "data" / "flameskimmer_JWST_CMD.csv",
    )
    args = parser.parse_args()

    tracks = build_tracks(args.evolution_archive, args.photometry_archive)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    header = ",".join(tracks.dtype.names or ())
    np.savetxt(
        args.output,
        np.column_stack([tracks[name] for name in tracks.dtype.names or ()]),
        delimiter=",",
        header=header,
        comments="",
        fmt=("%.0f", "%.6f", "%.3f", "%.5f", "%.6f", "%.6f", "%.6f"),
    )
    print(f"Wrote {tracks.size} Flame Skimmer CMD points to {args.output}")


if __name__ == "__main__":
    main()
