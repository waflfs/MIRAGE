"""Build jittles/data/cosmos3d_reference.npz from the research catalog.

This script is intentionally separate from the app: the 11.9 GB FITS catalog is
read once, while the deployed app uses a compact, pickle-free reference pack.
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

from astropy.io import fits
import numpy as np


FILTERS = ("f115w", "f150w", "f200w", "f277w", "f356w", "f444w")
APERTURES = (0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0)
CORE_DIAMETERS = (0.05, 0.1, 0.2, 0.3)
COMPACTNESS_BOUNDS = {
    "f115w": (0.50, 0.75),
    "f150w": (0.50, 0.75),
    "f200w": (0.50, 0.70),
    "f277w": (0.30, 0.45),
    "f356w": (0.27, 0.37),
    "f444w": (0.23, 0.30),
}
SURVEY_AREA_ARCMIN2 = {
    "f115w": 1932.0,
    "f150w": 1930.0,
    "f200w": 1282.0,
    "f277w": 1930.0,
    "f356w": 1241.0,
    "f444w": 1936.0,
}


def ab_magnitude_from_microjy(flux_microjy: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return -2.5 * np.log10((np.asarray(flux_microjy, dtype=float) * 1e-6) / 3631.0)


def core_tag(core: float) -> str:
    return f"{core:g}".replace(".", "p")


def build_pack(catalog_path: Path, types_path: Path, ee_path: Path, output_path: Path) -> None:
    object_types = np.load(types_path, mmap_mode="r")
    nonstellar = np.asarray(object_types != 1, dtype=bool)
    galaxies = np.asarray(object_types == 0, dtype=bool)
    nonstellar_catalog_indices = np.flatnonzero(nonstellar)
    nonstellar_positions = np.full(nonstellar.shape, -1, dtype=np.int32)
    nonstellar_positions[nonstellar_catalog_indices] = np.arange(
        nonstellar_catalog_indices.size, dtype=np.int32
    )
    output: dict[str, np.ndarray] = {
        "schema_version": np.array("1.2"),
        "catalog_name": np.array("COSMOS-3D v1.2 merged"),
    }

    with ee_path.open("rb") as handle:
        encircled_energy = pickle.load(handle)

    cmd_values: dict[str, np.ndarray] = {}
    with fits.open(catalog_path, memmap=True, mode="readonly") as hdul:
        catalog = hdul[1].data
        if len(catalog) != len(object_types):
            raise ValueError(
                f"Catalog/type lengths differ: {len(catalog):,} != {len(object_types):,}"
            )

        for filter_name in FILTERS:
            print(f"Processing {filter_name.upper()}...", flush=True)
            auto_flux = np.asarray(catalog[f"f_auto_{filter_name}"], dtype=float)
            auto_error = np.asarray(catalog[f"e_auto_{filter_name}"], dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"):
                snr = auto_flux / auto_error
            total_magnitude = ab_magnitude_from_microjy(auto_flux)
            output[f"constraint_total_{filter_name}_ab"] = total_magnitude[nonstellar].astype(np.float32)

            aperture_flux = np.asarray(catalog[f"f_aper_nat_{filter_name}"], dtype=float)
            aperture_error = np.asarray(catalog[f"e_aper_nat_{filter_name}"], dtype=float)
            small_flux = aperture_flux[:, APERTURES.index(0.1)]
            large_flux = aperture_flux[:, APERTURES.index(0.3)]
            small_error = aperture_error[:, APERTURES.index(0.1)]
            large_error = aperture_error[:, APERTURES.index(0.3)]
            with np.errstate(divide="ignore", invalid="ignore"):
                compactness = small_flux / large_flux
                compactness_error = np.abs(compactness) * np.sqrt(
                    (large_error / large_flux) ** 2 + (small_error / small_flux) ** 2
                )
            low, high = COMPACTNESS_BOUNDS[filter_name]
            overlaps_bounds = (
                ((compactness > low) & (compactness < high))
                | ((compactness + compactness_error > low) & (compactness - compactness_error < high))
            )
            point_source = overlaps_bounds & (snr > 5)
            point_catalog_indices = np.flatnonzero(nonstellar & point_source & np.isfinite(total_magnitude))
            point_positions = nonstellar_positions[point_catalog_indices]
            point_order = np.argsort(total_magnitude[point_catalog_indices])
            point_magnitude = total_magnitude[point_catalog_indices][point_order].astype(np.float32)
            output[f"point_{filter_name}"] = point_magnitude
            output[f"constraint_point_index_{filter_name}"] = point_positions[point_order].astype(np.int32)

            for core in CORE_DIAMETERS:
                core_flux = aperture_flux[:, APERTURES.index(core)]
                core_magnitude_all = ab_magnitude_from_microjy(core_flux)[nonstellar]
                core_positions = np.flatnonzero(np.isfinite(core_magnitude_all))
                core_order = np.argsort(core_magnitude_all[core_positions])
                core_magnitude = core_magnitude_all[core_positions][core_order].astype(np.float32)
                output[f"core_{filter_name}_{core_tag(core)}"] = core_magnitude
                output[f"constraint_core_index_{filter_name}_{core_tag(core)}"] = (
                    core_positions[core_order].astype(np.int32)
                )
                output[f"ee_{filter_name}_{core_tag(core)}"] = np.array(
                    float(encircled_energy[f"{filter_name}_{core:g}"]), dtype=np.float64
                )

            output[f"area_{filter_name}"] = np.array(SURVEY_AREA_ARCMIN2[filter_name])
            if filter_name in ("f200w", "f444w"):
                cmd_values[filter_name] = total_magnitude.copy()

            del aperture_flux, aperture_error, total_magnitude, auto_flux, auto_error

    # The notebook's CMD uses galaxies, F444W S/N > 5, and Vega offsets.
    # Reopen only the F444W error/flux columns to reproduce that S/N selection.
    with fits.open(catalog_path, memmap=True, mode="readonly") as hdul:
        catalog = hdul[1].data
        f444 = np.asarray(catalog["f_auto_f444w"], dtype=float)
        e444 = np.asarray(catalog["e_auto_f444w"], dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            snr444 = f444 / e444
    color_vega = (cmd_values["f200w"] - 1.71) - (cmd_values["f444w"] - 3.25)
    magnitude_vega = cmd_values["f444w"] - 3.25
    output["constraint_cmd_color_vega"] = color_vega[nonstellar].astype(np.float32)
    output["constraint_cmd_magnitude_vega"] = magnitude_vega[nonstellar].astype(np.float32)
    cmd_mask = galaxies & (snr444 > 5) & np.isfinite(color_vega) & np.isfinite(magnitude_vega)
    # Match the binned-box CMD in cmd_stuff block 11 exactly.
    color_edges = np.arange(-1.0, 5.0 + 0.25, 0.25)
    magnitude_edges = np.arange(13.0, 27.0 + 0.5, 0.5)
    counts, _, _ = np.histogram2d(
        color_vega[cmd_mask], magnitude_vega[cmd_mask], bins=[color_edges, magnitude_edges]
    )
    sparse_color_index = np.digitize(color_vega[cmd_mask], color_edges) - 1
    sparse_magnitude_index = np.digitize(magnitude_vega[cmd_mask], magnitude_edges) - 1
    valid_bins = (
        (sparse_color_index >= 0)
        & (sparse_color_index < counts.shape[0])
        & (sparse_magnitude_index >= 0)
        & (sparse_magnitude_index < counts.shape[1])
    )
    sparse_counts = np.zeros(color_vega[cmd_mask].shape, dtype=float)
    sparse_counts[valid_bins] = counts[
        sparse_color_index[valid_bins], sparse_magnitude_index[valid_bins]
    ]
    sparse_mask = valid_bins & (sparse_counts < 5)
    output["cmd_color_edges"] = color_edges.astype(np.float32)
    output["cmd_magnitude_edges"] = magnitude_edges.astype(np.float32)
    output["cmd_counts"] = counts.astype(np.int32)
    output["cmd_area_arcmin2"] = np.array(1282.0)
    output["cmd_sample_size"] = np.array(int(np.count_nonzero(cmd_mask)))
    output["cmd_sparse_color"] = color_vega[cmd_mask][sparse_mask].astype(np.float32)
    output["cmd_sparse_magnitude"] = magnitude_vega[cmd_mask][sparse_mask].astype(np.float32)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, **output)
    print(f"Wrote {output_path.resolve()} ({output_path.stat().st_size / 1024**2:.1f} MiB)")


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, default=repo_root / "catalogs" / "cosmos_v1.2_merged.fits")
    parser.add_argument("--types", type=Path, default=repo_root / "catalogs" / "type_in_3d.npy")
    parser.add_argument("--ee", type=Path, default=repo_root / "EE_dict.pkl")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "data" / "cosmos3d_reference.npz",
    )
    args = parser.parse_args()
    build_pack(args.catalog, args.types, args.ee, args.output)


if __name__ == "__main__":
    main()
