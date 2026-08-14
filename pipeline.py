"""Scientific calculations for the JWST extragalactic contamination tool.

This module implements the aperture-photometry method described in
Forney et al. (2026 draft), not the legacy surface-brightness method in jit/.
"""

from __future__ import annotations

from dataclasses import dataclass
from io import BytesIO, StringIO
from pathlib import Path
from typing import BinaryIO, Literal

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import LogFormatterMathtext
import numpy as np
from scipy.integrate import cumulative_trapezoid


FILTERS = ("f115w", "f150w", "f200w", "f277w", "f356w", "f444w")
CORE_DIAMETERS = (0.05, 0.1, 0.2, 0.3)
ContrastKind = Literal["auto", "contrast", "delta_mag"]
CARTER_F444W_FIDUCIAL_APPARENT_MAG_AB = 9.76
CARTER_F444W_FIDUCIAL_ABSOLUTE_MAG_AB = 4.57
CARTER_F444W_FIDUCIAL_DISTANCE_PC = 109.0


@dataclass(frozen=True)
class ContaminationResult:
    radius_arcsec: np.ndarray
    delta_mag: np.ndarray
    limiting_magnitude: np.ndarray
    surface_density_arcmin2: np.ndarray
    expected_number: np.ndarray
    probability: np.ndarray
    region_cutoff_arcsec: float
    filter_name: str
    core_diameter_arcsec: float
    expected_number_color_magnitude: np.ndarray | None = None
    probability_color_magnitude: np.ndarray | None = None
    cmd_color_range: tuple[float | None, float | None] | None = None
    cmd_magnitude_range: tuple[float | None, float | None] | None = None
    contrast_less: bool = False
    plot_min_arcsec: float | None = None


def _asset_key(filter_name: str, core_diameter: float) -> str:
    core_tag = f"{core_diameter:g}".replace(".", "p")
    return f"core_{filter_name}_{core_tag}"


def _ee_key(filter_name: str, core_diameter: float) -> str:
    core_tag = f"{core_diameter:g}".replace(".", "p")
    return f"ee_{filter_name}_{core_tag}"


def load_reference_pack(path: str | Path) -> dict[str, np.ndarray]:
    """Load the compact, pickle-free COSMOS-3D reference pack."""
    with np.load(path, allow_pickle=False) as archive:
        return {name: archive[name] for name in archive.files}


def _read_bytes(source: str | Path | bytes | bytearray | BinaryIO) -> bytes:
    if hasattr(source, "read"):
        raw = source.read()
        if isinstance(raw, str):
            return raw.encode("utf-8")
        return raw
    if isinstance(source, (bytes, bytearray)):
        return bytes(source)
    return Path(source).read_bytes()


def load_contrast_curve(
    source: str | Path | bytes | bytearray | BinaryIO,
    kind: ContrastKind = "auto",
) -> tuple[np.ndarray, np.ndarray]:
    """Read a two-column contrast curve and return separation and delta-mag.

    The parser accepts comma- or whitespace-delimited files, comments, and a
    single header row. In ``auto`` mode, values with a median below 0.1 are
    interpreted as flux ratios; otherwise they are interpreted as delta-mag.
    """
    raw = _read_bytes(source)
    text = raw.decode("utf-8-sig", errors="replace")
    if not text.strip():
        raise ValueError("The contrast-curve file is empty.")

    first_data_line = ""
    for line in text.splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            first_data_line = stripped
            break
    delimiter = "," if "," in first_data_line else None
    data = np.genfromtxt(
        StringIO(text),
        comments="#",
        delimiter=delimiter,
        invalid_raise=False,
    )
    if data.ndim == 1:
        if data.size != 2:
            raise ValueError("The contrast curve must contain exactly two columns.")
        data = data.reshape(1, 2)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("The contrast curve must contain two numeric columns.")

    separation = np.asarray(data[:, 0], dtype=float)
    values = np.asarray(data[:, 1], dtype=float)
    finite = np.isfinite(separation) & np.isfinite(values) & (separation >= 0)
    separation = separation[finite]
    values = values[finite]
    if separation.size < 2:
        raise ValueError("The contrast curve needs at least two valid rows.")

    order = np.argsort(separation)
    separation = separation[order]
    values = values[order]
    unique_separation, inverse = np.unique(separation, return_inverse=True)
    if unique_separation.size != separation.size:
        sums = np.bincount(inverse, weights=values)
        counts = np.bincount(inverse)
        values = sums / counts
        separation = unique_separation

    inferred_kind = kind
    if kind == "auto":
        inferred_kind = "contrast" if np.nanmedian(values) < 0.1 else "delta_mag"
    if inferred_kind == "contrast":
        if np.any(values <= 0):
            raise ValueError("Flux-ratio contrast values must be positive.")
        values = -2.5 * np.log10(values)
    elif inferred_kind != "delta_mag":
        raise ValueError(f"Unknown contrast kind: {kind}")

    return separation, values


def rescale_carter_f444w_curve(
    fiducial_delta_mag: np.ndarray,
    host_star_magnitude_ab: float,
    system_distance_pc: float,
    background_limit_mag: float | None = None,
) -> tuple[np.ndarray, float, float, float, float]:
    """Rescale the Carter F444W curve with a selectable background floor.

    The Carter magnitude-limit curve is decomposed in flux space into a
    host-dependent residual and a fixed background term.  The residual is
    scaled by the target host's apparent-flux offset, expressed using its
    absolute magnitude and distance.  The returned delta-magnitude curve can
    be supplied directly to :func:`calculate_contamination` with the same
    target host magnitude.

    The original Carter floor is used for decomposition. ``background_limit_mag``
    replaces only the floor added back during recomposition; when omitted, the
    original Carter floor is retained.

    Returns ``(rescaled_delta_mag, target_absolute_mag, host_offset_mag,
    fiducial_background_limit_mag, applied_background_limit_mag)``.
    """
    fiducial_delta_mag = np.asarray(fiducial_delta_mag, dtype=float)
    if fiducial_delta_mag.ndim != 1 or fiducial_delta_mag.size < 2:
        raise ValueError("The Carter curve must contain at least two delta-magnitude values.")
    if not np.all(np.isfinite(fiducial_delta_mag)):
        raise ValueError("The Carter curve contains non-finite delta-magnitude values.")
    if not np.isfinite(host_star_magnitude_ab):
        raise ValueError("Host-star magnitude must be finite.")
    if not np.isfinite(system_distance_pc) or system_distance_pc <= 0:
        raise ValueError("System distance must be positive and finite.")

    target_absolute_mag = host_star_magnitude_ab - 5.0 * np.log10(
        system_distance_pc / 10.0
    )
    host_offset_from_absolute_mag = (
        target_absolute_mag - CARTER_F444W_FIDUCIAL_ABSOLUTE_MAG_AB
        + 5.0 * np.log10(system_distance_pc / CARTER_F444W_FIDUCIAL_DISTANCE_PC)
    )
    # The quoted Carter absolute magnitude is rounded. Calibrate the
    # absolute-magnitude/distance relation to the documented m_F444W=9.76
    # fiducial pair so it exactly reproduces the supplied digitized curve.
    fiducial_zero_offset = (
        CARTER_F444W_FIDUCIAL_ABSOLUTE_MAG_AB
        + 5.0 * np.log10(CARTER_F444W_FIDUCIAL_DISTANCE_PC / 10.0)
        - CARTER_F444W_FIDUCIAL_APPARENT_MAG_AB
    )
    host_offset_mag = host_offset_from_absolute_mag + fiducial_zero_offset

    fiducial_limit_mag = (
        CARTER_F444W_FIDUCIAL_APPARENT_MAG_AB + fiducial_delta_mag
    )
    fiducial_background_limit_mag = float(np.max(fiducial_limit_mag))
    if background_limit_mag is None:
        applied_background_limit_mag = fiducial_background_limit_mag
    else:
        applied_background_limit_mag = float(background_limit_mag)
        if not np.isfinite(applied_background_limit_mag):
            raise ValueError("Background limiting magnitude must be finite.")
    fiducial_background_flux_limit = 10.0 ** (
        -0.4 * fiducial_background_limit_mag
    )
    applied_background_flux_limit = 10.0 ** (
        -0.4 * applied_background_limit_mag
    )
    fiducial_total_flux_limit = 10.0 ** (-0.4 * fiducial_limit_mag)
    fiducial_pure_flux_limit = np.clip(
        fiducial_total_flux_limit - fiducial_background_flux_limit,
        0.0,
        None,
    )
    shifted_pure_flux_limit = fiducial_pure_flux_limit * 10.0 ** (
        -0.4 * host_offset_mag
    )
    shifted_total_flux_limit = (
        shifted_pure_flux_limit + applied_background_flux_limit
    )
    shifted_limit_mag = -2.5 * np.log10(shifted_total_flux_limit)
    rescaled_delta_mag = shifted_limit_mag - host_star_magnitude_ab
    return (
        rescaled_delta_mag,
        float(target_absolute_mag),
        float(host_offset_mag),
        fiducial_background_limit_mag,
        applied_background_limit_mag,
    )


def prepare_radius_grid(
    separation: np.ndarray,
    delta_mag: np.ndarray,
    max_radius_arcsec: float,
    region_cutoff_arcsec: float,
    samples: int = 900,
    min_radius_arcsec: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate the curve and hold its final value beyond its last point."""
    separation = np.asarray(separation, dtype=float)
    delta_mag = np.asarray(delta_mag, dtype=float)
    if max_radius_arcsec <= separation[0]:
        raise ValueError("Maximum radius must exceed the contrast curve's inner working angle.")
    if region_cutoff_arcsec < 0 or region_cutoff_arcsec > max_radius_arcsec:
        raise ValueError("Region cutoff must lie between 0 and the maximum radius.")
    requested_min = separation[0] if min_radius_arcsec is None else float(min_radius_arcsec)
    if not np.isfinite(requested_min) or requested_min < 0:
        raise ValueError("Minimum separation must be finite and non-negative.")
    plot_min = max(separation[0], requested_min)
    if plot_min >= max_radius_arcsec:
        raise ValueError("Minimum separation must be less than the maximum separation.")

    radius = np.linspace(separation[0], max_radius_arcsec, samples)
    radius = np.unique(
        np.concatenate([radius, separation, [region_cutoff_arcsec, plot_min]])
    )
    radius = radius[(radius >= separation[0]) & (radius <= max_radius_arcsec)]
    gridded_delta_mag = np.interp(
        radius,
        separation,
        delta_mag,
        left=delta_mag[0],
        right=delta_mag[-1],
    )
    return radius, gridded_delta_mag


def integrate_surface_density(
    radius_arcsec: np.ndarray,
    surface_density_arcmin2: np.ndarray,
) -> np.ndarray:
    """Integrate a radial surface-density profile into expected counts."""
    radius_arcmin = np.asarray(radius_arcsec, dtype=float) / 60.0
    density = np.asarray(surface_density_arcmin2, dtype=float)
    return cumulative_trapezoid(2.0 * np.pi * radius_arcmin * density, radius_arcmin, initial=0.0)


def poisson_probability(expected_number: np.ndarray) -> np.ndarray:
    return -np.expm1(-np.asarray(expected_number, dtype=float))


def _normalise_range(
    values: tuple[float | None, float | None] | None,
    name: str,
) -> tuple[float | None, float | None] | None:
    if values is None:
        return None
    if len(values) != 2:
        raise ValueError(f"{name} must contain two numbers or blank bounds.")
    lower_raw, upper_raw = values
    lower = None if lower_raw is None else float(lower_raw)
    upper = None if upper_raw is None else float(upper_raw)
    if (lower is not None and not np.isfinite(lower)) or (
        upper is not None and not np.isfinite(upper)
    ):
        raise ValueError(f"{name} must contain finite numbers or blank bounds.")
    if lower is None and upper is None:
        return None
    if lower is not None and upper is not None and lower > upper:
        raise ValueError(f"{name} minimum cannot exceed its maximum.")
    return lower, upper


def _constrained_counts(
    sorted_magnitudes: np.ndarray,
    source_indices: np.ndarray,
    constraint_mask: np.ndarray,
    magnitude_limit: np.ndarray,
) -> np.ndarray:
    """Count constrained sources below each detection-magnitude limit."""
    selected_sorted = np.asarray(constraint_mask[source_indices], dtype=bool)
    cumulative = np.cumsum(selected_sorted, dtype=np.int64)
    positions = np.searchsorted(sorted_magnitudes, magnitude_limit, side="left")
    counts = np.zeros(positions.shape, dtype=np.int64)
    valid = positions > 0
    counts[valid] = cumulative[positions[valid] - 1]
    return counts


def _cmd_constraint_mask(
    reference_pack: dict[str, np.ndarray],
    filter_name: str,
    color_range: tuple[float | None, float | None] | None,
    magnitude_range: tuple[float | None, float | None] | None,
) -> np.ndarray | None:
    """Return the per-source mask for any active color/magnitude bounds."""
    if color_range is None and magnitude_range is None:
        return None
    required_keys = ("constraint_cmd_color_vega", f"constraint_total_{filter_name}_ab")
    missing = [key for key in required_keys if key not in reference_pack]
    if missing:
        raise ValueError(
            "The reference pack lacks CMD source mappings. Rebuild it with build_reference_pack.py."
        )
    cmd_color = np.asarray(reference_pack[required_keys[0]], dtype=float)
    cmd_magnitude = np.asarray(reference_pack[required_keys[1]], dtype=float)
    if cmd_color.shape != cmd_magnitude.shape:
        raise ValueError("The CMD constraint arrays have inconsistent lengths.")

    constraint_mask = np.ones(cmd_color.shape, dtype=bool)
    if color_range is not None:
        constraint_mask &= np.isfinite(cmd_color)
        if color_range[0] is not None:
            constraint_mask &= cmd_color >= color_range[0]
        if color_range[1] is not None:
            constraint_mask &= cmd_color <= color_range[1]
    if magnitude_range is not None:
        constraint_mask &= np.isfinite(cmd_magnitude)
        if magnitude_range[0] is not None:
            constraint_mask &= cmd_magnitude >= magnitude_range[0]
        if magnitude_range[1] is not None:
            constraint_mask &= cmd_magnitude <= magnitude_range[1]
    return constraint_mask


def calculate_contamination(
    reference_pack: dict[str, np.ndarray],
    contrast_separation_arcsec: np.ndarray,
    contrast_delta_mag: np.ndarray,
    filter_name: str,
    host_star_magnitude_ab: float,
    region_cutoff_arcsec: float,
    core_diameter_arcsec: float,
    max_radius_arcsec: float = 10.0,
    cmd_color_range: tuple[float | None, float | None] | None = None,
    cmd_magnitude_range: tuple[float | None, float | None] | None = None,
    min_radius_arcsec: float | None = None,
) -> ContaminationResult:
    """Calculate contamination probability as a function of separation.

    If CMD color and magnitude ranges are supplied, a second result curve is
    computed using only catalog objects inside both ranges, in addition to the
    contrast/aperture detection limit.
    """
    filter_name = filter_name.lower()
    if filter_name not in FILTERS:
        raise ValueError(f"Unsupported filter: {filter_name}")
    if core_diameter_arcsec not in CORE_DIAMETERS:
        raise ValueError(f"Unsupported core diameter: {core_diameter_arcsec:g} arcsec")
    cmd_color_range = _normalise_range(cmd_color_range, "CMD color range")
    cmd_magnitude_range = _normalise_range(cmd_magnitude_range, "CMD magnitude range")

    radius, delta_mag = prepare_radius_grid(
        contrast_separation_arcsec,
        contrast_delta_mag,
        max_radius_arcsec,
        region_cutoff_arcsec,
        min_radius_arcsec=min_radius_arcsec,
    )
    requested_min = (
        radius[0]
        if min_radius_arcsec is None
        else max(float(np.asarray(contrast_separation_arcsec, dtype=float).min()), float(min_radius_arcsec))
    )

    core_magnitudes = np.asarray(
        reference_pack[_asset_key(filter_name, core_diameter_arcsec)], dtype=float
    )
    point_magnitudes = np.asarray(reference_pack[f"point_{filter_name}"], dtype=float)
    survey_area = float(np.asarray(reference_pack[f"area_{filter_name}"]).item())
    encircled_energy = float(
        np.asarray(reference_pack[_ee_key(filter_name, core_diameter_arcsec)]).item()
    )
    if not 0 < encircled_energy <= 1:
        raise ValueError("Encircled energy must be between zero and one.")

    aperture_limit = host_star_magnitude_ab + delta_mag - 2.5 * np.log10(encircled_energy)
    total_limit = host_star_magnitude_ab + delta_mag

    counts_inside = np.searchsorted(core_magnitudes, aperture_limit, side="left")
    counts_outside = np.searchsorted(point_magnitudes, total_limit, side="left")
    counts = np.where(radius < region_cutoff_arcsec, counts_inside, counts_outside)
    density = counts.astype(float) / survey_area
    expected = integrate_surface_density(radius, density)
    probability = poisson_probability(expected)
    limiting_magnitude = np.where(radius < region_cutoff_arcsec, aperture_limit, total_limit)

    expected_color_magnitude = None
    probability_color_magnitude = None
    constraint_mask = _cmd_constraint_mask(
        reference_pack,
        filter_name,
        cmd_color_range,
        cmd_magnitude_range,
    )
    if constraint_mask is not None:
        required_index_keys = (
            f"constraint_core_index_{filter_name}_{f'{core_diameter_arcsec:g}'.replace('.', 'p')}",
            f"constraint_point_index_{filter_name}",
        )
        missing = [key for key in required_index_keys if key not in reference_pack]
        if missing:
            raise ValueError(
                "The reference pack lacks CMD source mappings. Rebuild it with build_reference_pack.py."
            )
        core_indices = np.asarray(reference_pack[required_index_keys[0]], dtype=np.int64)
        point_indices = np.asarray(reference_pack[required_index_keys[1]], dtype=np.int64)
        constrained_counts_inside = _constrained_counts(
            core_magnitudes,
            core_indices,
            constraint_mask,
            aperture_limit,
        )
        constrained_counts_outside = _constrained_counts(
            point_magnitudes,
            point_indices,
            constraint_mask,
            total_limit,
        )
        constrained_counts = np.where(
            radius < region_cutoff_arcsec,
            constrained_counts_inside,
            constrained_counts_outside,
        )
        constrained_density = constrained_counts.astype(float) / survey_area
        expected_color_magnitude = integrate_surface_density(radius, constrained_density)
        probability_color_magnitude = poisson_probability(expected_color_magnitude)

    return ContaminationResult(
        radius_arcsec=radius,
        delta_mag=delta_mag,
        limiting_magnitude=limiting_magnitude,
        surface_density_arcmin2=density,
        expected_number=expected,
        probability=probability,
        region_cutoff_arcsec=region_cutoff_arcsec,
        filter_name=filter_name,
        core_diameter_arcsec=core_diameter_arcsec,
        expected_number_color_magnitude=(
            None if expected_color_magnitude is None else expected_color_magnitude
        ),
        probability_color_magnitude=(
            None if probability_color_magnitude is None else probability_color_magnitude
        ),
        cmd_color_range=cmd_color_range,
        cmd_magnitude_range=cmd_magnitude_range,
        plot_min_arcsec=requested_min,
    )


def calculate_contrastless_contamination(
    reference_pack: dict[str, np.ndarray],
    filter_name: str,
    max_radius_arcsec: float = 10.0,
    cmd_color_range: tuple[float | None, float | None] | None = None,
    cmd_magnitude_range: tuple[float | None, float | None] | None = None,
    samples: int = 900,
    min_radius_arcsec: float | None = None,
) -> ContaminationResult:
    """Calculate an entirely background-limited point-source contamination curve."""
    filter_name = filter_name.lower()
    if filter_name not in FILTERS:
        raise ValueError(f"Unsupported filter: {filter_name}")
    if max_radius_arcsec <= 0:
        raise ValueError("Maximum radius must be positive.")
    requested_min = 0.0 if min_radius_arcsec is None else float(min_radius_arcsec)
    if not np.isfinite(requested_min) or requested_min < 0:
        raise ValueError("Minimum separation must be finite and non-negative.")
    if requested_min >= max_radius_arcsec:
        raise ValueError("Minimum separation must be less than the maximum separation.")
    if samples < 2:
        raise ValueError("At least two radial samples are required.")
    cmd_color_range = _normalise_range(cmd_color_range, "CMD color range")
    cmd_magnitude_range = _normalise_range(cmd_magnitude_range, "CMD magnitude range")

    radius = np.unique(
        np.concatenate([np.linspace(0.0, max_radius_arcsec, samples), [requested_min]])
    )
    point_magnitudes = np.asarray(reference_pack[f"point_{filter_name}"], dtype=float)
    survey_area = float(np.asarray(reference_pack[f"area_{filter_name}"]).item())
    density_value = point_magnitudes.size / survey_area
    density = np.full(radius.shape, density_value, dtype=float)
    expected = integrate_surface_density(radius, density)
    probability = poisson_probability(expected)

    expected_color_magnitude = None
    probability_color_magnitude = None
    constraint_mask = _cmd_constraint_mask(
        reference_pack,
        filter_name,
        cmd_color_range,
        cmd_magnitude_range,
    )
    if constraint_mask is not None:
        point_index_key = f"constraint_point_index_{filter_name}"
        if point_index_key not in reference_pack:
            raise ValueError(
                "The reference pack lacks CMD source mappings. Rebuild it with build_reference_pack.py."
            )
        point_indices = np.asarray(reference_pack[point_index_key], dtype=np.int64)
        constrained_count = int(np.count_nonzero(constraint_mask[point_indices]))
        constrained_density = np.full(
            radius.shape,
            constrained_count / survey_area,
            dtype=float,
        )
        expected_color_magnitude = integrate_surface_density(radius, constrained_density)
        probability_color_magnitude = poisson_probability(expected_color_magnitude)

    return ContaminationResult(
        radius_arcsec=radius,
        delta_mag=np.full(radius.shape, np.nan),
        limiting_magnitude=np.full(radius.shape, np.nan),
        surface_density_arcmin2=density,
        expected_number=expected,
        probability=probability,
        region_cutoff_arcsec=0.0,
        filter_name=filter_name,
        core_diameter_arcsec=np.nan,
        expected_number_color_magnitude=(
            None if expected_color_magnitude is None else expected_color_magnitude
        ),
        probability_color_magnitude=(
            None if probability_color_magnitude is None else probability_color_magnitude
        ),
        cmd_color_range=cmd_color_range,
        cmd_magnitude_range=cmd_magnitude_range,
        contrast_less=True,
        plot_min_arcsec=requested_min,
    )


def probability_figure(result: ContaminationResult):
    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=150, layout="constrained")
    ax.plot(
        result.radius_arcsec,
        result.probability,
        color="#2463eb",
        linewidth=2.4,
        label=("All point-source contaminants" if result.contrast_less else "All detectable contaminants"),
    )
    if result.probability_color_magnitude is not None:
        ax.plot(
            result.radius_arcsec,
            result.probability_color_magnitude,
            color="#d62728",
            linewidth=2.1,
            label="With color-magnitude constraints",
        )
    plot_min = (
        result.radius_arcsec[0]
        if result.plot_min_arcsec is None
        else float(result.plot_min_arcsec)
    )
    if (
        not result.contrast_less
        and plot_min <= result.region_cutoff_arcsec <= result.radius_arcsec[-1]
    ):
        ax.axvline(
            result.region_cutoff_arcsec,
            color="#e85d3f",
            linestyle="--",
            linewidth=1.5,
            label="Region cutoff",
        )
    ax.set(
        xlabel="Separation (arcsec)",
        ylabel="Probability of at least one contaminant",
    )
    ax.set_yscale("log")
    ax.set_ylim(bottom=1e-4)
    ax.set_xlim(plot_min, result.radius_arcsec[-1])
    if result.contrast_less:
        ax.text(
            0.5,
            0.98,
            "Background limited",
            transform=ax.transAxes,
            ha="center",
            va="top",
            color="#555555",
            fontsize=10,
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 2.0},
            zorder=5,
        )
    else:
        left_edge = plot_min
        right_edge = float(result.radius_arcsec[-1])
        cutoff = float(result.region_cutoff_arcsec)
        if cutoff > left_edge:
            ax.text(
                left_edge + 0.35 * (cutoff - left_edge),
                0.98,
                "Contrast limited",
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                color="#555555",
                fontsize=10,
                bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 2.0},
                zorder=5,
            )
        if cutoff < right_edge:
            ax.text(
                cutoff + 0.72 * (right_edge - cutoff),
                0.98,
                "Background limited",
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                color="#555555",
                fontsize=10,
                bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 2.0},
                zorder=5,
            )
    ax.grid(alpha=0.22)
    ax.legend(frameon=False, loc="lower right")
    return fig


def load_planet_tracks(track_source: str | Path, system_distance_pc: float) -> dict[float, np.ndarray]:
    """Load planetary-model CMD tracks at the user's system distance."""
    if system_distance_pc <= 0:
        raise ValueError("System distance must be positive.")

    distance_modulus = 5.0 * np.log10(system_distance_pc / 10.0)
    ages_to_plot = {10.0, 100.0, 1000.0}
    tracks: dict[float, np.ndarray] = {}

    source_path = Path(track_source)
    if source_path.suffix.lower() == ".npz":
        with np.load(source_path, allow_pickle=False) as archive:
            reference_keys: dict[float, tuple[float, str]] = {}
            for key in archive.files:
                if not key.endswith("_data") or "myr_" not in key or "pcs_data" not in key:
                    continue
                age_text, distance_text = key.removesuffix("_data").split("myr_")
                age = float(age_text)
                distance = float(distance_text.removesuffix("pcs"))
                if age not in ages_to_plot:
                    continue
                if age not in reference_keys or distance < reference_keys[age][0]:
                    reference_keys[age] = (distance, key)
            for age, (reference_distance, key) in reference_keys.items():
                cached = np.asarray(archive[key], dtype=float)
                data = np.column_stack(
                    [
                        cached[:, 0],
                        cached[:, 1] - 5.0 * np.log10(reference_distance / 10.0) + distance_modulus,
                        cached[:, 2],
                        cached[:, 3],
                    ]
                )
                data = data[np.isfinite(data).all(axis=1)]
                tracks[age] = data[np.argsort(data[:, 3])]
        return tracks

    table = np.genfromtxt(source_path, delimiter=",", names=True, encoding="utf-8")
    for age in np.unique(table["age_myr"]):
        if float(age) not in ages_to_plot:
            continue
        mask = table["age_myr"] == age
        data = np.column_stack(
            [
                table["F200W_minus_F444W"][mask],
                table["absolute_F444W"][mask] + distance_modulus,
                table["teff_K"][mask],
                table["mass_mjup"][mask],
            ]
        )
        data = data[np.isfinite(data).all(axis=1)]
        tracks[float(age)] = data[np.argsort(data[:, 3])]
    return tracks


def cmd_figure(
    reference_pack: dict[str, np.ndarray],
    planet_tracks: dict[float, np.ndarray],
    system_distance_pc: float,
):
    """Plot the binned COSMOS-3D CMD with age-styled planetary-model tracks."""
    x_edges = np.asarray(reference_pack["cmd_color_edges"], dtype=float)
    y_edges = np.asarray(reference_pack["cmd_magnitude_edges"], dtype=float)
    counts = np.asarray(reference_pack["cmd_counts"], dtype=float)
    survey_area = float(
        np.asarray(reference_pack.get("cmd_area_arcmin2", reference_pack["area_f200w"])).item()
    )
    rate = counts / survey_area
    rate_masked = np.ma.masked_less_equal(rate, 0.0)
    positive_rates = rate[rate > 0]
    if positive_rates.size == 0:
        raise ValueError("The CMD reference population has no populated bins.")
    base_colors = plt.cm.viridis(np.linspace(0.0, 1.0, 256))
    base_colors[0] = [1.0, 1.0, 1.0, 1.0]
    cmap = LinearSegmentedColormap.from_list("white_viridis", base_colors)
    cmap = cmap.with_extremes(bad="white")

    fig, ax = plt.subplots(figsize=(7.4, 5.5), dpi=150, layout="constrained")
    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        rate_masked.T,
        cmap=cmap,
        norm=LogNorm(vmin=float(positive_rates.min()), vmax=float(positive_rates.max())),
        shading="flat",
        rasterized=True,
    )
    colorbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    colorbar.set_label(r"Objects arcmin$^{-2}$ per bin")
    colorbar.ax.yaxis.set_major_formatter(LogFormatterMathtext())
    colorbar.update_ticks()

    age_styles = {10.0: "-", 100.0: "--", 1000.0: ":"}
    track_color = "C0"
    for age, data in sorted(planet_tracks.items()):
        if data.size == 0:
            continue
        data = data[data[:, 2] <= 2400.0]
        if data.size == 0:
            continue
        ax.plot(
            data[:, 0],
            data[:, 1],
            color=track_color,
            linestyle=age_styles.get(age, "-"),
            linewidth=1.7,
            zorder=2,
        )
        temperature = data[:, 2]
        order = np.argsort(temperature)
        temperature = temperature[order]
        color_track = data[order, 0]
        magnitude_track = data[order, 1]
        temperature, unique_index = np.unique(temperature, return_index=True)
        color_track = color_track[unique_index]
        magnitude_track = magnitude_track[unique_index]
        for requested_temperature in (2000.0, 1200.0, 800.0):
            if requested_temperature < temperature.min() or requested_temperature > temperature.max():
                continue
            color_at_temperature = np.interp(requested_temperature, temperature, color_track)
            magnitude_at_temperature = np.interp(requested_temperature, temperature, magnitude_track)
            ax.scatter(
                color_at_temperature,
                magnitude_at_temperature,
                s=30,
                facecolors="white",
                edgecolors="k",
                linewidths=1.2,
                zorder=3,
            )
            if age == 10.0:
                ax.text(
                    color_at_temperature + 0.08,
                    magnitude_at_temperature - 0.05,
                    f"{requested_temperature:g} K",
                    fontsize=12,
                    color="k",
                    zorder=4,
                )

    age_handles = [
        Line2D([0], [0], color="k", lw=1.7, linestyle=age_styles.get(age, "-"), label=f"{age:g} Myr")
        for age in sorted(planet_tracks)
    ]
    ax.legend( 
        handles=age_handles,
        loc="upper left",
        title="Age",
        frameon=True,
    )
    ax.text(
        0.98,
        0.98,
        f"Distance: {system_distance_pc:g} pc",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=14,
        bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none", "pad": 2.5},
        zorder=5,
    )

    ax.set(
        xlabel="F200W-F444W",
        ylabel=r"$m_{\mathrm{F444W}}$",
        xlim=(-1.0, 4.5),
        ylim=(10.0, 27.0),
    )
    ax.invert_yaxis()
    ax.grid(alpha=0.15)
    return fig


def result_csv_bytes(result: ContaminationResult) -> bytes:
    output = BytesIO()
    if result.contrast_less:
        columns = [
            result.radius_arcsec,
            result.surface_density_arcmin2,
            result.expected_number,
            result.probability,
        ]
        header_fields = [
            "radius_arcsec",
            "surface_density_arcmin2",
            "expected_number",
            "probability",
        ]
    else:
        columns = [
            result.radius_arcsec,
            result.delta_mag,
            result.limiting_magnitude,
            result.surface_density_arcmin2,
            result.expected_number,
            result.probability,
        ]
        header_fields = [
            "radius_arcsec",
            "delta_mag",
            "limiting_ab_mag",
            "surface_density_arcmin2",
            "expected_number",
            "probability",
        ]
    if (
        result.expected_number_color_magnitude is not None
        and result.probability_color_magnitude is not None
    ):
        columns.extend([result.expected_number_color_magnitude, result.probability_color_magnitude])
        header_fields.extend(
            ["expected_number_color_magnitude_constraint", "probability_color_magnitude_constraint"]
        )
    table = np.column_stack(columns)
    header = ",".join(header_fields)
    np.savetxt(output, table, delimiter=",", header=header, comments="")
    return output.getvalue()


def figure_png_bytes(figure, dpi: int = 300) -> bytes:
    """Serialize a Matplotlib figure for a high-resolution download."""
    output = BytesIO()
    figure.savefig(output, format="png", dpi=dpi, bbox_inches="tight")
    return output.getvalue()
