from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from pipeline import (
    CARTER_F444W_FIDUCIAL_APPARENT_MAG_AB,
    CORE_DIAMETERS,
    FILTERS,
    calculate_contrastless_contamination,
    calculate_contamination,
    cmd_figure,
    figure_png_bytes,
    load_contrast_curve,
    load_planet_tracks,
    load_reference_pack,
    probability_figure,
    rescale_carter_f444w_curve,
    result_csv_bytes,
)


APP_DIR = Path(__file__).resolve().parent
REFERENCE_PACK = APP_DIR / "data" / "cosmos3d_reference.npz"
DEFAULT_CARTER_CURVE = APP_DIR / "defaults" / "carter_f444w.csv"
PLANET_MODELS = {
    "Sonora Bobcat": {
        "path": APP_DIR / "data" / "sonora_bobcat_JWST_CMD_arrays.npz",
        "details": "solar-metallicity equilibrium chemistry",
    },
    "ATMO 2020": {
        "path": APP_DIR / "data" / "atmo_theoretical_JWST_CMD.csv",
        "details": "equilibrium chemistry",
    },
    "Sonora Flame Skimmer": {
        "path": APP_DIR / "data" / "flameskimmer_JWST_CMD.csv",
        "details": "equilibrium chemistry, [M/H] = 0.0, and C/O = 0.458 (solar)",
    },
}


def _parse_optional_bound(raw_value: str, label: str) -> float | None:
    value = raw_value.strip()
    if not value:
        return None
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(f"{label} must be a number or left blank.") from exc
    if not np.isfinite(parsed):
        raise ValueError(f"{label} must be finite or left blank.")
    return parsed


def _optional_range(lower_raw: str, upper_raw: str, label: str):
    lower = _parse_optional_bound(lower_raw, f"{label} minimum")
    upper = _parse_optional_bound(upper_raw, f"{label} maximum")
    if lower is None and upper is None:
        return None
    return lower, upper


st.set_page_config(page_title="MIRAGE", layout="wide")
st.title("MIRAGE")
st.caption(
    "Modeling Interloper Risk from Alignments with Galaxies in Exoplanet searches "
    "| Extragalactic contamination probabilities for JWST direct imaging"
)


@st.cache_resource(show_spinner="Loading the COSMOS-3D reference population...")
def get_reference_pack():
    return load_reference_pack(REFERENCE_PACK)


@st.cache_data
def get_carter_fiducial_background_limit() -> float:
    _, delta_mag = load_contrast_curve(DEFAULT_CARTER_CURVE)
    return float(CARTER_F444W_FIDUCIAL_APPARENT_MAG_AB + np.max(delta_mag))


st.sidebar.markdown(
    "[Read the MIRAGE README on GitHub](https://github.com/waflfs/MIRAGE)"
)
st.sidebar.header("Observation")
contrast_less = st.sidebar.checkbox(
    "Contrast-less",
    value=False,
    help=(
        "Treat the entire field as background limited and use the full selected-filter "
        "COSMOS-3D point-source density. Changes take effect when you press Calculate."
    ),
)
use_carter_curve = False
if not contrast_less:
    use_carter_curve = st.sidebar.checkbox(
        "Use fiducial Carter F444W curve",
        value=False,
        help=(
            "Use the packaged Carter et al. (2023) F444W curve. Its host-dependent "
            "component is rescaled for this system, while you set the background limit."
        ),
    )
with st.sidebar.form("inputs"):
    if use_carter_curve:
        filter_name = "f444w"
        st.caption("Contrast-curve filter: **F444W** (fiducial Carter curve)")
    else:
        filter_name = st.selectbox(
            "Filter" if contrast_less else "Contrast-curve filter",
            FILTERS,
            index=5,
            format_func=str.upper,
            key="filter_name",
        )
    if contrast_less:
        uploaded_curve = None
        contrast_kind = "auto"
        carter_background_limit = None
        host_magnitude = np.nan
        region_cutoff = 0.0
        core_diameter = 0.2
        st.caption("Contrast curve, host magnitude, region cutoff, and core aperture are not used.")
    elif use_carter_curve:
        uploaded_curve = None
        contrast_kind = "auto"
        carter_background_limit = st.number_input(
            "Background limiting magnitude (F444W, AB)",
            min_value=0.0,
            value=get_carter_fiducial_background_limit(),
            step=0.05,
            format="%.3f",
            help=(
                "This replaces the background floor added back after decomposing the "
                "fiducial Carter curve in flux space."
            ),
        )
        host_magnitude = st.number_input(
            "Host-star magnitude (AB)", value=9.76, step=0.05, format="%.3f"
        )
    else:
        uploaded_curve = st.file_uploader(
            "Contrast curve",
            type=("csv", "txt", "dat"),
            help="Two columns: separation in arcsec, then flux ratio or delta magnitude.",
        )
        contrast_kind = st.selectbox(
            "Second-column units",
            ("auto", "contrast", "delta_mag"),
            format_func=lambda value: {
                "auto": "Detect automatically",
                "contrast": "Flux ratio",
                "delta_mag": "Delta magnitude",
            }[value],
        )
        host_magnitude = st.number_input(
            "Host-star magnitude (AB)", value=9.76, step=0.05, format="%.3f"
        )
    with st.expander("Color-magnitude constraints"):
        apply_cmd_constraint = st.checkbox(
            "Apply color-magnitude constraints",
            value=False,
            help="Add a second probability curve using only catalog objects inside both ranges.",
        )
        st.caption("Ranges are inclusive. Color is Vega F200W-F444W; magnitude is AB in the selected filter.")
        st.markdown("**Color range (F200W-F444W, Vega)**")
        color_min_col, color_dash_col, color_max_col = st.columns([1, 0.15, 1])
        cmd_color_min = color_min_col.text_input(
            "Color minimum",
            value="",
            placeholder="min (optional)",
            label_visibility="collapsed",
        )
        color_dash_col.markdown("<div style='text-align:center; padding-top:8px'>-</div>", unsafe_allow_html=True)
        cmd_color_max = color_max_col.text_input(
            "Color maximum",
            value="",
            placeholder="max (optional)",
            label_visibility="collapsed",
        )
        st.markdown("**Magnitude range (selected filter, AB)**")
        magnitude_min_col, magnitude_dash_col, magnitude_max_col = st.columns([1, 0.15, 1])
        cmd_magnitude_min = magnitude_min_col.text_input(
            "Magnitude minimum",
            value="",
            placeholder="min (optional)",
            label_visibility="collapsed",
        )
        magnitude_dash_col.markdown("<div style='text-align:center; padding-top:8px'>-</div>", unsafe_allow_html=True)
        cmd_magnitude_max = magnitude_max_col.text_input(
            "Magnitude maximum",
            value="",
            placeholder="max (optional)",
            label_visibility="collapsed",
        )
    if not contrast_less:
        region_cutoff = st.number_input(
            "Region cutoff (arcsec)", min_value=0.0, value=2.0, step=0.1
        )
        core_diameter = st.selectbox(
            "Core-aperture diameter (arcsec)",
            CORE_DIAMETERS,
            index=2,
            format_func=lambda value: f'{value:g}"',
        )
    system_distance = st.number_input(
        "System distance (pc)", min_value=0.1, value=50.0, step=1.0
    )
    planet_model = st.selectbox(
        "Planetary model",
        tuple(PLANET_MODELS),
        help=(
            "Select the evolutionary/atmospheric model overlaid on the CMD. "
            "The contamination probability calculation is unchanged."
        ),
    )
    separation_interest = st.number_input(
        "Separation of interest (arcsec)", min_value=0.0, value=1.0, step=0.1
    )
    with st.expander("Plot range"):
        plot_min_col, plot_dash_col, plot_max_col = st.columns([1, 0.15, 1])
        min_radius = plot_min_col.number_input(
            "Minimum separation (arcsec)", min_value=0.0, value=0.0, step=0.5
        )
        plot_dash_col.markdown(
            "<div style='text-align:center; padding-top:8px'>-</div>",
            unsafe_allow_html=True,
        )
        max_radius = plot_max_col.number_input(
            "Maximum separation (arcsec)", min_value=0.1, value=10.0, step=0.5
        )
    submitted = st.form_submit_button("Calculate", type="primary", width="stretch")


calculation_state = st.session_state.get("calculation_state")
if not submitted and calculation_state is None:
    st.info("Choose the observation inputs and select **Calculate**.")
    if contrast_less:
        st.markdown(
            "Contrast-less mode treats the full field as background limited and uses "
            "the selected filter's complete COSMOS-3D point-source population."
        )
    else:
        st.markdown(
            "Inside the region cutoff, the tool conservatively counts all nonstellar sources "
            "detectable in the selected core aperture. Outside it, the tool counts only "
            "nonstellar sources brighter than the contrast-curve detection limit and "
            "consistent with the empirical point-source compactness bounds."
        )
    st.stop()

if submitted:
    uses_carter_rescaling = False
    carter_target_absolute_mag = None
    carter_host_offset_mag = None
    carter_background_limit_mag = None
    carter_fiducial_background_limit_mag = None
    if not REFERENCE_PACK.exists():
        st.error(f"Missing reference pack: {REFERENCE_PACK}. Run build_reference_pack.py first.")
        st.stop()
    planet_model_path = PLANET_MODELS[planet_model]["path"]
    if not planet_model_path.exists():
        st.error(f"Missing {planet_model} CMD tracks: {planet_model_path}.")
        st.stop()

    if contrast_less:
        curve_source = None
        st.info(
            f"Contrast-less mode: using the full {filter_name.upper()} point-source surface density."
        )
    else:
        if use_carter_curve:
            curve_source = DEFAULT_CARTER_CURVE
            uses_carter_rescaling = True
        else:
            if uploaded_curve is None:
                st.error("Upload a contrast curve or select the fiducial Carter F444W curve.")
                st.stop()
            curve_source = uploaded_curve

    try:
        reference_pack = get_reference_pack()
        cmd_color_range = (
            _optional_range(cmd_color_min, cmd_color_max, "CMD color")
            if apply_cmd_constraint
            else None
        )
        cmd_magnitude_range = (
            _optional_range(cmd_magnitude_min, cmd_magnitude_max, "CMD magnitude")
            if apply_cmd_constraint
            else None
        )
        if contrast_less:
            result = calculate_contrastless_contamination(
                reference_pack=reference_pack,
                filter_name=filter_name,
                max_radius_arcsec=float(max_radius),
                min_radius_arcsec=float(min_radius),
                cmd_color_range=cmd_color_range,
                cmd_magnitude_range=cmd_magnitude_range,
            )
        else:
            separation, delta_mag = load_contrast_curve(curve_source, contrast_kind)
            if uses_carter_rescaling:
                (
                    delta_mag,
                    carter_target_absolute_mag,
                    carter_host_offset_mag,
                    carter_fiducial_background_limit_mag,
                    carter_background_limit_mag,
                ) = rescale_carter_f444w_curve(
                    delta_mag,
                    host_star_magnitude_ab=float(host_magnitude),
                    system_distance_pc=float(system_distance),
                    background_limit_mag=float(carter_background_limit),
                )
                st.info(
                    "Using the Carter et al. (2023) F444W curve, rescaled in flux space "
                    f"for $M_{{F444W}}={carter_target_absolute_mag:.2f}$ AB at "
                    f"{float(system_distance):g} pc (host offset "
                    f"{carter_host_offset_mag:+.2f} mag); background limit = "
                    f"{carter_background_limit_mag:.2f} AB."
                )
            result = calculate_contamination(
                reference_pack=reference_pack,
                contrast_separation_arcsec=separation,
                contrast_delta_mag=delta_mag,
                filter_name=filter_name,
                host_star_magnitude_ab=float(host_magnitude),
                region_cutoff_arcsec=float(region_cutoff),
                core_diameter_arcsec=float(core_diameter),
                max_radius_arcsec=float(max_radius),
                min_radius_arcsec=float(min_radius),
                cmd_color_range=cmd_color_range,
                cmd_magnitude_range=cmd_magnitude_range,
            )
        planet_tracks = load_planet_tracks(planet_model_path, float(system_distance))
    except Exception as exc:
        st.error(f"Calculation failed: {exc}")
        st.stop()

    calculation_state = {
        "result": result,
        "reference_pack": reference_pack,
        "planet_tracks": planet_tracks,
        "planet_model": planet_model,
        "system_distance": float(system_distance),
        "min_radius": float(min_radius),
        "separation_interest": float(separation_interest),
        "max_radius": float(max_radius),
        "filter_name": filter_name,
        "contrast_less": bool(contrast_less),
        "uses_carter_rescaling": uses_carter_rescaling,
        "carter_target_absolute_mag": carter_target_absolute_mag,
        "carter_host_offset_mag": carter_host_offset_mag,
        "carter_background_limit_mag": carter_background_limit_mag,
        "carter_fiducial_background_limit_mag": carter_fiducial_background_limit_mag,
    }
    st.session_state["calculation_state"] = calculation_state
else:
    # Mode/input changes rerun the app, but should not discard the last successful output.
    result = calculation_state["result"]
    reference_pack = calculation_state["reference_pack"]
    planet_tracks = calculation_state["planet_tracks"]
    planet_model = calculation_state.get("planet_model", "Sonora Bobcat")
    system_distance = calculation_state["system_distance"]
    min_radius = calculation_state["min_radius"]
    separation_interest = calculation_state["separation_interest"]
    max_radius = calculation_state["max_radius"]
    filter_name = calculation_state["filter_name"]
    contrast_less = calculation_state["contrast_less"]
    uses_carter_rescaling = calculation_state.get("uses_carter_rescaling", False)
    carter_target_absolute_mag = calculation_state.get("carter_target_absolute_mag")
    carter_host_offset_mag = calculation_state.get("carter_host_offset_mag")
    carter_background_limit_mag = calculation_state.get("carter_background_limit_mag")
    carter_fiducial_background_limit_mag = calculation_state.get(
        "carter_fiducial_background_limit_mag",
        carter_background_limit_mag,
    )

if calculation_state["contrast_less"]:
    curve_source = None
else:
    curve_source = DEFAULT_CARTER_CURVE

if not result.radius_arcsec[0] <= separation_interest <= result.radius_arcsec[-1]:
    st.error(
        f"Separation of interest must be between {result.radius_arcsec[0]:.3f} and "
        f"{result.radius_arcsec[-1]:.3f} arcsec for this calculation."
    )
    st.stop()

interest_probability = float(np.interp(separation_interest, result.radius_arcsec, result.probability))
interest_expected = float(np.interp(separation_interest, result.radius_arcsec, result.expected_number))
if (
    result.probability_color_magnitude is not None
    and result.expected_number_color_magnitude is not None
):
    interest_probability_color_magnitude = float(
        np.interp(
            separation_interest,
            result.radius_arcsec,
            result.probability_color_magnitude,
        )
    )
    interest_expected_color_magnitude = float(
        np.interp(
            separation_interest,
            result.radius_arcsec,
            result.expected_number_color_magnitude,
        )
    )

metric_1, metric_2 = st.columns(2)
metric_1.metric(
    f"P(at least 1) at {separation_interest:g} arcsec",
    f"{interest_probability:.2%}",
)
metric_2.metric(
    f"Expected contaminants within {separation_interest:g} arcsec",
    f"{interest_expected:.4f}",
)
metric_3, metric_4 = st.columns(2)
metric_3.metric(f"P(at least 1) within {max_radius:g} arcsec", f"{result.probability[-1]:.2%}")
metric_4.metric(
    f"Expected contaminants within {max_radius:g} arcsec",
    f"{result.expected_number[-1]:.4f}",
)
if (
    result.probability_color_magnitude is not None
    and result.expected_number_color_magnitude is not None
):
    constraint_metric_1, constraint_metric_2 = st.columns(2)
    constraint_metric_1.metric(
        f"P(at least 1) with color-magnitude constraints at {separation_interest:g} arcsec",
        f"{interest_probability_color_magnitude:.2%}",
    )
    constraint_metric_2.metric(
        f"Expected contaminants with color-magnitude constraints within {separation_interest:g} arcsec",
        f"{interest_expected_color_magnitude:.4f}",
    )

left, right = st.columns(2, gap="large")
with left:
    probability_plot = probability_figure(result)
    st.pyplot(probability_plot, width="stretch")
    if contrast_less:
        probability_caption = (
            f"Cumulative probability using the full {filter_name.upper()} COSMOS-3D "
            "point-source surface density. The entire field is background limited."
        )
    else:
        probability_caption = (
            "Cumulative probability of at least one extragalactic contaminant within "
            "each separation. The dashed line marks the transition between the "
            "contrast-limited and background-limited counting regimes."
        )
        if uses_carter_rescaling:
            probability_caption += (
                " The fiducial Carter F444W curve is rescaled for the entered host "
                "magnitude and distance while retaining its fixed background floor."
            )
    if result.probability_color_magnitude is not None:
        probability_caption += (
            " The red curve shows the probability with color-magnitude constraints."
        )
    st.caption(probability_caption)
    st.download_button(
        "Download probability table",
        result_csv_bytes(result),
        file_name=f"{'contrastless' if contrast_less else 'contamination'}_{filter_name}_probability.csv",
        mime="text/csv",
        width="stretch",
        on_click="ignore",
    )
    st.download_button(
        "Download probability plot (300 DPI)",
        figure_png_bytes(probability_plot, dpi=300),
        file_name=f"{'contrastless' if contrast_less else 'contamination'}_{filter_name}_probability.png",
        mime="image/png",
        width="stretch",
        on_click="ignore",
    )
    plt.close(probability_plot)
with right:
    cmd_plot = cmd_figure(reference_pack, planet_tracks, float(system_distance))
    st.pyplot(cmd_plot, width="stretch")
    st.caption(
        "Binned colors show the COSMOS-3D extragalactic population as objects per "
        "square arcminute per CMD bin; empty bins are white. Colors and magnitudes "
        "are in the Vega system. Solid, dashed, and dotted lines are "
        f"{planet_model} tracks at 10, 100, and 1000 Myr, respectively "
        f"({PLANET_MODELS[planet_model]['details']}); white markers and labels identify "
        "800 K, 1200 K, and 2000 K on the 10 Myr track when covered by the model. "
        "All planetary tracks are shifted to the entered system distance."
    )
    st.download_button(
        "Download CMD plot (300 DPI)",
        figure_png_bytes(cmd_plot, dpi=300),
        file_name=f"cosmos3d_cmd_{planet_model.lower().replace(' ', '_')}.png",
        mime="image/png",
        width="stretch",
        on_click="ignore",
    )
    plt.close(cmd_plot)

with st.expander("Method and assumptions"):
    if contrast_less:
        st.markdown(
            f"""
- Contrast-less mode ignores the contrast curve, host magnitude, region cutoff, and core aperture.
- The entire field is background limited and uses all COSMOS-3D sources classified as point-like in {filter_name.upper()}.
- The point-source count is divided by the {filter_name.upper()} survey area to obtain a constant sky density.
- When enabled, color-magnitude constraints use the Vega F200W-F444W color range and the AB magnitude range in {filter_name.upper()}.
- Background sources are assumed to follow a spatial Poisson distribution.
"""
        )
    else:
        if uses_carter_rescaling:
            carter_method = f"""
- The packaged Carter F444W curve is decomposed in flux space using its fiducial background floor ({carter_fiducial_background_limit_mag:.2f} AB).
- The host-dependent residual is scaled using the entered host's inferred absolute magnitude ({carter_target_absolute_mag:.2f} AB) and distance ({system_distance:g} pc), equivalent to a {carter_host_offset_mag:+.2f} mag apparent-host offset from HIP 65426.
- The recomposed curve uses the selected background limiting magnitude ({carter_background_limit_mag:.2f} AB).
"""
        else:
            carter_method = """
- Uploaded contrast curves are used as supplied; their host magnitude and distance are assumed to be already appropriate for that observation.
"""
        st.markdown(
            f"""
- Contrast limits are combined with the entered **AB** host star magnitude.
- For separations below {region_cutoff:g} arcsec, the limiting core magnitude is
  `host star magnitude + delta magnitude - 2.5 log10(encircled energy)`.
- At and beyond {region_cutoff:g} arcsec, the total-magnitude limit is `host star magnitude + delta magnitude`.
- Source densities come from the COSMOS-3D v1.2 merged catalog after removal of matched stars.
- When enabled, color-magnitude constraints use the Vega F200W-F444W color range and the AB magnitude range in the selected contrast-curve filter, combined with the contrast/aperture detection limit.
- The final contrast value is held constant beyond the uploaded curve's outermost point.
- Background sources are assumed to follow a spatial Poisson distribution.
{carter_method}
"""
        )
