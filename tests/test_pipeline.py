from __future__ import annotations

from io import BytesIO
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pipeline import (
    CARTER_F444W_FIDUCIAL_APPARENT_MAG_AB,
    CARTER_F444W_FIDUCIAL_DISTANCE_PC,
    calculate_contrastless_contamination,
    calculate_contamination,
    cmd_figure,
    figure_png_bytes,
    integrate_surface_density,
    load_contrast_curve,
    load_reference_pack,
    load_planet_tracks,
    poisson_probability,
    probability_figure,
    rescale_carter_f444w_curve,
)


ROOT = Path(__file__).resolve().parent.parent


def test_load_contrast_ratio_and_delta_mag() -> None:
    source = BytesIO(b"0.2,0.01\n1.0,0.001\n")
    radius, delta = load_contrast_curve(source, "auto")
    np.testing.assert_allclose(radius, [0.2, 1.0])
    np.testing.assert_allclose(delta, [5.0, 7.5])


def test_integrate_constant_density() -> None:
    radius = np.linspace(0.2, 10.0, 5000)
    density = np.full_like(radius, 2.5)
    expected = integrate_surface_density(radius, density)
    analytic = np.pi * 2.5 * (10.0**2 - 0.2**2) / 3600.0
    np.testing.assert_allclose(expected[-1], analytic, rtol=2e-6)
    np.testing.assert_allclose(poisson_probability([0, 1]), [0, 1 - np.exp(-1)])


def test_calculation_switches_catalog_population_at_cutoff() -> None:
    pack = {
        "core_f444w_0p2": np.array([20.0, 21.0, 22.0], dtype=np.float32),
        "point_f444w": np.array([20.0], dtype=np.float32),
        "area_f444w": np.array(100.0),
        "ee_f444w_0p2": np.array(0.5),
    }
    result = calculate_contamination(
        pack,
        contrast_separation_arcsec=np.array([0.2, 2.0]),
        contrast_delta_mag=np.array([12.0, 12.0]),
        filter_name="f444w",
        host_star_magnitude_ab=9.0,
        region_cutoff_arcsec=1.0,
        core_diameter_arcsec=0.2,
        max_radius_arcsec=2.0,
    )
    inside = result.radius_arcsec < 1.0
    outside = result.radius_arcsec >= 1.0
    assert np.all(result.surface_density_arcmin2[inside] == 0.02)
    assert np.all(result.surface_density_arcmin2[outside] == 0.01)
    assert np.all(np.diff(result.probability) >= 0)


def test_contrastless_applies_background_limit_before_cmd_constraints() -> None:
    pack = {
        "point_f444w": np.array([18.0, 20.0, 22.0], dtype=np.float32),
        "area_f444w": np.array(100.0),
        "constraint_cmd_color_vega": np.array([0.0, 1.0, 2.0], dtype=np.float32),
        "constraint_total_f444w_ab": np.array([18.0, 20.0, 22.0], dtype=np.float32),
        "constraint_point_index_f444w": np.array([0, 1, 2], dtype=np.int32),
    }
    result = calculate_contrastless_contamination(
        pack,
        filter_name="f444w",
        background_limiting_magnitude_ab=21.0,
        max_radius_arcsec=2.0,
        cmd_magnitude_range=(None, 19.0),
    )
    assert result.contrast_less
    np.testing.assert_allclose(result.surface_density_arcmin2, 0.02)
    analytic = np.pi * 0.02 * 2.0**2 / 3600.0
    np.testing.assert_allclose(result.expected_number[-1], analytic, rtol=2e-6)
    assert result.expected_number_color_magnitude is not None
    np.testing.assert_allclose(result.expected_number_color_magnitude[-1], analytic / 2.0, rtol=2e-6)


def test_color_magnitude_constraint_adds_a_reduced_probability_curve() -> None:
    pack = {
        "core_f444w_0p2": np.array([18.0, 20.0, 22.0], dtype=np.float32),
        "point_f444w": np.array([18.0, 20.0], dtype=np.float32),
        "area_f444w": np.array(100.0),
        "ee_f444w_0p2": np.array(0.5),
        "constraint_cmd_color_vega": np.array([0.0, 1.0, 2.0], dtype=np.float32),
        "constraint_total_f444w_ab": np.array([18.0, 20.0, 22.0], dtype=np.float32),
        "constraint_core_index_f444w_0p2": np.array([0, 1, 2], dtype=np.int32),
        "constraint_point_index_f444w": np.array([0, 1], dtype=np.int32),
    }
    result = calculate_contamination(
        pack,
        contrast_separation_arcsec=np.array([0.2, 2.0]),
        contrast_delta_mag=np.array([12.0, 12.0]),
        filter_name="f444w",
        host_star_magnitude_ab=9.0,
        region_cutoff_arcsec=1.0,
        core_diameter_arcsec=0.2,
        max_radius_arcsec=2.0,
        cmd_color_range=(0.5, 2.5),
        cmd_magnitude_range=(19.0, 23.0),
    )
    assert result.probability_color_magnitude is not None
    assert result.expected_number_color_magnitude is not None
    assert result.cmd_color_range == (0.5, 2.5)
    assert result.cmd_magnitude_range == (19.0, 23.0)
    assert np.all(result.probability_color_magnitude <= result.probability)
    assert result.probability_color_magnitude[-1] < result.probability[-1]

    one_sided = calculate_contamination(
        pack,
        contrast_separation_arcsec=np.array([0.2, 2.0]),
        contrast_delta_mag=np.array([12.0, 12.0]),
        filter_name="f444w",
        host_star_magnitude_ab=9.0,
        region_cutoff_arcsec=1.0,
        core_diameter_arcsec=0.2,
        max_radius_arcsec=2.0,
        cmd_color_range=None,
        cmd_magnitude_range=(19.0, None),
    )
    assert one_sided.probability_color_magnitude is not None

    upper_only = calculate_contamination(
        pack,
        contrast_separation_arcsec=np.array([0.2, 2.0]),
        contrast_delta_mag=np.array([12.0, 12.0]),
        filter_name="f444w",
        host_star_magnitude_ab=9.0,
        region_cutoff_arcsec=1.0,
        core_diameter_arcsec=0.2,
        max_radius_arcsec=2.0,
        cmd_color_range=(None, 1.5),
        cmd_magnitude_range=(None, 24.0),
    )
    assert upper_only.probability_color_magnitude is not None


def test_probability_figure_labels_regimes_and_optional_curve() -> None:
    pack = {
        "core_f444w_0p2": np.array([18.0, 20.0], dtype=np.float32),
        "point_f444w": np.array([18.0], dtype=np.float32),
        "area_f444w": np.array(100.0),
        "ee_f444w_0p2": np.array(0.5),
        "constraint_cmd_color_vega": np.array([0.0, 1.0], dtype=np.float32),
        "constraint_total_f444w_ab": np.array([18.0, 20.0], dtype=np.float32),
        "constraint_core_index_f444w_0p2": np.array([0, 1], dtype=np.int32),
        "constraint_point_index_f444w": np.array([0], dtype=np.int32),
    }
    result = calculate_contamination(
        pack,
        np.array([0.2, 2.0]),
        np.array([12.0, 12.0]),
        "f444w",
        9.0,
        1.0,
        0.2,
        max_radius_arcsec=2.0,
        cmd_color_range=(0.5, 1.5),
        cmd_magnitude_range=(19.0, 21.0),
    )
    figure = probability_figure(result)
    labels = [text.get_text() for text in figure.axes[0].texts]
    assert "Contrast limited" in labels
    assert "Background limited" in labels
    assert len(figure.axes[0].lines) == 3  # two curves plus the cutoff marker
    plt.close(figure)


def test_plot_minimum_changes_view_without_truncating_calculation() -> None:
    pack = {
        "core_f444w_0p2": np.array([20.0, 21.0], dtype=np.float32),
        "point_f444w": np.array([20.0], dtype=np.float32),
        "area_f444w": np.array(100.0),
        "ee_f444w_0p2": np.array(0.5),
    }
    result = calculate_contamination(
        pack,
        np.array([0.2, 2.0]),
        np.array([12.0, 12.0]),
        "f444w",
        9.0,
        1.0,
        0.2,
        max_radius_arcsec=2.0,
        min_radius_arcsec=1.5,
    )
    assert result.plot_min_arcsec == 1.5
    assert result.radius_arcsec[0] == 0.2
    figure = probability_figure(result)
    np.testing.assert_allclose(figure.axes[0].get_xlim(), [1.5, 2.0])
    plt.close(figure)


def test_carter_rescaling_preserves_fiducial_and_background_floor() -> None:
    fiducial_delta = np.array([8.0, 10.0, 12.0], dtype=float)
    fiducial_mag = CARTER_F444W_FIDUCIAL_APPARENT_MAG_AB
    fiducial_distance = CARTER_F444W_FIDUCIAL_DISTANCE_PC
    (
        rescaled_delta,
        target_absolute,
        host_offset,
        fiducial_background_limit,
        background_limit,
    ) = (
        rescale_carter_f444w_curve(fiducial_delta, fiducial_mag, fiducial_distance)
    )
    np.testing.assert_allclose(rescaled_delta, fiducial_delta)
    np.testing.assert_allclose(
        target_absolute, fiducial_mag - 5.0 * np.log10(fiducial_distance / 10.0)
    )
    np.testing.assert_allclose(host_offset, 0.0, atol=1e-12)

    (
        shifted_delta,
        _,
        _,
        shifted_fiducial_background_limit,
        shifted_background_limit,
    ) = rescale_carter_f444w_curve(fiducial_delta, fiducial_mag + 2.0, fiducial_distance)
    shifted_limit = fiducial_mag + 2.0 + shifted_delta
    assert shifted_limit[0] > fiducial_mag + fiducial_delta[0]
    np.testing.assert_allclose(shifted_limit[-1], shifted_background_limit)
    np.testing.assert_allclose(shifted_background_limit, background_limit)
    np.testing.assert_allclose(fiducial_background_limit, background_limit)
    np.testing.assert_allclose(
        shifted_fiducial_background_limit,
        fiducial_background_limit,
    )

    custom_background = 25.0
    custom_delta, _, _, decomposition_floor, applied_floor = rescale_carter_f444w_curve(
        fiducial_delta,
        fiducial_mag,
        fiducial_distance,
        background_limit_mag=custom_background,
    )
    np.testing.assert_allclose(decomposition_floor, fiducial_background_limit)
    np.testing.assert_allclose(applied_floor, custom_background)
    np.testing.assert_allclose(np.max(fiducial_mag + custom_delta), custom_background)


def test_carter_hip65426_regression() -> None:
    pack = load_reference_pack(ROOT / "data" / "cosmos3d_reference.npz")
    radius, delta_mag = load_contrast_curve(ROOT / "defaults" / "carter_f444w.csv")
    result = calculate_contamination(
        pack,
        radius,
        delta_mag,
        filter_name="f444w",
        host_star_magnitude_ab=CARTER_F444W_FIDUCIAL_APPARENT_MAG_AB,
        region_cutoff_arcsec=2.0,
        core_diameter_arcsec=0.2,
        max_radius_arcsec=10.0,
    )
    assert 0.08 < result.probability[-1] < 0.09
    assert np.all(np.diff(result.probability) >= 0)


def test_all_planet_models_render_with_binned_colorbar() -> None:
    pack = load_reference_pack(ROOT / "data" / "cosmos3d_reference.npz")
    model_paths = (
        ROOT / "data" / "sonora_bobcat_JWST_CMD_arrays.npz",
        ROOT / "data" / "atmo_theoretical_JWST_CMD.csv",
        ROOT / "data" / "flameskimmer_JWST_CMD.csv",
    )
    for model_path in model_paths:
        tracks = load_planet_tracks(model_path, 50.0)
        assert sorted(tracks) == [10.0, 100.0, 1000.0]
        assert all(track.size for track in tracks.values())
        figure = cmd_figure(pack, tracks, 50.0)
        figure.canvas.draw()
        assert len(figure.axes) == 2
        assert figure.axes[0].get_title() == ""
        assert "Objects" in figure.axes[1].get_ylabel()
        assert figure_png_bytes(figure, dpi=300).startswith(b"\x89PNG")
        plt.close(figure)


def test_atmo_tracks_cover_all_displayed_temperature_markers() -> None:
    """Guard against accidentally truncating the older ATMO track by mass."""
    tracks = load_planet_tracks(ROOT / "data" / "atmo_theoretical_JWST_CMD.csv", 50.0)
    for age in (10.0, 100.0, 1000.0):
        shown = tracks[age][tracks[age][:, 2] <= 2400.0]
        assert shown[:, 2].min() <= 800.0
        assert shown[:, 2].max() >= 2000.0
