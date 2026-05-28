import numpy as np
import pandas as pd
import pytest

from gem.metabolism import (
    DEFAULT_BODY_TEMP_C,
    DEFAULT_TAXA_PARAMS,
    body_temp,
    calculate_biomass_metabolism,
    calculate_massspecific_metabolism,
    metabolism_parameters,
)


# --- body_temp ---------------------------------------------------------------
def test_body_temp_scalar_bird_and_mammal():
    # Software: Lookup of temperature constants for endotherm taxa.
    assert body_temp("endotherm-bird") == 41.5
    assert body_temp("endotherm-mammal") == 36.5


def test_body_temp_list_returns_array():
    # Software: Vectorization over a list of taxa returns a numpy array.
    out = body_temp(["endotherm-bird", "endotherm-mammal", "endotherm-bird"])
    assert isinstance(out, np.ndarray)
    np.testing.assert_array_equal(out, np.array([41.5, 36.5, 41.5]))


def test_body_temp_ectotherm_returns_nan():
    # Ecology: Ectotherms have no fixed body temperature; they match ambient temperature (represented as NaN).
    assert np.isnan(body_temp("ectotherm"))


def test_body_temp_mixed_list_has_nan_for_ectotherm():
    # Ecology: In a mixed list, endotherms have fixed temperatures while ectotherms are NaN.
    out = body_temp(["endotherm-bird", "ectotherm", "endotherm-mammal"])
    assert out[0] == 41.5
    assert np.isnan(out[1])
    assert out[2] == 36.5


def test_body_temp_unknown_taxa_raises():
    # Software: Raises error on unrecognized taxa name.
    with pytest.raises(ValueError, match="Unknown taxa"):
        body_temp("fish")


# --- metabolism_parameters ---------------------------------------------------
def test_metabolism_parameters_single_string_species():
    # Software: Retrieves metabolism coefficients (c_int, b) for a single species from default table.
    df = metabolism_parameters(species="robin", taxa="endotherm-bird")
    assert list(df.columns) == ["species", "taxa", "c_int", "b"]
    assert len(df) == 1
    assert df.loc[0, "species"] == "robin"
    assert df.loc[0, "c_int"] == 19.53
    assert df.loc[0, "b"] == 0.73


def test_metabolism_parameters_list_with_per_species_taxa():
    # Software: Retrieves per-species parameters matching each species to its taxa.
    df = metabolism_parameters(
        species=["robin", "vole", "frog"],
        taxa=["endotherm-bird", "endotherm-mammal", "ectotherm"],
    )
    assert len(df) == 3
    np.testing.assert_array_equal(df["c_int"].values, [19.53, 19.53, 17.4])
    np.testing.assert_array_equal(df["b"].values, [0.73, 0.73, 0.84])


def test_metabolism_parameters_single_taxa_broadcasts_to_all_species():
    # Software: Broadcasting: a single taxa applies to all species in the list.
    df = metabolism_parameters(
        species=["bird_a", "bird_b"], taxa="endotherm-bird",
    )
    assert (df["taxa"] == "endotherm-bird").all()


def test_metabolism_parameters_custom_params_overrides_taxa():
    # Software: Custom parameters override default table lookups.
    df = metabolism_parameters(
        species=["robin", "tenrec"],
        taxa=["endotherm-bird", None],
        custom_params={"tenrec": {"c_int": 18.0, "b": 0.80, "taxa": "tenrec"}},
        interactive=False,
    )
    assert df.loc[1, "taxa"] == "tenrec"
    assert df.loc[1, "c_int"] == 18.0
    assert df.loc[1, "b"] == 0.80


def test_metabolism_parameters_unknown_species_non_interactive_raises():
    # Software: Raises error when species is not found and interactive mode is off.
    with pytest.raises(ValueError, match="tenrec"):
        metabolism_parameters(
            species=["tenrec"], taxa=None, interactive=False,
        )


def test_metabolism_parameters_mismatched_lengths_raise():
    # Software: Raises error when species and taxa lists have different lengths.
    with pytest.raises(ValueError, match="length"):
        metabolism_parameters(
            species=["a", "b"], taxa=["endotherm-bird"], interactive=False,
        )


def test_metabolism_parameters_invalid_taxa_label_raises():
    # Software: Raises error on unrecognized taxa name.
    with pytest.raises(ValueError, match="Unknown taxa"):
        metabolism_parameters(
            species=["robin"], taxa=["fish"], interactive=False,
        )


def test_metabolism_parameters_interactive_prompt_uses_input(monkeypatch):
    # Software: In interactive mode, user input selects a taxa from the menu.
    monkeypatch.setattr("builtins.input", lambda *a, **kw: "2")
    df = metabolism_parameters(species=["tenrec"], taxa=None, interactive=True)
    assert df.loc[0, "taxa"] == "endotherm-bird"
    assert df.loc[0, "c_int"] == 19.53


def test_metabolism_parameters_interactive_custom_choice(monkeypatch):
    # Software: In interactive mode, user can enter custom c_int and b parameters.
    answers = iter(["4", "18.5", "0.81"])
    monkeypatch.setattr("builtins.input", lambda *a, **kw: next(answers))
    df = metabolism_parameters(species=["tenrec"], taxa=None, interactive=True)
    assert df.loc[0, "taxa"] == "custom"
    assert df.loc[0, "c_int"] == 18.5
    assert df.loc[0, "b"] == 0.81


# --- calculate_massspecific_metabolism ---------------------------------------
def test_calculate_massspecific_metabolism_matches_hand_formula():
    # Ecology: Validates the Boltzmann-Arrhenius metabolic scaling formula (mass, temperature, coefficients).
    mass = np.array([100.0])
    temp_C = np.array([41.5])
    c_int = np.array([19.53])
    b = np.array([0.73])

    E_a, k_B = 0.63, 8.617e-5
    T_K = 41.5 + 273.15
    expected = 3.0 * np.exp(19.53) * 100.0 ** (0.73 - 1.0) * np.exp(-E_a / (k_B * T_K))

    rate = calculate_massspecific_metabolism(mass, temp_C, c_int, b)
    np.testing.assert_allclose(rate, [expected], rtol=1e-12)


def test_calculate_massspecific_metabolism_vectorises_over_species():
    # Software: Vectorizes across multiple species with different masses, temperatures, and parameters.
    mass   = np.array([10.0, 100.0, 1000.0])
    temp_C = np.array([41.5, 36.5, 20.0])
    c_int  = np.array([19.53, 19.53, 17.4])
    b      = np.array([0.73, 0.73, 0.84])
    rate = calculate_massspecific_metabolism(mass, temp_C, c_int, b)
    assert rate.shape == (3,)
    assert np.all(rate > 0)


def test_calculate_massspecific_metabolism_runs_on_grid_shape():
    # Software: Runs on full spatial grid shape (X, Y, Species), required for ecosystem simulations.
    shape = (4, 5, 3)
    mass   = np.full(shape, 50.0)
    temp_C = np.full(shape, 36.5)
    c_int  = np.full(shape, 19.53)
    b      = np.full(shape, 0.73)
    rate = calculate_massspecific_metabolism(mass, temp_C, c_int, b)
    assert rate.shape == shape


def test_calculate_massspecific_metabolism_shape_mismatch_raises():
    # Software: Raises error when input arrays have incompatible shapes.
    with pytest.raises(AssertionError):
        calculate_massspecific_metabolism(
            mass_g=np.array([10.0, 20.0]),
            body_temp_C=np.array([36.5]),
            c_int=np.array([19.53, 19.53]),
            b=np.array([0.73, 0.73]),
        )


def test_calculate_massspecific_metabolism_nan_filled_from_ambient_scalar():
    # Ecology: Ectotherms (NaN body_temp) use ambient temperature; endotherms use their fixed temperature.
    mass   = np.array([20.0, 20.0])
    temp_C = np.array([41.5, np.nan])
    c_int  = np.array([19.53, 17.4])
    b      = np.array([0.73, 0.84])

    rate = calculate_massspecific_metabolism(
        mass, temp_C, c_int, b, ambient_temp_C=15.0,
    )

    expected = calculate_massspecific_metabolism(
        mass, np.array([41.5, 15.0]), c_int, b,
    )
    np.testing.assert_allclose(rate, expected, rtol=1e-12)


def test_calculate_massspecific_metabolism_nan_filled_from_ambient_array():
    # Ecology: Per-cell ambient temperature for ectotherms (allows spatially varying ectotherm metabolism).
    temp_C = np.array([41.5, np.nan, np.nan])
    ambient = np.array([0.0, 10.0, 25.0])
    mass   = np.array([20.0, 20.0, 20.0])
    c_int  = np.array([19.53, 17.4, 17.4])
    b      = np.array([0.73, 0.84, 0.84])

    rate = calculate_massspecific_metabolism(
        mass, temp_C, c_int, b, ambient_temp_C=ambient,
    )

    expected_temps = np.array([41.5, 10.0, 25.0])
    expected = calculate_massspecific_metabolism(mass, expected_temps, c_int, b)
    np.testing.assert_allclose(rate, expected, rtol=1e-12)


def test_calculate_massspecific_metabolism_nan_without_ambient_raises():
    # Software: Raises error if ectotherm (NaN temp) is encountered without ambient_temp_C provided.
    with pytest.raises(ValueError, match="ambient_temp_C"):
        calculate_massspecific_metabolism(
            mass_g=np.array([20.0]),
            body_temp_C=np.array([np.nan]),
            c_int=np.array([17.4]),
            b=np.array([0.84]),
        )


def test_calculate_massspecific_metabolism_ignores_ambient_without_nan():
    # Software: When all temperatures are defined (no NaN), ambient_temp_C is ignored.
    mass   = np.array([20.0])
    temp_C = np.array([41.5])
    c_int  = np.array([19.53])
    b      = np.array([0.73])
    with_ambient = calculate_massspecific_metabolism(
        mass, temp_C, c_int, b, ambient_temp_C=-999.0,
    )
    without = calculate_massspecific_metabolism(mass, temp_C, c_int, b)
    np.testing.assert_allclose(with_ambient, without, rtol=1e-12)


def test_calculate_massspecific_metabolism_fmr_multiplier_one_gives_resting():
    # Ecology: Default (fmr_multiplier=3) is 3× resting metabolic rate; fmr_multiplier=1 gives resting rate.
    mass = np.array([100.0])
    temp_C = np.array([36.5])
    c_int = np.array([19.53])
    b = np.array([0.73])
    field = calculate_massspecific_metabolism(mass, temp_C, c_int, b)
    resting = calculate_massspecific_metabolism(
        mass, temp_C, c_int, b, fmr_multiplier=1.0,
    )
    np.testing.assert_allclose(field, 3.0 * resting, rtol=1e-12)


# --- calculate_biomass_metabolism --------------------------------------------
def test_calculate_biomass_metabolism_matches_hand_formula():
    # Ecology: Converts mass-specific metabolic rate to biomass loss per time step.
    B0 = np.array([1000.0])
    fmr = calculate_massspecific_metabolism(
        mass_g=np.array([100.0]),
        body_temp_C=np.array([41.5]),
        c_int=np.array([19.53]),
        b=np.array([0.73]),
    )
    delta = calculate_biomass_metabolism(B0, fmr, dt=1.0)
    expected = -(1000.0 * fmr[0] * 86400.0) / 7000.0
    np.testing.assert_allclose(delta, [expected], rtol=1e-12)


def test_calculate_biomass_metabolism_is_negative():
    # Ecology: Metabolism causes biomass loss (negative delta).
    B0 = np.array([500.0, 50.0])
    fmr = np.array([0.02, 0.03])
    delta = calculate_biomass_metabolism(B0, fmr, dt=1.0)
    assert np.all(delta < 0)


def test_calculate_biomass_metabolism_scales_linearly_with_dt():
    # Ecology: Metabolic loss scales proportionally with time step (dt).
    B0 = np.array([100.0])
    fmr = np.array([0.02])
    one_day = calculate_biomass_metabolism(B0, fmr, dt=1.0)
    two_days = calculate_biomass_metabolism(B0, fmr, dt=2.0)
    np.testing.assert_allclose(two_days, 2.0 * one_day, rtol=1e-12)


def test_calculate_biomass_metabolism_runs_on_grid_shape():
    # Software: Runs on full spatial grid shape (X, Y, Species).
    shape = (4, 5, 3)
    B0 = np.full(shape, 100.0)
    fmr = np.full(shape, 0.02)
    delta = calculate_biomass_metabolism(B0, fmr, dt=1.0)
    assert delta.shape == shape
    assert np.all(delta < 0)


def test_calculate_biomass_metabolism_shape_mismatch_raises():
    # Software: Raises error when input arrays have incompatible shapes.
    with pytest.raises(AssertionError):
        calculate_biomass_metabolism(
            initial_biomass_g=np.array([100.0, 200.0]),
            mass_specific_metabolic_rate=np.array([0.02]),
            dt=1.0,
        )


def test_calculate_biomass_metabolism_zero_biomass_zero_delta():
    # Ecology: Zero biomass incurs zero metabolic loss.
    delta = calculate_biomass_metabolism(
        initial_biomass_g=np.array([0.0, 0.0]),
        mass_specific_metabolic_rate=np.array([0.02, 0.05]),
        dt=1.0,
    )
    np.testing.assert_array_equal(delta, [0.0, 0.0])


def test_calculate_biomass_metabolism_energy_density_override():
    # Ecology: Energy density of biomass affects metabolic loss (higher energy density → lower loss for same metabolic rate).
    B0 = np.array([100.0])
    fmr = np.array([0.02])
    default = calculate_biomass_metabolism(B0, fmr, dt=1.0)
    doubled = calculate_biomass_metabolism(
        B0, fmr, dt=1.0, energy_density_J_per_g=14000.0,
    )
    np.testing.assert_allclose(doubled, default / 2.0, rtol=1e-12)


# --- end-to-end glue ---------------------------------------------------------
def test_end_to_end_endotherm_pipeline():
    # Integration: Full workflow for endotherms: parameters → body temps → mass-specific rates.
    df = metabolism_parameters(
        species=["robin", "vole"],
        taxa=["endotherm-bird", "endotherm-mammal"],
    )
    temps_C = body_temp(df["taxa"].tolist())
    mass_g = np.array([20.0, 30.0])
    rate = calculate_massspecific_metabolism(
        mass_g=mass_g,
        body_temp_C=temps_C,
        c_int=df["c_int"].values,
        b=df["b"].values,
    )
    assert rate.shape == (2,)
    assert np.all(rate > 0)


def test_end_to_end_mixed_pipeline_with_ectotherm():
    # Integration + Ecology: Full workflow for mixed endotherm/ectotherm community; ectotherms use ambient temperature.
    df = metabolism_parameters(
        species=["robin", "frog", "vole"],
        taxa=["endotherm-bird", "ectotherm", "endotherm-mammal"],
    )
    temps_C = body_temp(df["taxa"].tolist())
    # Endotherm slots are populated, ectotherm slot is NaN.
    assert temps_C[0] == 41.5
    assert np.isnan(temps_C[1])
    assert temps_C[2] == 36.5

    rate = calculate_massspecific_metabolism(
        mass_g=np.array([20.0, 5.0, 30.0]),
        body_temp_C=temps_C,
        c_int=df["c_int"].values,
        b=df["b"].values,
        ambient_temp_C=15.0,
    )
    assert rate.shape == (3,)
    assert np.all(np.isfinite(rate))
    assert np.all(rate > 0)

    initial_biomass_g = np.array([2000.0, 500.0, 1500.0])
    delta = calculate_biomass_metabolism(initial_biomass_g, rate, dt=1.0)
    assert delta.shape == (3,)
    assert np.all(delta < 0)


# --- module-level constants are wired up correctly ---------------------------
def test_default_tables_have_expected_keys():
    # Software: Module initialization—check that default parameter and temperature tables are populated.
    assert set(DEFAULT_TAXA_PARAMS) == {"ectotherm", "endotherm-bird", "endotherm-mammal"}
    assert set(DEFAULT_BODY_TEMP_C) == {"endotherm-bird", "endotherm-mammal"}
