#!/usr/bin/env python3
"""Replay calibration from 2013 and run the 13 SMC/vaccine scenarios to 2030.

Edit the calibrated parameter set in inputs/smc_vaccine_scenarios.csv and
the schedule, coverage, and efficacy assumptions in
inputs/vaccine_dose_coverage_dummy.csv before submitting. The latter is a
generic SimpleVaccine proxy, not an R21 model.
"""

import os
import pathlib
from functools import partial

import numpy as np
import pandas as pd

from emod_api import schema_to_class as s2c
from emod_api.demographics.DemographicsTemplates import CrudeRate
from emod_api.interventions import common, utils as intervention_utils
import emod_api.campaign as campaign_api
import emodpy_malaria.demographics.MalariaDemographics as Demographics
import emodpy_malaria.interventions.bednet as itn
import emodpy_malaria.interventions.drug_campaign as dc
import emodpy_malaria.interventions.treatment_seeking as cm
from emodpy.emod_task import EMODTask
from emodpy_malaria.interventions.common import add_triggered_campaign_delay_event
from emodpy_malaria.interventions.scale_larval_habitats import (
    add_scale_larval_habitats,
)
from emodpy_malaria.reporters.builtin import add_malaria_summary_report
from idmtools.builders import SimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment

import manifest
from utils_slurm import build_burnin_df
from project_future_climate import WEATHER_DIR, build_projection


ROOT = pathlib.Path(__file__).resolve().parent
SCENARIO_FILE = ROOT / "inputs" / "smc_vaccine_scenarios.csv"
VACCINE_COVERAGE_FILE = ROOT / "inputs" / "vaccine_dose_coverage_dummy.csv"
SMC_FILE = ROOT / "inputs" / "SMC_intervention_northern_region.csv"
ITN_FILE = ROOT / "inputs" / "itn_irs_interventions.csv"
CASE_MANAGEMENT_FILE = ROOT / "inputs" / "case_management.csv"

# This must be the calibrated burn-in experiment that produced the 2013 states.
BURNIN_EXP_ID = "31f646d9-8002-4593-9530-2ae743835b6b"
SERIALIZE_YEARS = 50
START_YEAR = 2013
FORECAST_START_YEAR = 2026
FORECAST_END_YEAR = 2027
SIMULATION_END_YEAR = FORECAST_END_YEAR + 1
START_TIME = SERIALIZE_YEARS * 365 + 1
SIMULATION_DURATION = (SIMULATION_END_YEAR - START_YEAR) * 365
NODE_ID = 1
VECTOR_SPECIES = ["gambiae", "arabiensis", "funestus"]

VACCINE_RECEIPT_EVENTS = [
    "Received_Vaccine_Dose_1",
    "Received_Vaccine_Dose_2",
    "Received_Vaccine_Dose_3",
    "Received_Vaccine_Dose_4",
]


def set_param_fn(config):
    import emodpy_malaria.malaria_config as malaria_config

    config = malaria_config.set_team_defaults(config, manifest)
    malaria_config.add_species(config, manifest, VECTOR_SPECIES)
    config.parameters.Climate_Model = "CLIMATE_BY_DATA"
    config.parameters.Enable_Default_Reporting = 1
    config.parameters.Air_Temperature_Filename = os.path.join(
        "climate", "air_temperature.bin"
    )
    config.parameters.Land_Temperature_Filename = os.path.join(
        "climate", "land_temperature.bin"
    )
    config.parameters.Rainfall_Filename = os.path.join("climate", "rainfall.bin")
    config.parameters.Relative_Humidity_Filename = os.path.join(
        "climate", "relative_humidity.bin"
    )
    config.parameters.Simulation_Duration = SIMULATION_DURATION
    config.parameters.Run_Number = 0
    config.parameters.Start_Time = START_TIME
    config.parameters.Serialized_Population_Reading_Type = "READ"
    config.parameters.Serialization_Mask_Node_Read = 0
    config.parameters.Custom_Individual_Events = [
        "Received_Treatment",
        "Received_ITN",
        "Received_SMC",
        "Received_Campaign_Drugs",
        *VACCINE_RECEIPT_EVENTS,
    ]
    return config


def build_camp():
    campaign_api.set_schema(manifest.schema_file)
    campaign_api.schema_path = manifest.schema_file
    return campaign_api


def build_demog():
    demog = Demographics.from_template_node(
        lat=9.462781, lon=-0.285867, pop=1000, name="tamale"
    )
    national_cbr_1963 = 47.71
    northern_region_fr_1993 = 7.39
    tfr_1993 = 5.5
    national_cdr_1963 = 19.23
    northern_region_cbr_1963 = (
        national_cbr_1963 * (northern_region_fr_1993 / tfr_1993)
    )
    demog.SetBirthRate(CrudeRate(northern_region_cbr_1963))
    demog.SetMortalityRate(CrudeRate(national_cdr_1963))
    demog.SetEquilibriumAgeDistFromBirthAndMortRates(
        CrudeRate(northern_region_cbr_1963), CrudeRate(national_cdr_1963)
    )
    return demog


def simulation_day(year, day_of_year):
    return START_TIME + (int(year) - START_YEAR) * 365 + int(day_of_year)


def load_scenarios():
    scenarios = pd.read_csv(SCENARIO_FILE)
    required = {
        "scenario_id",
        "stage",
        "smc_future",
        "vaccine_future",
        "smc_coverage_mode",
        "smc_coverage",
        "vaccine_completion_coverage",
        "x_Temporary_Larval_Habitat",
        "scale_constant",
        "scale_temp_rain",
        "scale_water_veg",
    }
    missing = required.difference(scenarios.columns)
    if missing:
        raise ValueError(f"{SCENARIO_FILE} is missing: {sorted(missing)}")
    if scenarios["scenario_id"].duplicated().any():
        raise ValueError("Scenario IDs in the scenario CSV must be unique.")
    if set(scenarios["scenario_id"]) != {
        "S0", "S1", "S2", "S3",
        "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9",
    }:
        raise ValueError("The scenario CSV must define all S0-S3 and C1-C9 rows.")
    boolean_values = {"true": True, "1": True, "yes": True,
                      "false": False, "0": False, "no": False}
    for column in ("smc_future", "vaccine_future"):
        normalized = scenarios[column].astype(str).str.strip().str.lower()
        invalid = ~normalized.isin(boolean_values)
        if invalid.any():
            raise ValueError(
                f"{SCENARIO_FILE} has invalid boolean values in {column}: "
                f"{sorted(normalized.loc[invalid].unique())}"
            )
        scenarios[column] = normalized.map(boolean_values)
    numeric = [
        "smc_coverage",
        "vaccine_completion_coverage",
        "x_Temporary_Larval_Habitat",
        "scale_constant",
        "scale_temp_rain",
        "scale_water_veg",
    ]
    scenarios[numeric] = scenarios[numeric].apply(
        pd.to_numeric, errors="raise"
    )
    if not scenarios["vaccine_completion_coverage"].between(0, 1).all():
        raise ValueError("Vaccine completion coverage must be between 0 and 1.")
    if not scenarios["smc_coverage"].between(0, 1).all():
        raise ValueError("SMC coverage must be between 0 and 1.")
    if not scenarios["vaccine_completion_coverage"].between(0, 1).all():
        raise ValueError("Vaccine completion coverage must be between 0 and 1.")
    if not scenarios["x_Temporary_Larval_Habitat"].gt(0).all():
        raise ValueError("Larval habitat and scaling parameters must be positive.")
    for column in ("scale_constant", "scale_temp_rain", "scale_water_veg"):
        if not scenarios[column].gt(0).all():
            raise ValueError(f"{column} must be positive.")
    if not scenarios["smc_coverage_mode"].isin(
        {"scenario", "projected_2025"}
    ).all():
        raise ValueError(
            "SMC coverage mode must be 'scenario' or 'projected_2025'."
        )
    parameter_columns = [
        "x_Temporary_Larval_Habitat",
        "scale_constant",
        "scale_temp_rain",
        "scale_water_veg",
    ]
    if len(scenarios[parameter_columns].drop_duplicates()) != 1:
        raise ValueError(
            "All scenarios must use the same calibrated parameter set. "
            "Edit the four parameter columns consistently in the scenario CSV."
        )
    return scenarios


def load_vaccine_dose_coverages(completion_coverage):
    dose_table = pd.read_csv(VACCINE_COVERAGE_FILE)
    required = {
        "completion_coverage",
        "dose",
        "dose_age_months",
        "conditional_coverage",
        "product",
        "implementation",
        "vaccine_type",
        "vaccine_take",
        "vaccine_initial_effect",
        "vaccine_box_duration",
        "vaccine_decay_time_constant",
        "efficacy_is_multiplicative",
    }
    missing = required.difference(dose_table.columns)
    if missing:
        raise ValueError(
            f"{VACCINE_COVERAGE_FILE} is missing: {sorted(missing)}"
        )
    dose_table["completion_coverage"] = pd.to_numeric(
        dose_table["completion_coverage"], errors="raise"
    )
    dose_table["dose"] = pd.to_numeric(dose_table["dose"], errors="raise").astype(int)
    dose_table["dose_age_months"] = pd.to_numeric(
        dose_table["dose_age_months"], errors="raise"
    ).astype(int)
    dose_table["conditional_coverage"] = pd.to_numeric(
        dose_table["conditional_coverage"], errors="raise"
    )
    model_columns = [
        "product",
        "implementation",
        "vaccine_type",
        "vaccine_take",
        "vaccine_initial_effect",
        "vaccine_box_duration",
        "vaccine_decay_time_constant",
        "efficacy_is_multiplicative",
    ]
    for column in model_columns[3:]:
        dose_table[column] = pd.to_numeric(dose_table[column], errors="raise")
    selected = dose_table.loc[
        np.isclose(dose_table["completion_coverage"], completion_coverage)
    ].sort_values("dose")
    if selected["dose"].tolist() != [1, 2, 3, 4]:
        raise ValueError(
            "Vaccine dose CSV must provide exactly one row for each of doses 1-4 "
            f"at completion coverage {completion_coverage:.2f}."
        )
    if not selected["conditional_coverage"].between(0, 1).all():
        raise ValueError("Conditional vaccine dose coverage must be between 0 and 1.")
    if selected["dose_age_months"].tolist() != [6, 7, 9, 18]:
        raise ValueError("The primary routine schedule must be at 6, 7, 9 and 18 months.")
    if len(selected[model_columns].drop_duplicates()) != 1:
        raise ValueError(
            "All four doses must use one vaccine product and a consistent "
            "EMOD efficacy/waning parameterization."
        )
    if not selected["vaccine_take"].between(0, 1).all() or not selected[
        "vaccine_initial_effect"
    ].between(0, 1).all():
        raise ValueError("Vaccine take and initial effect must be between 0 and 1.")
    if not selected["efficacy_is_multiplicative"].isin({0, 1}).all():
        raise ValueError("efficacy_is_multiplicative must be 0 or 1.")
    implied_completion = selected["conditional_coverage"].prod()
    if not np.isclose(implied_completion, completion_coverage, atol=1e-5):
        raise ValueError(
            f"Dose conditional coverages imply {implied_completion:.4f} "
            f"completion, not {completion_coverage:.4f}. Update "
            f"{VACCINE_COVERAGE_FILE}."
        )
    return selected.to_dict(orient="records")


def load_smc_schedule():
    smc = pd.read_csv(SMC_FILE)
    required = {"region", "year", "start_day", "coverage"}
    missing = required.difference(smc.columns)
    if missing:
        raise ValueError(f"{SMC_FILE} is missing: {sorted(missing)}")
    smc = smc.loc[smc["region"].astype(str).str.lower().eq("northern")].copy()
    smc["year"] = pd.to_numeric(smc["year"], errors="raise").astype(int)
    smc["start_day"] = pd.to_numeric(smc["start_day"], errors="raise").astype(int)
    smc["coverage"] = pd.to_numeric(smc["coverage"], errors="raise")
    if not smc["coverage"].between(0, 1).all():
        raise ValueError("SMC schedule coverage must be between 0 and 1.")
    if smc.empty or int(smc["year"].max()) != 2025:
        raise ValueError(
            "Northern-region SMC schedule must include observed 2025 rounds "
            "to project their timing into 2026-2030."
        )
    latest_year = int(smc["year"].max())
    latest_rounds = smc.loc[smc["year"].eq(latest_year)].copy()
    if latest_rounds.empty:
        raise ValueError("No SMC rounds found for the latest schedule year.")
    return smc, latest_year, latest_rounds


def projected_smc_rounds(smc, latest_year, latest_rounds):
    rows = []
    for row in smc.itertuples(index=False):
        day_of_year = int(row.start_day) % 365
        rows.append(
            {
                "year": int(row.year),
                "day_of_year": day_of_year,
                "coverage": float(row.coverage),
                "projected": False,
            }
        )
    for year in range(max(latest_year + 1, FORECAST_START_YEAR), FORECAST_END_YEAR + 1):
        for row in latest_rounds.itertuples(index=False):
            rows.append(
                {
                    "year": year,
                    "day_of_year": int(row.start_day) % 365,
                    "coverage": float(row.coverage),
                    "projected": True,
                }
            )
    return rows


def extend_latest_schedule_through_forecast(
    schedule, forecast_start_year, forecast_end_year, schedule_name
):
    if schedule.empty:
        raise ValueError(f"No pickup entries found for {schedule_name}.")
    if "year" not in schedule.columns:
        raise ValueError(f"{schedule_name} schedule must contain a 'year' column.")

    schedule = schedule.copy()
    schedule["year"] = pd.to_numeric(schedule["year"], errors="raise").astype(int)
    latest_year = int(schedule["year"].max())
    latest_schedule = schedule.loc[schedule["year"].eq(latest_year)].copy()

    forecast_schedules = []
    for year in range(forecast_start_year, forecast_end_year + 1):
        annual_schedule = latest_schedule.copy()
        annual_schedule["year"] = year
        forecast_schedules.append(annual_schedule)

    historical_schedule = schedule.loc[
        schedule["year"].lt(forecast_start_year)
    ].copy()
    extended = pd.concat(
        [historical_schedule, *forecast_schedules],
        ignore_index=True,
    )
    return extended, latest_year


def add_smc_campaign(campaign, scenario, smc_rounds):
    for row in smc_rounds:
        if row["projected"]:
            if not scenario["smc_future"]:
                continue
            coverage = (
                float(scenario["smc_coverage"])
                if scenario["smc_coverage_mode"] == "scenario"
                else row["coverage"]
            )
        else:
            coverage = row["coverage"]
        dc.add_drug_campaign(
            campaign,
            campaign_type="MDA",
            drug_code="SPA",
            start_days=[simulation_day(row["year"], row["day_of_year"])],
            repetitions=1,
            coverage=coverage,
            target_group={"agemin": 0.25, "agemax": 5},
            node_ids=[NODE_ID],
            receiving_drugs_event_name="Received_SMC",
        )


def simple_vaccine_with_receipt_event(campaign, event_name, dose_parameters):
    intervention = s2c.get_class_with_defaults(
        "SimpleVaccine", campaign.schema_path
    )
    intervention.Vaccine_Type = dose_parameters["vaccine_type"]
    intervention.Vaccine_Take = float(dose_parameters["vaccine_take"])
    intervention.Efficacy_Is_Multiplicative = int(
        dose_parameters["efficacy_is_multiplicative"]
    )
    intervention.Waning_Config = intervention_utils.get_waning_from_params(
        schema_path=campaign.schema_path,
        initial=float(dose_parameters["vaccine_initial_effect"]),
        box_duration=int(dose_parameters["vaccine_box_duration"]),
        decay_time_constant=float(
            dose_parameters["vaccine_decay_time_constant"]
        ),
    )
    return [
        intervention,
        common.BroadcastEvent(campaign, Event_Trigger=event_name),
    ]


def add_routine_vaccine(campaign, scenario, vaccine_start_day):
    if not scenario["vaccine_future"]:
        return
    doses = load_vaccine_dose_coverages(
        float(scenario["vaccine_completion_coverage"])
    )
    trigger_events = ["Births", *VACCINE_RECEIPT_EVENTS[:3]]
    age_days = [
        round(int(row["dose_age_months"]) * 365 / 12)
        for row in doses
    ]
    for index, row in enumerate(doses):
        delay = (
            age_days[index]
            if index == 0
            else age_days[index] - age_days[index - 1]
        )
        nested_intervention = simple_vaccine_with_receipt_event(
            campaign,
            VACCINE_RECEIPT_EVENTS[index],
            row,
        )
        add_triggered_campaign_delay_event(
            campaign=campaign,
            start_day=vaccine_start_day,
            trigger_condition_list=[trigger_events[index]],
            listening_duration=-1,
            delay_period_constant=delay,
            demographic_coverage=float(row["conditional_coverage"]),
            node_ids=[NODE_ID],
            target_age_min=0,
            target_age_max=5,
            individual_intervention=nested_intervention,
        )


def build_scenario_campaign(
    scenario,
    burnin_row,
    smc_rounds,
    pickup_itns,
    case_management,
):
    campaign = build_camp()
    habitat_scales = {
        "NodeID": [NODE_ID],
        "WATER_VEGETATION.gambiae": [float(burnin_row["scale_water_veg"])],
        "CONSTANT.arabiensis": [float(burnin_row["scale_constant"])],
        "TEMPORARY_RAINFALL.arabiensis": [
            float(burnin_row["scale_temp_rain"])
        ],
        "CONSTANT.funestus": [float(burnin_row["scale_constant"])],
        "TEMPORARY_RAINFALL.funestus": [
            float(burnin_row["scale_temp_rain"])
        ],
    }
    add_scale_larval_habitats(
        campaign,
        df=pd.DataFrame(habitat_scales),
        start_day=START_TIME + 1,
        repetitions=-1,
    )

    for row in pickup_itns.drop_duplicates(
        subset=["year", "start_day_of_year", "coverage"]
    ).itertuples(index=False):
        itn.add_itn_scheduled(
            campaign,
            start_day=simulation_day(
                int(row.year), int(row.start_day_of_year) % 365
            ),
            demographic_coverage=float(row.coverage),
            repetitions=1,
            node_ids=[NODE_ID],
            receiving_itn_broadcast_event="Received_ITN",
        )

    for row in case_management.itertuples(index=False):
        cm_start = simulation_day(
            int(row.year), int(row.start_day_of_year) % 365
        )
        drug = str(row.drug).replace("–", "-").split("-")
        drug = [
            "Chloroquine" if name.strip() == "Quinine" else name.strip()
            for name in drug
        ]
        coverage = float(row.u5_effective_coverage)
        cm.add_treatment_seeking(
            campaign,
            start_day=cm_start,
            drug=drug,
            targets=[
                {
                    "trigger": "NewClinicalCase",
                    "coverage": coverage,
                    "agemin": 0,
                    "agemax": 5,
                    "rate": 0.3,
                },
                {
                    "trigger": "NewClinicalCase",
                    "coverage": coverage * 0.75,
                    "agemin": 5,
                    "agemax": 115,
                    "rate": 0.3,
                },
                {
                    "trigger": "NewSevereCase",
                    "coverage": min(coverage * 1.15, 1),
                    "agemin": 0,
                    "agemax": 115,
                    "rate": 0.5,
                },
            ],
        )

    add_smc_campaign(campaign, scenario, smc_rounds)
    add_routine_vaccine(
        campaign,
        scenario,
        simulation_day(FORECAST_START_YEAR, 0),
    )
    return campaign


def configure_scenario_simulation(
    simulation,
    combination_index,
    burnin_rows,
    scenarios,
    smc_rounds,
    pickup_itns,
    case_management,
):
    scenario_count = len(scenarios)
    burnin_index, scenario_index = divmod(int(combination_index), scenario_count)
    burnin_row = burnin_rows.iloc[burnin_index]
    scenario = scenarios.iloc[scenario_index]
    simulation.task.set_parameter(
        "Serialized_Population_Filenames",
        burnin_row["Serialized_Population_Filenames"],
    )
    simulation.task.set_parameter(
        "Serialized_Population_Path",
        os.path.join(burnin_row["serialized_file_path"], "output"),
    )
    simulation.task.set_parameter("Run_Number", int(burnin_row["Run_Number"]))
    simulation.task.set_parameter(
        "x_Temporary_Larval_Habitat",
        float(burnin_row["x_Temporary_Larval_Habitat"]),
    )
    simulation.task.create_campaign_from_callback(
        partial(
            build_scenario_campaign,
            scenario=scenario.to_dict(),
            burnin_row=burnin_row.to_dict(),
            smc_rounds=smc_rounds,
            pickup_itns=pickup_itns,
            case_management=case_management,
        )
    )
    return {
        "Run_Number": int(burnin_row["Run_Number"]),
        "Scenario": str(scenario["scenario_id"]),
        "Stage": str(scenario["stage"]),
        "SMC_Future": bool(scenario["smc_future"]),
        "Vaccine_Future": bool(scenario["vaccine_future"]),
        "SMC_Coverage": float(scenario["smc_coverage"]),
        "Vaccine_Completion_Coverage": float(
            scenario["vaccine_completion_coverage"]
        ),
        "x_Temporary_Larval_Habitat": float(
            burnin_row["x_Temporary_Larval_Habitat"]
        ),
        "scale_constant": float(burnin_row["scale_constant"]),
        "scale_temp_rain": float(burnin_row["scale_temp_rain"]),
        "scale_water_veg": float(burnin_row["scale_water_veg"]),
    }


def general_sim():
    if not BURNIN_EXP_ID:
        raise ValueError("Set BURNIN_EXP_ID to a completed calibrated burn-in experiment.")
    build_projection()
    scenarios = load_scenarios()
    smc_data, latest_smc_year, latest_smc_rounds = load_smc_schedule()
    smc_rounds = projected_smc_rounds(
        smc_data, latest_smc_year, latest_smc_rounds
    )
    itn_data = pd.read_csv(ITN_FILE)
    pickup_itns = itn_data.loc[
        (itn_data["phase"] == "pickup") & (itn_data["intervention"] == "itn")
    ].copy()
    pickup_itns, latest_itn_year = extend_latest_schedule_through_forecast(
        pickup_itns,
        FORECAST_START_YEAR,
        FORECAST_END_YEAR,
        "pickup ITN",
    )
    case_management = pd.read_csv(CASE_MANAGEMENT_FILE)
    case_management = case_management.loc[
        case_management["phase"].eq("pickup")
    ].copy()
    case_management, latest_cm_year = extend_latest_schedule_through_forecast(
        case_management,
        FORECAST_START_YEAR,
        FORECAST_END_YEAR,
        "pickup case-management",
    )
    print(
        f"Repeating the latest pickup ITN schedule ({latest_itn_year}) and "
        f"case-management schedule ({latest_cm_year}) in each forecast year "
        f"({FORECAST_START_YEAR}-{FORECAST_END_YEAR}) for every scenario."
    )

    selected_parameters = scenarios.iloc[0][
        [
            "x_Temporary_Larval_Habitat",
            "scale_constant",
            "scale_temp_rain",
            "scale_water_veg",
        ]
    ]
    platform = Platform(
        manifest.default_platform,
        job_directory=manifest.job_directory,
        docker_image=manifest.plat_image,
        max_running_jobs=manifest.max_running_jobs,
        max_workers=12,
        max_job=12,
    )
    burnin_df = build_burnin_df(
        BURNIN_EXP_ID, platform, SERIALIZE_YEARS * 365
    )
    numeric_parameters = list(selected_parameters.index)
    for column in numeric_parameters:
        burnin_df[column] = pd.to_numeric(burnin_df[column], errors="raise")
    selected_mask = np.ones(len(burnin_df), dtype=bool)
    for column, value in selected_parameters.items():
        selected_mask &= np.isclose(burnin_df[column], float(value))
    burnin_rows = burnin_df.loc[selected_mask].reset_index(drop=True)
    if burnin_rows.empty:
        available = burnin_df[numeric_parameters].drop_duplicates().to_dict(
            orient="records"
        )
        raise ValueError(
            "No burn-in snapshots match the calibrated parameters in "
            f"{SCENARIO_FILE}. Available parameter sets: {available}"
        )
    if burnin_rows["Run_Number"].duplicated().any():
        raise ValueError(
            "Multiple burn-in snapshots share a Run_Number for the selected "
            "parameter set; inspect the burn-in experiment tags."
        )

    task = EMODTask.from_default2(
        config_path="config.json",
        eradication_path=manifest.eradication_path,
        campaign_builder=build_camp,
        schema_path=manifest.schema_file,
        param_custom_cb=set_param_fn,
        ep4_custom_cb=None,
        demog_builder=build_demog,
        plugin_report=None,
    )
    task.common_assets.add_directory(str(WEATHER_DIR), relative_path="climate")

    builder = SimulationBuilder()
    builder.add_sweep_definition(
        partial(
            configure_scenario_simulation,
            burnin_rows=burnin_rows,
            scenarios=scenarios,
            smc_rounds=smc_rounds,
            pickup_itns=pickup_itns,
            case_management=case_management,
        ),
        range(len(burnin_rows) * len(scenarios)),
    )
    for year in range(FORECAST_START_YEAR, FORECAST_END_YEAR + 1):
        start_day = simulation_day(year, 0)
        add_malaria_summary_report(
            task,
            manifest,
            start_day=start_day,
            end_day=start_day + 365,
            reporting_interval=30,
            age_bins=[0.25, 0.5, 5, 115],
            max_number_reports=13,
            pretty_format=True,
            filename_suffix=f"Scenario_{year}",
        )

    experiment_name = (
        f"{manifest.user}_SMC_Vaccine_Forecast_{FORECAST_START_YEAR}_{FORECAST_END_YEAR}"
    )
    experiment = Experiment.from_builder(builder, task, name=experiment_name)
    print(
        f"Submitting {len(scenarios)} scenarios × {len(burnin_rows)} matched "
        f"Run_Number replicates ({len(scenarios) * len(burnin_rows)} simulations)."
    )
    print(
        "Vaccine is an explicitly generic acquisition-blocking proxy, not a "
        "validated R21 parameterization. Its toy efficacy inputs are in "
        f"{VACCINE_COVERAGE_FILE}."
    )
    experiment.run(wait_until_done=False, platform=platform)
    print(
        f"Submitted experiment {experiment.uid}; analyze it with "
        f"analyze_smc_vaccine_scenarios.py --experiment-id {experiment.uid}."
    )


if __name__ == "__main__":
    import emod_malaria.bootstrap as dtk

    dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    os.chmod(manifest.eradication_path, 0o755)
    general_sim()
