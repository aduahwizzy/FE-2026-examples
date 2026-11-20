import os
import pathlib
from functools import partial

import pandas as pd
from emod_api.demographics.DemographicsTemplates import CrudeRate
from idmtools.builders import SimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment
from utils_slurm import build_burnin_df

from emodpy.emod_task import EMODTask
import emod_api.campaign as camp
import emodpy_malaria.demographics.MalariaDemographics as Demographics
import emod_api.demographics.PreDefinedDistributions as Distributions
import emodpy_malaria.interventions.bednet as itn
import emodpy_malaria.interventions.drug_campaign as dc
import emodpy_malaria.interventions.treatment_seeking as cm
from emodpy_malaria.reporters.builtin import *
from emodpy_malaria.interventions.scale_larval_habitats import (
    add_scale_larval_habitats,
)

import manifest


HERE = os.path.dirname(os.path.abspath(__file__))
# BURNIN_EXP_ID = "c696ab6d-5158-4df3-a5d1-b2f2042e0b5e"  # Set this to the experiment ID of the completed burn-in.
# BURNIN_EXP_ID = "c5da823a-7eac-4cf2-ac2b-8dc21f102212"  # Set this to the experiment ID of the completed burn-in.
# BURNIN_EXP_ID = "519d6a11-3f85-43cf-9434-ec6ed1a946f0"  # Set this to the experiment ID of the completed burn-in.
# BURNIN_EXP_ID = "d273fc32-5650-45a4-b2a2-4580b069ddc2"  # Set this to the experiment ID of the completed burn-in.
# BURNIN_EXP_ID = "2e204e86-d11e-4706-9884-fa50df3d1888"  # Set this to the experiment ID of the completed burn-in.
BURNIN_EXP_ID = "31f646d9-8002-4593-9530-2ae743835b6b"  # Set this to the experiment ID of the completed burn-in.
SERIALIZE_YEARS = 50
PICKUP_START_YEAR = 2013
PICKUP_YEARS = 13
PICKUP_START_TIME = SERIALIZE_YEARS * 365 + 1
VECTOR_SPECIES = ["gambiae", "arabiensis", "funestus"]


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

    config.parameters.Simulation_Duration = PICKUP_YEARS * 365
    config.parameters.Run_Number = 0
    config.parameters.Serialized_Population_Reading_Type = "READ"
    config.parameters.Serialization_Mask_Node_Read = 0
    config.parameters.Start_Time = PICKUP_START_TIME
    config.parameters.Custom_Individual_Events = [
        "Received_Treatment",
        "Received_ITN",
        "Received_SMC",
        "Received_Campaign_Drugs",
    ]
    return config


def build_camp():
    camp.set_schema(manifest.schema_file)
    camp.schema_path = manifest.schema_file
    return camp


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
    # demog.SetEquilibriumVitalDynamics()
    # age_distribution = Distributions.AgeDistribution_SSAfrica
    # demog.SetAgeDistribution(age_distribution)
    return demog


def pickup_day(year, source_day):
    """Convert a source's cumulative day to a day relative to pickup year."""
    return (int(year) - PICKUP_START_YEAR) * 365 + int(source_day) % 365


def pickup_absolute_day(year, source_day):
    """Convert an intervention date to the resumed simulation's day axis."""
    return PICKUP_START_TIME + pickup_day(year, source_day)


def read_pickup_interventions():
    itn_irs = pd.read_csv(
        os.path.join(HERE, "inputs", "itn_irs_interventions.csv")
    )
    pickup_itns = itn_irs[
        (itn_irs["phase"] == "pickup") & (itn_irs["intervention"] == "itn")
    ]
    pickup_irs = itn_irs[
        (itn_irs["phase"] == "pickup") & (itn_irs["intervention"] == "irs")
    ]

    smc = pd.read_csv(os.path.join(HERE, "inputs", "SMC_intervention_northern_region.csv"))
    case_management = pd.read_csv(os.path.join(HERE, "inputs", "case_management.csv"))
    case_management = case_management[case_management["phase"] == "pickup"]
    # pickup_smc = smc[smc["phase"] == "pickup"]
    # print(f"Read {len(pickup_itns)} pickup ITN records, {len(pickup_irs)} pickup IRS records, and {len(smc)} SMC records.")
    return pickup_itns, pickup_irs, smc, case_management


def build_pickup_campaign(
    scale_constant,
    scale_temp_rain,
    scale_water_veg,
):
    campaign = build_camp()
    scales_by_habitat = {
        "CONSTANT": float(scale_constant),
        "TEMPORARY_RAINFALL": float(scale_temp_rain),
        "WATER_VEGETATION": float(scale_water_veg),
    }
    species_habitats = {
        "gambiae": ("WATER_VEGETATION",),
        "arabiensis": ("CONSTANT", "TEMPORARY_RAINFALL"),
        "funestus": ("CONSTANT", "TEMPORARY_RAINFALL"),
    }
    habitat_scales = {"NodeID": [1]}
    for species, habitats in species_habitats.items():
        for habitat_type in habitats:
            habitat_scales[f"{habitat_type}.{species}"] = [
                scales_by_habitat[habitat_type]
            ]
    add_scale_larval_habitats(
        campaign,
        df=pd.DataFrame(habitat_scales),
        start_day=PICKUP_START_TIME + 1,
        repetitions=-1,
    )

    pickup_itns, pickup_irs, pickup_smc, case_management = read_pickup_interventions()
    for _, row in pickup_itns.iterrows():
        itn.add_itn_scheduled(
            campaign,
            start_day=pickup_absolute_day(row["year"], row["start_day_of_year"]),
            demographic_coverage=float(row["coverage"]),
            repetitions=1,
            node_ids=[1],
            # receiving_itn_broadcast_event="Received_ITN",
        )

    # for i in range(len(pickup_smc)):
    #     row = pickup_smc.iloc[i]
    #     smc_day = pickup_absolute_day(row["year"], row["start_day"])
    #     coverage = float(row['coverage'])
    #     dc.add_drug_campaign(
    #         campaign,
    #         campaign_type="MDA",
    #         drug_code="SPA",
    #         start_days=[smc_day],
    #         repetitions=1,
    #         coverage=coverage,
    #         target_group={"agemin": 0.25, "agemax": 5},
    #         node_ids=[1],
    #         # receiving_drugs_event_name="Received_SMC",
    #     )

    for i in range(len(case_management)):
        row = case_management.iloc[i]
        cm_start_day = pickup_absolute_day(
            row["year"], row["start_day_of_year"]
        )
        coverage = float(row["u5_effective_coverage"])
        drug = str(row['drug'])
        if('–' in drug):
            drug = drug.split('–')
        else:
            drug = drug.replace('Quinine', 'Chloroquine')
            drug = [drug]
        cm.add_treatment_seeking(
        campaign,
        start_day=cm_start_day,
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
        # broadcast_event_name="Received_Treatment",
    )

    return campaign


def attach_burnin_state_and_settings(simulation, burnin_df, row_index):
    row = burnin_df.iloc[int(row_index)]
    simulation.task.set_parameter(
        "Serialized_Population_Filenames",
        row["Serialized_Population_Filenames"],
    )
    simulation.task.set_parameter(
        "Serialized_Population_Path",
        os.path.join(row["serialized_file_path"], "output"),
    )
    simulation.task.set_parameter("Run_Number", int(row["Run_Number"]))
    simulation.task.set_parameter(
        "x_Temporary_Larval_Habitat",
        float(row["x_Temporary_Larval_Habitat"]),
    )
    simulation.task.create_campaign_from_callback(
        partial(
            build_pickup_campaign,
            scale_constant=float(row["scale_constant"]),
            scale_temp_rain=float(row["scale_temp_rain"]),
            scale_water_veg=float(row["scale_water_veg"]),
        )
    )
    simulation.task.set_parameter(
        "Custom_Individual_Events",
        [
            "Received_Treatment",
            "Received_ITN",
            "Received_SMC",
            "Received_Campaign_Drugs",
        ],
    )
    return {
        "Run_Number": int(row["Run_Number"]),
        "x_Temporary_Larval_Habitat": float(row["x_Temporary_Larval_Habitat"]),
        "scale_constant": float(row["scale_constant"]),
        "scale_temp_rain": float(row["scale_temp_rain"]),
        "scale_water_veg": float(row["scale_water_veg"]),
    }


def general_sim(selected_platform):
    if not BURNIN_EXP_ID:
        raise ValueError(
            "Set BURNIN_EXP_ID to the experiment ID of your completed burn-in."
        )

    platform = Platform(
        selected_platform,
        job_directory=manifest.job_directory,
        docker_image=manifest.plat_image,
        max_running_jobs=manifest.max_running_jobs,
        plat_image=manifest.plat_image,
        max_workers=12, 
        max_job=12
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
    task.common_assets.add_directory(
        os.path.join(HERE, "climate", "emod_weather_northern_region"),
        relative_path="climate",
    )

    burnin_df = build_burnin_df(
        BURNIN_EXP_ID, platform, SERIALIZE_YEARS * 365
    )
    pickup_itns, pickup_irs, pickup_smc, case_management = read_pickup_interventions()
    if not pickup_irs.empty:
        print(
            "IRS records are present in itn_irs_interventions.csv but are not "
            "scheduled: the file does not specify insecticide/effect parameters. "
            "Add an IRS intervention with project-specific efficacy and waning "
            "before treating these as modeled IRS campaigns."
        )
    required_columns = {
        "Run_Number",
        "x_Temporary_Larval_Habitat",
        "scale_constant",
        "scale_temp_rain",
        "scale_water_veg",
        "Serialized_Population_Filenames",
        "serialized_file_path",
    }
    missing_columns = required_columns.difference(burnin_df.columns)
    if missing_columns:
        raise ValueError(
            "The burn-in experiment is missing required sweep tags/serialization "
            f"metadata: {', '.join(sorted(missing_columns))}"
        )

    builder = SimulationBuilder()
    builder.add_sweep_definition(
        partial(attach_burnin_state_and_settings, burnin_df=burnin_df),
        range(len(burnin_df)),
    )

    add_event_recorder(
        task,
        event_list=["HappyBirthday", "Births"],
        start_day=PICKUP_START_TIME + 1,
        end_day=PICKUP_START_TIME + PICKUP_YEARS * 365,
        node_ids=[1],
        min_age_years=0,
        max_age_years=100,
    )
    for year in range(PICKUP_YEARS):
        calendar_year = PICKUP_START_YEAR + year
        start_day = PICKUP_START_TIME + year * 365
        add_malaria_summary_report(
            task,
            manifest,
            start_day=start_day,
            end_day=start_day + 365,
            reporting_interval=30,
            age_bins=[0.25, 5, 115],
            max_number_reports=13,
            pretty_format=True,
            filename_suffix=f"Monthly_U5_{calendar_year}",
        )

    user = os.getlogin()
    experiment = Experiment.from_builder(
        builder,
        task,
        name=f"{user}_FE_example_pickup_2014_2025",
    )
    print(f"Configured {len(burnin_df)} pickup simulations, one per burn-in state.")
    experiment.run(wait_until_done=False, platform=platform)
    print(
        f"Submitted pickup experiment {experiment.uid}; "
        "simulations are running asynchronously."
    )
    return experiment


if __name__ == "__main__":
    import emod_malaria.bootstrap as dtk

    dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    os.chmod(manifest.eradication_path, 0o755)
    general_sim(manifest.default_platform)
