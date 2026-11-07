import pathlib
import os
import numpy as np
import pandas as pd
from functools import partial
from itertools import product
 
#idmtools   
from idmtools.assets import Asset, AssetCollection
from idmtools.builders import SimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment
from emod_api.demographics.DemographicsTemplates import CrudeRate


#emodpy
from emodpy.emod_task import EMODTask
import emod_api.config.default_from_schema_no_validation as dfs
import emod_api.campaign as camp

#emodpy-malaria
import emodpy_malaria.demographics.MalariaDemographics as Demographics
import emod_api.demographics.PreDefinedDistributions as Distributions
import emodpy_malaria.interventions.treatment_seeking as cm
import emodpy_malaria.interventions.bednet as itn
import emodpy_malaria.interventions.drug_campaign as dc
from emodpy_malaria.reporters.builtin import *
from idmtools.builders import SimulationBuilder
from emodpy_malaria.interventions.scale_larval_habitats import (
    add_scale_larval_habitats,
)

import manifest

sim_years = 50
num_seeds = 1
sim_start_year = 1963
serialize_years = 50
TFR_1993 = 5.5
NORTHERN_REGION_FR_1993 = 7.39
NATIONAL_CBR_1963 = 47.71
NATIONAL_CDR_1963 = 19.23
NORTHERN_REGION_CBR_1963 = NATIONAL_CBR_1963*(NORTHERN_REGION_FR_1993/TFR_1993)
NORTHERN_REGION_CDR_1963 = NATIONAL_CDR_1963*1.0

VECTOR_SPECIES = ["gambiae", "arabiensis", "funestus"]

# Habitat types defined by the standard emodpy-malaria presets used by
# conf.add_species(...). ScaleLarvalHabitat can only scale habitat types that
# are configured for a species; it cannot create a habitat type that is absent.
#
# In the standard presets:
#   gambiae     -> WATER_VEGETATION
#   arabiensis  -> CONSTANT, TEMPORARY_RAINFALL
#   funestus    -> CONSTANT, TEMPORARY_RAINFALL
# Keep this mapping consistent if you customize species habitats in set_param_fn.
SPECIES_HABITATS = {
    "gambiae": ("WATER_VEGETATION",),
    "arabiensis": ("CONSTANT", "TEMPORARY_RAINFALL"),
    "funestus": ("CONSTANT", "TEMPORARY_RAINFALL"),
}

# Keep habitat multipliers strictly positive in this sweep. In particular,
# WATER_VEGETATION is the only configured habitat for gambiae in SPECIES_HABITATS;
# a zero multiplier can make its total larval capacity zero and may trigger a
# floating-point exception in some Eradication builds during egg laying.
# To test near-zero habitat, use a small positive value (for example 0.01) and
# confirm the behavior in a short test run.
scale_water_veg_values = list(np.logspace(-0.5, 0.7, 2))+[0.01]
scale_temp_rain_values = list(np.logspace(0.0, 1.5, 3))
scale_constant_values = [1.0]
# scale_temp_rain_values = [0.5, 1.0, 2.0]
# scale_water_veg_values = [0.5, 1.0, 2.0]

if any(value <= 0 or not np.isfinite(value)
       for value in scale_water_veg_values + scale_temp_rain_values + scale_constant_values):
    raise ValueError("Habitat scaling candidates must be finite and strictly positive.")

# All combinations of the three habitat-specific scaling factors.
# Current values give 1 * 5 * 4 = 20 habitat configurations.
habitat_scale_combinations = list(
    product(
        scale_constant_values,
        scale_temp_rain_values,
        scale_water_veg_values,
    )
)

# Preserve the original global habitat sweep, which scales ALL habitat types.
x_temporary_larval_habitat_values = np.logspace(0.1, 0.9, 3)

def build_interventions_data():
    itns = pd.read_csv('/home/anaphase21/emod-tutorials/FE-2026-examples/inputs/itn_irs_interventions.csv')
    itns = itns[(itns['intervention']=='itn') & (itns['phase']=='burnin')]
    case_management = pd.read_csv('/home/anaphase21/emod-tutorials/FE-2026-examples/inputs/case_management.csv')
    cm_cols = ['start_day_of_year', 'drug', 'u5_coverage', 'proportion_with_fever', 'u5_effective_coverage', 'phase']
    case_management = case_management[cm_cols]
    case_management = case_management[case_management['phase']=='burnin']

    itn_starts_coverages = {}
    cm_starts_coverages_drugs = {}
    
    for index, row in itns.iterrows():
        itn_starts_coverages[row['start_day_of_year']] = row['coverage']
    for i in range(len(case_management)):
        row = case_management.iloc[i]
        drug = row['drug']
        coverage = row['u5_effective_coverage']
        cm_starts_coverages_drugs[int(row['start_day_of_year'])] = (drug, float(coverage))
    return itn_starts_coverages, cm_starts_coverages_drugs

itn_coverages, cm_coverages = build_interventions_data()

def build_larval_habitat_scale_df(habitat_index=None):
    """Build habitat scaling factors only for habitat types configured per species.

    The standard presets do not give every species every habitat type. Omitting
    absent species-habitat pairs avoids Eradication's ValidateSpeciesAndHabitats
    exception. The factors are applied to node 1 only.
    """
    if habitat_index is None:
        scale_constant, scale_temp_rain, scale_water_veg = (1.0, 1.0, 1.0)
    else:
        scale_constant, scale_temp_rain, scale_water_veg = (
            habitat_scale_combinations[int(habitat_index)]
        )

    scales_by_habitat = {
        "CONSTANT": float(scale_constant),
        "TEMPORARY_RAINFALL": float(scale_temp_rain),
        "WATER_VEGETATION": float(scale_water_veg),
    }

    scales = {"NodeID": [1]}
    for species in VECTOR_SPECIES:
        if species not in SPECIES_HABITATS:
            raise ValueError(
                f"No habitat mapping is defined for species {species!r}. "
                "Update SPECIES_HABITATS to match its configured Habitats."
            )
        for habitat_type in SPECIES_HABITATS[species]:
            scales[f"{habitat_type}.{species}"] = [
                scales_by_habitat[habitat_type]
            ]

    return pd.DataFrame(scales)


def set_param_fn(config):
    """
    This function is a callback that is passed to emod-api.config to set config parameters, including the malaria defaults.
    """
    import emodpy_malaria.malaria_config as conf
    config = conf.set_team_defaults(config, manifest)
    conf.add_species(config, manifest, VECTOR_SPECIES)
    config.parameters.Climate_Model = "CLIMATE_BY_DATA"

    # Defaults to 0 (off) in the installed Eradication build's schema - turn on so
    # InsetChart.json gets written to each simulation's output folder.
    config.parameters.Enable_Default_Reporting = 1
    config.parameters.Air_Temperature_Filename = os.path.join('climate',
        'air_temperature.bin')
    config.parameters.Land_Temperature_Filename = os.path.join('climate',
        'land_temperature.bin')
    config.parameters.Rainfall_Filename = os.path.join('climate',
        'rainfall.bin')
    config.parameters.Relative_Humidity_Filename = os.path.join('climate',
        'relative_humidity.bin')

    config.parameters.Simulation_Duration = serialize_years*365
    config.parameters.Run_Number = 0

    config.parameters.Serialized_Population_Writing_Type = "TIMESTEP"
    config.parameters.Serialization_Time_Steps = [365 * serialize_years]
    config.parameters.Serialization_Mask_Node_Write = 0
    config.parameters.Serialization_Precision = "REDUCED"
    return config

def get_smc_start_days_and_coverages():
    smc_interventions = pd.read_csv('/home/anaphase21/emod-tutorials/FE-2026-examples/inputs/interventions_SMC.csv')
    smc_start_days = smc_interventions['start_day'].tolist()
    smc_coverages = smc_interventions['coverage'].tolist()
    smc_start_days_coverages = {}
    for smc_start_day, smc_coverage in zip(smc_start_days, smc_coverages):
        smc_start_days_coverages[smc_start_day] = smc_coverage
    return smc_start_days_coverages

def build_camp(habitat_index=None, scale_start_day=1):
    """
    This function builds a campaign input file for the DTK using emod_api.
    """

    camp.set_schema(manifest.schema_file)
    camp.schema_path = manifest.schema_file

    # Apply the selected habitat-type multipliers throughout the simulation.
    # repetitions=-1 follows the approach used in FE_project_template.
    habitat_scales = build_larval_habitat_scale_df(habitat_index)
    add_scale_larval_habitats(
        camp,
        df=habitat_scales,
        start_day=int(scale_start_day),
        repetitions=-1,
    )

    for cm_start_day, drug_coverage in cm_coverages.items():
        drug = drug_coverage[0]
        if('–' in drug):
            drug = drug.split('–')
        else:
            drug = drug.replace('Quinine', 'Chloroquine')
            drug = [drug]

        coverage = drug_coverage[1]
        cm.add_treatment_seeking(
            camp,
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

    for itn_start_day, itn_coverage in itn_coverages.items():
        itn.add_itn_scheduled(
            camp,
            start_day=itn_start_day,
            demographic_coverage=itn_coverage,
            repetitions=1,
            # timesteps_between_repetitions=365,
            node_ids=[1],
            # receiving_itn_broadcast_event="Received_ITN",
        )
    # itn.add_itn_scheduled(
    #     camp,
    #     start_day=ITN_START_DAY,
    #     demographic_coverage=ITN_COVERAGE,
    #     repetitions=1,
    #     # timesteps_between_repetitions=365,
    #     node_ids=[1],
    #     receiving_itn_broadcast_event="Received_ITN",
    # )
    
    # smc_start_days = [
    #     year * 365 + SMC_FIRST_ROUND_DAY + round_index * SMC_DAYS_BETWEEN_ROUNDS
    #     for year in range(sim_years)
    #     for round_index in range(SMC_ROUNDS_PER_YEAR)
    # ]
    # smc_start_days_and_coverages = get_smc_start_days_and_coverages()
    # i = 0
    # for smc_start_day, smc_coverage in smc_start_days_and_coverages.items():
    #     if(smc_start_day > serialize_years * 365):
    #         smc_start_day = i+15
    #     dc.add_drug_campaign(
    #         camp,
    #         campaign_type="MDA",
    #         drug_code="SPA",
    #         node_ids=[1],
    #         start_days=[smc_start_day],
    #         repetitions=1,
    #         coverage=smc_coverage,
    #         target_group={"agemin": 0.25, "agemax": 5},
    #         receiving_drugs_event_name="Received_SMC",
    #     )
    #     i += 1

    return camp


def build_demog():
    """
    This function builds a demographics input file for the DTK using emod_api.
    """

    demog = Demographics.from_template_node(lat=9.462781, lon=-0.285867, pop=1000, name="tamale")
    # demog.SetBirthRate(CrudeRate(NORTHERN_REGION_CBR_1963))
    # demog.SetMortalityRate(CrudeRate(NORTHERN_REGION_CDR_1963))
    # demog.SetEquilibriumAgeDistFromBirthAndMortRates(CrudeRate(NORTHERN_REGION_CBR_1963), CrudeRate(NORTHERN_REGION_CDR_1963))
    demog.SetEquilibriumVitalDynamics()

    age_distribution = Distributions.AgeDistribution_SSAfrica
    demog.SetAgeDistribution(age_distribution)
    return demog

def set_param(simulation, param, value):
    """
    Set specific parameter value
    Args:
        simulation: idmtools Simulation
        param: parameter
        value: new value
    Returns:
        dict
    """
    return simulation.task.set_parameter(param, value)

def update_habitat_scaling(simulation, habitat_index, scale_start_day=1):
    """Set a simulation's habitat-scaling campaign and record sweep metadata."""
    habitat_index = int(habitat_index)
    scale_start_day = int(scale_start_day)

    campaign_callback = partial(
        build_camp,
        habitat_index=habitat_index,
        scale_start_day=scale_start_day,
    )
    simulation.task.create_campaign_from_callback(campaign_callback)

    scale_constant, scale_temp_rain, scale_water_veg = (
        habitat_scale_combinations[habitat_index]
    )
    # Returned values become simulation metadata/tags for later analysis.
    return {
        "scale_constant": float(scale_constant),
        "scale_temp_rain": float(scale_temp_rain),
        "scale_water_veg": float(scale_water_veg),
    }


def general_sim(selected_platform):
    """
    This function is designed to be a parameterized version of the sequence of things we do 
    every time we run an emod experiment. 
    """

    # Set platform and associated values, such as the maximum number of jobs to run at one time.
    # All cluster-specific settings live in manifest.py so users only edit that file.
    platform = Platform(selected_platform, job_directory=manifest.job_directory, docker_image=manifest.plat_image,
                            partition=manifest.partition, time=manifest.sim_time,
                            modules=[manifest.singularity_module],
                            max_running_jobs=manifest.max_running_jobs,
                            plat_image=manifest.plat_image)

    # create EMODTask 
    print("Creating EMODTask (from files)...")

    
    task = EMODTask.from_default2(
        config_path="config.json",
        eradication_path=manifest.eradication_path,
        campaign_builder=build_camp,
        schema_path=manifest.schema_file,
        param_custom_cb=set_param_fn,
        ep4_custom_cb=None,
        demog_builder=build_demog,
        plugin_report=None
    )
    
    
    # set the singularity image to be used when running this experiment
    # task.set_sif(manifest.SIF_PATH, platform)
    task.common_assets.add_directory(os.path.join(os.path.dirname(__file__),
        "climate", "emod_weather_northern_region"), relative_path="climate")

    builder = SimulationBuilder()

    # Sweep over stochastic realizations, global habitat scaling, and
    # habitat-type-specific scaling combinations. These sweeps form a Cartesian
    # product, so every combination is evaluated with each Run_Number.
    builder.add_sweep_definition(
        partial(set_param, param='Run_Number'),
        range(num_seeds),
    )
    builder.add_sweep_definition(
        partial(set_param, param='x_Temporary_Larval_Habitat'),
        x_temporary_larval_habitat_values,
    )
    builder.add_multiple_parameter_sweep_definition(
        update_habitat_scaling,
        {
            "scale_start_day": [1],
            "habitat_index": list(range(len(habitat_scale_combinations))),
        },
    )

    # Configure reports before constructing the experiment so they are included
    # in the task copied into each simulation.
    add_event_recorder(
        task,
        event_list=[
            "HappyBirthday", "Births"
        ],
        start_day=1,
        end_day=serialize_years * 365,
        node_ids=[1],
        min_age_years=0,
        max_age_years=100,
    )
    for year in range(serialize_years):
        start_day = 365 * year
        sim_year = sim_start_year + year
        add_malaria_summary_report(
            task,
            manifest,
            start_day=start_day,
            end_day=365 + year * 365,
            reporting_interval=30,
            age_bins=[0.25, 5, 115],
            max_number_reports=13,
            pretty_format=True,
            filename_suffix=f'Monthly_U5_{sim_year}',
        )

    user = os.getlogin()
    experiment = Experiment.from_builder(
        builder,
        task,
        name=f'{user}_FE_example_burnin{serialize_years}_habitat_sweep',
    )

    expected_runs = (
        num_seeds
        * len(x_temporary_larval_habitat_values)
        * len(habitat_scale_combinations)
    )
    print(
        f"Configured {expected_runs} simulations: {num_seeds} seeds x "
        f"{len(x_temporary_larval_habitat_values)} global habitat values x "
        f"{len(habitat_scale_combinations)} habitat-type combinations."
    )

    # The last step is to call run() on the ExperimentManager to run the simulations.
    experiment.run(wait_until_done=True, platform=platform)


    # Check result
    if not experiment.succeeded:
        print(f"Experiment {experiment.uid} failed.\n")
        exit()

    print(f"Experiment {experiment.uid} succeeded.")



if __name__ == "__main__":
    import emod_malaria.bootstrap as dtk
    import pathlib

    dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    os.chmod(manifest.eradication_path, 0o755)

    selected_platform = manifest.default_platform
    general_sim(selected_platform)