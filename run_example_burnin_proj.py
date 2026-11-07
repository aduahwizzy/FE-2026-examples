import pathlib
import os
import numpy as np
import pandas as pd
from functools import \
    partial 
 
#idmtools   
from idmtools.assets import Asset, AssetCollection
from idmtools.builders import SimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment


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
from itertools import product
from emodpy_malaria.interventions.scale_larval_habitats import (
    add_scale_larval_habitats,
)

import manifest

sim_years = 50
num_seeds = 3
sim_start_year = 1948
serialize_years = 65

# Example intervention defaults; adjust these values for your campaign.
CM_COVERAGE_U5 = 0.8
ITN_START_DAY = 165
ITN_COVERAGE = 0.582
SMC_COVERAGE = 0.90
SMC_FIRST_ROUND_DAY = 180
SMC_ROUNDS_PER_YEAR = 4
SMC_DAYS_BETWEEN_ROUNDS = 30

VECTOR_SPECIES = ["gambiae", "arabiensis", "funestus"]

# Replace these with the candidate values you want to test.
scale_constant_values = [1.0]
scale_temp_rain_values = [0.5, 1.0, 2.0]
scale_water_veg_values = [0.5, 1.0, 2.0]

habitat_scale_combinations = list(
    product(
        scale_constant_values,
        scale_temp_rain_values,
        scale_water_veg_values,
    )
)

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

    cm.add_treatment_seeking(
        camp,
        start_day=1,
        drug=["Artemether", "Lumefantrine"],
        targets=[
            {
                "trigger": "NewClinicalCase",
                "coverage": CM_COVERAGE_U5,
                "agemin": 0,
                "agemax": 5,
                "rate": 0.3,
            },
            {
                "trigger": "NewClinicalCase",
                "coverage": CM_COVERAGE_U5 * 0.75,
                "agemin": 5,
                "agemax": 115,
                "rate": 0.3,
            },
            {
                "trigger": "NewSevereCase",
                "coverage": min(CM_COVERAGE_U5 * 1.15, 1),
                "agemin": 0,
                "agemax": 115,
                "rate": 0.5,
            },
        ],
        broadcast_event_name="Received_Treatment",
    )

    itn.add_itn_scheduled(
        camp,
        start_day=ITN_START_DAY,
        demographic_coverage=ITN_COVERAGE,
        repetitions=1,
        # timesteps_between_repetitions=365,
        node_ids=[1],
        receiving_itn_broadcast_event="Received_ITN",
    )

    add_scale_larval_habitats(
    camp,
    df=habitat_scales,
    start_day=int(scale_start_day),
    repetitions=-1,
)
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

    demog = Demographics.from_template_node(lat=1, lon=2, pop=1000, name="Navrongo")

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

    builder.add_sweep_definition(partial(set_param, param='Run_Number'), range(num_seeds))
    builder.add_sweep_definition(partial(set_param, param='x_Temporary_Larval_Habitat'), np.logspace(-0.5,1,3))
   ## reports are still located her

   # create experiment from builder
    user = os.getlogin()
    experiment = Experiment.from_builder(builder, task, name=f'{user}_FE_example_burnin50')

    # create experiment from builder
    add_event_recorder(task, event_list=[
                           "HappyBirthday", "Births", "Received_Treatment",
                           "Received_ITN", "Received_SMC"
                       ],
                       start_day=1, end_day=serialize_years*365, 
                       node_ids=[1], min_age_years=0,
                       max_age_years=100)
    for year in range(serialize_years):
        start_day = 0 + 365 * year
        sim_year = sim_start_year + year
        add_malaria_summary_report(task, manifest, start_day=start_day,
                                    end_day=365+year*365, reporting_interval=30,
                                    age_bins=[0.25, 5, 115],
                                    max_number_reports=13,
                                    pretty_format=True, 
                                    filename_suffix=f'Monthly_U5_{sim_year}')

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