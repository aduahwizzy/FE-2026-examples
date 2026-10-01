import pathlib
import os
import numpy as np
from functools import \
    partial 
 
#idmtools   
from idmtools.assets import Asset, AssetCollection  
from idmtools.builders import SimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment
from utils_slurm import build_burnin_df

#emodpy
from emodpy.emod_task import EMODTask
import emod_api.config.default_from_schema_no_validation as dfs
import emod_api.campaign as camp

#emodpy-malaria
import emodpy_malaria.demographics.MalariaDemographics as Demographics
import emod_api.demographics.PreDefinedDistributions as Distributions
from emodpy_malaria.reporters.builtin import *
from idmtools.builders import SimulationBuilder

import manifest

# sim_years = 10
num_seeds = 5
sim_start_year = 2040
serialize_years = 50
pickup_years=10
burnin_exp_id = '5243b85d-303f-4abb-8b16-bdc0261be1a5'

def set_param_fn(config):
    """
    This function is a callback that is passed to emod-api.config to set config parameters, including the malaria defaults.
    """
    import emodpy_malaria.malaria_config as conf
    config = conf.set_team_defaults(config, manifest)
    conf.add_species(config, manifest, ["gambiae", "arabiensis", "funestus"])
    config.parameters.Climate_Model = "CLIMATE_BY_DATA"

    # Defaults to 0 (off) in the installed Eradication build's schema - turn on so
    # InsetChart.json gets written to each simulation's output folder.
    config.parameters.Enable_Default_Reporting = 1
    config.parameters.Air_Temperature_Filename = os.path.join('climate',
        'example_air_temperature_daily.bin')
    config.parameters.Land_Temperature_Filename = os.path.join('climate',
        'example_air_temperature_daily.bin')
    config.parameters.Rainfall_Filename = os.path.join('climate',
        'example_rainfall_daily.bin')
    config.parameters.Relative_Humidity_Filename = os.path.join('climate', 
        'example_relative_humidity_daily.bin')

    config.parameters.Simulation_Duration = pickup_years*365
    config.parameters.Run_Number = 0

    config.parameters.Serialized_Population_Reading_Type = "READ"
    config.parameters.Serialization_Mask_Node_Read = 0
    config.parameters.Serialization_Time_Steps = [serialize_years*365]
    
    #Add calibrated larval habitat
    config.parameters.x_Temporary_Larval_Habitat = 0.98
    return config


def build_camp():
    """
    This function builds a campaign input file for the DTK using emod_api.
    """

    camp.set_schema(manifest.schema_file)
    
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

def update_serialize_parameters(simulation, df, x: int):

   path = df["serialized_file_path"][x]
   seed = int(df["Run_Number"][x])

   simulation.task.set_parameter("Serialized_Population_Filenames", df["Serialized_Population_Filenames"][x])
   simulation.task.set_parameter("Serialized_Population_Path", os.path.join(path, "output"))
   simulation.task.set_parameter("Run_Number", seed) #match pickup simulation run number to burnin simulation

   return {"Run_Number":seed}


def general_sim(selected_platform):
    """
    This function is designed to be a parameterized version of the sequence of things we do 
    every time we run an emod experiment. 
    """

    # Set platform and associated values, such as the maximum number of jobs to run at one time.
    # All cluster-specific settings live in manifest.py so users only edit that file.
    platform = Platform(selected_platform, job_directory=manifest.job_directory,
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
    task.common_assets.add_directory(os.path.join(manifest.input_dir,
        "example_weather", "out"), relative_path="climate")

    builder = SimulationBuilder()
    burnin_df = build_burnin_df(burnin_exp_id, platform, serialize_years*365)

    builder.add_sweep_definition(partial(update_serialize_parameters, df=burnin_df), range(len(burnin_df.index)))
    # builder.add_sweep_definition(partial(set_param, param='Run_Number'), range(num_seeds))
    builder.add_sweep_definition(partial(set_param, param='x_Temporary_Larval_Habitat'), np.logspace(-0.5,1,5))
   ## reports are still located here

   # create experiment from builder
    user = os.getlogin()
    experiment = Experiment.from_builder(builder, task, name=f'{user}_FE_example_pickup50')

    # create experiment from builder
    add_event_recorder(task, event_list=["HappyBirthday", "Births"],
                       start_day=1, end_day=pickup_years*365, 
                       node_ids=[1], min_age_years=0,
                       max_age_years=100)
    for year in range(pickup_years):
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