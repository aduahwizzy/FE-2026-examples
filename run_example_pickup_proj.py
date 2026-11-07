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
from utils_slurm import build_burnin_df  # same helper used by FE_project_template/site_setup.py

# Pickup configuration. Set BURNIN_EXPERIMENT_ID to the UID printed/stored by
# the successful burn-in experiment you want to continue from.
BURNIN_EXPERIMENT_ID = "523cbc34-5966-4503-9c88-0e5999f2e0b9"
burnin_years = 50                 # must match the burn-in serialization time
pickup_years = 10                  # duration of the continuation run
num_seeds = 3                      # new Run_Number values for each burn-in state
sim_start_year = 1963              # start year used for burn-in report naming
serialize_years = burnin_years
pickup_start_year = sim_start_year + burnin_years

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
scale_water_veg_values = list(np.logspace(-0.5, 0.7, 4))
scale_temp_rain_values = list(np.logspace(0, 2.7, 5))
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
x_temporary_larval_habitat_values = np.logspace(-0.5, 1, 4)

def build_interventions_data(phase="pickup"):
    itns = pd.read_csv('/home/anaphase21/emod-tutorials/FE-2026-examples/inputs/itn_irs_interventions.csv')
    itns = itns[(itns['intervention']=='itn') & (itns['phase']==phase)]
    case_management = pd.read_csv('/home/anaphase21/emod-tutorials/FE-2026-examples/inputs/case_management.csv')
    cm_cols = ['start_day_of_year', 'drug', 'u5_coverage', 'proportion_with_fever', 'u5_effective_coverage', 'phase']
    case_management = case_management[cm_cols]
    case_management = case_management[case_management['phase']==phase]

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

itn_coverages, cm_coverages = build_interventions_data(phase="pickup")

def build_larval_habitat_scale_df(
    habitat_index=None,
    scale_constant=None,
    scale_temp_rain=None,
    scale_water_veg=None,
):
    """Build habitat scaling factors only for habitat types configured per species.

    Pickup passes the exact scale values recorded in the matching burn-in row.
    This avoids relying on the pickup grid's ordering to recover a habitat index.
    The factors are applied to node 1 only.
    """
    supplied_scales = (scale_constant, scale_temp_rain, scale_water_veg)
    if all(value is not None for value in supplied_scales):
        scale_constant, scale_temp_rain, scale_water_veg = map(float, supplied_scales)
    elif habitat_index is not None:
        scale_constant, scale_temp_rain, scale_water_veg = (
            habitat_scale_combinations[int(habitat_index)]
        )
    else:
        scale_constant, scale_temp_rain, scale_water_veg = (1.0, 1.0, 1.0)

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
    """Configure EMOD to load the matching 50-year burn-in serialization."""
    import emodpy_malaria.malaria_config as conf

    config = conf.set_team_defaults(config, manifest)
    conf.add_species(config, manifest, VECTOR_SPECIES)

    # A serialized burn-in population can contain delayed interventions that
    # reference custom events (for example, Received_Treatment). Register these
    # before Eradication deserializes the population. ReportEventRecorder alone
    # does not register them.
    # This installed schema exposes Custom_Individual_Events (plural).
    # Do not assign Custom_Individual_Event (singular): that key is absent
    # from this schema and raises KeyError while the task is being created.
    # existing_custom_events = config.parameters.Custom_Individual_Events or []
    # if isinstance(existing_custom_events, str):
    #     existing_custom_events = [existing_custom_events]

    # custom_events = sorted(
    #     set(existing_custom_events)
    #     | {"Received_Treatment", "Received_ITN", "Received_SMC"}
    # )
    # config.parameters.Custom_Individual_Events = custom_events
    # print("Custom_Individual_Events in pickup config:", custom_events)

    config.parameters.Climate_Model = "CLIMATE_BY_DATA"

    config.parameters.Enable_Default_Reporting = 1
    config.parameters.Air_Temperature_Filename = os.path.join("climate", "air_temperature.bin")
    config.parameters.Land_Temperature_Filename = os.path.join("climate", "land_temperature.bin")
    config.parameters.Rainfall_Filename = os.path.join("climate", "rainfall.bin")
    config.parameters.Relative_Humidity_Filename = os.path.join("climate", "relative_humidity.bin")

    config.parameters.Run_Number = 0

    # These settings must match the burn-in's serialization at burnin_years * 365.
    config.parameters.Serialized_Population_Reading_Type = "READ"
    config.parameters.Serialization_Time_Steps = [365 * burnin_years]
    config.parameters.Serialization_Mask_Node_Read = 0
    config.parameters.Simulation_Duration = pickup_years * 365 + 30
    return config

def get_smc_start_days_and_coverages():
    smc_interventions = pd.read_csv('/home/anaphase21/emod-tutorials/FE-2026-examples/inputs/interventions_SMC.csv')
    smc_start_days = smc_interventions['start_day'].tolist()
    smc_coverages = smc_interventions['coverage'].tolist()
    smc_start_days_coverages = {}
    for smc_start_day, smc_coverage in zip(smc_start_days, smc_coverages):
        smc_start_days_coverages[smc_start_day] = smc_coverage
    return smc_start_days_coverages

def build_camp(
    habitat_index=None,
    scale_start_day=1,
    scale_constant=None,
    scale_temp_rain=None,
    scale_water_veg=None,
):
    """
    This function builds a campaign input file for the DTK using emod_api.
    """

    camp.set_schema(manifest.schema_file)
    camp.schema_path = manifest.schema_file

    # Apply the selected habitat-type multipliers throughout the simulation.
    # repetitions=-1 follows the approach used in FE_project_template.
    habitat_scales = build_larval_habitat_scale_df(
        habitat_index=habitat_index,
        scale_constant=scale_constant,
        scale_temp_rain=scale_temp_rain,
        scale_water_veg=scale_water_veg,
    )
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
        )

    for itn_start_day, itn_coverage in itn_coverages.items():
        itn.add_itn_scheduled(
            camp,
            start_day=itn_start_day,
            demographic_coverage=itn_coverage,
            repetitions=1,
            # timesteps_between_repetitions=365,
            node_ids=[1],
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

def update_pickup_from_burnin(simulation, burnin_index, burnin_df):
    """Load one burn-in state and restore its matching habitat parameters."""
    row = burnin_df.iloc[int(burnin_index)]

    x_tlh = float(row["x_Temporary_Larval_Habitat"])
    scale_constant = float(row["scale_constant"])
    scale_temp_rain = float(row["scale_temp_rain"])
    scale_water_veg = float(row["scale_water_veg"])

    # utils_slurm.build_burnin_df returns the run directory; serialized files
    # are in its output/ subdirectory. The helper also determines the expected
    # state filename from the burn-in timestep and core count.
    serialized_path = os.path.join(str(row["serialized_file_path"]), "output")
    serialized_filenames = row["Serialized_Population_Filenames"]
    if isinstance(serialized_filenames, str):
        serialized_filenames = [serialized_filenames]

    simulation.task.set_parameter("Serialized_Population_Filenames", serialized_filenames)
    simulation.task.set_parameter("Serialized_Population_Path", serialized_path)
    simulation.task.set_parameter("x_Temporary_Larval_Habitat", x_tlh)

    # Reapply precisely the habitat scales used to generate this serialization.
    # The values come from burn-in sweep tags, rather than from an assumed index.
    if min(scale_constant, scale_temp_rain, scale_water_veg, x_tlh) <= 0:
        raise ValueError(
            "A burn-in row has a non-positive habitat multiplier. Such a state is "
            "not used for pickup; rerun burn-in with positive habitat multipliers."
        )
    campaign_callback = partial(
        build_camp,
        scale_start_day=1,
        scale_constant=scale_constant,
        scale_temp_rain=scale_temp_rain,
        scale_water_veg=scale_water_veg,
    )
    simulation.task.create_campaign_from_callback(campaign_callback)

    # These returned values become pickup simulation tags, facilitating comparison
    # of each pickup run to the serialized burn-in state from which it started.
    tags = {
        "x_Temporary_Larval_Habitat": x_tlh,
        "scale_constant": scale_constant,
        "scale_temp_rain": scale_temp_rain,
        "scale_water_veg": scale_water_veg,
    }
    if "simid" in burnin_df.columns:
        tags["burnin_simid"] = str(row["simid"])
    return tags

def general_sim(selected_platform):
    """Run pickup simulations, one for each burn-in serialization and pickup seed."""
    platform = Platform(
        selected_platform,
        job_directory=manifest.job_directory,
        docker_image=manifest.plat_image,
        partition=manifest.partition,
        time=manifest.sim_time,
        modules=[manifest.singularity_module],
        max_running_jobs=manifest.max_running_jobs,
        plat_image=manifest.plat_image,
    )

    if BURNIN_EXPERIMENT_ID == "REPLACE_WITH_SUCCESSFUL_BURNIN_EXPERIMENT_UID":
        raise ValueError(
            "Set BURNIN_EXPERIMENT_ID to the UID of the successful burn-in experiment."
        )

    print(f"Retrieving burn-in serializations from experiment {BURNIN_EXPERIMENT_ID}...")
    burnin_df = build_burnin_df(
        BURNIN_EXPERIMENT_ID, platform, burnin_years * 365
    )
    if burnin_df.empty:
        raise RuntimeError(
            "No burn-in simulations were returned. Check the experiment UID, platform, "
            "and whether the burn-in experiment is accessible from this platform."
        )

    required_columns = {
        "serialized_file_path",
        "Serialized_Population_Filenames",
        "x_Temporary_Larval_Habitat",
        "scale_constant",
        "scale_temp_rain",
        "scale_water_veg",
    }
    missing_columns = required_columns.difference(burnin_df.columns)
    if missing_columns:
        raise KeyError(
            "The burn-in experiment is missing required sweep tags/columns: "
            f"{sorted(missing_columns)}. Confirm the burn-in was submitted with the "
            "habitat-sweep callback that returns all three scale_* tags."
        )

    # Exclude non-positive habitat multiplier rows (e.g. from an earlier run that
    # included scale_water_veg=0). Failed rows with no valid serialized state must
    # also be excluded; see the run notes below before submitting at scale.
    scale_columns = [
        "x_Temporary_Larval_Habitat",
        "scale_constant",
        "scale_temp_rain",
        "scale_water_veg",
    ]
    valid_scales = np.ones(len(burnin_df), dtype=bool)
    for column in scale_columns:
        values = pd.to_numeric(burnin_df[column], errors="coerce").to_numpy(dtype=float)
        valid_scales &= np.isfinite(values) & (values > 0)
    skipped = int((~valid_scales).sum())
    if skipped:
        print(f"Excluding {skipped} burn-in rows with non-positive/invalid scale values.")
        burnin_df = burnin_df.loc[valid_scales].reset_index(drop=True)
    if burnin_df.empty:
        raise RuntimeError("No valid positive-scale burn-in rows remain for pickup.")

    print(
        f"Found {len(burnin_df)} candidate burn-in rows. Pickup will create "
        f"{len(burnin_df)} x {num_seeds} = {len(burnin_df) * num_seeds} simulations, "
        "provided each row has a successful serialization file."
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

    # set_param_fn has already registered the custom events in the task config.
    # Avoid setting schema keys again here, especially the unsupported singular
    # spelling Custom_Individual_Event.

    # Keep the climate, demographics, species configuration, and baseline inputs
    # consistent with the burn-in. Change this directory only if pickup is meant
    # to use a different climate series.
    task.common_assets.add_directory(
        os.path.join(
            os.path.dirname(__file__), "climate", "emod_weather_northern_region"
        ),
        relative_path="climate",
    )

    builder = SimulationBuilder()
    builder.add_sweep_definition(
        partial(update_pickup_from_burnin, burnin_df=burnin_df),
        range(len(burnin_df)),
    )
    # builder.add_sweep_definition(
    #     partial(set_param, param="Run_Number"),
    #     range(num_seeds),
    # )

    add_event_recorder(
        task,
        event_list=["HappyBirthday", "Births"],
        start_day=1,
        end_day=pickup_years * 365 + 30,
        node_ids=[1],
        min_age_years=0,
        max_age_years=100,
    )

    # Produce the same monthly under-five reports as the burn-in script, but
    # over the pickup period only.
    for year in range(pickup_years):
        start_day = 365 * year
        end_day = 365 * (year + 1)
        if year == pickup_years - 1:
            end_day += 30  # optionally include January of the following year
        sim_year = pickup_start_year + year
        add_malaria_summary_report(
            task,
            manifest,
            start_day=start_day,
            end_day=end_day,
            reporting_interval=30,
            age_bins=[0.25, 5, 115],
            max_number_reports=13,
            pretty_format=True,
            filename_suffix=f"Monthly_U5_{sim_year}",
        )

    user = os.getlogin()
    experiment = Experiment.from_builder(
        builder,
        task,
        name=f"{user}_FE_example_pickup{pickup_years}_from_{BURNIN_EXPERIMENT_ID[:8]}",
    )

    experiment.run(wait_until_done=True, platform=platform)

    if not experiment.succeeded:
        print(f"Pickup experiment {experiment.uid} failed.")
        raise RuntimeError(f"Pickup experiment {experiment.uid} failed.")

    print(f"Pickup experiment {experiment.uid} succeeded.")
    print(f"Pickup experiment ID: {experiment.uid}")



if __name__ == "__main__":
    import emod_malaria.bootstrap as dtk
    import pathlib

    dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    os.chmod(manifest.eradication_path, 0o755)

    selected_platform = manifest.default_platform
    general_sim(selected_platform)