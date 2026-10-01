import os
import datetime
import pandas as pd
import numpy as np
import sys
import re
import random
from idmtools.entities import IAnalyzer	
from idmtools.entities.simulation import Simulation
import manifest
from pathlib import Path
import glob

## For plotting
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates


expts = {
        #'week2_weather' : '2c090358-cb7b-44e5-a2fd-842a6c23a5b7'
        'anaphase21_FE_example_sim_pick' : '33aa9989-b39e-424e-8607-ab62e97a1ddd'
    }

class InsetChartAnalyzer(IAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, sweep_variables=None, channels=None, working_dir=".", start_year=0):
        super(InsetChartAnalyzer, self).__init__(working_dir=working_dir, filenames=["output/InsetChart.json"])
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.inset_channels = channels or ['Statistical Population', 'New Clinical Cases', 'Blood Smear Parasite Prevalence',
                                           'Infectious Vectors']
        self.expt_name = expt_name
        self.start_year = start_year

    def map(self, data, simulation: Simulation):
        simdata = pd.DataFrame({x: data[self.filenames[0]]['Channels'][x]['Data'] for x in self.inset_channels})
        simdata['Time'] = simdata.index
        simdata['Day'] = simdata['Time'] % 365
        simdata['Year'] = simdata['Time'].apply(lambda x: int(x / 365) + self.start_year)
        simdata['date'] = simdata.apply(
            lambda x: datetime.date(int(x['Year']), 1, 1) + datetime.timedelta(int(x['Day']) - 1), axis=1)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
            elif sweep_var == 'Run_Number' :
                simdata[sweep_var] = 0
        return simdata

    def reduce(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(os.path.join(self.working_dir, self.expt_name)):
            os.mkdir(os.path.join(self.working_dir, self.expt_name))

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(os.path.join(self.working_dir, self.expt_name, 'All_Age_InsetChart.csv'), index=False)

class MonthlyPfPRAnalyzer(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir='./', start_year=0,
                 burnin=None, filenames=None, filter_exists=False):

        super(MonthlyPfPRAnalyzer, self).__init__(working_dir=working_dir,
                                                   filenames=filenames#"output/MalariaSummaryReport_Monthly_U5_2040.json"]
                                                   )
     
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.start_year = start_year
        self.burnin = burnin
        self.filter_exists = filter_exists
        # self.filenames = [f"output/MalariaSummaryReport_Monthly_U5_204{index}.json" for index in range (0, 10)]
        # print(f"self.filenames: {self.filenames}")

# Define the directory
        # dir_path = Path(simulation.get_path() + '/output/')

        # 1. Match all text files in the directory
        # malaria_report_files = list(dir_path.glob('*.json'))
        # self.filenames = [str(file) for file in malaria_report_files]
        # print(self.filenames)


    # def filter(self, simulation: Simulation):
    #     # print("#################"+simulation.get_path())
    #     base_path = "home/anaphase21/emod-tutorials/FE-2026-examples/experiments/"
    #     self.filenames = ["output/MalariaSummaryReport_Monthly_U5_2040.json"]
    #     expt_id = expts[self.expt_name]
    #     if self.filter_exists:
    #         simulation_path = os.path.join(base_path, "e_"+self.expt_name+"_"+expt_id, simulation.id, "output", self.filenames[0])
    #         # print("@@@@@@@@@@@@@@@@"+simulation_path)
    #         file = os.path.join(simulation_path)
    #         print("++++++++++++++++"+file)
    #         return os.path.exists(file)
    #     else:
    #         return True

    def filter(self, simulation: Simulation):
        if self.filter_exists:
            file = os.path.join(simulation.get_path(), self.filenames[0])
            return os.path.exists(file)
        else:
            return True
    
    def map(self, data, simulation: Simulation):
        adf = pd.DataFrame()
        fname = self.filenames[0]
        age_bins = data[self.filenames[0]]['Metadata']['Age Bins']
      
        for age in range(len(age_bins)):
            d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin'][:-1]
            pfpr = [x[age] for x in d]
          
            d = data[fname]['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin'][:-1]
            clinical_cases = [x[age] for x in d]
         
            d = data[fname]['DataByTimeAndAgeBins']['Annual Severe Incidence by Age Bin'][:-1]
            severe_cases = [x[age] for x in d]

            d = data[fname]['DataByTimeAndAgeBins']['New Infections by Age Bin'][:-1]
            new_infections = [x[age] for x in d]

            # d = data[fname]['DataByTimeAndAgeBins']['Annual Mild Anemia by Age Bin'][:-1]
            # mild_anemia = [x[age] for x in d]

            # d = data[fname]['DataByTimeAndAgeBins']['Pf Gametocyte Prevalence by Age Bin'][:-1]
            # gamete_prev = [x[age] for x in d]
            
            # d = data[fname]['DataByTimeAndAgeBins']['Mean Log Parasite Density by Age Bin'][:-1]
            # dens = [x[age] for x in d]            

            d = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin'][:-1]
            pop = [x[age] for x in d]
            # print(data[fname]['DataByTimeAndAgeBins'].keys())

            # d = data[fname]['DataByTimeAndAgeBins']['coverage by Age Bin'][:-1]
            # coverage = [x[age] for x in d]

            # d = data[fname]['DataByTimeAndAgeBins']['cm_start'][:-1]
            # cm_start = [x[age] for x in d]

            simdata = pd.DataFrame({'month': range(1, len(pfpr)+1),
                                    'PfPR': pfpr,
                                    'Cases': clinical_cases,
                                    'Severe_cases': severe_cases,
                                    'New_infections': new_infections,
                                    # 'Anemia': mild_anemia,
                                    # 'Mean_density': dens,
                                    # 'Gametocyte_prevalence': gamete_prev,
                                    'Pop': pop
                                    })
                       
            simdata['agebin'] = age_bins[age]


            adf = pd.concat([adf, simdata])

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                try:
                    adf[sweep_var] = simulation.tags[sweep_var]
                except:
                    adf[sweep_var] = '-'.join([str(x) for x in simulation.tags[sweep_var]])
        
        adf['Year'] = None
        adf['Year'] = np.floor((adf['month']-1)/12)+self.start_year
        
        return adf

    def reduce(self, all_data):

        selected = [data for sim, data in all_data.items()]
        print(len(selected))
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return

        if not os.path.exists(os.path.join(self.working_dir, self.expt_name)):
            os.mkdir(os.path.join(self.working_dir, self.expt_name))

        print(f'\nSaving outputs to: {os.path.join(self.working_dir, self.expt_name)}')

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, f'PfPR_ClinicalIncidence_Monthly_{self.start_year}.csv')),
                   index=False)
        
if __name__ == "__main__":

    from idmtools.analysis.analyze_manager import AnalyzeManager
    from idmtools.core import ItemType
    from idmtools.core.platform_factory import Platform

    start_year = 2000

    jdir = manifest.job_directory
    wdir = os.path.join(jdir, 'my_outputs')
    
    if not os.path.exists(wdir):
        os.mkdir(wdir)

    serialize_years = 50  # Same as in run_example_burnin.py
    step = 'pickup'

    sweep_variables = ['Run_Number'] 

    # set desired InsetChart channels to analyze and plot
    # channels_inset_chart = ['Statistical Population', 'True Prevalence', 'New Clinical Cases','Infectious Vectors','Rainfall','Air Temperature']
    channels_inset_chart = ['Statistical Population', 'New Clinical Cases', 
                        'Adult Vectors', 'Infected']
    sweep_variables = ['Run_Number']

    if step == 'pickup':
        sweep_variables = ['Run_Number']
    
    with Platform('Container',job_directory=jdir, docker_image=manifest.plat_image) as platform:

        for expt_name, exp_id in expts.items():
            analyzers_burnin = [InsetChartAnalyzer(expt_name=expt_name,
                                    channels=channels_inset_chart,
                                    start_year=start_year - serialize_years,
                                    sweep_variables=sweep_variables,
                                    working_dir=wdir),
                                ]

            # analyzers_pickup = [InsetChartAnalyzer(expt_name=expt_name,
            #                         channels=channels_inset_chart,
            #                         start_year=start_year,
            #                         sweep_variables=sweep_variables,
            #                         working_dir=wdir),
            #                     MonthlyPfPRAnalyzer(expt_name=expt_name,
            #                         start_year=start_year,
            #                         sweep_variables=sweep_variables,
            #                         working_dir=wdir)
            #                     ]

        if step == 'burnin':
            am = AnalyzeManager(configuration={}, ids=[(exp_id, ItemType.EXPERIMENT)],
                            analyzers=analyzers_burnin, partial_analyze_ok=True)
            am.analyze()

        elif step == 'pickup':
            for i in range(0, 10):
                malaria_summary_report_dataset_name = f'output/MalariaSummaryReport_Monthly_U5_{start_year+i}.json'
                analyzers_pickup = [MonthlyPfPRAnalyzer(expt_name=expt_name,
                                    start_year=start_year+i,
                                    sweep_variables=sweep_variables,
                                    working_dir=wdir, filenames=[malaria_summary_report_dataset_name])
                                ]
                am = AnalyzeManager(configuration={}, ids=[(exp_id, ItemType.EXPERIMENT)],
                                analyzers=analyzers_pickup, partial_analyze_ok=True)
                am.analyze()

            analyzers_pickup = [InsetChartAnalyzer(expt_name=expt_name,
                                    channels=channels_inset_chart,
                                    start_year=2000,
                                    sweep_variables=sweep_variables,
                                    working_dir=wdir)]
            am = AnalyzeManager(configuration={}, ids=[(exp_id, ItemType.EXPERIMENT)],
                                            analyzers=analyzers_pickup, partial_analyze_ok=True)
            am.analyze()
        else:
            print('Please define step, options are burnin or pickup')
            
            
 
    # read in analyzed InsetChart data
    expt_name=list(expts.keys())[0]
    inset_chart_dataset_name = 'All_Age_InsetChart.csv'
    df = pd.read_csv(os.path.join(wdir, expt_name, inset_chart_dataset_name))
    df['date'] = pd.to_datetime(df['date'])
    df = df.groupby(['date'] + sweep_variables)[channels_inset_chart].agg(np.mean).reset_index()

    # make InsetChart plot
    fig1 = plt.figure('InsetChart', figsize=(12, 6))
    fig1.subplots_adjust(hspace=0.5, left=0.08, right=0.97)
    fig1.suptitle(f'Analyzer: InsetChartAnalyzer')
    axes = [fig1.add_subplot(2, 2, x + 1) for x in range(4)]
    for ch, channel in enumerate(channels_inset_chart):
        ax = axes[ch]
        
        ax.tick_params(axis='x', labelrotation=45)

        for p, pdf in df.groupby(sweep_variables):
            ax.plot(pdf['date'], pdf[channel], '-', linewidth=0.8, label=p)
        ax.set_title(channel)
        ax.set_ylabel(channel)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=12))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    if len(sweep_variables) > 0:
        axes[-1].legend(title=', '.join(sweep_variables))
    plt.tight_layout()
    fig1.savefig(os.path.join(wdir, expt_name, 'InsetChart.png'))
    plt.clf()
    
    # expt_name=list(expts.keys())[0]
    channels_inset_chart = ['PfPR', 'Cases', 'Severe_cases', 'New_infections', 'Pop']
    malaria_report_files = os.listdir(os.path.join(wdir, expt_name))
    malaria_summary_report_dataset_names = [file for file in malaria_report_files if file.endswith('.csv') and 'PfPR_ClinicalIncidence' in file]
    # malaria_summary_report_dataset_names = [f'PfPR_ClinicalIncidence_Monthly_{2040+i}.csv' for i in range(0, 10)]
    for i, malaria_summary_report_dataset_name in enumerate(malaria_summary_report_dataset_names):
        df = pd.read_csv(os.path.join(wdir, expt_name, malaria_summary_report_dataset_names[i]))
        channels_inset_chart = ['PfPR', 'Cases', 'Severe_cases', 'New_infections', 'Pop']
        df = df.groupby(['month'] + sweep_variables)[channels_inset_chart].agg(np.mean).reset_index()

        # make InsetChart plot
        fig1 = plt.figure('InsetChart', figsize=(12, 6))
        fig1.subplots_adjust(hspace=0.5, left=0.08, right=0.97)
        fig1.suptitle(f'{malaria_summary_report_dataset_name.split(".")[0].replace("40", "00")}')
        axes = [fig1.add_subplot(2, 3, x + 1) for x in range(5)]
        for ch, channel in enumerate(channels_inset_chart):
            ax = axes[ch]
            
            ax.tick_params(axis='x', labelrotation=0)

            for p, pdf in df.groupby(sweep_variables):
                ax.plot(pdf['month'], pdf[channel], '-', linewidth=0.8, label=p)
            ax.set_title(channel)
            ax.set_ylabel(channel)
            # ax.xaxis.set_major_locator(mdates.MonthLocator(interval=12))
            # ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
            ax.set_xticks(range(1, 13))
            ax.set_xlim(1, 12)
        if len(sweep_variables) > 0:
            axes[-1].legend(title=', '.join(sweep_variables))
        plt.tight_layout()
        fig1.savefig(os.path.join(wdir, expt_name, f'Malaria_Summary_InsetChart_{malaria_summary_report_dataset_name.split(".")[0]}.png'))
        plt.clf()