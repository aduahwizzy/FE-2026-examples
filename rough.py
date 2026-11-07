import pandas as pd
import numpy as np

smc_data = pd.read_csv(
    "/home/anaphase21/emod-tutorials/FE-2026-examples/inputs/SMC_Northern_Region.csv",
    thousands=",",
)
# print(len(smc_data.columns.tolist()))

renamed_cols = ['region', 'district', 'year', 'smc_target_round_1', 'smc_registered_round_1', 'smc_covered_round_1', 'smc_target_round_2', 'smc_registered_round_2', 'smc_covered_round_2', 'smc_target_round_3', 'smc_registered_round_3', 'smc_covered_round_3', 'smc_target_round_4', 'smc_registered_round_4', 'smc_covered_round_4', 'smc_target_round_5', 'smc_registered_round_5', 'smc_covered_round_5']
smc_data.columns = renamed_cols
northern_region_data = smc_data[smc_data['region'] == 'Northern']
northern_region_data = northern_region_data.fillna(0)
# print(smc_data.head(50))
# print(northern_region_data.head())
numeric = ['smc_target_round_1', 'smc_registered_round_1', 'smc_covered_round_1', 'smc_target_round_2', 'smc_registered_round_2', 'smc_covered_round_2', 'smc_target_round_3', 'smc_registered_round_3', 'smc_covered_round_3', 'smc_target_round_4', 'smc_registered_round_4', 'smc_covered_round_4', 'smc_target_round_5', 'smc_registered_round_5', 'smc_covered_round_5']
# northern_region_data[numeric] = northern_region_data[numeric].apply(pd.to_numeric, errors='coerce')
# for col in numeric:
#     smc_data[col] = smc_data[col].astype('Int64')
data = []
years = northern_region_data['year'].unique().tolist()
for year in years:
    year_data = northern_region_data[northern_region_data['year'] == year]
    row = year_data.iloc[0]
    year_data = [str(row['region']), year]+year_data[numeric].sum(axis=0).to_numpy().tolist()
    # year_data = year_data.groupby('year').sum().reset_index()
    data.append(year_data)
northern_region_data = pd.DataFrame(data, columns=['region', 'year'] + numeric)
# northern_region_data = pd.concat(data, ignore_index=True)
northern_region_data['coverage_round_1'] = northern_region_data['smc_covered_round_1'] / northern_region_data['smc_registered_round_1']
northern_region_data['coverage_round_2'] = northern_region_data['smc_covered_round_2'] / northern_region_data['smc_registered_round_2']
northern_region_data['coverage_round_3'] = northern_region_data['smc_covered_round_3'] / northern_region_data['smc_registered_round_3']
northern_region_data['coverage_round_4'] = northern_region_data['smc_covered_round_4'] / northern_region_data['smc_registered_round_4']
northern_region_data['coverage_round_5'] = northern_region_data['smc_covered_round_5'] / northern_region_data['smc_registered_round_5']
northern_region_data.to_csv("/home/anaphase21/emod-tutorials/FE-2026-examples/inputs/SMC_Northern_Region_Cleaned.csv", index=False)
northern_region_data = northern_region_data.fillna(0)
# print(northern_region_data.iloc[0])
# print(northern_region_data.head())
def get_smc_start_days_and_coverages():
    smc_interventions = pd.read_csv('/home/anaphase21/emod-tutorials/FE-2026-examples/inputs/interventions_SMC.csv')
    smc_start_days = smc_interventions['start_day'].tolist()
    smc_coverages = smc_interventions['coverage'].tolist()
    smc_start_days_coverages = {}
    for smc_start_day, smc_coverage in zip(smc_start_days, smc_coverages):
        smc_start_days_coverages[smc_start_day] = smc_coverage
    return smc_start_days_coverages

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
itn, cm = build_interventions_data()
a = cm[21482][0].split('–')
print(a)