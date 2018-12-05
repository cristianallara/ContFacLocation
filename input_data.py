# generate the data depending on the problem structure
import pandas as pd
import os


def read_data(datafolder):
    curPath = os.path.join(os.path.abspath(os.path.curdir), datafolder)

    # suppliers
    supplier_location = pd.read_csv(os.path.join(curPath, 'suppliers_location.csv'), index_col=0)
    suppliers = list(supplier_location.index)
    xi = {i: supplier_location.loc[i].values[0] for i in suppliers}
    yi = {i: supplier_location.loc[i].values[1] for i in suppliers}
    raw_material = pd.read_csv(os.path.join(curPath, 'raw_material_cost.csv'), index_col=0)

    # time periods
    time_periods = raw_material.columns.values.tolist()

    # markets
    markets_location = pd.read_csv(os.path.join(curPath, 'markets_location.csv'), index_col=0)
    markets = list(markets_location.index)
    xj = {j: markets_location.loc[j].values[0] for j in markets}
    yj = {j: markets_location.loc[j].values[1] for j in markets}

    # facilities
    facilities_data = pd.read_csv(os.path.join(curPath, 'facilities_data.csv'), index_col=0)
    centr_facilities = ['cf' + str(n) for n in range(1, int(facilities_data.loc['centr']['number']) + 1)]
    distr_facilities = ['df' + str(n) for n in range(1, int(facilities_data.loc['distr']['number']) + 1)]
    facility_types = {'distr': distr_facilities, 'centr': centr_facilities}
    facilities = list(distr_facilities + centr_facilities)
    cv = {k: facilities_data.loc[n]['cv'] for n in facility_types.keys() for k in facilities if k in facility_types[n]}
    mc = {k: facilities_data.loc[n]['mc'] for n in facility_types.keys() for k in facilities if k in facility_types[n]}

    # availability and demand
    availability = pd.read_csv(os.path.join(curPath, 'availability.csv'), index_col=0)
    a = {(i, t): availability[t][i] for i in suppliers for t in time_periods}
    demand = pd.read_csv(os.path.join(curPath, 'demand.csv'), index_col=0)
    d = {(j, t): demand[t][j] for j in markets for t in time_periods}

    # cost of raw material
    RM = {(i, t): raw_material[t][i] for i in suppliers for t in time_periods}

    # facility costs
    FIC = {(k, t): facilities_data.loc[n]['FIC'] for n in facility_types.keys() for k in facilities if
           k in facility_types[n] for t in time_periods}
    VIC = {(k, t): facilities_data.loc[n]['VIC'] for n in facility_types.keys() for k in facilities if
           k in facility_types[n] for t in time_periods}
    FOC = {(k, t): facilities_data.loc[n]['FOC'] for n in facility_types.keys() for k in facilities if
           k in facility_types[n] for t in time_periods}
    VOC = {(k, t): facilities_data.loc[n]['VOC'] for n in facility_types.keys() for k in facilities if
           k in facility_types[n] for t in time_periods}

    # transportation costs
    trans_costs = pd.read_csv(os.path.join(curPath, 'coeff_trans_costs.csv'), index_col=0, header=None)
    ft1 = {(i, k, t): trans_costs.loc['ft1'].values[0] for i in suppliers for k in facilities for t in time_periods}
    ft2 = {(k, j, t): trans_costs.loc['ft2'].values[0] for k in facilities for j in markets for t in time_periods}
    vt1 = {(i, k, t): trans_costs.loc['vt1'].values[0] for i in suppliers for k in facilities for t in time_periods}
    vt2 = {(k, j, t): trans_costs.loc['vt2'].values[0] for k in facilities for j in markets for t in time_periods}

    # interest rate
    interest_rate = 0.01
    interest_factor = {t: 1/(1+interest_rate)**(time_periods.index(t)) for t in time_periods}

    data = [suppliers, xi, yi, time_periods, markets, xj, yj, centr_facilities, distr_facilities, facilities, cv, mc,
            a, d, RM, FIC, VIC, FOC, VOC, ft1, ft2, vt1, vt2, interest_factor]

    return data

