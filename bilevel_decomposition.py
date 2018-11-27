from pyomo.environ import *
from input_data import read_data
from full_space_minlp import create_multiperiod_minlp
from full_space_mip import create_multiperiod_mip
from grid_discretization import discretize_space
import time

start_time = time.time()

# ########################################## User-defined parameters #############################################

# bilevel decomposition
max_iter_bilevel = 20
opt_tol_bilevel = 0.005         # optimality tolerance for the nested decomposition

# MILP
time_limit_mip = 300             # time limit in seconds for mip_LB
opt_tol_mip = 0.0

# grid
dist_min = 0.0                  # arbitrary
p_x = 5                         # start grid with p_x X p_y partitions
p_y = 5
n_x = 5                         # adds p_x + n_x and p_y + n_y in each iteration of the bilevel decomposition
n_y = 5

# ################################################################################################################

# Input data
data = read_data()

# Create minlp model
minlp = create_multiperiod_minlp(data)
# nlpsolver = SolverFactory('gams')
# nlpsolver.solve(minlp,
#                 tee=True,
#                 add_options=['option reslim=3600; option optcr = 0.005;'],
#                 # keepfiles=True,
#                 solver='baron',
#                 load_solutions=True
#                 )


# Tightest rectangle that includes all points
x_min = min(min(minlp.suppl_x[i] for i in minlp.suppl), min(minlp.mkt_x[j] for j in minlp.mkt))
x_max = max(max(minlp.suppl_x[i] for i in minlp.suppl), max(minlp.mkt_x[j] for j in minlp.mkt))
y_min = min(min(minlp.suppl_y[i] for i in minlp.suppl), min(minlp.mkt_y[j] for j in minlp.mkt))
y_max = max(max(minlp.suppl_y[i] for i in minlp.suppl), max(minlp.mkt_y[j] for j in minlp.mkt))

mip_sol = {}
nlp_sol = {}
opt_gap = {}

iter_list = range(1, max_iter_bilevel+1)
for iter_ in iter_list:

    mip = create_multiperiod_mip(data, p_x, p_y)
    dist_supp, dist_mkt, xp_mid, yp_mid, length_x, length_y = discretize_space(mip, x_min, x_max, y_min,
                                                                                       y_max, p_x, p_y, dist_min)

    for i in mip.suppl:
        for p in mip.part:
            # print(i, p, dist_supp[i, p])
            mip.dist_supp[i, p] = dist_supp[i, p]
    for j in mip.mkt:
        for p in mip.part:
            # print(j, p, dist_mkt[j, p])
            mip.dist_mkt[j, p] = dist_mkt[j, p]

    # Solve MIP for LB
    mipsolver = SolverFactory('gurobi')
    mipsolver.options['mipgap'] = opt_tol_mip
    mipsolver.options['timelimit'] = time_limit_mip
    mipsolver.options['threads'] = 6
    mipsolver.solve(mip)  # , tee=True)

    mip_sol[iter_] = mip.obj()
    LB = max(mip_sol[idx] for idx in iter_list if idx <= iter_)

    w_p = {}
    z_supp2fac_p = {}
    z_fac2mkt_fx = {}
    for t in mip.t:
        for p in mip.part:
            for k in mip.fac:
                w_p[k, p, t] = round(mip.w[k, p, t].value)
                if w_p[k, p, t] != 0:
                    print (k, p, t, w_p[k, p, t])
                for i in mip.suppl:
                    z_supp2fac_p[i, k, p, t] = round(mip.z_supp2fac[i, k, p, t].value)
                for j in mip.mkt:
                    z_fac2mkt_fx[k, j, p, t] = round(mip.z_fac2mkt[k, j, p, t].value)

    # Fix discrete solutions
    w_fx = {(k, t): sum(w_p[k, p, t] for p in mip.part) for k in mip.fac for t in mip.t}
    for t in mip.t:
        for k in minlp.fac:
            minlp.w[k, t].fix(w_fx[k, t])
            # print('w', k, t, w_fx[k, t])

    z_supp2fac_fx = {(i, k, t): sum(z_supp2fac_p[i, k, p, t] for p in mip.part) for i in mip.suppl
                     for k in mip.fac for t in mip.t}
    for t in mip.t:
        for i in minlp.suppl:
            for k in minlp.fac:
                minlp.z_supp2fac[i, k, t].fix(z_supp2fac_fx[i, k, t])

    z_fac2mkt_fx = {(k, j, t): sum(z_fac2mkt_fx[k, j, p, t] for p in mip.part) for j in mip.mkt
                    for k in mip.fac for t in mip.t}
    for t in mip.t:
        for k in minlp.fac:
            for j in minlp.mkt:
                minlp.z_fac2mkt[k, j, t].fix(z_fac2mkt_fx[k, j, t])

    # Update bounds for (fac_x, fac_y)
    fac_x_max = {}
    fac_x_min = {}
    fac_y_max = {}
    fac_y_min = {}
    for k in mip.fac:
        if w_fx[k, t] == 0:
            fac_x_max[k] = 0
            fac_x_min[k] = 0
            fac_y_max[k] = 0
            fac_y_min[k] = 0
        else:
            for p in mip.part:
                if w_p[k, p, t] == 1.0:
                    fac_x_max[k] = xp_mid[p] + length_x/2
                    fac_x_min[k] = xp_mid[p] - length_x/2
                    fac_y_max[k] = yp_mid[p] + length_y/2
                    fac_y_min[k] = yp_mid[p] - length_y/2
    for k in mip.fac:
        minlp.x_max[k] = fac_x_max[k]
        minlp.x_min[k] = fac_x_min[k]
        minlp.y_max[k] = fac_y_max[k]
        minlp.y_min[k] = fac_y_min[k]
        # print(k, 'x_max: ', fac_x_max[k], 'x_min: ', fac_x_min[k])
        # print(k, 'y_max: ', fac_y_max[k], 'y_min: ', fac_y_min[k])

    # Solve NLP for UB
    nlpsolver = SolverFactory('gams')
    nlpsolver.solve(minlp,
                    # tee=True,
                    add_options=['option reslim=300; option optcr = 0.0;'],
                    # keepfiles=True,
                    solver='baron',
                    load_solutions=True
                    )
    nlp_sol[iter_] = minlp.obj()
    minlp.fac_x.pprint()
    minlp.fac_y.pprint()
    print('solution of iteration', iter_, ': ', nlp_sol)
    UB = min(nlp_sol[idx] for idx in iter_list if idx <= iter_)
    opt_gap[iter_] = (UB - LB) / UB
    print('LB: ', LB)
    print('UB: ', UB)
    print('opt gap: ', opt_gap[iter_])

    if opt_gap[iter_] <= opt_tol_bilevel or iter_ == max_iter_bilevel:
        fac_x = {k: minlp.fac_x[k].value}
        fac_y = {k: minlp.fac_y[k].value}
        cost = UB
        elapsed_time = time.time() - start_time
        print("Solution Time (s)", elapsed_time)
        break
    else:
        p_x += n_x
        p_y += n_y





