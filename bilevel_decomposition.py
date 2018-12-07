from pyomo.environ import *
from input_data import read_data
from full_space_minlp import create_multiperiod_minlp
from full_space_mip import create_multiperiod_mip
from uniform_grid import discretize_space
from nonuniform_grid import refine_gride
import time

start_time = time.time()

# ########################################## User-defined parameters #############################################

# case study from folder
datafolder = 'data1_singleplant'

# bilevel decomposition
max_iter_bilevel = 2
opt_tol_bilevel = 0.01           # optimality tolerance for the nested decomposition

# MILP
time_limit_mip = 3600           # time limit in seconds for mip_LB
opt_tol_mip = 0.0

# grid
dist_min = 0.5                  # arbitrary
p_x = 2                        # start grid with p_x*p_y partitions
p_y = 2
n_x = 2                         # adds p_x + n_x and p_y + n_y in each iteration of the bilevel decomposition
n_y = 2

# ################################################################################################################

# Input data
data = read_data(datafolder)

# Create minlp model
minlp = create_multiperiod_minlp(data, dist_min)
# nlpsolver = SolverFactory('gams')
# nlpsolver.solve(minlp,
#                 tee=True,
#                 add_options=['option reslim=3600; option optcr = 0.0;'],
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

# Star algorithm with uniform grid p_x * p_y

select_part = {}
for p in range(1, p_x*p_y + 1):
    select_part[p] = 0

mip = create_multiperiod_mip(data, p_x*p_y, select_part, if_complement=False)
dist_supp, dist_mkt, max_dist_supp, max_dist_mkt, xp_min, xp_max, yp_min, yp_max, mapping \
    = discretize_space(mip, x_min, x_max, y_min, y_max, p_x, p_y)

for iter_ in iter_list:

    # ################################### Master problem ###################################
    for i in mip.suppl:
        for p in mip.part:
            if dist_supp[i, p] >= dist_min:
                mip.dist_supp[i, p] = dist_supp[i, p]
            else:
                mip.dist_supp[i, p] = dist_min
            if max_dist_supp[i, p] <= dist_min:
                for k in mip.fac:
                    for t in mip.t:
                        mip.w[k, p, t].fix(0.0)
    for j in mip.mkt:
        for p in mip.part:
            if dist_mkt[j, p] >= dist_min:
                mip.dist_mkt[j, p] = dist_mkt[j, p]
            else:
                mip.dist_mkt[j, p] = dist_min
            if max_dist_mkt[j, p] <= dist_min:
                for k in mip.fac:
                    for t in mip.t:
                        mip.w[k, p, t].fix(0.0)

    # Solve MIP for LB
    mipsolver = SolverFactory('gurobi')
    mipsolver.options['mipgap'] = opt_tol_mip
    mipsolver.options['timelimit'] = time_limit_mip
    mipsolver.options['threads'] = 6
    mipsolver.solve(mip, tee=True)

    mip_sol[iter_] = mip.obj()
    LB = mip_sol[iter_]

    w_p = {}
    b_p = {}
    z_supp2fac_p = {}
    z_fac2mkt_fx = {}
    for t in mip.t:
        for p in mip.part:
            for k in mip.fac:
                w_p[k, p, t] = round(mip.w[k, p, t].value)
                b_p[k, p, t] = round(mip.b[k, p, t].value)
                if w_p[k, p, t] != 0:
                    print((k, p, t), w_p[k, p, t])
                for i in mip.suppl:
                    z_supp2fac_p[i, k, p, t] = round(mip.z_supp2fac[i, k, p, t].value)
                for j in mip.mkt:
                    z_fac2mkt_fx[k, j, p, t] = round(mip.z_fac2mkt[k, j, p, t].value)

    # ##################################### Subproblem #########################################

    # Fix discrete solutions
    w_fx = {(k, t): sum(w_p[k, p, t] for p in mip.part) for k in mip.fac for t in mip.t}
    for t in mip.t:
        for k in minlp.fac:
            minlp.w[k, t].fix(w_fx[k, t])
            print('w', k, t, w_fx[k, t])

    b_fx = {(k, t): sum(b_p[k, p, t] for p in mip.part) for k in mip.fac for t in mip.t}
    for t in mip.t:
        for k in minlp.fac:
            minlp.b[k, t].fix(b_fx[k, t])
            print('b', k, t, b_fx[k, t])

    z_supp2fac_fx = {(i, k, t): sum(z_supp2fac_p[i, k, p, t] for p in mip.part) for i in mip.suppl
                     for k in mip.fac for t in mip.t}
    for t in mip.t:
        for i in minlp.suppl:
            for k in minlp.fac:
                minlp.z_supp2fac[i, k, t].fix(z_supp2fac_fx[i, k, t])
                print('z_supp', i, k, t, z_supp2fac_fx[i, k, t])

    z_fac2mkt_fx = {(k, j, t): sum(z_fac2mkt_fx[k, j, p, t] for p in mip.part) for j in mip.mkt
                    for k in mip.fac for t in mip.t}
    for t in mip.t:
        for k in minlp.fac:
            for j in minlp.mkt:
                minlp.z_fac2mkt[k, j, t].fix(z_fac2mkt_fx[k, j, t])
                print('z_mkt', k, j, t, z_fac2mkt_fx[k, j, t])

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
                    fac_x_max[k] = xp_max[p]
                    fac_x_min[k] = xp_min[p]
                    fac_y_max[k] = yp_max[p]
                    fac_y_min[k] = yp_min[p]
    for k in mip.fac:
        minlp.x_max[k] = fac_x_max[k]
        minlp.x_min[k] = fac_x_min[k]
        minlp.y_max[k] = fac_y_max[k]
        minlp.y_min[k] = fac_y_min[k]
        print(k, 'x_max: ', fac_x_max[k], 'x_min: ', fac_x_min[k])
        print(k, 'y_max: ', fac_y_max[k], 'y_min: ', fac_y_min[k])

    # Solve NLP for UB
    nlpsolver = SolverFactory('gams')
    nlpsolver.solve(minlp,
                    tee=True,
                    add_options=['option reslim=300; option optcr = 0.0;'],
                    keepfiles=True,
                    solver='baron',
                    load_solutions=True,
                    symbolic_solver_labels=True
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

    # ################################### Optimality check and Pruning ###################################

    if opt_gap[iter_] <= opt_tol_bilevel or iter_ == max_iter_bilevel:
        fac_x = {k: minlp.fac_x[k].value}
        fac_y = {k: minlp.fac_y[k].value}
        cost = UB
        elapsed_time = time.time() - start_time
        print("Solution Time (s)", elapsed_time)
        break
    else:
        # Choose active partition
        select_part = {}
        for p in mip.part:
            if sum(w_p[k, p, t] for k in mip.fac for t in mip.t) == 0:
                select_part[p] = 0
            else:
                select_part[p] = 1
                for k in mip.fac:
                    if sum(w_p[k, p, t] for t in mip.t) >= 1:
                        if xp_min[p] == minlp.fac_x[k].value:
                            for pp in mip.part:
                                if xp_max[pp] == xp_min[p] and yp_min[pp] <= minlp.fac_y[k].value <= yp_max[pp]:
                                    select_part[pp] = 1
                        if xp_max[p] == minlp.fac_x[k].value:
                            for pp in mip.part:
                                if xp_min[pp] == xp_max[p] and yp_min[pp] <= minlp.fac_y[k].value <= yp_max[pp]:
                                    select_part[pp] = 1
                        if yp_min[p] == minlp.fac_y[k].value:
                            for pp in mip.part:
                                if yp_max[pp] == yp_min[p] and xp_min[pp] <= minlp.fac_x[k].value <= xp_max[pp]:
                                    select_part[pp] = 1
                        if yp_max[p] == minlp.fac_y[k].value:
                            for pp in mip.part:
                                if yp_min[pp] == yp_max[p] and xp_min[pp] <= minlp.fac_x[k].value <= xp_max[pp]:
                                    select_part[pp] = 1
        print(select_part)

        # Solve MILP to check LB if at least one on the non-selected partitions has to be picked
        n_part = len(mip.part)
        mip_comp = create_multiperiod_mip(data, n_part, select_part, if_complement=True)
        for i in mip_comp.suppl:
            for p in mip_comp.part:
                if dist_supp[i, p] >= dist_min:
                    mip_comp.dist_supp[i, p] = dist_supp[i, p]
                else:
                    mip_comp.dist_supp[i, p] = dist_min
                if max_dist_supp[i, p] <= dist_min:
                    for k in mip_comp.fac:
                        for t in mip_comp.t:
                            mip_comp.w[k, p, t].fix(0.0)
        for j in mip_comp.mkt:
            for p in mip_comp.part:
                if dist_mkt[j, p] >= dist_min:
                    mip_comp.dist_mkt[j, p] = dist_mkt[j, p]
                else:
                    mip_comp.dist_mkt[j, p] = dist_min
                if max_dist_mkt[j, p] <= dist_min:
                    for k in mip_comp.fac:
                        for t in mip_comp.t:
                            mip_comp.w[k, p, t].fix(0.0)

        # Solve MILP with added constraint for LB^c
        mipcompsolver = SolverFactory('gurobi')
        mipcompsolver.options['mipgap'] = opt_tol_mip
        mipcompsolver.options['timelimit'] = 100
        mipcompsolver.options['threads'] = 6
        mipcompsolver.solve(mip_comp)
        LB_comp = mip_comp.obj()
        print('LB^c: ', LB_comp)
        if LB_comp > UB:
            prune = True
        else:
            prune = False
        print('Would it be pruned: ', prune)
        # TODO: add pruning routine

        dist_supp, dist_mkt, max_dist_supp, max_dist_mkt, xp_min, xp_max, yp_min, yp_max, mapping \
            = refine_gride(mip, xp_min, xp_max, yp_min, yp_max, n_x, n_y, select_part, dist_supp, dist_mkt,
                           max_dist_supp, max_dist_mkt, mapping)
        n_part = len(mapping)
        mip = create_multiperiod_mip(data, n_part, {}, if_complement=False)


