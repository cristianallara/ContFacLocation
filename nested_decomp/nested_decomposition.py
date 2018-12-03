from pyomo.environ import *
from input_data import read_data
from nested_decomp.minlp_formulation import create_minlp
from nested_decomp.mip_LB import create_mip
from nested_decomp.lagrangean_relaxation_mip import create_LR
from uniform_grid import discretize_space
from nested_decomp.forward import forward_pass
from nested_decomp.backward import backward_pass
import time

from full_space_mip import create_multiperiod_mip

start_time = time.time()

# ########################################## User-defined parameters #############################################

# nested decomposition
max_iter_nested = 1000
opt_tol_nested = 0.01           # optimality tolerance for the nested decomposition
benders_cut = True
integer_cut = False
strength_benders_cut = False

# bilevel decomposition
max_iter_bilevel = 1
opt_tol_bilevel = 0.005         # optimality tolerance for the nested decomposition

# MILP
time_limit_mip = 100             # time limit in seconds for mip_LB
opt_tol_mip = 0.0001

# grid
dist_min = 0.5                  # arbitrary
p_x = 2                         # start grid with p_x X p_y partitions
p_y = 2
n_x = 5                         # adds p_x + n_x and p_y + n_y in each iteration of the bilevel decomposition
n_y = 5

# ################################################################################################################


# Input data
data = read_data()

# Create minlp model
minlp = create_minlp(data)

# Tightest rectangle that includes all points
x_min = min(min(minlp.suppl_x[i] for i in minlp.suppl), min(minlp.mkt_x[j] for j in minlp.mkt))
x_max = max(max(minlp.suppl_x[i] for i in minlp.suppl), max(minlp.mkt_x[j] for j in minlp.mkt))
y_min = min(min(minlp.suppl_y[i] for i in minlp.suppl), min(minlp.mkt_y[j] for j in minlp.mkt))
y_max = max(max(minlp.suppl_y[i] for i in minlp.suppl), max(minlp.mkt_y[j] for j in minlp.mkt))

# #################################### Full-space mip #######################################################
full_mip = create_multiperiod_mip(data, p_x, p_y)
dist_supp, dist_mkt, xp_mid, yp_mid, length_x, length_y = discretize_space(full_mip, x_min, x_max, y_min,
                                                                           y_max, p_x, p_y, dist_min)
for i in full_mip.suppl:
    for p in full_mip.part:
        full_mip.dist_supp[i, p] = dist_supp[i, p]
for j in full_mip.mkt:
    for p in full_mip.part:
        full_mip.dist_mkt[j, p] = dist_mkt[j, p]
opt = SolverFactory('gurobi')
opt.options['mipgap'] = opt_tol_mip
opt.options['timelimit'] = time_limit_mip
opt.options['threads'] = 6
opt.options['FeasibilityTol'] = 1e-9
opt.solve(full_mip, tee=True)
elapsed_time = time.time() - start_time
print('Full space MIP time (s):', elapsed_time)
start_time = time.time()
full_mip.pprint()

# ###########################################################################################################

# create MIP relaxation for every time period
mip = create_mip(data, p_x, p_y)

# create Lagrangean relaxation for every time period
lr = create_LR(data, p_x, p_y)

# define grid discretization
dist_supp, dist_mkt, xp_mid, yp_mid, length_x, length_y = discretize_space(mip, x_min, x_max, y_min,
                                                                           y_max, p_x, p_y, dist_min)

# Retrieve duals for cut generation
for t in mip.t:
    mip.Bl[t].dual = Suffix(direction=Suffix.IMPORT)

for i in mip.suppl:
    for p in mip.part:
        mip.dist_supp[i, p] = dist_supp[i, p]
for j in mip.mkt:
    for p in mip.part:
        mip.dist_mkt[j, p] = dist_mkt[j, p]


mip.iters = RangeSet(max_iter_nested)
w_par_iter = {}
cost_fp_iter = {}
tot_cost_forward = {}
nested_UB = {}
mltp_iter = {}
cost_bp_iter = {}
nested_LB = {}
nested_gap = {}
cost_LR_iter = {}
cost_LP_iter = {}
z_supp2fac_par = {}
z_fac2mkt_par = {}

for it in range(1, max_iter_bilevel + 1):

    for iter_ in mip.iters:
        print("iteration: ", iter_)

        # Forward Pass
        for t in minlp.t:
            print(t)

            w, cost_forward = forward_pass(mip.Bl[t], mip.Bl[t].w, opt_tol_mip, time_limit_mip)

            for k in mip.fac:
                for p in mip.part:
                    mip.w_prev_par[k, p, t] = w[k, p]
                    if strength_benders_cut:
                        lr.w_prev_par[k, p, t] = w[k, p]
                    # if w[k, p] != 0:
                    #     print(k, p, w[k, p])
                    w_par_iter[k, p, t, iter_] = w[k, p]

                    # store other discrete variables for bilevel decomposition
                    for i in mip.suppl:
                        z_supp2fac_par[i, k, p, iter_] = mip.Bl[t].z_supp2fac[i, k, p].value
                    for j in mip.mkt:
                        z_fac2mkt_par[k, j, p, iter_] = mip.Bl[t].z_fac2mkt[k, j, p].value

            cost_fp_iter[t, iter_] = cost_forward
            print('cost', cost_fp_iter[t, iter_])

        tot_cost_forward[iter_] = sum(cost_fp_iter[t, iter_] for t in mip.t)
        nested_UB = min(tot_cost_forward[kk] for kk in mip.iters if kk <= iter_)
        print('nested decomp UB: ', nested_UB)
        elapsed_time = time.time() - start_time
        print('elapsed time (s):', elapsed_time)
        print('------------------------------------------------')

        # Backward Pass
        for t in reversed(list(mip.t)):
            print(t)

            cost_backward, cost_LR, cost_LP, mltp = backward_pass(mip.t.ord(t), mip.Bl[t], mip.Bl[t].equal,
                                                                  time_limit_mip, lr.Bl[t], benders_cut, integer_cut,
                                                                  strength_benders_cut)
            cost_bp_iter[t, iter_] = cost_backward

            if mip.t.ord(t) != 1:
                if strength_benders_cut:
                    cost_LR_iter[t, iter_] = cost_LR
                if benders_cut:
                    cost_LP_iter[t, iter_] = cost_LP

                t_prev = mip.t[mip.t.ord(t) - 1]
                for k in mip.fac:
                    for p in mip.part:
                        mltp_iter[k, p, t, iter_] = mltp[k, p]

                # add Benders cut for cost-to-go function approximation
                if benders_cut:
                    mip.Bl[t_prev].fut_cost.add(expr=(mip.Bl[t_prev].alphafut >= cost_LP_iter[t, iter_]
                                                      + sum(mltp_iter[k, p, t, iter_] *
                                                            (w_par_iter[k, p, t_prev, iter_] -
                                                             mip.Bl[t_prev].w[k, p])
                                                            for k in mip.fac for p in mip.part)))
                    if strength_benders_cut:
                        lr.Bl[t_prev].fut_cost.add(expr=(lr.Bl[t_prev].alphafut >= cost_LP_iter[t, iter_]
                                                         + sum(mltp_iter[k, p, t, iter_] *
                                                               (w_par_iter[k, p, t_prev, iter_] -
                                                                lr.Bl[t_prev].w[k, p])
                                                               for k in lr.fac for p in lr.part)))

                # add integer optimality cut to cost-to-go function approximation
                if integer_cut:
                    mip.Bl[t_prev].fut_cost.add(expr=(mip.Bl[t_prev].alphafut >= cost_bp_iter[t, iter_]
                                                      * (sum((w_par_iter[k, p, t_prev, iter_] - 1) *
                                                             mip.Bl[t_prev].w[k, p] for k in mip.fac for p in mip.part)
                                                         + sum((mip.Bl[t_prev].w[k, p] - 1) *
                                                               w_par_iter[k, p, t_prev, iter_]
                                                               for k in mip.fac for p in mip.part))
                                                      + cost_bp_iter[t, iter_]))
                    if strength_benders_cut:
                        lr.Bl[t_prev].fut_cost.add(expr=(lr.Bl[t_prev].alphafut >= cost_bp_iter[t, iter_]
                                                         * (sum((w_par_iter[k, p, t_prev, iter_] - 1) *
                                                                lr.Bl[t_prev].w[k, p] for k in lr.fac for p in lr.part)
                                                            + sum((lr.Bl[t_prev].w[k, p] - 1) *
                                                                  w_par_iter[k, p, t_prev, iter_]
                                                                  for k in lr.fac for p in lr.part))
                                                         + cost_bp_iter[t, iter_]))

                # add Strengthened Benders cut for cost-to-go function approximation
                if strength_benders_cut:
                    mip.Bl[t_prev].fut_cost.add(expr=(mip.Bl[t_prev].alphafut >= cost_LR_iter[t, iter_]
                                                      + sum(mltp_iter[k, p, t, iter_] *
                                                            (w_par_iter[k, p, t_prev, iter_] -
                                                             mip.Bl[t_prev].w[k, p])
                                                            for k in mip.fac for p in mip.part)))
                    lr.Bl[t_prev].fut_cost.add(expr=(lr.Bl[t_prev].alphafut >= cost_LR_iter[t, iter_]
                                                     + sum(mltp_iter[k, p, t, iter_] *
                                                           (w_par_iter[k, p, t_prev, iter_] -
                                                            lr.Bl[t_prev].w[k, p]) for k in lr.fac for p in lr.part)))

            print('cost', cost_bp_iter[t, iter_])
            print('alphafut', mip.Bl[t].alphafut.value)

        # Compute lower bound
        nested_LB[iter_] = cost_bp_iter[mip.t[1], iter_]
        print('nested decomp LB: ', nested_LB[iter_])
        # Compute optimality gap
        nested_gap[iter_] = (nested_UB - nested_LB[iter_]) / nested_UB
        print('nested opt gap: ', nested_gap[iter_])
        print('------------------------------------------------')

        if nested_gap[iter_] <= opt_tol_nested:
            last_iter = iter_
            print("nested decomposition converged")
            break

    elapsed_time = time.time() - start_time
    print("Solution Time (s)", elapsed_time)

    # w_p = {}
    # z_supp2fac_p = {}
    # z_fac2mkt_p = {}
    # for t in mip.t:
    #     for p in mip.part:
    #         for k in mip.fac:
    #             w_p[k, p, t] = round(w_par_iter[k, p, t, last_iter])
    #             for i in mip.suppl:
    #                 z_supp2fac_p[i, k, p] = round(z_supp2fac_par[i, k, p, last_iter])
    #             for j in mip.mkt:
    #                 z_fac2mkt_p[k, j, p] = round(z_fac2mkt_par[k, j, p, last_iter])
    #
    # # Fix discrete solutions
    # w_fx = {(k, t): sum(w_p[k, p, t] for p in mip.part) for k in mip.fac for t in mip.t}
    # for t in mip.t:
    #     for k in mip.fac:
    #         minlp.Bl[t].w[k].fix(w_fx[k, t])
    #         # print('w', k, w_fx[k])
    #
    # z_supp2fac_fx = {(i, k): sum(z_supp2fac_p[i, k, p] for p in mip.part) for i in mip.suppl
    #                  for k in mip.fac}
    # for t in mip.t:
    #     for i in mip.suppl:
    #         for k in minlp.fac:
    #             minlp.Bl[t].z_supp2fac[i, k].fix(z_supp2fac_fx[i, k])
    #
    # z_fac2mkt_fx = {(k, j): sum(z_fac2mkt_p[k, j, p] for p in mip.part) for j in mip.mkt
    #                 for k in mip.fac}
    # for t in mip.t:
    #     for k in mip.fac:
    #         for j in minlp.mkt:
    #             minlp.Bl[t].z_fac2mkt[k, j].fix(z_fac2mkt_fx[k, j])
    #
    # # Update bounds for (fac_x, fac_y)
    # fac_x_max = {}
    # fac_x_min = {}
    # fac_y_max = {}
    # fac_y_min = {}
    # for k in mip.fac:
    #     if w_fx[k] == 0:
    #         fac_x_max[k] = 0
    #         fac_x_min[k] = 0
    #         fac_y_max[k] = 0
    #         fac_y_min[k] = 0
    #     else:
    #         for p in mip.part:
    #             if w_p[k, p] == 1.0:
    #                 fac_x_max[k] = xp_mid[p] + length_x / 2
    #                 fac_x_min[k] = xp_mid[p] - length_x / 2
    #                 fac_y_max[k] = yp_mid[p] + length_y / 2
    #                 fac_y_min[k] = yp_mid[p] - length_y / 2
    # for k in mip.fac:
    #     minlp.x_max[k] = fac_x_max[k]
    #     minlp.x_min[k] = fac_x_min[k]
    #     minlp.y_max[k] = fac_y_max[k]
    #     minlp.y_min[k] = fac_y_min[k]
        # print(k, 'x_max: ', fac_x_max[k], 'x_min: ', fac_x_min[k])
    #
    # # Solve NLP for UB
    # nlpsolver = SolverFactory('gams')
    # nlpsolver.solve(minlp.Bl[t],
    #                 # tee=True,
    #                 add_options=['option reslim=300; option optcr = 0.0;'],
    #                 # keepfiles=True,
    #                 solver='baron',
    #                 load_solutions=True
    #                 )
    # nlp_sol[iter_] = minlp.Bl[t].obj()
    # print('solution of iteration', iter_, ': ', nlp_sol)
    # UB = min(nlp_sol[idx] for idx in iter_list if idx <= iter_)
    # opt_gap[iter_] = (UB - LB) / UB
    # print('LB: ', LB)
    # print('UB: ', UB)
    # print('opt gap: ', opt_gap[iter_])
    #
    # if opt_gap[iter_] <= bilevel_opt_tol or iter_ == max_iter:
    #     fac_x = {k: minlp.Bl[t].fac_x[k].value}
    #     fac_y = {k: minlp.Bl[t].fac_y[k].value}
    #     cost = UB
    #     break
    # else:
    #     p_x += n_x
    #     p_y += n_y







