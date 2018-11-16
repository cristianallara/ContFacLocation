from pyomo.environ import *
from input_data import read_data
from minlp_formulation import create_minlp
from mip_LB import create_mip
from grid_discretization import discretize_space
from forward import forward_pass
from backward import backward_pass
import time

start_time = time.time()

# User-defined parameters:
max_iter_nested = 1
opt_tol_nested = 0.01           # optimality tolerance for the nested decomposition
max_iter_bilevel = 1
opt_tol_bilevel = 0.005         # optimality tolerance for the nested decomposition
time_limit_mip = 40             # time limit in seconds for mip_LB
opt_tol_mip = 0.0001
dist_min = 0.5                  # arbitrary
p_x = 15                        # start grid with p_x X p_y partitions
p_y = 15
n_x = 5                         # adds p_x + n_x and p_y + n_y in each iteration of the bilevel decomposition
n_y = 5

# Input data
data = read_data()

# Create model
minlp = create_minlp(data)

# Tightest rectangle that includes all points
x_min = min(min(minlp.suppl_x[i] for i in minlp.suppl), min(minlp.mkt_x[j] for j in minlp.mkt))
x_max = max(max(minlp.suppl_x[i] for i in minlp.suppl), max(minlp.mkt_x[j] for j in minlp.mkt))
y_min = min(min(minlp.suppl_y[i] for i in minlp.suppl), min(minlp.mkt_y[j] for j in minlp.mkt))
y_max = max(max(minlp.suppl_y[i] for i in minlp.suppl), max(minlp.mkt_y[j] for j in minlp.mkt))

# create MIP relaxation for every time period
mip = create_mip(data, p_x, p_y)
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
w_prev_par_iter = {}
cost_fp_iter = {}
tot_cost_forward = {}
nested_UB = {}
mltp_iter = {}
cost_bp_iter = {}
nested_LB = {}
nested_gap = {}

for iter_ in mip.iters:

    # Forward Pass
    for t in minlp.t:
        print("Time period", t)

        w, cost_forward = forward_pass(mip.Bl[t], mip.Bl[t].w, opt_tol_mip, time_limit_mip)

        for k in mip.fac:
            for p in mip.part:
                mip.w_prev_par[k, p, t] = w[k, p]
                if w[k, p] != 0:
                    print(k, p, w[k, p])
                w_prev_par_iter[k, p, t, iter_] = w[k, p]
        cost_fp_iter[t, iter_] = cost_forward
        print('cost', cost_fp_iter[t, iter_])

    tot_cost_forward[iter_] = sum(cost_fp_iter[t, iter_] for t in mip.t)
    nested_UB = min(tot_cost_forward[kk] for kk in mip.iters if kk <= iter_)
    print('nested decomp UB: ', nested_UB)
    elapsed_time = time.time() - start_time
    print('elapsed time (s):', elapsed_time)

    # Backward Pass

    for t in reversed(list(mip.t)):
        print("Time period", t)

        mltp, cost_backward = backward_pass(mip.t.ord(t), mip.Bl[t], mip.Bl[t].equal, time_limit_mip)
        cost_bp_iter[t, iter_] = cost_backward
        print('cost', cost_bp_iter[t, iter_])
        print('alphafut', mip.Bl[t].alphafut.value)  # TODO: why alpha=0 in t1? Debug!!

        if mip.t.ord(t) != 1:
            t_prev = mip.t[mip.t.ord(t) - 1]
            for k in mip.fac:
                for p in mip.part:
                    mltp_iter[k, p, t, iter_] = mltp[k, p]

            # add Benders cut for cost to go function approximation
            mip.Bl[t_prev].fut_cost.add(expr=(mip.Bl[t_prev].alphafut >= cost_bp_iter[t, iter_]
                                              + sum(mltp_iter[k, p, t, iter_] *
                                                    (w_prev_par_iter[k, p, t_prev, iter_] - mip.Bl[t_prev].w[k, p])
                                                    for k in mip.fac for p in mip.part)))

    # Compute lower bound
    nested_LB[iter_] = cost_bp_iter[mip.t[1], iter_]
    print('nested decomp LB: ', nested_LB[iter_])
    # Compute optimality gap
    nested_gap[iter_] = (nested_UB - nested_LB[iter_]) / nested_UB
    print('nested opt gap: ', nested_gap[iter_])




