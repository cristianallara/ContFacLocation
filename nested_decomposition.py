from input_data import read_data
from minlp_formulation import create_minlp
from forward import forward_pass
import time

start_time = time.time()

# User-defined parameters:
max_iter_nested = 1
opt_tol_nested = 0.01           # optimality tolerance for the nested decomposition
max_iter_bilevel = 1
opt_tol_bilevel = 0.005         # optimality tolerance for the nested decomposition
bilevel_single_period = False   # Boolean for the option of solving the bilevel inside each block of the forward pass
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

iters = range(1, max_iter_nested+1)

for iter_ in iters:

    # Forward Pass
    for t in minlp.t:
        print("Time period", t)

        w, cost = forward_pass(t, minlp, data, x_min, x_max, y_min, y_max, p_x, p_y, n_x, n_y, dist_min,
                                 bilevel_single_period, opt_tol_mip, time_limit_mip, max_iter_bilevel, opt_tol_bilevel)

        elapsed_time = time.time() - start_time
        print('elapsed time (s):', elapsed_time)