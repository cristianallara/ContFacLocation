from input_data import read_data
from minlp_formulation import create_minlp
from forward import forward_pass

# User-defined parameters:
max_iter = 1
opt_tol = 1  # optimality tolerance in %
dist_min = 0.5 # arbitrary
p_x = 2 # start grid with p_x X p_y partitions
p_y = 2

# Input data
data = read_data()

# Create model
minlp = create_minlp(data)

# Tightest rectangle that includes all points
x_min = min(min(minlp.suppl_x[i] for i in minlp.suppl), min(minlp.mkt_x[j] for j in minlp.mkt))
x_max = max(max(minlp.suppl_x[i] for i in minlp.suppl), max(minlp.mkt_x[j] for j in minlp.mkt))
y_min = min(min(minlp.suppl_y[i] for i in minlp.suppl), min(minlp.mkt_y[j] for j in minlp.mkt))
y_max = max(max(minlp.suppl_y[i] for i in minlp.suppl), max(minlp.mkt_y[j] for j in minlp.mkt))

iters = range(1, max_iter+1)

for iter_ in iters:

    # Forward Pass
    for t in minlp.t:
        print("Time period", t)

        # state_vars_minlp = [minlp.Bl[t].w, minlp.Bl[t].fac_x, minlp.Bl[t].fac_y]
        state_vars_mip = [minlp.Bl[t].w]

        [w], cost = forward_pass(t, minlp, data, state_vars_mip, x_min, x_max, y_min, y_max, p_x, p_y, dist_min,
                                 opt_tol=0.0, time_limit=40, max_iter=1)
