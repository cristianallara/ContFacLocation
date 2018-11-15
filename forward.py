__author__ = "Cristiana L. Lara"

from pyomo.environ import *
from bilevel_decomposition import bilevel_for_single_period


def forward_pass(t, minlp, data, state_vars, x_min, x_max, y_min, y_max,  p_x, p_y, dist_min, opt_tol, time_limit,
                 max_iter, bilevel_opt_tol):

    # Solve bilevel decomposition for single period
    cost = bilevel_for_single_period(t, minlp, data, x_min, x_max, y_min, y_max, p_x, p_y, dist_min, opt_tol, time_limit,
                                     max_iter, bilevel_opt_tol)

    state_vars_val = [{} for i in range(len(state_vars))]

    # Fix the state variable as parameter for next t
    for i in range(len(state_vars)):
        for key in state_vars[i].keys():
            state_vars_val[i][key] = state_vars[i][key].value

    # Store obj value to compute UB
    # cost = bl.obj() - bl.alphafut.value

    return state_vars_val, cost
