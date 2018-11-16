__author__ = "Cristiana L. Lara"

from pyomo.environ import *
from bilevel_decomposition import bilevel_for_single_period
from mip_LB import create_mip
from grid_discretization import discretize_space


def forward_pass(t, minlp, data, x_min, x_max, y_min, y_max, p_x, p_y, n_x, n_y, dist_min, bilevel_single_period,
                 opt_tol_mip, time_limit_mip, max_iter_bilevel, opt_tol_bilevel):

    if bilevel_single_period:
        # TODO: check if this works
        # Solve bilevel decomposition for single period
        cost_nlp, w_fx, fac_x, fac_y = bilevel_for_single_period(t, minlp, data, x_min, x_max, y_min, y_max, p_x, p_y,
                                                                 n_x, n_y, dist_min, opt_tol_mip, time_limit_mip,
                                                                 max_iter_bilevel, opt_tol_bilevel)

        # Fix the state variable as parameter for next t
        state_vars_val = [w_fx, fac_x, fac_y]

        # Store obj value to compute UB
        cost = cost_nlp # TODO: add the - bl.alphafut.value term

    else:
        # Solve MIP relaxation for single period
        mip = create_mip(data, p_x, p_y)
        dist_supp, dist_mkt, xp_mid, yp_mid, length_x, length_y = discretize_space(mip, x_min, x_max, y_min,
                                                                                   y_max, p_x, p_y, dist_min)

        for i in mip.suppl:
            for p in mip.part:
                mip.dist_supp[i, p] = dist_supp[i, p]
        for j in mip.mkt:
            for p in mip.part:
                mip.dist_mkt[j, p] = dist_mkt[j, p]

        # Solve MIP for LB
        mipsolver = SolverFactory('gurobi')
        mipsolver.options['mipgap'] = opt_tol_mip
        mipsolver.options['timelimit'] = time_limit_mip
        mipsolver.options['threads'] = 6
        mipsolver.solve(mip.Bl[t]) # , tee=True)

        # Fix the state variable as parameter for next t
        w_fx = {}
        for k in mip.fac:
            for p in mip.part:
                w_fx[k, p] = mip.Bl[t].w[k, p].value
                print(k, p, w_fx[k, p])

    # Store obj value to compute UB
    cost = mip.Bl[t].obj() - mip.Bl[t].alphafut.value
    print(t, 'cost', cost)

    return w_fx, cost
