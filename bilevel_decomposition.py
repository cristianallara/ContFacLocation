from pyomo.environ import *
from mip_LB import create_mip
from grid_discretization import discretize_space


def bilevel_for_single_period(t, minlp, data, x_min, x_max, y_min, y_max,  p_x, p_y, dist_min, opt_tol, time_limit,
                              max_iter):

    for iter_ in range(1, max_iter+1):

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
        mipsolver.options['mipgap'] = opt_tol
        mipsolver.options['timelimit'] = time_limit
        mipsolver.options['threads'] = 6
        mipsolver.solve(mip.Bl[t], tee=True)
        LB = {iter_: mip.Bl[t].obj()}

        for k in mip.fac:
            for p in mip.part:
                print(k, p, mip.Bl[t].w[k, p].value)

        # Fix discrete solutions
        w_fx = {k: sum(mip.Bl[t].w[k, p].value for p in mip.part) for k in mip.fac}
        z_supp2fac_fx = {(i, k): sum(mip.Bl[t].z_supp2fac[i, k, p].value for p in mip.part) for i in mip.suppl
                         for k in mip.fac}
        z_fac2mkt_fx = {(k, j): sum(mip.Bl[t].z_fac2mkt[k, j, p].value for p in mip.part) for j in mip.mkt
                        for k in mip.fac}
        for k in mip.fac:
            print('w:', k, w_fx[k])
            for i in mip.suppl:
                print('z_supp2fac:', i, k, z_supp2fac_fx[i, k])
            for j in mip.mkt:
                print('z_fac2mkt:', k, j, z_fac2mkt_fx[k, j])

        # Update bounds for (fac_x, fac_y)
        fac_x_max = {}
        fac_x_min = {}
        fac_y_max = {}
        fac_y_min = {}
        for k in mip.fac:
            if w_fx[k] == 0:
                fac_x_max[k] = 0
                fac_x_min[k] = 0
                fac_y_max[k] = 0
                fac_y_min[k] = 0
            else:
                for p in mip.part:
                    if mip.Bl[t].w[k, p].value == 1:
                        fac_x_max[k] = xp_mid[p] + length_x/2
                        fac_x_min[k] = xp_mid[p] - length_x/2
                        fac_y_max[k] = yp_mid[p] + length_y/2
                        fac_y_min[k] = yp_mid[p] - length_y/2
        for k in mip.fac:
            print(k, fac_x_max[k])

        cost = 0

    return cost



