from pyomo.environ import *
from mip_LB import create_mip
from grid_discretization import discretize_space


def bilevel_for_single_period(t, minlp, data, x_min, x_max, y_min, y_max,  p_x, p_y, dist_min, opt_tol, time_limit,
                              max_iter, bilevel_opt_tol):
    mip_sol = {}
    nlp_sol = {}
    opt_gap = {}

    iter_list = range(1, max_iter+1)
    for iter_ in iter_list:

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
        mipsolver.solve(mip.Bl[t])#, tee=True)

        mip_sol[iter_] = mip.Bl[t].obj()
        LB = max(mip_sol[idx] for idx in iter_list if idx <= iter_)

        w_p = {}
        z_supp2fac_p = {}
        z_fac2mkt_fx = {}
        for p in mip.part:
            for k in mip.fac:
                w_p[k, p] = round(mip.Bl[t].w[k, p].value)
                for i in mip.suppl:
                    z_supp2fac_p[i, k, p] = round(mip.Bl[t].z_supp2fac[i, k, p].value)
                for j in mip.mkt:
                    z_fac2mkt_fx[k, j, p] = round(mip.Bl[t].z_fac2mkt[k, j, p].value)

        # Fix discrete solutions
        w_fx = {k: sum(w_p[k, p] for p in mip.part) for k in mip.fac}
        for k in minlp.fac:
            minlp.Bl[t].w[k].fix(w_fx[k])
            # print('w', k, w_fx[k])

        z_supp2fac_fx = {(i, k): sum(z_supp2fac_p[i, k, p] for p in mip.part) for i in mip.suppl
                         for k in mip.fac}
        for i in minlp.suppl:
            for k in minlp.fac:
                minlp.Bl[t].z_supp2fac[i, k].fix(z_supp2fac_fx[i, k])

        z_fac2mkt_fx = {(k, j): sum(z_fac2mkt_fx[k, j, p] for p in mip.part) for j in mip.mkt
                        for k in mip.fac}
        for k in minlp.fac:
            for j in minlp.mkt:
                minlp.Bl[t].z_fac2mkt[k, j].fix(z_fac2mkt_fx[k, j])

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
                    if w_p[k, p] == 1.0:
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

        # Solve MIP for LB
        nlpsolver = SolverFactory('gams')
        nlpsolver.solve(minlp.Bl[t],
                        # tee=True,
                        add_options=['option reslim=300; option optcr = 0.0;'],
                        # keepfiles=True,
                        solver='baron',
                        load_solutions=True
                        )
        nlp_sol[iter_] = minlp.Bl[t].obj()
        UB = min(nlp_sol[idx] for idx in iter_list if idx <= iter_)
        opt_gap[iter_] = (UB - LB) / UB
        print('LB: ', LB)
        print('UB: ', UB)
        print('opt gap: ', opt_gap[iter_])

        if opt_gap[iter_] <= bilevel_opt_tol or iter_ == max_iter:
            cost = UB
            break
        else:
            p_x += 10
            p_y += 10

    return cost



