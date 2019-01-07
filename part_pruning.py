from full_space_mip import create_multiperiod_mip
from pyomo.environ import *


def prune_partitions(mip, UB, data, active_part, n_part, b_part, pruned_fac, dist_supp, max_dist_supp, dist_mkt,
                     max_dist_mkt, dist_min):
    prune_part = {}
    for p in mip.part:
        print(p)
        print(p, active_part[p])
        prune_part[p] = False  # initialization
        if active_part[p] == 0:
            mip_comp = create_multiperiod_mip(data, n_part, b_part, None, False, p, True, pruned_fac, [], False)

            for i in mip_comp.suppl:
                for pp in mip_comp.part:
                    if dist_supp[i, pp] >= dist_min:
                        mip_comp.dist_supp[i, pp] = dist_supp[i, pp]
                    else:
                        mip_comp.dist_supp[i, pp] = dist_min
                    if max_dist_supp[i, pp] <= dist_min:
                        for k in mip_comp.fac:
                            for t in mip_comp.t:
                                mip_comp.w[k, pp, t].fix(0.0)
            for j in mip_comp.mkt:
                for pp in mip_comp.part:
                    if dist_mkt[j, pp] >= dist_min:
                        mip_comp.dist_mkt[j, pp] = dist_mkt[j, pp]
                    else:
                        mip_comp.dist_mkt[j, pp] = dist_min
                    if max_dist_mkt[j, pp] <= dist_min:
                        for k in mip_comp.fac:
                            for t in mip_comp.t:
                                mip_comp.w[k, pp, t].fix(0.0)

            for k in pruned_fac:
                for pp in mip.part:
                    for t in mip.t:
                        mip.w[k, pp, t].fix(0.0)
                        mip.b[k, pp, t].fix(0.0)
                        mip.f_fac[k, pp, t].fix(0.0)
                        for i in mip.suppl:
                            mip.f_supp2fac[i, k, pp, t].fix(0.0)
                            mip.z_supp2fac[i, k, pp, t].fix(0.0)
                        for j in mip.mkt:
                            mip.f_fac2mkt[k, j, pp, t].fix(0.0)
                            mip.z_fac2mkt[k, j, pp, t].fix(0.0)

            # Solve MILP with added constraint for LB^c
            mipcompsolver = SolverFactory('gurobi')
            # mipcompsolver.options['mipgap'] = opt_tol_mip
            mipcompsolver.options['timelimit'] = 10
            mipcompsolver.options['threads'] = 6
            results = mipcompsolver.solve(mip_comp)  # , tee=True)
            best_bound = results.problem.lower_bound
            LB_comp = best_bound
            print('LB^c: ', LB_comp)
            if LB_comp > UB:
                prune_part[p] = True
            else:
                prune_part[p] = False
            print('Prune partition ', p, '? ', prune_part[p])

    return prune_part