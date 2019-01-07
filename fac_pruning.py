from full_space_mip import create_multiperiod_mip
from pyomo.environ import *


def prune_facilities(mip, UB, data, n_part, b_part, pruned_fac, dist_supp, max_dist_supp, dist_mkt, max_dist_mkt,
                     dist_min):
    prune = {}

    for type_fac in [mip.centr, mip.distr]:
        for k in type_fac:
            prune[k] = False  # initialize as false
            if sum(b_part[k, p] for p in mip.part) == 0:
                # Generate MILP restricting facility k to be built in the complement set of partitions
                mip_comp = create_multiperiod_mip(data, n_part, b_part, k, True, None, False, pruned_fac, [], False)

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
                mipcompsolver.options['mipgap'] = 0.0001
                mipcompsolver.options['timelimit'] = 40
                mipcompsolver.options['relax_integrality'] = 1
                mipcompsolver.options['threads'] = 6
                results = mipcompsolver.solve(mip_comp)  # , tee=True)
                best_bound = results.problem.lower_bound
                LB_comp = best_bound
                print('LB^c: ', LB_comp)
                if LB_comp > UB:
                    prune[k] = True
                    if sum(b_part[k, p] for p in mip.part) == 0:
                        pruned_fac_type = [k_ for k_ in type_fac if type_fac.ord(k_) >= type_fac.ord(k)]
                        pruned_fac.extend(pruned_fac_type)
                        print('pruned facilities:', pruned_fac)
                        break
                else:
                    prune[k] = False
                print('Prune complementary partitions for facility ', k, '? ', prune[k])

    return pruned_fac
