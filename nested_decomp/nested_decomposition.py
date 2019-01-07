from pyomo.environ import *
from nested_decomp.mip_LB import create_mip
from nested_decomp.lagrangean_relaxation_mip import create_LR
from uniform_grid import discretize_space
from nested_decomp.forward import forward_pass
from nested_decomp.backward import backward_pass
import time


def run_nested_decomposition(data, x_min, x_max, y_min, y_max, p_x, p_y, full_mip, start_time):

    # create MIP relaxation for every time period
    mip = create_mip(data, p_x, p_y, 100)

    # # create Lagrangean relaxation for every time period
    # lr = create_LR(data, p_x, p_y)

    # define grid discretization
    dist_supp, dist_mkt, max_dist_supp, max_dist_mkt, xp_min, xp_max, yp_min, yp_max, mapping = \
        discretize_space(mip, x_min, x_max, y_min, y_max, p_x, p_y)

    # Retrieve duals for cut generation
    for t in mip.t:
        mip.Bl[t].dual = Suffix(direction=Suffix.IMPORT)

    for i in mip.suppl:
        for p in mip.part:
            mip.dist_supp[i, p] = dist_supp[i, p]
    for j in mip.mkt:
        for p in mip.part:
            mip.dist_mkt[j, p] = dist_mkt[j, p]

    w_par_iter = {}
    cost_fp_iter = {}
    tot_cost_forward = {}
    mltp_iter = {}
    cost_bp_iter = {}
    nested_LB = {}
    nested_gap = {}
    cost_LR_iter = {}
    cost_LP_iter = {}
    z_supp2fac_par = {}
    z_fac2mkt_par = {}

    for iter_ in mip.iters:
        print("iteration: ", iter_)

        # Forward Pass
        for t in full_mip.t:
            print(t)

            w, cost_forward = forward_pass(mip.Bl[t], mip.Bl[t].w, 0.00001, 300)

            for k in mip.fac:
                for p in mip.part:
                    mip.w_prev_par[k, p, t] = w[k, p]
                    # if strength_benders_cut:
                    #     lr.w_prev_par[k, p, t] = w[k, p]
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
                                                                  30, None, True, True, False)
            cost_bp_iter[t, iter_] = cost_backward

            if mip.t.ord(t) != 1:
                # if strength_benders_cut:
                #     cost_LR_iter[t, iter_] = cost_LR
                # if benders_cut:
                cost_LP_iter[t, iter_] = cost_LP

                t_prev = mip.t[mip.t.ord(t) - 1]
                for k in mip.fac:
                    for p in mip.part:
                        mltp_iter[k, p, t, iter_] = mltp[k, p]

                # add Benders cut for cost-to-go function approximation
                # if benders_cut:
                mip.Bl[t_prev].fut_cost.add(expr=(mip.Bl[t_prev].alphafut >= cost_LP_iter[t, iter_]
                                                  + sum(mltp_iter[k, p, t, iter_] *
                                                        (w_par_iter[k, p, t_prev, iter_] -
                                                         mip.Bl[t_prev].w[k, p])
                                                        for k in mip.fac for p in mip.part)))
                    # if strength_benders_cut:
                    #     lr.Bl[t_prev].fut_cost.add(expr=(lr.Bl[t_prev].alphafut >= cost_LP_iter[t, iter_]
                    #                                      + sum(mltp_iter[k, p, t, iter_] *
                    #                                            (w_par_iter[k, p, t_prev, iter_] -
                    #                                             lr.Bl[t_prev].w[k, p])
                    #                                            for k in lr.fac for p in lr.part)))

                # add integer optimality cut to cost-to-go function approximation
                # if integer_cut:
                mip.Bl[t_prev].fut_cost.add(expr=(mip.Bl[t_prev].alphafut >= cost_bp_iter[t, iter_]
                                                  * (sum((w_par_iter[k, p, t_prev, iter_] - 1) *
                                                         mip.Bl[t_prev].w[k, p] for k in mip.fac for p in mip.part)
                                                     + sum((mip.Bl[t_prev].w[k, p] - 1) *
                                                           w_par_iter[k, p, t_prev, iter_]
                                                           for k in mip.fac for p in mip.part))
                                                  + cost_bp_iter[t, iter_]))
                    # if strength_benders_cut:
                    #     lr.Bl[t_prev].fut_cost.add(expr=(lr.Bl[t_prev].alphafut >= cost_bp_iter[t, iter_]
                    #                                      * (sum((w_par_iter[k, p, t_prev, iter_] - 1) *
                    #                                             lr.Bl[t_prev].w[k, p] for k in lr.fac for p in lr.part)
                    #                                         + sum((lr.Bl[t_prev].w[k, p] - 1) *
                    #                                               w_par_iter[k, p, t_prev, iter_]
                    #                                               for k in lr.fac for p in lr.part))
                    #                                      + cost_bp_iter[t, iter_]))

                # add Strengthened Benders cut for cost-to-go function approximation
                # if strength_benders_cut:
                #     mip.Bl[t_prev].fut_cost.add(expr=(mip.Bl[t_prev].alphafut >= cost_LR_iter[t, iter_]
                #                                       + sum(mltp_iter[k, p, t, iter_] *
                #                                             (w_par_iter[k, p, t_prev, iter_] -
                #                                              mip.Bl[t_prev].w[k, p])
                #                                             for k in mip.fac for p in mip.part)))
                #     lr.Bl[t_prev].fut_cost.add(expr=(lr.Bl[t_prev].alphafut >= cost_LR_iter[t, iter_]
                #                                      + sum(mltp_iter[k, p, t, iter_] *
                #                                            (w_par_iter[k, p, t_prev, iter_] -
                #                                             lr.Bl[t_prev].w[k, p]) for k in lr.fac for p in lr.part)))

            print('cost', cost_bp_iter[t, iter_])
            print('alphafut', mip.Bl[t].alphafut.value)

        # Compute lower bound
        nested_LB[iter_] = cost_bp_iter[mip.t[1], iter_]
        print('nested decomp LB: ', nested_LB[iter_])
        # Compute optimality gap
        nested_gap[iter_] = (nested_UB - nested_LB[iter_]) / nested_UB
        print('nested opt gap: ', nested_gap[iter_])
        print('------------------------------------------------')

        if nested_gap[iter_] <= 0.0001:
            last_iter = iter_
            print("nested decomposition converged")
            break

    elapsed_time = time.time() - start_time
    print("Solution Time (s)", elapsed_time)







