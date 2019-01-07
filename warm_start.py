

def warm_start_MIP(mip, minlp, xp_min, xp_max, yp_min, yp_max, w_fx, b_fx, z_supp2fac_fx, z_fac2mkt_fx, p_y, n_y,
                   pruned_fac, n_part, iter_):

    # warm-start
    w_p_new = {}
    b_p_new = {}
    z_supp2fac_p_new = {}
    z_fac2mkt_p_new = {}
    for k in mip.fac:
        if k not in pruned_fac:
            for t in mip.t:
                # initialize w
                if w_fx[k, t] == 1:
                    for p in range(1, n_part + 1):
                        if xp_min[p] < minlp.fac_x[k].value < xp_max[p] and yp_min[p] < minlp.fac_y[k].value < yp_max[p]:
                            w_p_new[k, p, t] = 1
                        elif xp_min[p] <= minlp.fac_x[k].value <= xp_max[p] \
                                and yp_min[p] <= minlp.fac_y[k].value <= yp_max[p]:
                            if minlp.fac_x[k].value == xp_min[p] and yp_min[p] < minlp.fac_y[k].value < yp_max[p]:
                                p_adj = p - (p_y * iter_ * n_y)
                            elif minlp.fac_x[k].value == xp_max[p] and yp_min[p] < minlp.fac_y[k].value < yp_max[p]:
                                p_adj = p + (p_y * iter_ * n_y)
                            elif minlp.fac_y[k].value == yp_min[p] and xp_min[p] < minlp.fac_x[k].value < xp_max[p]:
                                p_adj = p - 1
                            elif minlp.fac_y[k].value == yp_max[p] and xp_min[p] < minlp.fac_x[k].value < xp_max[p]:
                                p_adj = p + 1
                            elif minlp.fac_y[k].value == yp_min[p] and minlp.fac_x[k].value == xp_min[p]:
                                p_adj = p - (p_y * iter_ * n_y) - 1
                            elif minlp.fac_y[k].value == yp_max[p] and minlp.fac_x[k].value == xp_min[p]:
                                p_adj = p - (p_y * iter_ * n_y) + 1
                            elif minlp.fac_y[k].value == yp_min[p] and minlp.fac_x[k].value == xp_max[p]:
                                p_adj = p + (p_y * iter_ * n_y) - 1
                            else:
                                p_adj = p + (p_y * iter_ * n_y) + 1
                            w_p_new[k, p_adj, t] = 1
                            w_p_new[k, p, t] = 0
                        else:
                            w_p_new[k, p, t] = 0
                else:
                    for p in range(1, n_part + 1):
                        w_p_new[k, p, t] = 0

                # initialize b
                if b_fx[k, t] == 1:
                    for p in range(1, n_part + 1):
                        if xp_min[p] < minlp.fac_x[k].value < xp_max[p] and yp_min[p] < minlp.fac_y[k].value < \
                                yp_max[p]:
                            b_p_new[k, p, t] = 1
                        elif xp_min[p] <= minlp.fac_x[k].value <= xp_max[p] \
                                and yp_min[p] <= minlp.fac_y[k].value <= yp_max[p]:
                            if minlp.fac_x[k].value == xp_min[p] and yp_min[p] < minlp.fac_y[k].value < yp_max[p]:
                                p_adj = p - (p_y * iter_*  n_y)
                            elif minlp.fac_x[k].value == xp_max[p] and yp_min[p] < minlp.fac_y[k].value < yp_max[p]:
                                p_adj = p + (p_y * iter_*  n_y)
                            elif minlp.fac_y[k].value == yp_min[p] and xp_min[p] < minlp.fac_x[k].value < xp_max[p]:
                                p_adj = p - 1
                            elif minlp.fac_y[k].value == yp_max[p] and xp_min[p] < minlp.fac_x[k].value < xp_max[p]:
                                p_adj = p + 1
                            elif minlp.fac_y[k].value == yp_min[p] and minlp.fac_x[k].value == xp_min[p]:
                                p_adj = p - (p_y * iter_*  n_y) - 1
                            elif minlp.fac_y[k].value == yp_max[p] and minlp.fac_x[k].value == xp_min[p]:
                                p_adj = p - (p_y * iter_*  n_y) + 1
                            elif minlp.fac_y[k].value == yp_min[p] and minlp.fac_x[k].value == xp_max[p]:
                                p_adj = p + (p_y * iter_*  n_y) - 1
                            else:
                                p_adj = p + (p_y * iter_*  n_y) + 1
                            b_p_new[k, p_adj, t] = 1
                            b_p_new[k, p, t] = 0
                        else:
                            b_p_new[k, p, t] = 0
                else:
                    for p in range(1, n_part + 1):
                        b_p_new[k, p, t] = 0

                # initialize z_supp2fac
                for i in mip.suppl:
                    if z_supp2fac_fx[i, k, t] == 1:
                        for p in range(1, n_part + 1):
                            if xp_min[p] < minlp.fac_x[k].value < xp_max[p] and yp_min[p] < minlp.fac_y[k].value < \
                                    yp_max[p]:
                                z_supp2fac_p_new[i, k, p, t] = 1
                            elif xp_min[p] <= minlp.fac_x[k].value <= xp_max[p] \
                                 and yp_min[p] <= minlp.fac_y[k].value <= yp_max[p]:
                                if minlp.fac_x[k].value == xp_min[p] and yp_min[p] < minlp.fac_y[k].value < yp_max[p]:
                                    p_adj = p - (p_y * iter_*  n_y)
                                elif minlp.fac_x[k].value == xp_max[p] and yp_min[p] < minlp.fac_y[k].value < yp_max[p]:
                                    p_adj = p + (p_y * iter_*  n_y)
                                elif minlp.fac_y[k].value == yp_min[p] and xp_min[p] < minlp.fac_x[k].value < xp_max[p]:
                                    p_adj = p - 1
                                elif minlp.fac_y[k].value == yp_max[p] and xp_min[p] < minlp.fac_x[k].value < xp_max[p]:
                                    p_adj = p + 1
                                elif minlp.fac_y[k].value == yp_min[p] and minlp.fac_x[k].value == xp_min[p]:
                                    p_adj = p - (p_y * iter_*  n_y) - 1
                                elif minlp.fac_y[k].value == yp_max[p] and minlp.fac_x[k].value == xp_min[p]:
                                    p_adj = p - (p_y * iter_*  n_y) + 1
                                elif minlp.fac_y[k].value == yp_min[p] and minlp.fac_x[k].value == xp_max[p]:
                                    p_adj = p + (p_y * iter_*  n_y) - 1
                                else:
                                    p_adj = p + (p_y * iter_*  n_y) + 1
                                z_supp2fac_p_new[i, k, p_adj, t] = 1
                                z_supp2fac_p_new[i, k, p, t] = 0
                            else:
                                z_supp2fac_p_new[i, k, p, t] = 0
                    else:
                        z_supp2fac_p_new[i, k, p, t] = 0

                # initialize z_fac2mkt
                for j in mip.mkt:
                    if z_fac2mkt_fx[k, j, t] == 1:
                        for p in range(1, n_part + 1):
                            if xp_min[p] < minlp.fac_x[k].value < xp_max[p] and yp_min[p] < minlp.fac_y[
                                k].value < \
                                    yp_max[p]:
                                z_fac2mkt_p_new[k, j, p, t] = 1
                            elif xp_min[p] <= minlp.fac_x[k].value <= xp_max[p] \
                                 and yp_min[p] <= minlp.fac_y[k].value <= yp_max[p]:
                                if minlp.fac_x[k].value == xp_min[p] and yp_min[p] < minlp.fac_y[k].value < \
                                        yp_max[p]:
                                    p_adj = p - (p_y * iter_*  n_y)
                                elif minlp.fac_x[k].value == xp_max[p] and yp_min[p] < minlp.fac_y[k].value < \
                                        yp_max[p]:
                                    p_adj = p + (p_y * iter_*  n_y)
                                elif minlp.fac_y[k].value == yp_min[p] and xp_min[p] < minlp.fac_x[k].value < \
                                        xp_max[p]:
                                    p_adj = p - 1
                                elif minlp.fac_y[k].value == yp_max[p] and xp_min[p] < minlp.fac_x[k].value < \
                                        xp_max[p]:
                                    p_adj = p + 1
                                elif minlp.fac_y[k].value == yp_min[p] and minlp.fac_x[k].value == xp_min[p]:
                                    p_adj = p - (p_y * iter_*  n_y) - 1
                                elif minlp.fac_y[k].value == yp_max[p] and minlp.fac_x[k].value == xp_min[p]:
                                    p_adj = p - (p_y * iter_*  n_y) + 1
                                elif minlp.fac_y[k].value == yp_min[p] and minlp.fac_x[k].value == xp_max[p]:
                                    p_adj = p + (p_y * iter_*  n_y) - 1
                                else:
                                    p_adj = p + (p_y * iter_*  n_y) + 1
                                z_fac2mkt_p_new[k, j, p_adj, t] = 1
                                z_fac2mkt_p_new[k, j, p, t] = 0
                            else:
                                z_fac2mkt_p_new[k, j, p, t] = 0
                    else:
                        for p in range(1, n_part + 1):
                            z_fac2mkt_p_new[k, j, p, t] = 0

    print('b_p_new')
    for k in mip.fac:
        if k not in pruned_fac:
            for p in range(1, n_part + 1):
                for t in mip.t:
                    if b_p_new[k, p, t] == 1:
                        print(k, p, t, b_p_new[k, p, t])
    print('w_p_new')
    for k in mip.fac:
        if k not in pruned_fac:
            for p in range(1, n_part + 1):
                for t in mip.t:
                    if w_p_new[k, p, t] == 1:
                        print(k, p, t, w_p_new[k, p, t])

    warm_start = [w_p_new, b_p_new, z_supp2fac_p_new, z_fac2mkt_p_new]
    return warm_start
