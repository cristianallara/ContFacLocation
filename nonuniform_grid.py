from pyomo.environ import *
import copy


def refine_gride(mip, xp_min, xp_max, yp_min, yp_max, n_x, n_y, mapping, iter_, p_y):

    # add new indexes after repartition
    list_idx_part = []
    list_old_p = []
    for p in mip.part:
        for pp in list(range(1, 1 + n_y * n_x)):
            x_idx = ceil(pp / n_x)
            if pp <= n_y:
                y_idx = pp
            else:
                y_idx = pp - n_y * (ceil((pp - n_y) / n_y))
            new_idx = copy.deepcopy(mapping[p])
            new_idx.extend([x_idx, y_idx])
            list_idx_part.append(new_idx)
            list_old_p.append(p)
    # print(list_old_p)
    # print(list_idx_part)

    # order them by x_idx then by y_idx
    for i in range(len(list_idx_part)):
        for i_ in range(i):
            for j in range(len(list_idx_part[i])):
                if j % 2 == 0:
                    if list_idx_part[i_][j] < list_idx_part[i][j]:
                        break
                    elif list_idx_part[i_][j] > list_idx_part[i][j]:
                        list_idx_part[i_], list_idx_part[i] = list_idx_part[i], list_idx_part[i_]
                        list_old_p[i_], list_old_p[i] = list_old_p[i], list_old_p[i_]

    for i in range(len(list_idx_part)):
        if ceil((i+1)/ (p_y*iter_*n_y)) % 2 == 0:
            for i_ in range(i):
                idx_diff = 0
                for j in range(len(list_idx_part[i_])):
                    if j % 2 == 0 and list_idx_part[i_][j] != list_idx_part[i][j]:
                        idx_diff += 1
                        break
                if idx_diff == 0:
                    for l in range(len(list_idx_part[i_])):
                        if j % 2 == 1:
                            if list_idx_part[i_][l] < list_idx_part[i][l]:
                                break
                            elif list_idx_part[i_][l] > list_idx_part[i][l]:
                                list_idx_part[i_], list_idx_part[i] = list_idx_part[i], list_idx_part[i_]
                                list_old_p[i_], list_old_p[i] = list_old_p[i], list_old_p[i_]
                                break
    # print(list_old_p)
    # print(list_idx_part)

    # Calculate dimensions of regions on the grid and min distance to fix points
    dist_supp_new = {}
    dist_mkt_new = {}
    dx_suppl_new = {}
    dy_suppl_new = {}
    dx_mkt_new = {}
    dy_mkt_new = {}
    max_dist_supp_new = {}
    max_dist_mkt_new = {}
    xp_min_new = {}
    xp_max_new = {}
    yp_min_new = {}
    yp_max_new = {}
    xp_mid_new = {}
    yp_mid_new = {}
    for p_idx in range(len(list_idx_part)):
        p_old_idx = list_old_p[p_idx]
        # if region not refined
        length_x_new = (xp_max[p_old_idx] - xp_min[p_old_idx]) / n_x
        length_y_new = (yp_max[p_old_idx] - yp_min[p_old_idx]) / n_y
        x_idx = list_idx_part[p_idx][-2]
        y_idx = list_idx_part[p_idx][-1]
        xp_min_new[p_idx+1] = xp_min[p_old_idx] + length_x_new * (x_idx - 1)
        xp_max_new[p_idx+1] = xp_min[p_old_idx] + length_x_new * x_idx
        yp_min_new[p_idx+1] = yp_min[p_old_idx] + length_y_new * (y_idx - 1)
        yp_max_new[p_idx+1] = yp_min[p_old_idx] + length_y_new * y_idx
        xp_mid_new[p_idx+1] = xp_min_new[p_idx+1] + (length_x_new / 2)
        yp_mid_new[p_idx+1] = yp_min_new[p_idx+1] + (length_y_new / 2)
        for i in mip.suppl:
            dx_suppl_new[i, p_idx+1] = max(abs(mip.suppl_x[i] - xp_mid_new[p_idx+1]) - length_x_new / 2, 0)
            dy_suppl_new[i, p_idx+1] = max(abs(mip.suppl_y[i] - yp_mid_new[p_idx+1]) - length_y_new / 2, 0)
            dist_supp_new[i, p_idx+1] = sqrt(dx_suppl_new[i, p_idx+1] ** 2 + dy_suppl_new[i, p_idx+1] ** 2)
            max_dist_supp_new[i, p_idx+1] = dist_supp_new[i, p_idx+1] + sqrt(length_x_new ** 2 + length_y_new ** 2)
        for j in mip.mkt:
            dx_mkt_new[j, p_idx+1] = max(abs(mip.mkt_x[j] - xp_mid_new[p_idx+1]) - length_x_new / 2, 0)
            dy_mkt_new[j, p_idx+1] = max(abs(mip.mkt_y[j] - yp_mid_new[p_idx+1]) - length_y_new / 2, 0)
            dist_mkt_new[j, p_idx+1] = sqrt(dx_mkt_new[j, p_idx+1] ** 2 + dy_mkt_new[j, p_idx+1] ** 2)
            max_dist_mkt_new[j, p_idx+1] = dist_mkt_new[j, p_idx+1] + sqrt(length_x_new ** 2 + length_y_new ** 2)

    # update index mapping
    mapping_new = {}
    for p_idx in range(len(list_idx_part)):
        mapping_new[p_idx + 1] = list_idx_part[p_idx]
    # print(mapping_new)

    return dist_supp_new, dist_mkt_new, max_dist_supp_new, max_dist_mkt_new, xp_min_new, xp_max_new, yp_min_new, \
           yp_max_new, mapping_new, list_old_p
