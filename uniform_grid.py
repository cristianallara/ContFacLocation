from pyomo.environ import *


def discretize_space(mip, x_min, x_max, y_min, y_max, p_x, p_y):

    # Define grid parameters
    length_x = (x_max - x_min) / p_x
    length_y = (y_max - y_min) / p_y

    xp_mid = {p: length_x * ((ceil(mip.part.ord(p) / p_x)) - 1) + length_x / 2 for p in mip.part}
    yp_mid = {}
    for p in mip.part:
        if p <= p_y:
            yp_mid[p] = length_y * (mip.part.ord(p) - 1) + length_y / 2
        else:
            yp_mid[p] = length_y * (mip.part.ord(p) - p_y * (ceil(mip.part.ord(p) / p_y) - 1) - 1) + length_x / 2

    xp_min = {p: length_x * ((ceil(mip.part.ord(p) / p_x)) - 1) for p in mip.part}
    yp_min = {}
    for p in mip.part:
        if p <= p_y:
            yp_min[p] = length_y * (mip.part.ord(p) - 1)
        else:
            yp_min[p] = length_y * (mip.part.ord(p) - p_y * (ceil(mip.part.ord(p) / p_y) - 1) - 1)
    xp_max = {p: length_x * (ceil(mip.part.ord(p) / p_x)) for p in mip.part}
    yp_max = {}
    for p in mip.part:
        if p <= p_y:
            yp_max[p] = length_y * (mip.part.ord(p))
        else:
            yp_max[p] = length_y * (mip.part.ord(p) - p_y * (ceil(mip.part.ord(p) / p_y) - 1))

    dx_suppl = {(i, p): max(abs(mip.suppl_x[i] - xp_mid[p]) - length_x / 2, 0) for i in mip.suppl for p in mip.part}
    dy_suppl = {(i, p): max(abs(mip.suppl_y[i] - yp_mid[p]) - length_y / 2, 0) for i in mip.suppl for p in mip.part}
    dx_mkt = {(j, p): max(abs(mip.mkt_x[j] - xp_mid[p]) - length_x / 2, 0) for j in mip.mkt for p in mip.part}
    dy_mkt = {(j, p): max(abs(mip.mkt_y[j] - yp_mid[p]) - length_y / 2, 0) for j in mip.mkt for p in mip.part}

    dist_supp = {}
    max_dist_supp = {}
    for i in mip.suppl:
        for p in mip.part:
            dist_supp[i, p] = sqrt(dx_suppl[i, p] ** 2 + dy_suppl[i, p] ** 2)
            max_dist_supp[i, p] = dist_supp[i, p] + sqrt(length_x**2 + length_y**2)
            # print('min:', i, p, dist_supp[i, p])
            # print('max:', i, p, max_dist_supp[i, p])
    dist_mkt = {}
    max_dist_mkt = {}
    for j in mip.mkt:
        for p in mip.part:
            dist_mkt[j, p] = sqrt(dx_mkt[j, p] ** 2 + dy_mkt[j, p] ** 2)
            max_dist_mkt[j, p] = dist_mkt[j, p] + sqrt(length_x**2 + length_y**2)
            # print('min:', j, p, dist_mkt[j, p])
            # print('max:', j, p, max_dist_mkt[j, p])

    # index mapping
    mapping = {}
    for p in mip.part:
        x_idx = ceil(p / p_x)
        if p <= p_y:
            y_idx = p
        else:
            y_idx = p - p_y * (ceil((p - p_y) / p_y))
        mapping[p] = [x_idx, y_idx]
    # print(mapping)

    return dist_supp, dist_mkt, max_dist_supp, max_dist_mkt, xp_min, xp_max, yp_min, yp_max, mapping
