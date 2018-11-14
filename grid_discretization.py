from pyomo.environ import *


def discretize_space(mip, x_min, x_max, y_min, y_max, p_x, p_y, dist_min):
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
    for i in mip.suppl:
        for p in mip.part:
            if mip.suppl_x[i] < xp_min[p] or mip.suppl_x[i] > xp_max[p] or mip.suppl_y[i] < yp_min[p] or mip.suppl_y[
                i] > \
                    yp_max[p]:
                dist_supp[i, p] = sqrt(dx_suppl[i, p] ** 2 + dy_suppl[i, p] ** 2)
            else:
                dist_supp[i, p] = dist_min
    dist_mkt = {}
    for j in mip.mkt:
        for p in mip.part:
            if mip.mkt_x[j] < xp_min[p] or mip.mkt_x[j] > xp_max[p] or mip.mkt_y[j] < yp_min[p] or mip.mkt_y[j] > \
                    yp_max[p]:
                dist_mkt[j, p] = sqrt(dx_mkt[j, p] ** 2 + dy_mkt[j, p] ** 2)
            else:
                dist_mkt[j, p] = dist_min

    return dist_supp, dist_mkt, xp_mid, yp_mid, length_x, length_y