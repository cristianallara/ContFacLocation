from pyomo.environ import *


def create_multiperiod_minlp(data, dist_min, pruned_fac):

    suppliers, xi, yi, time_periods, markets, xj, yj, centr_facilities, distr_facilities, facilities, cv, mc, a, d, \
    RM, FIC, VIC, FOC, VOC, ft1, ft2, vt1, vt2, interest_factor = data

    m = ConcreteModel()

    # Set declarations
    m.suppl = Set(initialize=suppliers, ordered=True)
    m.fac = Set(initialize=facilities, ordered=True)
    m.distr = Set(within=m.fac, initialize=distr_facilities, ordered=True)
    m.centr = Set(within=m.fac, initialize=centr_facilities, ordered=True)
    m.mkt = Set(initialize=markets, ordered=True)
    m.t = Set(initialize=time_periods, ordered=True)

    # Parameters
    m.suppl_x = Param(m.suppl, initialize=xi)
    m.suppl_y = Param(m.suppl, initialize=yi)
    m.mkt_x = Param(m.mkt, initialize=xj)
    m.mkt_y = Param(m.mkt, initialize=yj)
    m.conv = Param(m.fac, initialize=cv)
    m.cap = Param(m.fac, initialize=mc)
    m.avail = Param(m.suppl, m.t, initialize=a)
    m.dd = Param(m.mkt, m.t, initialize=d)
    m.RM = Param(m.suppl, m.t, initialize=RM)
    m.FIC = Param(m.fac, m.t, initialize=FIC)
    m.VIC = Param(m.fac, m.t, initialize=VIC)
    m.FOC = Param(m.fac, m.t, initialize=FOC)
    m.VOC = Param(m.fac, m.t, initialize=VOC)
    m.FTC_supp2fac = Param(m.suppl, m.fac, m.t, initialize=ft1)
    m.VTC_supp2fac = Param(m.suppl, m.fac, m.t, initialize=vt1)
    m.FTC_fac2mkt = Param(m.fac, m.mkt, m.t, initialize=ft2)
    m.VTC_fac2mkt = Param(m.fac, m.mkt, m.t, initialize=vt2)
    m.i_factor = Param(m.t, initialize=interest_factor)

    # Bounds
    x_min = min(min(m.suppl_x[i] for i in m.suppl), min(m.mkt_x[j] for j in m.mkt))
    x_max = max(max(m.suppl_x[i] for i in m.suppl), max(m.mkt_x[j] for j in m.mkt))
    y_min = min(min(m.suppl_y[i] for i in m.suppl), min(m.mkt_y[j] for j in m.mkt))
    y_max = max(max(m.suppl_y[i] for i in m.suppl), max(m.mkt_y[j] for j in m.mkt))
    m.x_min = Param(m.fac, initialize={k: x_min for k in m.fac}, mutable=True)
    m.x_max = Param(m.fac, initialize={k: x_max for k in m.fac}, mutable=True)
    m.y_min = Param(m.fac, initialize={k: y_min for k in m.fac}, mutable=True)
    m.y_max = Param(m.fac, initialize={k: y_max for k in m.fac}, mutable=True)

    dist_max = sqrt((x_max - x_min)**2 + (y_max-y_min)**2)

    # Variables declaration
    m.f_supp2fac = Var(m.suppl, m.fac, m.t, domain=NonNegativeReals)
    m.f_fac2mkt = Var(m.fac, m.mkt, m.t, domain=NonNegativeReals)
    m.f_fac = Var(m.fac, m.t, domain=NonNegativeReals)
    m.fac_x = Var(m.fac, domain=NonNegativeReals, bounds=(x_min, x_max))
    m.fac_y = Var(m.fac, domain=NonNegativeReals, bounds=(y_min, y_max))
    m.dist_supp2fac = Var(m.suppl, m.fac, domain=NonNegativeReals, bounds=(dist_min, dist_max))
    m.dist_fac2mkt = Var(m.fac, m.mkt, domain=NonNegativeReals, bounds=(dist_min, dist_max))
    m.cost_supp2fac = Var(m.t, domain=NonNegativeReals)
    m.cost_fac2mkt = Var(m.t, domain=NonNegativeReals)
    m.inv_cost = Var(m.t, domain=NonNegativeReals)
    m.op_cost = Var(m.t, domain=NonNegativeReals)
    m.w = Var(m.fac, m.t, domain=Binary) # binary for when facility is in place and operating
    m.w_prev = Var(m.fac, m.t, domain=Binary)
    m.b = Var(m.fac, m.t, domain=Binary) # binary for building facility k in time period t
    m.z_supp2fac = Var(m.suppl, m.fac, m.t, domain=Binary)
    m.z_fac2mkt = Var(m.fac, m.mkt, m.t, domain=Binary)

    # Objective function
    def obj_rule(m):
        return sum(m.i_factor[t] * (m.inv_cost[t] + m.op_cost[t] + m.cost_supp2fac[t] + m.cost_fac2mkt[t]) for t in m.t)
    m.obj = Objective(rule=obj_rule, sense=minimize)

    # Constraints
    def investment_cost(m, t):
        return m.inv_cost[t] == sum((m.FIC[k, t] + m.VIC[k, t] * m.cap[k]) * m.b[k, t] for k in m.fac
                                    if k not in pruned_fac)
    m.investment_cost = Constraint(m.t, rule=investment_cost)

    def operating_cost(m, t):
        return m.op_cost[t] == sum(m.FOC[k, t]*m.w[k, t] + m.VOC[k, t] * m.f_fac[k, t] for k in m.fac
                                   if k not in pruned_fac)
    m.operating_cost = Constraint(m.t, rule=operating_cost)

    def transp_cost_supp2fac(m, t):
        return m.cost_supp2fac[t] == sum(m.RM[i, t]*m.f_supp2fac[i, k, t] + m.FTC_supp2fac[i, k, t]*m.z_supp2fac[i, k, t]
                                       + m.VTC_supp2fac[i, k, t]*m.dist_supp2fac[i, k]*m.f_supp2fac[i, k, t]
                                       for i in m.suppl for k in m.fac if k not in pruned_fac)
    m.transp_cost_supp2fac = Constraint(m.t, rule=transp_cost_supp2fac)

    def transp_cost_fac2mkt(m, t):
        return m.cost_fac2mkt[t] == sum(m.FTC_fac2mkt[k, j, t] * m.z_fac2mkt[k, j, t] + m.VTC_fac2mkt[k, j, t]
                                      * m.dist_fac2mkt[k, j] * m.f_fac2mkt[k, j, t] for k in m.fac if k not in pruned_fac
                                        for j in m.mkt)
    m.transp_cost_fac2mkt = Constraint(m.t, rule=transp_cost_fac2mkt)

    def massmal1(m, k, t):
        if k not in pruned_fac:
            return m.f_fac[k, t] == sum(m.f_supp2fac[i, k, t] * m.conv[k] for i in m.suppl)
        return Constraint.Skip
    m.massmal1 = Constraint(m.fac, m.t, rule=massmal1)

    def massmal2(m, k, t):
        if k not in pruned_fac:
            return m.f_fac[k, t] == sum(m.f_fac2mkt[k, j, t] for j in m.mkt)
        return Constraint.Skip
    m.massmal2 = Constraint(m.fac, m.t, rule=massmal2)

    def demand(m, j, t):
        return sum(m.f_fac2mkt[k, j, t] for k in m.fac if k not in pruned_fac) == m.dd[j, t]
    m.demand = Constraint(m.mkt, m.t, rule=demand)

    def availability(m, i, t):
        return sum(m.f_supp2fac[i, k, t] for k in m.fac if k not in pruned_fac) <= m.avail[i, t]
    m.availability = Constraint(m.suppl, m.t, rule=availability)

    def dist_supp2fac_rule(m, i, k):
        if k not in pruned_fac:
            return m.dist_supp2fac[i, k] >= sqrt((m.suppl_x[i] - m.fac_x[k])**2 + (m.suppl_y[i] - m.fac_y[k])**2)
        return Constraint.Skip
    m.dist_supp2fac_rule = Constraint(m.suppl, m.fac, rule=dist_supp2fac_rule)

    def dist_fac2mkt_rule(m, k, j):
        if k not in pruned_fac:
            return m.dist_fac2mkt[k, j] >= sqrt((m.mkt_x[j] - m.fac_x[k])**2 + (m.mkt_y[j] - m.fac_y[k])**2)
        return Constraint.Skip
    m.dist_fac2mkt_rule = Constraint(m.fac, m.mkt, rule=dist_fac2mkt_rule)

    def x_bigM_min(m, k):
        if k not in pruned_fac:
            return m.fac_x[k] >= m.x_min[k] * sum(m.b[k, t] for t in m.t)
        return Constraint.Skip
    m.x_bigM_min = Constraint(m.fac, rule=x_bigM_min)

    def x_bigM_max(m, k):
        if k not in pruned_fac:
            return m.fac_x[k] <= m.x_max[k] * sum(m.b[k, t] for t in m.t)
        return Constraint.Skip
    m.x_bigM_max = Constraint(m.fac, rule=x_bigM_max)

    def y_bigM_min(m, k):
        if k not in pruned_fac:
            return m.fac_y[k] >= m.y_min[k] * sum(m.b[k, t] for t in m.t)
        return Constraint.Skip
    m.y_bigM_min = Constraint(m.fac, rule=y_bigM_min)

    def y_bigM_max(m, k):
        if k not in pruned_fac:
            return m.fac_y[k] <= m.y_max[k] * sum(m.b[k, t] for t in m.t)
        return Constraint.Skip
    m.y_bigM_max = Constraint(m.fac, rule=y_bigM_max)

    def supp2facmigM(m, i, k, t):
        if k not in pruned_fac:
            return m.f_supp2fac[i, k, t] <= m.avail[i, t]*m.z_supp2fac[i, k, t]
        return Constraint.Skip
    m.supp2facmigM = Constraint(m.suppl, m.fac, m.t, rule=supp2facmigM)

    def fac2mktmigM(m, k, j, t):
        if k not in pruned_fac:
            return m.f_fac2mkt[k, j, t] <= m.dd[j, t]*m.z_fac2mkt[k, j, t]
        return Constraint.Skip
    m.fac2mktmigM = Constraint(m.fac, m.mkt, m.t, rule=fac2mktmigM)

    def facmigM(m, k, t):
        if k not in pruned_fac:
            return m.f_fac[k, t] <= m.cap[k]*m.w[k, t]
        return Constraint.Skip
    m.facmigM = Constraint(m.fac, m.t, rule=facmigM)

    def logic_1(m, i, k, t):
        if k not in pruned_fac:
            return m.w[k, t] >= m.z_supp2fac[i, k, t]
        return Constraint.Skip
    m.logic_1 = Constraint(m.suppl, m.fac, m.t, rule=logic_1)

    def logic_2(m, k, t):
        if k not in pruned_fac:
            return sum(m.z_supp2fac[i, k, t] for i in m.suppl) >= m.w[k, t]
        return Constraint.Skip
    m.logic_2 = Constraint(m.fac, m.t, rule=logic_2)

    def logic_3(m, k, j, t):
        if k not in pruned_fac:
            return m.w[k, t] >= m.z_fac2mkt[k, j, t]
        return Constraint.Skip
    m.logic_3 = Constraint(m.fac, m.mkt, m.t, rule=logic_3)

    def logic_4(m, k, t):
        if k not in pruned_fac:
            return sum(m.z_fac2mkt[k, j, t] for j in m.mkt) >= m.w[k, t]
        return Constraint.Skip
    m.logic_4 = Constraint(m.fac, m.t, rule=logic_4)

    def logic_5(m, k, t):
        if k not in pruned_fac:
            if m.t.ord(t) == 1:
                return m.w[k, t] == m.b[k, t]
            else:
                return m.w[k, t] == m.w[k, m.t[m.t.ord(t) - 1]] + m.b[k, t]
        return Constraint.Skip
    m.logic_5 = Constraint(m.fac, m.t, rule=logic_5)

    def logic_6(m, k):
        if k not in pruned_fac:
            return sum(m.b[k, t] for t in m.t) <= 1
        return Constraint.Skip
    m.logic_6 = Constraint(m.fac, rule=logic_6)

    def sym_1(m, l, u):
        if l not in pruned_fac and u not in pruned_fac:
            if m.distr.ord(l) < m.distr.ord(u):
                return m.fac_x[l] >= m.fac_x[u]
            return Constraint.Skip
        return Constraint.Skip
    m.sym_1 = Constraint(m.distr, m.distr, rule=sym_1)

    def sym_2(m, n, v):
        if n not in pruned_fac and v not in pruned_fac:
            if m.centr.ord(n) < m.centr.ord(v):
                return m.fac_x[n] >= m.fac_x[v]
            return Constraint.Skip
        return Constraint.Skip
    m.sym_2 = Constraint(m.centr, m.centr, rule=sym_2)

    def sym_3(m, l, u, t):
        if l not in pruned_fac and u not in pruned_fac:
            if m.distr.ord(l) < m.distr.ord(u):
                return m.w[l, t] >= m.w[u, t]
            return Constraint.Skip
        return Constraint.Skip
    m.sym_3 = Constraint(m.distr, m.distr, m.t, rule=sym_3)

    def sym_4(m, n, v, t):
        if n not in pruned_fac and v not in pruned_fac:
            if m.centr.ord(n) < m.centr.ord(v):
                return m.w[n, t] >= m.w[v, t]
            return Constraint.Skip
        return Constraint.Skip
    m.sym_4 = Constraint(m.centr, m.centr, m.t, rule=sym_4)

    return m