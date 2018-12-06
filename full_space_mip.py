from pyomo.environ import *


def create_multiperiod_mip(data, n_part, select_part, if_complement):

    suppliers, xi, yi, time_periods, markets, xj, yj, centr_facilities, distr_facilities, facilities, cv, mc, a, d, \
    RM, FIC, VIC, FOC, VOC, ft1, ft2, vt1, vt2, interest_factor = data

    m = ConcreteModel()

    # Set declarations
    m.part = Set(initialize=RangeSet(n_part), ordered=True)
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
    m.dist_supp = Param(m.suppl, m.part, default=0, mutable=True)
    m.dist_mkt = Param(m.mkt, m.part, default=0, mutable=True)

    # Variables declaration
    m.f_supp2fac = Var(m.suppl, m.fac, m.part, m.t, within=NonNegativeReals)
    m.f_fac2mkt = Var(m.fac, m.mkt, m.part, m.t, within=NonNegativeReals)
    m.f_fac = Var(m.fac, m.part, m.t, within=NonNegativeReals)
    m.cost_supp2fac = Var(m.t, within=NonNegativeReals)
    m.cost_fac2mkt = Var(m.t, within=NonNegativeReals)
    m.inv_cost = Var(m.t, within=NonNegativeReals)
    m.op_cost = Var(m.t, within=NonNegativeReals)
    m.w = Var(m.fac, m.part, m.t, within=Binary)  # binary for when facility is in place and operating
    m.b = Var(m.fac, m.part, m.t, within=Binary)  # binary for building facility k in time period t
    m.z_supp2fac = Var(m.suppl, m.fac, m.part, m.t, within=Binary)
    m.z_fac2mkt = Var(m.fac, m.mkt, m.part, m.t, within=Binary)

    # Objective function
    def obj_rule(m):
        return sum(m.i_factor[t] * (m.inv_cost[t] + m.op_cost[t] + m.cost_supp2fac[t] + m.cost_fac2mkt[t]) for t in m.t)
    m.obj = Objective(rule=obj_rule, sense=minimize)

    # Constraints
    def investment_cost(m, t):
        return m.inv_cost[t] == sum((m.FIC[k, t] + m.VIC[k, t]*m.cap[k])*m.b[k, p, t] for k in m.fac for p in m.part)
    m.investment_cost = Constraint(m.t, rule=investment_cost)

    def operating_cost(m, t):
        return m.op_cost[t] == sum(m.FOC[k, t]*m.w[k, p, t] + m.VOC[k, t] * m.f_fac[k, p, t] for k in m.fac
                                   for p in m.part)
    m.operating_cost = Constraint(m.t, rule=operating_cost)

    def transp_cost_supp2fac(m, t):
        return m.cost_supp2fac[t] == sum(m.RM[i, t]*m.f_supp2fac[i, k, p, t] + m.FTC_supp2fac[i, k, t]
                                         * m.z_supp2fac[i, k, p, t] + m.VTC_supp2fac[i, k, t]*m.dist_supp[i, p]
                                         * m.f_supp2fac[i, k, p, t] for i in m.suppl for k in m.fac for p in m.part)
    m.transp_cost_supp2fac = Constraint(m.t, rule=transp_cost_supp2fac)

    def transp_cost_fac2mkt(m, t):
        return m.cost_fac2mkt[t] == sum(m.FTC_fac2mkt[k, j, t] * m.z_fac2mkt[k, j, p, t] + m.VTC_fac2mkt[k, j, t]
                                        * m.dist_mkt[j, p] * m.f_fac2mkt[k, j, p, t] for k in m.fac
                                        for j in m.mkt for p in m.part)
    m.transp_cost_fac2mkt = Constraint(m.t, rule=transp_cost_fac2mkt)

    def mass_bal1(m, k, p, t):
        return m.f_fac[k, p, t] == sum(m.f_supp2fac[i, k, p, t] * m.conv[k] for i in m.suppl)
    m.mass_bal1 = Constraint(m.fac, m.part, m.t, rule=mass_bal1)

    def mass_bal2(m, k, p, t):
        return m.f_fac[k, p, t] == sum(m.f_fac2mkt[k, j, p, t] for j in m.mkt)
    m.mass_bal2 = Constraint(m.fac, m.part, m.t, rule=mass_bal2)

    def demand(m, j, t):
        return sum(m.f_fac2mkt[k, j, p, t] for k in m.fac for p in m.part) == m.dd[j, t]
    m.demand = Constraint(m.mkt, m.t, rule=demand)

    def availability(m, i, t):
        return sum(m.f_supp2fac[i, k, p, t] for k in m.fac for p in m.part) <= m.avail[i, t]
    m.availability = Constraint(m.suppl, m.t, rule=availability)

    def supp2fac_bigM(m, i, k, p, t):
        return m.f_supp2fac[i, k, p, t] <= m.avail[i, t]*m.z_supp2fac[i, k, p, t]
    m.supp2fac_bigM = Constraint(m.suppl, m.fac, m.part, m.t, rule=supp2fac_bigM)

    def fac2mkt_bigM(m, k, j, p, t):
        return m.f_fac2mkt[k, j, p, t] <= m.dd[j, t]*m.z_fac2mkt[k, j, p, t]
    m.fac2mkt_bigM = Constraint(m.fac, m.mkt, m.part, m.t, rule=fac2mkt_bigM)

    def fac_bigM(m, k, p, t):
        return m.f_fac[k, p, t] <= m.cap[k]*m.w[k, p, t]
    m.fac_bigM = Constraint(m.fac, m.part, m.t, rule=fac_bigM)

    def logic_1(m, i, k, p, t):
        return m.w[k, p, t] >= m.z_supp2fac[i, k, p, t]
    m.logic_1 = Constraint(m.suppl, m.fac, m.part, m.t, rule=logic_1)

    def logic_2(m, k, p, t):
        return sum(m.z_supp2fac[i, k, p, t] for i in m.suppl) >= m.w[k, p, t]
    m.logic_2 = Constraint(m.fac, m.part, m.t, rule=logic_2)

    def logic_3(m, k, j, p, t):
        return m.w[k, p, t] >= m.z_fac2mkt[k, j, p, t]
    m.logic_3 = Constraint(m.fac, m.mkt, m.part, m.t, rule=logic_3)

    def logic_4(m, k, p, t):
        return sum(m.z_fac2mkt[k, j, p, t] for j in m.mkt) >= m.w[k, p, t]
    m.logic_4 = Constraint(m.fac, m.part, m.t, rule=logic_4)

    def logic_5(m, k, p, t):
        if m.t.ord(t) == 1:
            return m.w[k, p, t] == m.b[k, p, t]
        else:
            return m.w[k, p, t] == m.w[k, p, m.t[m.t.ord(t) - 1]] + m.b[k, p, t]
    m.logic_5 = Constraint(m.fac, m.part, m.t, rule=logic_5)

    def logic_6(m, k, t):
        return sum(m.w[k, p, t] for p in m.part) <= 1
    m.logic_6 = Constraint(m.fac, m.t, rule=logic_6)

    def logic_7(m, i, k, t):
        return sum(m.z_supp2fac[i, k, p, t] for p in m.part) <= 1
    m.logic_7 = Constraint(m.suppl, m.fac, m.t, rule=logic_7)

    def logic_8(m, k, j, t):
        return sum(m.z_fac2mkt[k, j, p, t] for p in m.part) <= 1
    m.logic_8 = Constraint(m.fac, m.mkt, m.t, rule=logic_8)

    def logic_9(m, k):
        return sum(m.b[k, p, t] for p in m.part for t in m.t) <= 1
    m.logic_9 = Constraint(m.fac, rule=logic_9)

    def sym_1(m, l, u, t):
        if m.distr.ord(l) < m.distr.ord(u):
            return sum(m.w[l, p, t] for p in m.part) >= sum(m.w[u, p, t] for p in m.part)
        return Constraint.Skip
    m.sym_1 = Constraint(m.distr, m.distr, m.t, rule=sym_1)

    def sym_2(m, n, v, t):
        if m.centr.ord(n) < m.centr.ord(v):
            return sum(m.w[n, p, t] for p in m.part ) >= sum(m.w[v, p, t] for p in m.part)
        return Constraint.Skip
    m.sym_2 = Constraint(m.centr, m.centr, m.t, rule=sym_2)

    def sym_3(m, l, u, p, t):
        if m.distr.ord(l) < m.distr.ord(u):
            return sum(m.w[l, pp, t] for pp in m.part if m.part.ord(pp) <= m.part.ord(p)) >= m.w[u, p, t]
        return Constraint.Skip
    m.sym_3 = Constraint(m.distr, m.distr, m.part, m.t, rule=sym_3)

    def sym_4(m, n, v, p, t):
        if m.centr.ord(n) < m.centr.ord(v):
            return sum(m.w[n, pp, t] for pp in m.part if m.part.ord(pp) <= m.part.ord(p)) >= m.w[v, p, t]
        return Constraint.Skip
    m.sym_4 = Constraint(m.centr, m.centr, m.part, m.t, rule=sym_4)

    def complement_set_partitions(m):
        if if_complement:
            return sum(m.b[k, p, t] for p in m.part if select_part[p] == 0 for t in m.t for k in m.fac) >= 1
        return Constraint.Skip
    m.complement_set_partitions = Constraint(rule=complement_set_partitions)

    return m
