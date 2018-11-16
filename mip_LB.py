from pyomo.environ import *


def create_mip(data, p_x, p_y):

    suppliers, xi, yi, time_periods, markets, xj, yj, centr_facilities, distr_facilities, facilities, cv, mc, a, d, \
    RM, FIC, VIC, FOC, VOC, ft1, ft2, vt1, vt2, interest_factor = data

    m = ConcreteModel()

    # Set declarations
    m.part = Set(initialize=RangeSet(p_x*p_y), ordered=True)
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

    #Decomposition Parameters
    m.w_prev_par = Param(m.fac, m.part, m.t, initialize=0, mutable=True)

    # Block of Equations per time period
    def planning_block_rule(b, t):

        # Variables declaration
        b.f_supp2fac = Var(m.suppl, m.fac, m.part, within=NonNegativeReals)
        b.f_fac2mkt = Var(m.fac, m.mkt, m.part, within=NonNegativeReals)
        b.f_fac = Var(m.fac, m.part, within=NonNegativeReals)
        b.cost_supp2fac = Var(within=NonNegativeReals)
        b.cost_fac2mkt = Var(within=NonNegativeReals)
        b.inv_cost = Var(within=NonNegativeReals)
        b.op_cost = Var(within=NonNegativeReals)
        b.w = Var(m.fac, m.part, within=Binary) # binary for when facility is in place and operating
        b.w_prev = Var(m.fac, m.part, within=Binary)
        b.b = Var(m.fac, m.part, within=Binary) # binary for building facility k in time period t
        b.z_supp2fac = Var(m.suppl, m.fac, m.part,  within=Binary)
        b.z_fac2mkt = Var(m.fac, m.mkt, m.part, within=Binary)
        b.alphafut = Var(within=NonNegativeReals)

        # Objective function
        def obj_rule(_b):
            return m.i_factor[t] * (_b.inv_cost + _b.op_cost + _b.cost_supp2fac + _b.cost_fac2mkt) + _b.alphafut
        b.obj = Objective(rule=obj_rule, sense=minimize)

        # Constraints
        def investment_cost(_b):
            return _b.inv_cost == sum((m.FIC[k, t] + m.VIC[k, t]*m.cap[k])*_b.b[k, p] for k in m.fac for p in m.part)
        b.investment_cost = Constraint(rule=investment_cost)

        def operating_cost(_b):
            return _b.op_cost == sum(m.FOC[k, t]*_b.w[k, p] + m.VOC[k, t] * _b.f_fac[k, p] for k in m.fac
                                     for p in m.part)
        b.operating_cost = Constraint(rule=operating_cost)

        def transp_cost_supp2fac(_b):
            return _b.cost_supp2fac == sum(m.RM[i, t]*_b.f_supp2fac[i, k, p] + m.FTC_supp2fac[i, k, t]
                                           *_b.z_supp2fac[i, k, p] + m.VTC_supp2fac[i, k, t]*m.dist_supp[i, p]
                                           *_b.f_supp2fac[i, k, p] for i in m.suppl for k in m.fac for p in m.part)
        b.transp_cost_supp2fac = Constraint(rule=transp_cost_supp2fac)

        def transp_cost_fac2mkt(_b):
            return _b.cost_fac2mkt == sum(m.FTC_fac2mkt[k, j, t] * _b.z_fac2mkt[k, j, p] + m.VTC_fac2mkt[k, j, t]
                                          * m.dist_mkt[j, p] * _b.f_fac2mkt[k, j, p] for k in m.fac
                                          for j in m.mkt for p in m.part)
        b.transp_cost_fac2mkt = Constraint(rule=transp_cost_fac2mkt)

        def mass_bal1(_b, k, p):
            return _b.f_fac[k, p] == sum(_b.f_supp2fac[i, k, p] * m.conv[k] for i in m.suppl)
        b.mass_bal1 = Constraint(m.fac, m.part, rule=mass_bal1)

        def mass_bal2(_b, k, p):
            return _b.f_fac[k, p] == sum(_b.f_fac2mkt[k, j, p] for j in m.mkt)
        b.mass_bal2 = Constraint(m.fac, m.part, rule=mass_bal2)

        def demand(_b, j):
            return sum(_b.f_fac2mkt[k, j, p] for k in m.fac for p in m.part) == m.dd[j, t]
        b.demand = Constraint(m.mkt, rule=demand)

        def availability(_b, i):
            return sum(_b.f_supp2fac[i, k, p] for k in m.fac for p in m.part) <= m.avail[i, t]
        b.availability = Constraint(m.suppl, rule=availability)

        def max_cap(_b, k, p):
            return _b.f_fac[k, p] <= m.cap[k]*_b.w[k, p]
        b.max_cap = Constraint(m.fac, m.part, rule=max_cap)

        def supp2fac_bigM(_b, i, k, p):
            return _b.f_supp2fac[i, k, p] <= m.avail[i, t]*_b.z_supp2fac[i, k, p]
        b.supp2fac_bigM = Constraint(m.suppl, m.fac, m.part, rule=supp2fac_bigM)

        def fac2mkt_bigM(_b, k, j, p):
            return _b.f_fac2mkt[k, j, p] <= m.dd[j, t]*_b.z_fac2mkt[k, j, p]
        b.fac2mkt_bigM = Constraint(m.fac, m.mkt, m.part, rule=fac2mkt_bigM)

        def fac_bigM(_b, k, p):
            return _b.f_fac[k, p] <= m.cap[k]*_b.w[k, p]
        b.fac_bigM = Constraint(m.fac, m.part, rule=fac_bigM)

        def logic_1(_b, i, k, p):
            return _b.w[k, p] >= _b.z_supp2fac[i, k, p]
        b.logic_1 = Constraint(m.suppl, m.fac, m.part, rule=logic_1)

        def logic_2(_b, k, p):
            return sum(_b.z_supp2fac[i, k, p] for i in m.suppl) >= _b.w[k, p]
        b.logic_2 = Constraint(m.fac, m.part, rule=logic_2)

        def logic_3(_b, k, j, p):
            return _b.w[k, p] >= _b.z_fac2mkt[k, j, p]
        b.logic_3 = Constraint(m.fac, m.mkt, m.part, rule=logic_3)

        def logic_4(_b, k, p):
            return sum(_b.z_fac2mkt[k, j, p] for j in m.mkt) >= _b.w[k, p]
        b.logic_4 = Constraint(m.fac, m.part, rule=logic_4)

        def logic_5(_b, k, p):
            return _b.w[k, p] == _b.w_prev[k, p] + _b.b[k, p]
        b.logic_5 = Constraint(m.fac, m.part, rule=logic_5)

        def logic_6(_b, k):
            return sum(_b.w[k, p] for p in m.part) <= 1
        b.logic_6 = Constraint(m.fac, rule=logic_6)

        def logic_7(_b, i, k):
            return sum(_b.z_supp2fac[i, k, p] for p in m.part) <= 1
        b.logic_7 = Constraint(m.suppl, m.fac, rule=logic_7)

        def logic_8(_b, k, j):
            return sum(_b.z_fac2mkt[k, j, p] for p in m.part) <= 1
        b.logic_8 = Constraint(m.fac, m.mkt, rule=logic_8)

        def logic_9(_b, k):
            return sum(_b.b[k, p] for p in m.part) <= 1
        b.logic_9 = Constraint(m.fac, rule=logic_9)

        def equal(_b, k, p):
            if m.t.ord(t) > 1:
                return _b.w_prev[k, p] == m.w_prev_par[k, p, m.t[m.t.ord(t)-1]]
            return _b.w_prev[k, p] == 0
        b.equal = Constraint(m.fac, m.part, rule=equal)

        def sym_1(_b, l, u):
            if m.distr.ord(l) < m.distr.ord(u):
                return sum(_b.w[l, p] + _b.w_prev[l, p] for p in m.part) >= sum(_b.w[u, p] for p in m.part)
            return Constraint.Skip
        b.sym_1 = Constraint(m.distr, m.distr, rule=sym_1)

        def sym_2(_b, n, v):
            if m.centr.ord(n) < m.centr.ord(v):
                return sum(_b.w[n, p] + _b.w_prev[n, p] for p in m.part) >= sum(_b.w[v, p] for p in m.part)
            return Constraint.Skip
        b.sym_2 = Constraint(m.centr, m.centr, rule=sym_2)

        def sym_3(_b, l, u, p):
            if m.distr.ord(l) < m.distr.ord(u):
                return sum(_b.w[l, pp] for pp in m.part if m.part.ord(pp) <= m.part.ord(p)) >= _b.w[u, p]
            return Constraint.Skip
        b.sym_3 = Constraint(m.distr, m.distr, m.part, rule=sym_3)

        def sym_4(_b, n, v, p):
            if m.centr.ord(n) < m.centr.ord(v):
                return sum(_b.w[n, pp] for pp in m.part if m.part.ord(pp) <= m.part.ord(p)) >= _b.w[v, p]
            return Constraint.Skip
        b.sym_4 = Constraint(m.centr, m.centr, m.part, rule=sym_4)

        b.fut_cost = ConstraintList()

    m.Bl = Block(m.t, rule=planning_block_rule)
    return m

