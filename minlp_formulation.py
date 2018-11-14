from pyomo.environ import *


def create_minlp(data):

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

    #Decomposition Parameters
    m.w_prev_par = Param(m.fac, m.t, initialize=0, mutable=True)

    # Bounds
    x_min = min(min(m.suppl_x[i] for i in m.suppl), min(m.mkt_x[j] for j in m.mkt))
    x_max = max(max(m.suppl_x[i] for i in m.suppl), max(m.mkt_x[j] for j in m.mkt))
    y_min = min(min(m.suppl_y[i] for i in m.suppl), min(m.mkt_y[j] for j in m.mkt))
    y_max = max(max(m.suppl_y[i] for i in m.suppl), max(m.mkt_y[j] for j in m.mkt))
    dist_max = sqrt((x_max-x_min)**2 + (y_max-y_min)**2)
    dist_min = 0.5 # arbitrary

    # Block of Equations per time period
    def planning_block_rule(b, t):

        # Variables declaration
        b.f_supp2fac = Var(m.suppl, m.fac, domain=NonNegativeReals)
        b.f_fac2mkt = Var(m.fac, m.mkt, domain=NonNegativeReals)
        b.f_fac = Var(m.fac, domain=NonNegativeReals)
        b.fac_x = Var(m.fac, domain=NonNegativeReals, bounds=(x_min, x_max))
        b.fac_y = Var(m.fac, domain=NonNegativeReals, bounds=(y_min, y_max))
        b.dist_supp2fac = Var(m.suppl, m.fac, domain=NonNegativeReals, bounds=(dist_min, dist_max))
        b.dist_fac2mkt = Var(m.fac, m.mkt, domain=NonNegativeReals, bounds=(dist_min, dist_max))
        b.cost_supp2fac = Var(domain=NonNegativeReals)
        b.cost_fac2mkt = Var(domain=NonNegativeReals)
        b.inv_cost = Var(domain=NonNegativeReals)
        b.op_cost = Var(domain=NonNegativeReals)
        b.w = Var(m.fac, domain=Binary) # binary for when facility is in place and operating
        b.w_prev = Var(m.fac, domain=Binary)
        b.b = Var(m.fac, domain=Binary) # binary for building facility k in time period t
        b.z_supp2fac = Var(m.suppl, m.fac, domain=Binary)
        b.z_fac2mkt = Var(m.fac, m.mkt, domain=Binary)

        # Objective function
        def obj_rule(_b):
            return m.i_factor[t] * (_b.inv_cost + _b.op_cost + _b.cost_supp2fac + _b.cost_fac2mkt)
        b.obj = Objective(rule=obj_rule, sense=minimize)

        # Constraints
        def investment_cost(_b):
            return _b.inv_cost == sum(m.FIC[k, t]*_b.b[k] + m.VIC[k, t] * _b.f_fac[k] for k in m.fac)
        b.investment_cost = Constraint(rule=investment_cost) # TODO: change VIC to depend on cap

        def operating_cost(_b):
            return _b.op_cost == sum(m.FOC[k, t]*_b.w[k] + m.VOC[k, t] * _b.f_fac[k] for k in m.fac)
        b.operating_cost = Constraint(rule=operating_cost)

        def transp_cost_supp2fac(_b):
            return _b.cost_supp2fac == sum(m.RM[i, t]*_b.f_supp2fac[i, k] + m.FTC_supp2fac[i, k, t]*_b.z_supp2fac[i, k]
                                           + m.VTC_supp2fac[i, k, t]*_b.dist_supp2fac[i, k]*_b.f_supp2fac[i, k]
                                           for i in m.suppl for k in m.fac)
        b.transp_cost_supp2fac = Constraint(rule=transp_cost_supp2fac)

        def transp_cost_fac2mkt(_b):
            return _b.cost_fac2mkt == sum(m.FTC_fac2mkt[k, j, t] * _b.z_fac2mkt[k, j] + m.VTC_fac2mkt[k, j, t]
                                          * _b.dist_fac2mkt[k, j] * _b.f_fac2mkt[k, j] for k in m.fac for j in m.mkt)
        b.transp_cost_fac2mkt = Constraint(rule=transp_cost_fac2mkt)

        def mass_bal1(_b, k):
            return _b.f_fac[k] == sum(_b.f_supp2fac[i, k] * m.conv[k] for i in m.suppl)
        b.mass_bal1 = Constraint(m.fac, rule=mass_bal1)

        def mass_bal2(_b, k):
            return _b.f_fac[k] == sum(_b.f_fac2mkt[k, j] for j in m.mkt)
        b.mass_bal2 = Constraint(m.fac, rule=mass_bal2)

        def demand(_b, j):
            return sum(_b.f_fac2mkt[k, j] for k in m.fac) == m.dd[j, t]
        b.demand = Constraint(m.mkt, rule=demand)

        def availability(_b, i):
            return sum(_b.f_supp2fac[i, k] for k in m.fac) <= m.avail[i, t]
        b.availability = Constraint(m.suppl, rule=availability)

        def max_cap(_b, k):
            return _b.f_fac[k] <= m.cap[k]*_b.w[k]
        b.max_cap = Constraint(m.fac, rule=max_cap)

        def dist_supp2fac_rule(_b, i, k):
            return _b.dist_supp2fac[i, k] >= sqrt((m.suppl_x[i] - _b.fac_x[k])**2 + (m.suppl_y[i] - _b.fac_y[k])**2)
        b.dist_supp2fac_rule = Constraint(m.suppl, m.fac, rule=dist_supp2fac_rule)

        def dist_fac2mkt_rule(_b, k, j):
            return _b.dist_fac2mkt[k, j] >= sqrt((m.mkt_x[j] - _b.fac_x[k])**2 + (m.mkt_y[j] - _b.fac_y[k])**2)
        b.dist_fac2mkt_rule = Constraint(m.fac, m.mkt, rule=dist_fac2mkt_rule)

        def x_bigM(_b, k):
            return _b.fac_x[k] <= x_max * _b.w[k]
        b.x_bigM = Constraint(m.fac, rule=x_bigM)

        def y_bigM(_b, k):
            return _b.fac_y[k] <= y_max * _b.w[k]
        b.y_bigM = Constraint(m.fac, rule=y_bigM)

        def supp2fac_bigM(_b, i, k):
            return _b.f_supp2fac[i, k] <= m.avail[i, t]*_b.z_supp2fac[i, k]
        b.supp2fac_bigM = Constraint(m.suppl, m.fac, rule=supp2fac_bigM)

        def fac2mkt_bigM(_b, k, j):
            return _b.f_fac2mkt[k, j] <= m.dd[j, t]*_b.z_fac2mkt[k, j]
        b.fac2mkt_bigM = Constraint(m.fac, m.mkt, rule=fac2mkt_bigM)

        def fac_bigM(_b, k):
            return _b.f_fac[k] <= m.cap[k]*_b.w[k]
        b.fac_bigM = Constraint(m.fac, rule=fac_bigM)

        def logic_1(_b, i, k):
            return _b.w[k] >= _b.z_supp2fac[i, k]
        b.logic_1 = Constraint(m.suppl, m.fac, rule=logic_1)

        def logic_2(_b, k):
            return sum(_b.z_supp2fac[i, k] for i in m.suppl) >= _b.w[k]
        b.logic_2 = Constraint(m.fac, rule=logic_2)

        def logic_3(_b, k, j):
            return _b.w[k] >= _b.z_fac2mkt[k, j]
        b.logic_3 = Constraint(m.fac, m.mkt, rule=logic_3)

        def logic_4(_b, k):
            return sum(_b.z_fac2mkt[k, j] for j in m.mkt) >= _b.w[k]
        b.logic_4 = Constraint(m.fac, rule=logic_4)

        def logic_5(_b, k):
            return _b.w[k] == _b.w_prev[k] + _b.b[k]
        b.logic_5 = Constraint(m.fac, rule=logic_5)

        def equal_1(_b, k):
            return _b.w_prev[k] == m.w_prev_par[k, t]
        b.equal_1 = Constraint(m.fac, rule=equal_1)

        def sym_1(_b, l, u):
            if m.distr.ord(l) < m.distr.ord(u):
                return _b.fac_x[l] >= _b.fac_x[u]
            return Constraint.Skip
        b.sym_1 = Constraint(m.distr, m.distr, rule=sym_1)

        def sym_2(_b, n, v):
            if m.centr.ord(n) < m.centr.ord(v):
                return _b.fac_x[n] >= _b.fac_x[v]
            return Constraint.Skip
        b.sym_2 = Constraint(m.centr, m.centr, rule=sym_2)

        def sym_3(_b, l, u):
            if m.distr.ord(l) < m.distr.ord(u):
                return _b.w[l] + _b.w_prev[l] >= _b.w[u]
            return Constraint.Skip
        b.sym_3 = Constraint(m.distr, m.distr, rule=sym_3)

        def sym_4(_b, n, v):
            if m.centr.ord(n) < m.centr.ord(v):
                return _b.w[n] + _b.w_prev[n] >= _b.w[v]
            return Constraint.Skip
        b.sym_4 = Constraint(m.centr, m.centr, rule=sym_4)

    m.Bl = Block(m.t, rule=planning_block_rule)
    return m

