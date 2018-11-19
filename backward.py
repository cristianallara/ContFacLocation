from pyomo.environ import *


def backward_pass(t_ord, bl, equality_constraints, time_limit_mip, bl_LR):

    # if t_ord != 4:
    #     bl.fut_cost.pprint()

    # Solve relaxed LP model
    lpsolver = SolverFactory('gurobi')
    if t_ord != 1:
        lpsolver.options['relax_integrality'] = 1
    lpsolver.options['timelimit'] = time_limit_mip
    lpsolver.options['threads'] = 6
    lpsolver.options['FeasibilityTol'] = 1e-9
    lpsolver.solve(bl)  # , tee=True)
    # bl.w.pprint()

    if t_ord != 1:
        cost_LP = bl.obj()
        # get LP multipliers
        mltp = {key: - bl.dual[equality_constraints[key]] for key in equality_constraints.keys()}

        # # Initialize multipliers for Lagrangean relaxation
        for key in equality_constraints.keys():
            bl_LR.mltp[key] = mltp[key]

        # # Solve Lagrangean relaxation
        lrsolver = SolverFactory('gurobi')
        lrsolver.options['timelimit'] = time_limit_mip
        lrsolver.options['threads'] = 6
        lrsolver.solve(bl_LR)  # , tee=True)
        lrsolver.options['FeasibilityTol'] = 1e-9
        cost_LR = bl_LR.obj()

        # Solve discrete model
        mipsolver = SolverFactory('gurobi')
        mipsolver.options['timelimit'] = time_limit_mip
        mipsolver.options['threads'] = 6
        mipsolver.options['FeasibilityTol'] = 1e-9
        mipsolver.solve(bl)  # , tee=True)
        cost = bl.obj()

    else:
        cost = bl.obj()
        cost_LR = 0
        cost_LP = 0
        mltp = {}

    return cost, cost_LR, cost_LP, mltp
