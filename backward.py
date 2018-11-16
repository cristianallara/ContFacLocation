from pyomo.environ import *


def backward_pass(t_ord, bl, equality_constraints, time_limit_mip):

    # Solve the model
    mipsolver = SolverFactory('gurobi')
    if t_ord != 1:
        mipsolver.options['relax_integrality'] = 1
    mipsolver.options['timelimit'] = time_limit_mip
    mipsolver.options['threads'] = 6
    mipsolver.solve(bl) #, tee=True)

    if t_ord != 1:
        mltp = {key: - bl.dual[equality_constraints[key]] for key in equality_constraints.keys()}
    else:
        mltp = {}
    # Get optimal value
    cost = bl.obj()

    return mltp, cost