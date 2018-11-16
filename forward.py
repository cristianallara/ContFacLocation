from pyomo.environ import *


def forward_pass(bl, state_vars, opt_tol_mip, time_limit_mip):

    # Solve MIP for LB
    mipsolver = SolverFactory('gurobi')
    mipsolver.options['mipgap'] = opt_tol_mip
    mipsolver.options['timelimit'] = time_limit_mip
    mipsolver.options['threads'] = 6
    mipsolver.solve(bl) # , tee=True)

    # Fix the state variable as parameter for next t

    state_vars_val = {key: state_vars[key].value for key in state_vars.keys()}

    # Store obj value to compute UB
    cost = bl.obj() - bl.alphafut.value

    return state_vars_val, cost

