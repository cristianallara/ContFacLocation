from pyomo.environ import *


def backward_pass(t_ord, bl, equality_constraints, time_limit_mip, bl_LR, benders_cut, integer_cut,
                  strength_benders_cut):

    if t_ord != 4:
        bl.fut_cost.pprint()

    if t_ord != 1:
        if benders_cut or strength_benders_cut:
            # Solve relaxed LP model
            lpsolver = SolverFactory('gurobi')
            lpsolver.options['relax_integrality'] = 1
            lpsolver.options['timelimit'] = time_limit_mip
            lpsolver.options['threads'] = 6
            lpsolver.options['FeasibilityTol'] = 1e-9
            lpsolver.solve(bl)  # , tee=True)
            cost_LP = bl.obj()
            cost = cost_LP
            cost_LR = 0
            bl.w.pprint()

            # get LP multipliers
            mltp = {key: - bl.dual[equality_constraints[key]] for key in equality_constraints.keys()}

        if strength_benders_cut:
            # Initialize multipliers for Lagrangean relaxation
            for key in equality_constraints.keys():
                bl_LR.mltp[key] = mltp[key]

            # Solve Lagrangean relaxation
            lrsolver = SolverFactory('gurobi')
            lrsolver.options['timelimit'] = time_limit_mip
            lrsolver.options['threads'] = 6
            lrsolver.solve(bl_LR)  # , tee=True)
            lrsolver.options['FeasibilityTol'] = 1e-9
            cost_LR = bl_LR.obj()
            cost = cost_LR

        if integer_cut:
            # Solve discrete model
            mipsolver = SolverFactory('gurobi')
            mipsolver.options['timelimit'] = time_limit_mip
            mipsolver.options['threads'] = 6
            mipsolver.options['FeasibilityTol'] = 1e-9
            mipsolver.solve(bl)  # , tee=True)
            cost = bl.obj()
            if not strength_benders_cut:
                cost_LR = 0
            if not benders_cut:
                cost_LP = 0
    else:
        mipsolver = SolverFactory('gurobi')
        mipsolver.options['timelimit'] = time_limit_mip
        mipsolver.options['threads'] = 6
        mipsolver.options['FeasibilityTol'] = 1e-9
        mipsolver.solve(bl)  # , tee=True)

        bl.w.pprint()
        cost = bl.obj()
        cost_LR = 0
        cost_LP = 0
        mltp = {}

    return cost, cost_LR, cost_LP, mltp
