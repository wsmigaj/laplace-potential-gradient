#!/usr/bin/env python

# Help Python find the bempp module
# (ADJUST THE PATH BELOW TO MATCH YOUR BEM++ INSTALLATION LOCATION)
bemppInstallPrefix = "/home/wojtek/Projects/BEM/bempp-2.0.2/bempp/install/"
import sys
sys.path.append(bemppInstallPrefix + "bempp/python")

from bempp.lib import *
from bempp.laplace_gradient import createLaplace3dSingleLayerPotentialOperatorGradient
import numpy as np

# Boundary conditions

def evalDirichletData(point):
    x, y, z = point
    r = np.sqrt(x**2 + y**2 + z**2)
    return 2 * x * z / r**5 - y / r**3

# Exact solution

def evalExactNeumannData(point):
    x, y, z = point
    r = np.sqrt(x**2 + y**2 + z**2)
    return -6 * x * z / r**6 + 2 * y / r**4

# Load mesh

grid = createGridFactory().importGmshGrid(
    "triangular", bemppInstallPrefix + "bempp/examples/meshes/sphere-h-0.1.msh")

# Create quadrature strategy

accuracyOptions = createAccuracyOptions()
# Increase by 2 the order of quadrature rule used to approximate
# integrals of regular functions on pairs on elements
accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2)
# Increase by 2 the order of quadrature rule used to approximate
# integrals of regular functions on single elements
accuracyOptions.singleRegular.setRelativeQuadratureOrder(2)
quadStrategy = createNumericalQuadratureStrategy(
    "float64", "float64", accuracyOptions)

# Create assembly context

assemblyOptions = createAssemblyOptions()
assemblyOptions.switchToAcaMode(createAcaOptions())
context = createContext(quadStrategy, assemblyOptions)

# Initialize spaces

pwiseConstants = createPiecewiseConstantScalarSpace(context, grid)
pwiseLinears = createPiecewiseLinearContinuousScalarSpace(context, grid)

# Construct elementary operators

slpOp = createLaplace3dSingleLayerBoundaryOperator(
    context, pwiseConstants, pwiseLinears, pwiseConstants)
dlpOp = createLaplace3dDoubleLayerBoundaryOperator(
    context, pwiseLinears, pwiseLinears, pwiseConstants)
idOp = createIdentityOperator(
    context, pwiseLinears, pwiseLinears, pwiseConstants)

# Form the left- and right-hand-side operators

lhsOp = slpOp
rhsOp = -0.5 * idOp + dlpOp

# Construct the grid function representing the (input) Dirichlet data

dirichletData = createGridFunction(
    context, pwiseLinears, pwiseLinears, evalDirichletData)

# Construct the right-hand-side grid function

rhs = rhsOp * dirichletData

# Initialize the solver

solver = createDefaultIterativeSolver(lhsOp)
solver.initializeSolver(defaultGmresParameterList(1e-5))

# Solve the equation

solution = solver.solve(rhs)
print solution.solverMessage()

# Extract the solution in the form of a grid function

solFun = solution.gridFunction()

# Create potential operators

slPotOp = createLaplace3dSingleLayerPotentialOperator(context)
slPotOpGradient = createLaplace3dSingleLayerPotentialOperatorGradient(context)

# Define points at which the potential operators should be evaluated.
# We want to compare the y derivative of the single-layer potential operator
# calculated numerically (via finite differences) with one calculated using
# the custom operator slPotOpGradient.

spacing = 0.1
y = np.arange(2 - spacing, 3 + spacing + 1e-8, spacing)
evaluationPoints = np.vstack((np.zeros_like(y), y, np.zeros_like(y)))
evaluationOptions = createEvaluationOptions()
field = slPotOp.evaluateAtPoints(solFun, evaluationPoints, evaluationOptions)

numericalYDerivative = (field[0,2:] - field[0,:-2]) / (2 * spacing)

y = y[1:-1]
evaluationPoints = np.vstack((np.zeros_like(y), y, np.zeros_like(y)))
gradient = slPotOpGradient.evaluateAtPoints(solFun, evaluationPoints, evaluationOptions)

# Note: gradient is an array of dimensions (3, number_of_points).
# Each column stores the three components of the gradient of the
# potential operator at a particular point.
analyticalYDerivative = gradient[1,:]

print "Numerical y derivative", "Analytical y derivative"
print np.column_stack((numericalYDerivative, analyticalYDerivative))
