module EconPDEs
using LinearAlgebra
import NLsolve: nlsolve
using DataStructures
using ForwardDiff
using Interpolations
#using BandedMatrices
#using ReverseDiff
#using DifferentialEquations
#using LSODA
#using SteadyStateDiffEq
##############################################################################
##
## Load files
##
##############################################################################
include("finiteschemesolve.jl")
include("pdesolve.jl")
include("utils.jl")


##############################################################################
##
## Exported methods and types 
##
##############################################################################
export finiteschemesolve,
pdesolve,
simulate,
OrderedDict
end