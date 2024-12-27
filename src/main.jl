using LinearAlgebra, Statistics, UnPack, Zygote, ArrayPadding, Porcupine
using Porcupine: keys, values, pairs
using Zygote: Buffer, @ignore_derivatives

include("ops.jl")
include("del.jl")

# using Pkg
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\beans\ArrayPadding.jl"
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\beans\Porcupine.jl;up"
