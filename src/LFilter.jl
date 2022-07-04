module LFilter

using BMesh, LMesh
using LinearAlgebra, SparseArrays

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
       @eval Base.Experimental.@optlevel 3
   end

include("viz.jl")
include("filter.jl")

export Filter
export x2proj
export dproj2dx!

# Precompile
precompile(Neighbours,Mesh2D,Float64))
precompile(Neighbours,Mesh3D,Float64))
precompile(Map_x_ρ,(Float64,Mesh2D))
precompile(Map_x_ρ,(Float64,Mesh3D))

end # module LFilter
