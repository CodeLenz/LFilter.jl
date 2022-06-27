# LFilter
Spatial filter with projection for Topology Optimization of continuum structures

```julia

using BMesh, LMesh, TMeshes
using LFEM 
using LFilter

# Create a 2D solid mesh
m = Simply_supported2D(10,10,:solid2D);

# Create a filter with radius 0.3m
filter = Filter(m,0.3)

# Create an input vector (zeros)
x = zeros(m.bmesh.ne)

# Put an 1.0 in element 55
x[55]=1.0

# Filter 
ρ = x2proj(x,filter)

# Create a gmsh output file
name = "example.pos"
Gmsh_init(name,m)

# Export original (x) and filtered (ρ) fields
Gmsh_element_scalar(m,x,name,"Original")
Gmsh_element_scalar(m,ρ,name,"Filtered")

# It is also possible to fix some elements during the filtering/optimization
# Lets set element 56 to 1E-2 and 100 to 1.0
m.options[:Fixed_elements]=[56 1E-3 ; 100 1.0]

# Create a filter with radius 0.3m
filter2 = Filter(m,0.3)

# Filter 
ρ2 = x2proj(x,filter2)

# Export the new field
Gmsh_element_scalar(m,ρ2,name,"Filtered and fixed")


```
