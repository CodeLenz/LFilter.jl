# LFilter
Spatial filter with projection for Topology Optimization of continuum structures

```julia

using BMesh, LMesh, LFilter

# Create a 2D solid mesh
m = Simply_supported2D(10,10,:solid2D)

# Create a filter with radius 0.3m
filter = Filter(m,0.3)

# Create an input vector
x = ones(m.bmesh.ne)

# Filter 
ρ = x2proj(x,filter)

```
