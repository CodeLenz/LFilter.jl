
# Neighbours and weigths
function Neighbours(mesh::Mesh,radius::Float64)

    # Basic assertion
    @assert radius>0.0 "Neighbours:: radius must be >0.0"
    
    # Alias
    ne = Get_ne(mesh)

    # Maximum number of neighbours in the mesh
    num_max = 0

    # Tables to store the neighbours and weigts
    viz = Vector{Int64}[]
    weigths = Vector{Float64}[]

     # Loop for each element
     for ele in mesh

         # Pega a posição do centróide do elemento central
         c_ele = Centroid(mesh,ele) 

         # Local Arrays
         viz_ele = Int64[]
         weigths_ele = Float64[]

         # Scan the mesh and store the neighbours of ele
         for j in mesh

             # Centroid of j
             c_j = Centroid(mesh,j) 

             # Distance
             d = norm(c_j.-c_ele)

             # If inside the radius...
             if d<=radius
                 push!(viz_ele,j)
                 push!(weigths_ele, 1.0 - d/radius)
             end #if

         end #j

        # Assert if the element has neighbours
        length(viz_ele)>1 || throw("Neighbours:: element $ele does not have neighbours")

        # Maximum number of neighbours until now
        num_max = max(num_max,length(viz_ele))

        # Store the neighbours for ele
        push!(viz,viz_ele)
        push!(weigths,weigths_ele)

     end #ele

     # Return 
     return viz, weigths

end

