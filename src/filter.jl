
#
# Main structure for filtering
#
mutable struct Filter

   # Radius
   radius::Float64

   # Spatial
   MAP::SparseMatrixCSC{Float64,Int64}

   # Positions with fixed elements
   posfix::Array{Int64}

   # Values of fixed elements
   valfix::Array{Float64}

   # Sharpness 
   β::Float64
   
   # Center
   η::Float64

   # Densidade mínima
   ρ_min::Float64 


   function Filter(mesh::Mesh,radius=0.0;β=1.0,ρ_min=1E-3,η=0.5) 

      # It just makes sense  for solid elements
      @assert contains(string(mesh.bmesh.etype),"solid") "Filter::only for meshes with :solid elements"
      
      # Basic assertions
      @assert radius>=0.0 "Filter:: radius must be >=0"
      @assert β>=1.0 "Filter:: β must be >= 1.0"
      @assert 0<ρ_min<1 "Filter:: ρ_min must be in (0,1)"
      @assert  0<η<1 "Filter:: η must be in (0,1)"
      
      # If radius is 0.0, find a minimum radius 
      if radius==0.0
         radius = max(mesh.bmesh.Lx/mesh.bmesh.nx, mesh.bmesh.Ly/mesh.bmesh.ny)*1.5
      end

      # Cria a matrix que mapeia as variáveis de projeto matemáticas
      # para o vetor de densidades relativas (FILTRO)
      # MAP é uma matriz esparsa
      MAP = Map_x_ρ(radius,mesh)

      # Force a 1<->1  mapping between ρ and x for each fixed element in the design
      # and generate two lists: one containing the fixed elementos and ohter with
      # the fixed values
      posfix = Int64[]
      valfix = Float64[]
      if haskey(mesh.options,:Fixed_elements)

         # Alias and number of fixed elements 
         fixed_elements =  mesh.options[:Fixed_elements]
         nfixed = size(fixed_elements,1)

         @inbounds for i=1:nfixed
             ele = Int(fixed_elements[i,1])
             val = fixed_elements[i,2]
             MAP[ele,:] .= 0.0
             MAP[ele,ele] = 1.0    
             push!(posfix,ele)
             push!(valfix,val)
         end

         # Remove zeros
         dropzeros!(MAP)

      end    
    
      new(radius,MAP,posfix,valfix,β,η,ρ_min)

   end

end



#
# Monta a matriz F que mapeia x para ρ, ou seja,
# filtro de densidades padrão
#
function Map_x_ρ(R::Float64,mesh::Mesh)

    # Obtem os vizinhos e pesos para o raio/malha
    viz, weigths = Neighbours(mesh,R)

    # Sparse matrix
    VI = Int64[] 
    VJ = Int64[]
    VV = Float64[]

    # For each element, put the neighbours in columns. Lets
    # divide by the sum of the weigths 
    @inbounds for ele in mesh

        # Sum of weigths
        somat = sum(weigths[ele])

        # Scan the neighbours
        @inbounds for par in zip(viz[ele],weigths[ele]/somat) 

            push!(VI,ele)
            push!(VJ,par[1])
            push!(VV,par[2])
             
        end # zip
    end # ele

    return sparse(VI,VJ,VV)
    
end


#
# Operador para a projeção tanh. Este operador
# não é linear (depende de ρ), portanto devemos 
# fazer na forma algoritmica
#
function Map_ρ_tanh(β::Float64,η::Float64,ρ::Vector{Float64})

    # Primeiro alocamos um vetor com a mesma dimensão 
    # de ρ
    ρ_proj = Array{Float64}(undef,length(ρ))

    # Cte
    cte = tanh(β*η)

    cte2 = (cte+tanh(β*(1-η)))

    cte2 !=0.0 || throw("Map_ρ_tanh:: check β and η")

    # Loop sobre os elementos, calculando a projeção
    @inbounds for ele in LinearIndices(ρ)

        ρ_proj[ele] = (tanh(β*(ρ[ele]-η))+cte)/cte2

    end

    return ρ_proj

end

#
# Operador de correção de derivada d ρ_proj / d ρ 
# que será utilizado na regra da cadeia para correção de derivadas
#
function dρ_projdρ(β::Float64,η::Float64,ρ::Vector{Float64})

    # Primeiro alocamos um vetor com a mesma dimensão 
    # de ρ
    operador = Array{Float64}(undef,length(ρ))

    # Cte
    cte =  (tanh(β*η)+tanh(β*(1-η)))

    # Assertion
    cte != 0.0 || throw("dρ_projdρ:: check  β and η")

    # Loop sobre os elementos, calculando a projeção
    @inbounds for ele in LinearIndices(ρ)

         operador[ele] = (β*sech(β*(ρ[ele]-η))^2)/cte

    end
    
    return Diagonal(operador)

end


#
# do the Map 
#
function x2ρ(x::Vector{Float64},filter::Filter)

    # Assert fixed elements
    x[filter.posfix] .= filter.valfix

    # Map
    ρ = filter.MAP*x
     
end


#
# do the proj -> this is the one to use 
#
function x2proj(x::Vector{Float64},filter::Filter)
   
    # Mapeamento para as densidades - Filtro
    ρ = x2ρ(x,filter)

    # Projeção heaviside
    ρ_proj = Map_ρ_tanh(filter.β,filter.η,ρ)

    # Correção do rho mínimo
    return filter.ρ_min .+ (1.0-filter.ρ_min).*ρ_proj


end

#
# Fix the sensitivities
#
function dproj2dx!(x::Vector{Float64},dproj::Vector{Float64},filter::Filter)

    # Mapeamento para as densidades - Filtro
    ρ = x2ρ(x,filter)
    
    # Operador de correção da projeção (matriz diagonal)
    R = dρ_projdρ(filter.β,filter.η,ρ) 

    # Corrige as derivadas devido ao filtro
    # Cuidado com o MAP'   
    dx = R*filter.MAP'*dproj

    # Corrige para o valor mínimo
    dx .= dx.*(1.0-filter.ρ_min)

    # Assert fixed elements
    dx[filter.posfix] .= 0.0

    dproj.= dx

end

