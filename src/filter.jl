
#
# Main structure for filtering
#
mutable struct Filter

   # Radius
   radius::Float64

   # Spatial
   MAP::SparseMatrixCSC

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

      # If radius is 0.0, find a minimum radius 
      if radius==zero(radius)
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

         for i=1:nfixed
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
function Map_x_ρ(R::T,mesh::Mesh) where T

    # Alias
    ne = mesh.bmesh.ne

    # Obtem os vizinhos e pesos para o raio/malha
    viz, weigths = Neighbours(mesh,R)

    # Sparse matrix
    VI = Int64[]
    VJ = Int64[]
    VV = T[]

    # For each element, put the neighbours in columns. Lets
    # divide by the sum of the weigths 
    for ele=1:ne

        # Sum of weigths
        somat = sum(weigths[ele])

        # Scan the neighbours
        for par in zip(viz[ele],weigths[ele]/somat) 

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
function Map_ρ_tanh(β::T,η::T,ρ::Array{T,1}) where T

    # Primeiro alocamos um vetor com a mesma dimensão 
    # de ρ
    ρ_proj = Array{T}(undef,length(ρ))

    # Cte
    cte = tanh(β*η)

    # Loop sobre os elementos, calculando a projeção
    for ele in LinearIndices(ρ)

        ρ_proj[ele] = (tanh(β*(ρ[ele]-η))+cte)/(cte+tanh(β*(1-η)))

    end

    return ρ_proj

end

#
# Operador de correção de derivada d ρ_proj / d ρ 
# que será utilizado na regra da cadeia para correção de derivadas
#
function dρ_projdρ(β::T,η::T,ρ::Array{T,1}) where T

    # Primeiro alocamos um vetor com a mesma dimensão 
    # de ρ
    operador = Array{T}(undef,length(ρ))

    # Loop sobre os elementos, calculando a projeção
    for ele in LinearIndices(ρ)

         operador[ele] = (β*sech(β*(ρ[ele]-η))^2)/(tanh(β*η)+tanh(β*(1-η)))

    end
    
    return Diagonal(operador)

end


#
# do the Map 
#
function x2ρ(x::Array{T},filter::Filter) where T

    # Assert fixed elements
    x[filter.posfix] .= filter.valfix

    # Map
    ρ = filter.MAP*x
     
end


#
# do the proj -> this is the one to use 
#
function x2proj(x::Array{T},filter::Filter) where T

    
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
function dproj2dx!(x::Array{T},dproj::Array{T},filter::Filter) where T


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

