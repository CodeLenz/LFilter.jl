@testset "Neighbours" begin

   # Mesh 1 x 1 with 10 divisions. 
   m = Simply_supported2D(10,10,:solid2D)

   # Distance between centroids is 0.1. If we set  
   # radius larger than sqrt(2)*0.1 we will cover 
   # a first order radius
   #   ......
   #
   #31       36 37  38     40  
   #21 22    26 27  28     30
   #11 12 .. 16 17  18     20
   #1  2                   10
   #
   # 10% larger than sqrt(2)*0.1
   radius = 1.1*sqrt(2.0)*(0.1)

   # Build the neighbourd list and weigths
   viz, weights = LFilter.Neighbours(m,radius)

   # Check sizes
   @test length(viz)==m.bmesh.ne
   @test length(weights)==m.bmesh.ne

   #               Check some elements

   #                First element

   # Neighbours
   @test all(viz[1].==[1; 2; 11; 12])

   # Weights
   distances = 0.1*[0; 1; 1; sqrt(2.0)]
   @test isapprox(weights[1],1.0.-(distances/radius))
   
   #                    Element 27

   # Neighbours
   @test all(viz[27].==[16; 17; 18; 26; 27; 28; 36; 37; 38])

   # Weights
   distances = 0.1*[sqrt(2.0); 1; sqrt(2.0); 1; 0; 1; sqrt(2.0); 1; sqrt(2.0)]
   @test isapprox(weights[27],1.0.-(distances/radius))
   
   # Check for no neighbourd
   @test_throws String LFilter.Neighbours(m,0.001)

   # Check for valid radius
   @test_throws AssertionError LFilter.Neighbours(m,0.0)
   @test_throws AssertionError LFilter.Neighbours(m,-1.0)



end