#    Function Layer_update
#
#------------------------------------------------------------------------
#     Description:  Updates Layers by todays event and shifting to one timestepahead. 
#
#     Author: Thomas Skaugen
#     Revised: 17.12.2019
#--------------------------------------------------------------------------

function LayerUpdate(ddist, outx, Layers, layerUH, nodaysvector, NoL)

# ddist: 1 dim array float
# Layers: 2 dim array float
# layersUH: 2 dim array float
# NoL: scalar, integer    
# nodaysvector: 1 dim array integer
# outx: scalar, float
        
  for j in 1 : NoL
    qlayer = zeros(nodaysvector[j])    
    qlayer .= ddist[j]*outx .* layerUH[j,1:nodaysvector[j]]    #finds response  in mm!!!!  for the actual layer, en vektor     
        
    if(nodaysvector[j] > 1)
      Layers[j,(1:(nodaysvector[j]-1))] .= Layers[j,2:nodaysvector[j]] .+ qlayer[2:nodaysvector[j]]# flytter the level of the matrix one timestep ahead
      Layers[j,nodaysvector[j]] = 0.0       
    end 

    if(nodaysvector[j] == 1)
      Layers[j,1:nodaysvector[j]] .= qlayer[1:nodaysvector[j]]
    end
  end

return Layers
end                
