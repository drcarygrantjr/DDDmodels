#    Function Layer_update
#
#------------------------------------------------------------------------
#     Description:  Updates Layers by todays event and shifting to one timestepahead. 
#
#     Author: Thomas Skaugen
#     Revised: 16.12.2019
#--------------------------------------------------------------------------

function BogLayerUpdate(outbog, BogLayers, UHBog, nodaysvector)

    qlayer = zeros(nodaysvector)
    qlayer .= outbog*UHBog    #finds response  in mm!!!!  for bog, a vector     
    if(nodaysvector > 1)    
      BogLayers[1:(nodaysvector-1)] .= BogLayers[2:nodaysvector] .+ qlayer[2:nodaysvector]# shifts the level of the matrix one timestep ahead
      BogLayers[nodaysvector] = 0.0
    end 
    if(nodaysvector == 1)
	  BogLayers[1:nodaysvector] .= 0.0 
      #BogLayers[1:nodaysvector] .= qlayer
    end

return BogLayers
end            


