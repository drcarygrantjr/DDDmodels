#    Function Layer_capacity_init
#
#------------------------------------------------------------------------
#     Description:  Calculates current capacity in groundwater layers from subsurface state. 
#
#     Author: Thomas Skaugen
#     Revised: 13.06.2023
#--------------------------------------------------------------------------

function LayerCapacityInit(Layers, nodaysvector, Magkap, NoL)

ddistx = zeros(Float64,NoL)
aktMag = zeros(Float64,NoL)

#Below are the states (in mm) for each saturation level
for j in reverse(1:NoL)
                                      #state due to from groundwater state NOT minus current timestep
  aktMag[j] = sum(Layers[j,1:nodaysvector[j]])

  if (aktMag[j] < Magkap[j])
   ddistx[j] = Magkap[j] - aktMag[j]
  end
        
end

return ddistx  # ddistx informs on current capacity for each level in mm.
end

