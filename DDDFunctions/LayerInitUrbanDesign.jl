#    Function LayerInitUrbanDesign
#
#------------------------------------------------------------------------
#     Description:  Initializes Layers with drawn Subsurface state. 
#
#     Author: Thomas Skaugen
#     Revised: 13.06.2023
#--------------------------------------------------------------------------

function LayerInitUrbanDesign(ddist, outx, Layers, layerUH, nodaysvector, NoL)

# ddist: 1 dim array float
# Layers: 2 dim array float
# layersUH: 2 dim array float
# NoL: scalar, integer    
# nodaysvector: 1 dim array integer
# outx: scalar, float
        
  for j in 1 : NoL
    qlayer = zeros(nodaysvector[j])    
    qlayer .= ddist[j]*outx .* layerUH[j,1:nodaysvector[j]]    #finds response  in mm!!!!  for the actual layer, en vektor     
    Layers[j,1:nodaysvector[j]] .= qlayer[1:nodaysvector[j]]    
  end

return Layers
end                
