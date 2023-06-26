# Function for distributing input to the differnt groundwater layers according to their individual
# level of saturation
#
# Author: Thomas Skaugen
# Revised: 16.12.2019
#################################################################################################
                                       
# ddist provides the weights to distribute outx to different layers
# The sum of ddsit[1:NoL] is always 1 an the higher ddist the more water is received by the layer
################################################################################################

function GrvInputDistribution(outx, NoL,ddistx, ddist)

redoutx = outx                       #outx is reduced sucessively by the layers, starting from the slowest, no NoL
ddist[1:NoL] .= 0.0

for j in reverse(1:NoL)                   # Remember NoL is the slowest layer, 1 is the fastest  
    
    if(redoutx > 0.0)    
       differ = ddistx[j]-redoutx
       if (differ < 0)                #i.e the deficit i layer j is les than input, input(outx (redoutx)) > deficit       
         ddist[j] = ddistx[j]/outx    # divide by outx informs us what frafction of outx needed for this layer
         redoutx = redoutx-ddistx[j] # must reduce  the input correspondingly 
       end
       if (differ >= 0)               # i.e. deficit in layer j is more than input
         if (j < NoL)
             ddist[j]  =  1.0 - sum(ddist[(j+1):NoL])
         end
         if (j == Int(NoL))
             ddist[j] = 1.0
         end
         if(j > 1) 
             ddist[(j-1):1] .= 0.0  # We allow for layer 1 (overland flow layer) receives the rest in case input > layer capacity
         end
         redoutx = 0.0
       end
    
    end
   end   

return ddist
end
