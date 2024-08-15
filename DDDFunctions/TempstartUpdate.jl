#    Function tempstartUpdate
#
#------------------------------------------------------------------------
#     Description:  Updates the temperature matrix for snowpacktemperature estimation
#
#     Author: Thomas Skaugen
#     Revised: 19.02.2024
#--------------------------------------------------------------------------

function TempstartUpdate(tempmatrix, temps, len)

# tempmatrix: timeseries of length 5 days 8length dep on temporal resolution) for 10 elevation zones
# temps: this timesteps temperatures

  dummy = tempmatrix      

      dummy[1:(len-1), 1:10] .= dummy[2:len, 1:10]# flytter the level of the matrix one timestep ahead
      dummy[len, 1:10] = temps[1:10]
      tempmatrix = dummy 

return tempmatrix
end                
