#    Function Temperaturevector
#
#------------------------------------------------------------------------
#     Description:  Calculates a temperature vector for snowpack temp estimation 
#
#     Author: Thomas Skaugen
#     Revised: 8.1.2020
#--------------------------------------------------------------------------


function TemperatureVector(tempmat,idim,len)
  STempvec = zeros(len)  
  STempvec[1:len] .= tempmat[reverse(1:len),idim]
  return STempvec
end

