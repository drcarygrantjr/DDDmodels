function NewSnowDensityEB(T)

New_Snow_Density = 0.05
T =  T * 9.0/5.0 + 32.0                     # from Celsius to Fahrenhet
  if(T > 0.0) 
      newdensity = New_Snow_Density + (T/100)^2
  else
      newdensity = New_Snow_Density 
  end
return newdensity
end
