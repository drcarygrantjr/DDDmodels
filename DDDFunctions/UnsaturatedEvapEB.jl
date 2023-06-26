#    Function UnsaturatedEvapEB
#
#C------------------------------------------------------------------------
#
#     Description: Calculates soil moisture (unsaturated) and evapotranspiration (Priestley-Taylor). Some evapotranspiration may evaporate from saturated layes
#     Author: Thomas Skaugen
#     Revised: 16.12.2019
#
#     Parameters:
#     ea       Actual evapotransp
#     ET       Potential evapotransiration (Priestly-Taylor)   
#     sm       Soil moisture content, mm
#     D       Saturation deficit (D>Z) mm
#     M       Total sudsurface water reservoir mm     
#--------------------------------------------------------------------------

function UnsaturatedEvapEB(toSoil, eatemp, sm, M, D, ET)
     
      outx = 0.0
      ea = 0.0    # actual evapotranspiration drawn from sm
      eaS = 0.0   # actual evapotranspiration to be drawn from S (Layers)

#-------------------------------------------------------------------------
#  Actual evapotranspiration
#-------------------------------------------------------------------------
      if (eatemp > -10.0) #Function of temp, areal mean, see Sælthun 1996, p. 9
        #ET is potential evapotranspiration which becomes actual due to deficit, Priestly Taylor
        saturation = (M-D+sm+toSoil)/M
        ea = min(ET, ET*(1-exp(-4*saturation)))
      else      
        ea = 0.0
      end   
      
#------------------------------------------------------------------------
#     Updating sm caused by evapotranspiration
#------------------------------------------------------------------------
      sm = sm-ea
      
      if (sm < 0) # prevents breaching the water balance, must adjust ea      
        ea = ea + sm #sm is negative
        eaS = -sm    #sm is negative and eaS (evapotranspiration to be drawn from S) becomes positive
        sm = 0
      end
      
return sm, ea, eaS
end
        
