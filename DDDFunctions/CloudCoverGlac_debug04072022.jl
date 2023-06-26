 

function CloudCover(regn,Timeresinsec)
 
 grense = 0.98 
 pgrense =  (Timeresinsec)/3600 # 1.0 # precip limit for  full Cloud cover
 pgrense2 = pgrense/24# precip limit for high/low Cloud cover
   
#Estimating Cloudcover 
 if regn > pgrense # 1 for 1 time funker bra, 3 for 3 timer?
   Cl = 1.0
   Wa = 1.0
 end      
    
 if (regn > pgrense2 && regn <=pgrense)
   Cl = rand(Uniform(grense,0.99))
   Wa = (Cl+1.667)/2.667
 end
    
 if (regn <= pgrense2)
   Cl = rand(Uniform(0.01,grense)) 
   Wa = (Cl+1.667)/2.667
 end       

 #if (Cl >= grense)
 #  Wa = (Cl+1.667)/2.667       # parametrization of Relative Humidity from Herrero and Polo(HESS, 2012) 
 #end     

 #if (Cl < grense)
 #   Wa = 0.8349                  # mean values from Filefjell
 #end       

 #if (Cl < (grense - 0.1))
 #   Wa = (Cl+1.167)/2.667  # parametrization of Relative Humidity from Herrero and Polo(HESS, 2012)
 #end

 if (Wa > 1.0)
    Wa = 1.0# 0.8349  #1.0
 end
 #println(Cl, " ", Wa)  

 return Cl, Wa

end
