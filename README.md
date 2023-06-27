# DDDmodels
## Hydrological models based on the principle of distance distribution dynamics (DDD)
Rainfall-runoff models based on Distance Distribution Dynamics (DDD) have been developed with the aim of reducing the number of free, calibrated model parameters to a minimum. Subsurface- and runoff dynamics are parameterised from catchment distance distributions and estimates of water celerity (subsurface flow, overland flow and streamflow), and not from calibration against runoff. A GIS analysis determines the shape of linear reservoirs from the distance distributions. Celerities, and hence the scale of the linear reservoirs, can be estimated from a recession analysis or from model calibration against recession events. The models for natural catchments (the DDD model) and for urban areas (the DDDUrban) keep track of the moisture input for 10 elevations zones of equal area, so the spatial distribution, accumulation and melt of snow is carried out independently for each elevation zone. Snow melt and evapotranspirations are estimated using an energybalance model, where each of the energy elements are estimated using proxy models. 

The DDD model for natural catchments has been used operationally at the Flood Forecasting Service at the Norwegian Water Resources and Energy Directorate since 2013. 

DDDurban is very much based on the DDD model, but the model calculates rainfall-runoff separately for permeable- and impermable areas. In addition, infiltration capacity is specified in the parameterfile. 

DDD-Design is based on DDDUrban but is a stochastic model intended used for obtaining design values. In this model subsurface saturation state, extreme precipitation, and temporal distribution of extreme precipitation are drawn from an Exponential distribution, a General Pareto distribution (GPD) and a Beta distribution respectively. This model simulates runoff from drawn states and can hence simulate short time series (of length twice the concentration time of the catchment) the sufficient number of times to estimate an extremevalue distribution of runoff.   

The files located at this repository should let you run the models for the sample catchments located in each model directory.  The models have a lot of common subroutines we can be downloaded from the dedicated DDDFunctions directory. There are main Juliascripts in Jupyter notebook format, "Run_DDD* " for each model (the * signifies which model) which loads the necessary subroutines. NOTE that you need, of course, to change all the paths in the "Run.." scripts. In the "Run.." scripts you can set various controls, so that you can run and calibrate the models (for DDD model and DDDUrban). Not all of these controls are fully operational since the development of the models are ongoing projects. Running the DDD model and DDDUrban (kal=0) and calibrating (kal=1) work fine. In each model directory there are "getting started" documents which may or may not be completely updated, but will, nevertheless, help you to get the model running. In addition, you will find example files of parameterfiles, and input- and outputfiles. 

If you do not want to download the entire repository, just a folder, this will do the job: https://download-directory.github.io/

###References
* Skaugen T. and C. Onof, 2014. A rainfall runoff model parameterized form GIS and runoff data. Hydrol. Process. 28, 4529-4542,DOI:10.1002/hyp.9968
* Skaugen, T., I. O. Peerebom and A. Nilsson, 2015. Use of a parsimonious rainfall-runoff model for predicting hydrological response in ungauged basins. Hydrol. Process. 29, 1999-2013, DOI:10.1002/hyp.10315
* Skaugen, T. and Z. Mengistu, 2016. Estimating catchment scale groundwater dynamics from recession analysis- enhanced constraining of hydrological models. Hydrol. Earth. Syst. Sci. 20, 4963-4981, doi: 10.5194/hess-20-4963-2016.
* Skaugen, T. and Weltzien, I. H., 2016. A model for the spatial distribution of snow water equivalent parameterised from the spatial variability of precipitation, The Cryosphere. 10, 1947-1963, doi:10.5194/tc-10_1947_2016
* Skaugen, T., H. Luijting, T. Saloranta, D. Vikhamar-Schuler and K. Müller, 2018. In search of operational snow model structures for the future - comparing four snowmodels for 17 catchments in Norway. Hydrology Research, 49.6, https://doi.org/10.2166/nh.2018.198
* Skaugen, T. D. Lawrence and R. Z. Ortega, 2020. A parameter parsimonious approach for catchment scale urban hydrology – Which processes are important? Journal of Hydrology X, 8, https://doi.org/10.1016/j.hydroa.2020.100060
* Skaugen T., A. E. Stavang, D. Lawrence and K. M. Møen, 2023. Catchment response times– understanding runoff dynamics from catchment distances and celerities. Hydrol. Sciences Journal, https://doi.org/10.1080/02626667.2023.2201449
