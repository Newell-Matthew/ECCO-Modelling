#    Function LayerEstimation
#
#------------------------------------------------------------------------
#     Description:  Calculates the total capacity of groundwater reservoir and for its Layers  
#
#     Author: Thomas Skaugen
#     Revised: 16.12.2019
#--------------------------------------------------------------------------

# Dupuit-Forchheimer Subsurface Capacities Redistribution
function RedistributeCapacitiesDupuitForchheimer(ssRes1, NoL, M, k)
    # Create a new array for redistributed capacities
    weighted_capacities = zeros(Float64, NoL)
	
	#Create an array for weights
	w = zeros(Float64, NoL)
	
	#Create weights for subsurface capacities
	w[NoL] = 1/k[NoL]
    w[4] = 1/k[4]
    w[3] = 1/k[3]
    w[2] = 1/k[2]
	
    # Apply Dupuit-Forchheimer distribution principles
    weighted_capacities[NoL] = M *(w[NoL]/sum(w[2:NoL]))
    weighted_capacities[4] = M *(w[4]/sum(w[2:NoL]))
    weighted_capacities[3] = M *(w[3]/sum(w[2:NoL]))
    weighted_capacities[2] = M *(w[2]/sum(w[2:NoL]))

    # Update original ssRes1 with the weghted capacities for subsurface layers
    ssRes1[2:NoL] .= weighted_capacities[2:NoL]

    return ssRes1
end




function LayerEstimation(GshInt,GscInt,Timeresinsec,maxDl,midDL, MAD, area2, NoL, gtcel, k)

# using Distributions
# include("F:\\HB\\HB-modellering\\DDDtestbenk\\DDD_Julia\\DDDFunctions\\SingleUH.jl")

 # GshInt: scalar float = 1.5
 # GscInt: scalar float = 0.22 
 # midDL: scalar float = 129.84
 # maxDl: scalar float = 538
 # MAD: scalar float = 2.35
 # Timeresinsec: scalar float = 86400
 # NoL: scalar integer = 5
 # area2: scalar float = 2300000 
 # gtcel: scalar float = 0.99
 
mLam = GshInt*GscInt
varLam = GshInt*(GscInt)^2                           #Yevjevich p.145
meanIntk = mLam*midDL/Timeresinsec                   #mean celerity estimated through Integrated Celerity
antBox = Int(trunc(maxDl/(meanIntk*Timeresinsec)))+1 #Temporal length UH_MAD
UH_MAD = zeros(Float64,antBox)
sRes = zeros(Float64,antBox) # saturation sum

#Unit hydrograph for MAD
UH_MAD = SingleUH(meanIntk,Timeresinsec, midDL, maxDl, 0)

StSt = (1000*MAD*Timeresinsec)/(area2)         # Steady state Input eq. output in mm
sRes[1] = 0
sRes[2:antBox] .= StSt.*UH_MAD[2:antBox]

for i in 3: antBox
  sRes[i:antBox] .= sRes[i:antBox] + StSt.*UH_MAD[i:antBox]
end

mRes = sum(sRes)
Fact = mLam/mRes
stdRes = (varLam/Fact^2)^0.5                   # see Haan p.51

GshRes = mRes^2/stdRes^2
GscRes = stdRes^2/mRes

MLev = [1/(NoL-1):1/(NoL-1):1.0;]              # (sequence)Quantiles  to calculate reservoir levels [0.1:0.1:0.9;]

MLev[NoL-1] = gtcel                            # quantile for start overland flow
Res_prob = zeros(Float64,(NoL-1))
Magkap = zeros(Float64,NoL)
g = Gamma(GshRes,GscRes) 
#calculates the reservoir levels associated with quantiles. Mean is GshRes*GscRes
Res_prob .= quantile.(g,MLev)

#Capasity of Layers
ssRes1 = zeros(Float64,NoL)
ssRes1[1] = 2000                               # capacity of overland flow level

for i in 2:(NoL-1)
  ssRes1[i] = Res_prob[NoL-i+1]-Res_prob[(NoL-i)]
end

ssRes1[NoL] = Res_prob[1]                     # capasity for the first slowest level         
                       
M = Res_prob[(NoL-1)]                           # Total groundwater reservoir

# Apply Dupuit-Forchheimer redistribution
ssRes1 = RedistributeCapacitiesDupuitForchheimer(ssRes1, NoL, M, k)

# Update Magkap with the new capacities
Magkap = ssRes1                                 # capasity for Layers

return Magkap, M
end
