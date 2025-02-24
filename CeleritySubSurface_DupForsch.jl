#    Function CeleritySubSurface
#
#------------------------------------------------------------------------
#     Description:  Calculates Sub surface celerities based on recession analysis (lower case Lambdas). 
#
#     Author: Thomas Skaugen
#     Revised: 17.11.2017
#--------------------------------------------------------------------------


function CeleritySubSurface(NoL, Gshape, Gscale, midDL, Timeresinsec) 

# using Distributions   
#  NoL = 5
#  Gshape = 2.1
#  Gscale = 0.0022 
#  midDL = 129.84 
#  Timeresinsec = 86400

a = 3
k = zeros(Float64,NoL)                   #celerity of subsurface (and overland) flow
probvec = zeros(Float64,NoL)            #all leves and overland flow level
meanVel = Gshape*Gscale*midDL/Timeresinsec

println("Hei ", a)
dp = 1/(NoL-1)                              #Overland flow level (nol=1], fixed celerity and extremely high capacity(2000 mm)

for i in reverse(1:(NoL-1))
    probvec[i+1] = i*dp - dp/2                #celerities are estimated at center of level, hence dp/2
end
probvec = 1 .- probvec
probvec[1] = 0.99                            #Quantile in celerity distribution for overland flow fixed at 0.99

g = Gamma(Gshape, Gscale)    
k[1:NoL] = quantile.(g,probvec[1:NoL])*midDL/Timeresinsec
# celerities k[5:3] is given the same celerity
# k[4] is  twice that of k[5:3]
# k[1] is the overland flow but assigned outside of this subroutine
##Dupuit Forschheimer 
k[NoL] = 4*meanVel/(3 + a)
k[4] = 4*meanVel/(3 + a)
k[3] = 4*meanVel/(3 + a)
k[2] = k[3]*a

return k                        #Celerities [m/s] 
end
