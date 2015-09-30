#########################
## 7 # ELASTIC NETWORK ##  -> @ELN <-
#########################
import math
import FUNC 

## ELASTIC NETWORK ##

# Only the decay function is defined here, the network 
# itself is set up through the Topology class

# The function to determine the decay scaling factor for the elastic network 
# force constant, based on the distance and the parameters provided.
# This function is very versatile and can be fitted to most commonly used 
# profiles, including a straight line (rate=0)
def decayFunction(distance,shift,rate,power):
    return math.exp(-rate*math.pow(distance-shift,power))

def rubberBands(atomList,lowerBound,upperBound,decayFactor,decayPower,forceConstant,minimumForce):
    out = []
    u2  = upperBound**2
    while len(atomList) > 3:
        bi,xi = atomList.pop(0)
        for bj,xj in atomList[2:]:
            # Mind the nm/A conversion -- This has to be standardized! Global use of nm?
            d2 = FUNC.distance2(xi,xj)/100
            
            if d2 < u2:
                dij  = math.sqrt(d2)
                fscl = decayFunction(dij,lowerBound,decayFactor,decayPower)
                if fscl*forceConstant > minimumForce:
                    out.append({"atoms":(bi,bj),"parameters": (dij,"RUBBER_FC*%f"%fscl)})
    return out



