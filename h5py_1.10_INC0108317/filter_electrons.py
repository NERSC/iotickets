import numpy as np
import pandas as pd
import math
#import line_profiler 
#from memory_profiler import profile

eta_range = [0.0, 0.8, 1.3, 2.0, 2.2, 2.3, 2.4]
eff_area = np.array([0.1752, 0.1862, 0.1411, 0.1534, 0.1903, 0.2243, 0.2687])

#@profile
def filter_electrons(df):
    erng=np.searchsorted(eta_range,df.eta)
    np.place(erng, erng>0, erng-1)
    iso = df.chHadIso + np.maximum(0.0, df.gammaIso+df.neuHadIso+(eff_area[erng]*df.rhoIso))
    dEtaIn = abs(df.dEtaIn)
    dPhiIn = abs(df.dPhiIn)
    eoverp = abs(1.0 - df.eoverp)
    d0 = abs(df.d0) 
    dz = abs(df.dz) 
    scEta = abs(df.scEta)
    f_scEta0 = scEta < 1.479
    f_scEta1 = scEta >= 1.479
    filter0 = df.isConv != 1
    filter1 = ((iso < 0.126*df.pt) | (dEtaIn < 0.01520) | (dPhiIn < 0.21600) | (df.sieie < 0.01140)) & f_scEta0
    filter2 = ((iso < 0.144*df.pt) | (dEtaIn < 0.01130) | (dPhiIn < 0.23700) | (df.sieie < 0.03520)) & f_scEta1
    filter3 = ((df.hovere < 0.18100) | (eoverp < (0.20700*df.ecalEnergy)) | (d0 < 0.05640) | (dz < 0.47200) | (df.nMissingHits <= 2)) & f_scEta0
    filter4 = ((df.hovere < 0.11600) | (eoverp < (0.17400*df.ecalEnergy)) | (d0 < 0.22200) | (dz < 0.92100) | (df.nMissingHits <= 3)) & f_scEta1
    f = (filter0 & (filter1 | filter2)) & (filter3 | filter4)
    return df[(abs(df.eta) < 2.5) & (df.pt >= 10) & f] # & (filter0 & (filter1 | filter2)) & (filter3 | filter4)] 


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: ", sys.argv[0], "input-dir group-name"
        sys.exit(-1)

    in_dir = sys.argv[1]
    grp_name = sys.argv[2]
