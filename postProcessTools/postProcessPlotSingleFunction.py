import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt

def hemeCollateFunction(path, stepSize, DX, lets):
    IOlets = DX*lets

    IOP = [[] for x in range(np.shape(IOlets)[0])]
    IOV = [[] for x in range(np.shape(IOlets)[0])]

    # steps gridX gridY gridZ velX velY velZ pressure

    print("starting to process", path)
    f = open(path, "r")
        
    stepCounter = 0

    outs=[]
    ioV = [[] for x in range(np.shape(IOlets)[0])]

    notStart = False

    for i in f:  
        i = np.float_(i.split())
        i = np.concatenate((i, np.linalg.norm(i[4:7])), axis=None)
        
        if len(i)==1:
            # gap in iteration output with space
            continue
        if int(i[0]) != stepCounter*stepSize:
            #print(stepCounter)
            stepCounter = stepCounter+1
            
            if notStart:       
                outs = np.array(outs)
                ioV = np.array(ioV) 

                for io in range(np.shape(IOlets)[0]):
                    IOP[io] = IOP[io]+[np.max(ioV[io], axis=0)[7], np.average(ioV[io], axis=0)[7], np.min(ioV[io], axis=0)[7]]
                    IOV[io] = IOV[io]+[np.max(ioV[io], axis=0)[8], np.average(ioV[io], axis=0)[8], np.min(ioV[io], axis=0)[8]]

            notStart=True
            
            outs=[]
            ioV = [[] for x in range(np.shape(IOlets)[0])]

        outs.append(i)

        dist = 1e10
        minIO = -1
        for io in range(np.shape(IOlets)[0]):
            if np.linalg.norm(i[1:4] - IOlets[io])<dist:
                dist = np.linalg.norm(i[1:4] - IOlets[io])
                minIO = io

        ioV[minIO].append(i)
          
    f.close()

    return np.array(IOP), np.array(IOV)
