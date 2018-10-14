"""
Calculate the probability distribution to be used for identifying the structures
Uses the steinhardt module to efficiently calculate q values.
Takes one argument - name of the input trajectory file.
By default calculates the q4,q6 space. However this can be easily modified by changing the variable qspace
written by srm
"""
from __future__ import print_function
import steinhardt as st
import numpy as np
import sys
import logging
import os
#import matplotlib.pyplot as plt
##############################################################
#nput to modify
##############################################################
qspace = [4,6]
cutoff = 0.0001
rcut   = 3.63
histomin = 0.0
histomax = 0.60
histobins = 1000
histoarea = ((histomax-histomin)/float(histobins))**2
#cutoff <= sum(histo(n)-histo(n-1))/(number of non zero grids)
##############################################################
logger = logging.getLogger(__name__)
handler = logging.FileHandler("runlog.log")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False  


clogger = logging.getLogger('log1')
handler = logging.FileHandler("convergence.dat")
formatter = logging.Formatter('%(message)s')
handler.setFormatter(formatter)
clogger.addHandler(handler)
clogger.setLevel(logging.DEBUG)
clogger.propagate = False  
#set up logger


#get the infile
infile = sys.argv[1]
natoms = int(sys.argv[2])

#checkfor dimension of qspace
if len(qspace)!=2:
    logger.error("only 2d qspace")
    raise SystemExit()

logger.info("infile name is %s"%infile)

#convert the data from infile to systems
#systems = st.traj_to_systems2(infile,natoms=1372,nsteps=100)
#systems = st.traj_to_systems(infile)
systems,files = st.traj_to_systems3(infile,natoms)

#some histo params
xaxis = np.linspace(histomin,histomax,histobins)

#now create an oldhisto
#cannot be normed at this point- we have to manually norm it later
oldhisto,edgex,edgey = np.histogram2d([],[],bins=(xaxis,xaxis))

#arrays for collecting the vals
aq4=[]
aq6=[]

#now start calculating one by one
converged=False
for i in range(len(systems)):
    systems[i].read_particle_file()
    systems[i].set_neighbordistance(rcut)
    systems[i].set_reqd_qs(qspace)
    systems[i].calculate_aq()
    os.remove(files[i])

    #collect the arrays
    #aq4 = np.concatenate((aq4,systems[i].gaqvals(qspace[0])))
    #aq6 = np.concatenate((aq6,systems[i].gaqvals(qspace[1])))
    aq4 = systems[i].gaqvals(qspace[0])
    aq6 = systems[i].gaqvals(qspace[1])
    #now remove te=he systems
    systems[i] = 0.0
    
    #now calculate the newhisto
    newhisto,edgex,edgey = np.histogram2d(aq4,aq6,bins=(xaxis,xaxis))
    newhisto = newhisto + oldhisto
    #norm both old and new histo
    oldhistosum = np.sum(oldhisto)
    newhistosum = np.sum(newhisto)
    #print(oldhistosum)
    #create the normed histos
    if oldhistosum>0:
        normedold = oldhisto/(float(oldhistosum)*histoarea)
    else:
        normedold = oldhisto
    if newhistosum>0:
        normednew = newhisto/(float(newhistosum)*histoarea)
    else:
        normednew = newhisto

    #print type(normednew)
    #find the difference between old and new
    diff = normednew-normedold
    diff = np.abs(diff)*histoarea
    diff = np.sum(diff)
    #diff = newhisto-oldhisto
    #print diff
    #print np.sum(diff)
    #print "old sum is: %f new sum is %f"%(np.sum(normedold),np.sum(normednew))
    #find number of non zero entries in new
    #nonzero = np.count_nonzero(normednew)

    #diffsum = np.sum(np.abs(diff))

    #diff_factor = diffsum/float(nonzero)

    #maxdiff = np.amax(np.abs(diff))
    #test = normednew*histoarea
    #test = np.sum(test)
    #print(test)
    #print the diff factor
    logger.info("diff after adding %d systems is %f "%(i+1,diff))
    print("diff after adding %d systems is %f "%(i+1,diff))
    clogger.info("%d %f "%(i+1,diff))
    #now update the infos
    oldhisto=newhisto

    if (diff<=cutoff):
        logger.info("histograms are converged after %d systems"%i+1)
        converged=True
        break

#now lets test this much
#once done, create histo and save the data
if not converged:
    logger.info("Did not achieve the required convergence cutoff of %f"%cutoff)
#save as np array format
outfile = infile+".histo.npy"
np.save(outfile,normednew)
logger.info("saved npy format as %s"%outfile)

#save in human readable format
#xmesh,ymesh = np.meshgrid(xaxis-1,xaxis-1)
#plt.contourf(xmesh,ymesh,newhisto)
summ = 0
outfile = infile+".histo.dat"
with open(outfile,'w') as fout:
    for i in range(histobins-1):
        for j in range(histobins-1):
            if normednew[i,j]!=0:
                summ+=normednew[i,j]
                fout.write("%f %f %f\n"%(edgex[i],edgey[j],float(normednew[i,j])))

logger.info("saved human readable histo as %s"%outfile)
logger.info("convergence data is logged in convergence.dat")
logger.info("added total of %f points, exiting.."%summ)
