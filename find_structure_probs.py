"""
Read in the structure histograms of different structure -BCC,
FCC,HCP and LQD.  
"""
from __future__ import print_function
import steinhardt as st
import numpy as np
import sys
import logging

##############################################################
#input to modify
##############################################################
bccfile = 'traj.bcc.histo.npy'
fccfile = 'traj.fcc.histo.npy'
hcpfile = 'traj.hcp.histo.npy'
lqdfile = 'traj.lqd.histo.npy'
qspace = [4,6]
cutoff = 0.00001
rcut   = 3.63
histomin = 0.0
histomax = 0.60
histobins = 1000
##############################################################

#loggers
logger = logging.getLogger(__name__)
handler = logging.FileHandler("runlog.log")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False  

clogger = logging.getLogger('log1')
handler = logging.FileHandler("structure.dat")
formatter = logging.Formatter('%(message)s')
handler.setFormatter(formatter)
clogger.addHandler(handler)
clogger.setLevel(logging.DEBUG)
clogger.propagate = False 

clogger2 = logging.getLogger('log2')
handler = logging.FileHandler("structure-non-normalised.dat")
formatter = logging.Formatter('%(message)s')
handler.setFormatter(formatter)
clogger2.addHandler(handler)
clogger2.setLevel(logging.DEBUG)
clogger2.propagate = False 

#get the infile
filename = sys.argv[1]


logger.info("infile name is %s"%filename)

#convert the data from infile to systems
#systems = st.traj_to_systems(infile)
#some histo params
xaxis = np.linspace(histomin,histomax,histobins)

#read in the histofiles of various structures
bcchisto = np.load(bccfile,mmap_mode='r')
fcchisto = np.load(fccfile,mmap_mode='r')
hcphisto = np.load(hcpfile,mmap_mode='r')
lqdhisto = np.load(lqdfile,mmap_mode='r')

#convert the data from infile to systems
systems = st.traj_to_systems(filename)
            
for count,sys in enumerate(systems):
    sys.set_neighbordistance(rcut)
    sys.set_reqd_qs(qspace)
    sys.calculate_aq()
    #get q4 and q6 arrays
    aq4=sys.gaqvals(qspace[0])
    aq6=sys.gaqvals(qspace[1])
    #get the bins of all the values to be binnerd
    #aqs = np.stack((aq4,aq6),axis=-1)
    #print type(aqs)
    #print aqs.shape
    #now find the bins of each values
    netbcc=0
    netfcc=0
    nethcp=0
    netlqd=0
    netudf=0

    for i in range(len(aq4)):
        aq4bin = np.digitize([aq4[i]],bins=xaxis,right=True)[0]-1
        aq6bin = np.digitize([aq6[i]],bins=xaxis,right=True)[0]-1

        bccprob = bcchisto[aq4bin,aq6bin]
        fccprob = fcchisto[aq4bin,aq6bin]
        hcpprob = hcphisto[aq4bin,aq6bin]
        lqdprob = lqdhisto[aq4bin,aq6bin]
        udfprob = 0.0
        #we need to correct the probabilities now
        probsum = bccprob+fccprob+hcpprob+lqdprob

        #account for zero or less probsum
        if probsum<=cutoff:
            udfprob=1.0-probsum
            probsum=1

        clogger2.info("%f %f %f %f %f"%(bccprob,fccprob,hcpprob,lqdprob,udfprob))
        bccprob/=probsum
        fccprob/=probsum
        hcpprob/=probsum
        lqdprob/=probsum
        #print the values
        #print(bccprob)
        clogger.info("%f %f %f %f %f"%(bccprob,fccprob,hcpprob,lqdprob,udfprob))

        netbcc+=bccprob
        netfcc+=fccprob
        nethcp+=hcpprob
        netlqd+=lqdprob
        netudf+=udfprob

    netbcc/=float(len(aq4))
    netfcc/=float(len(aq4))
    nethcp/=float(len(aq4))
    netlqd/=float(len(aq4))
    netudf/=float(len(aq4))

    print("%d %f %f %f %f %f"%(count+1,netbcc,netfcc,nethcp,netlqd,netudf))

    #newhisto,edgex,edgey = np.histogram2d(aq4,aq6,bins=(xaxis,xaxis))
    #now get all nonzero bins of this histo
    #newhistoint = (newhisto>0).astype(int)
    #print newhistoint
    #bccprob = newhistoint*bcchisto
    #fccprob = newhistoint*fcchisto
    #hcpprob = newhistoint*hcphisto
    #lqdprob = newhistoint*lqdhisto
    #finally print each prob and we are done?

