"""
Includes routines for creation of structural histograms in qspace
and identification of structures based on this.

"""
from __future__ import print_function
import steinhardt as st
import numpy as np
import sys
import logging
import os

def create_histograms(systems,qspace,histovars,cutoff,textfile=False,outfile=""):
    """
    Calculate the probability distributions for identifying the structures.

    Parameters
    ----------

    systems : array like
        araray of System objects with atomic positions and simulation box
        read in.
        All parameters required for calculation of q values need to be set
        beforehand. See-
        steinhardt.System.set_neighbordistance()
        steinhardt.System.set_reqd_qs([x,y]) : only 2 d qspace can be used
    
    qspace : array of dim 2
        The qspace on which the histogram is to be calculated. It should be 
        two dimensional.

    histovars : array of dim 3
        Array of three quantities required for calculation of histograms, namely,
        [histomin, histomax, histobins]. histomin is the minimum value of histogram,
        histomax is the maximum and histobins are the number of bins required.

    cutoff : float
        converge threshold for the histograms.

    textfile : bool
        If true, a human readable text file for the histograms is also written.
        Otherwise only a .npy file is written.

    outfile : string
        If provided, the output histogram would be written in this filename.
        Otherwise the same name as input file with .histo.npy / .histo.txt 
        extension is used.
    
    Returns
    -------

    converged : bool
        True if the histograms are converged, false otherwise

    pdiff : float
        The difference between the actual value and required convergence.

    """
    #some histo params
    xaxis = np.linspace(histovars[0],histovars[1],histovars[2])
    histoarea = ((histovars[1]-histovars[0])/float(histovars[2]))**2

    #now create an oldhisto
    #cannot be normed at this point- we have to manually norm it later
    oldhisto,edgex,edgey = np.histogram2d([],[],bins=(xaxis,xaxis))

    #arrays for collecting the vals
    aq4=[]
    aq6=[]


    converged=False
    for i in range(len(systems)):
        systems[i].set_reqd_qs(qspace)
        systems[i].calculate_aq()
        
        aq4 = systems[i].gaqvals(qspace[0])
        aq6 = systems[i].gaqvals(qspace[1])
        
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
        
        #now update the infos
        oldhisto=newhisto

        if (diff<=cutoff):
            converged=True
            break

    pdiff = diff - cutoff

    #save as np array format
    if outfile == "":
        outfile = "structure.histo.npy"
    np.save(outfile,normednew)
    
    if textfile:
        summ = 0
        outfile = "structure.histo.dat"
        with open(outfile,'w') as fout:
            for i in range(histobins-1):
                for j in range(histobins-1):
                    if normednew[i,j]!=0:
                        summ+=normednew[i,j]
                        fout.write("%f %f %f\n"%(edgex[i],edgey[j],float(normednew[i,j])))

    
    return converged, pdiff

    
def find_structure_probs(systems, qspace, histovars,cutoff,histofiles,peratom=False,normed=True, textfile=False, outfile="", assigntosystems=False):
    """
    Calculate the probability distributions for identifying the structures.

    Parameters
    ----------

    systems : array like
        araray of System objects with atomic positions and simulation box
        read in.
        All parameters required for calculation of q values need to be set
        beforehand. See-
        steinhardt.System.set_neighbordistance()
        steinhardt.System.set_reqd_qs([x,y]) : only 2 d qspace can be used

    qspace : array of dim 2
        The qspace on which the histogram is to be calculated. It should be 
        two dimensional.

    histovars : array of dim 3
        Array of three quantities required for calculation of histograms, namely,
        [histomin, histomax, histobins]. histomin is the minimum value of histogram,
        histomax is the maximum and histobins are the number of bins required.

    cutoff : float
        A value below which the histogram entry is considered to be zero. 

    histofiles : array like
        A list of file names which contain the histogram information in the .npy 
        format.

    peratom : bool
        If True, it will print out the per atom information for each time step
        False otherwise

    normed : bool
        If True the structural files will be written out with normalised probabilities
        False otherwise
    
    assigntosystems : bool
        If True the structure is assigned to individual atoms of the system and 
        systems are returned. The values of the structure would be from 1 to n, in order
        of the histogram files provided. If it does not ebelong to all structure, 
        a value of 0 will be given.
        This is not yet implemented : Should not use.

    """
    xaxis = np.linspace(histovars[0],histovars[1],histovars[2])

    #now open the .npy files in mmap mode
    probmaps = []
    for file in histofiles:
        probmaps.append(np.load(file))
    nstructs = len(probmaps)

    alltimedata = []
    alltimeperatomdata = []
    for count,sys in enumerate(systems):
        sys.set_reqd_qs(qspace)
        sys.calculate_aq()
        #get q4 and q6 arrays
        aq4=sys.gaqvals(qspace[0])
        aq6=sys.gaqvals(qspace[1])
        snapatomdata = []
        #get atom ids for this system
        atomids = [atom.gid() for atom in systems[count].gallatoms()]

        netstruct = np.zeros(nstructs+1)

        for i in range(len(aq4)):
            aq4bin = np.digitize([aq4[i]],bins=xaxis,right=True)[0]-1
            aq6bin = np.digitize([aq6[i]],bins=xaxis,right=True)[0]-1

            sprobs = np.zeros(nstructs+1)
            for j in range(nstructs):
                sprobs[j] =  probmaps[j][aq4bin, aq6bin]
            #now add zero for the unknown structrs
            sprobs[-1] = 0.0

            probsum = np.sum(sprobs)
            #account for zero probsum
            if probsum < cutoff:
                sprobs[-1] = 1.0 - probsum
                probsum = 1

            if normed:
                sprobs /= probsum

            #write out the per atom values here if reqd
            ###########################################
            ###########################################
            if peratom:
                sprobsmod = np.insert(sprobs,0,atomids[i])
                snapatomdata.append(sprobsmod)

            netstruct += sprobs

        alltimeperatomdata.append(snapatomdata)
        #now norm netstructs
        netstruct /= float(len(aq4))
        alltimedata.append(netstruct)

    if outfile == "":
        outfile = ".".join("structure",)
    np.save(outfile,alltimedata)

    if peratom:
        atomoutfile = ".".join([outfile,"atom"])
        np.save(atomoutfile,alltimeperatomdata)












