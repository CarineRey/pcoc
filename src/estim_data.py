#  estim_data.py
#  
#  Copyright 2017 Carine Rey <carine.rey@ens-lyon.fr>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import math
import os
import re
from Bio import AlignIO
from ete3 import Tree
import logging
logger = logging.getLogger("pcoc.estim_data")

DEV = False

def read_info(nf, n_sites):
    f=open(nf, "r")
    ll=[]
    i=1
    #logger.debug(nf)
    for l in f:
        if l[0]!="S":
            el = l.split()
            pos = int(el[0].replace("[","").replace("]",""))
            lnL = float(el[3])
            while pos != i:
                ll.append(None)
                i+=1
            ll.append(lnL)
            i+=1
    pos = n_sites
    while pos >= i:
        ll.append(None)
        i+=1

    f.close()
    return ll

def read_info_exception(f, n_sites, O, NorP, C1C1_C1C2):
    res = [None]*n_sites
    if os.path.exists(f):
        #logger.debug("%s exists", f)
        pass
    else:
        logger.debug("%s does not exist", f)
        pass
    try:
        res = read_info(f, n_sites)
    except:
        if DEV:
            logger.error("Pb in read_info(%s)", f)
        else:
            if O == "withOneChange" and C1C1_C1C2 == 11:
                pass
            elif O == "noOneChange" and C1C1_C1C2 == 12:
                pass
            else:
                logger.error("Pb in read_info(%s)", f)
                #logger.warning("f: %s, O: %s, NorP: %s, C1C1_C1C2: %s", f, O, NorP, C1C1_C1C2)
    #logger.debug("O: %s, k: %s, NorP: %s, res: %s", O, NorP, C1C1_C1C2, res)
    return res
    
def prob_ap(x,y):
    if x and y:
        return math.exp(x-y)/(math.exp(x-y)+1)
    else:
        return 0

def prob_ap_no_log(X,Y):
    if X and Y:
        return X/(X+Y)
    else:
        return 0

def likelihood_mean(x):
    if x == [None]:
        return None
    l = []
    for z in x:
        if z:
            l.append(math.exp(z))
    if l:
        return(float(sum(l))/len(l))
    else:
        return None

def likelihood_max(x):
    if x == [None]:
        return None
    return(max(x))
    

def sensspec(tab):
    #logger.debug ("P: %s ( TP %s + FN %s)", tab["TP"] + tab["FN"],  tab["TP"], tab["FN"])
    #logger.debug ("N: %s ( TN %s + FP %s)", tab["TN"] + tab["FP"],  tab["TN"], tab["FP"])
    if tab["TP"]+tab["FN"]==0:
      ss="NA"
    else:
      ss=tab["TP"]/(tab["TP"]+tab["FN"])
    
    if tab["TN"]+tab["FP"]==0:
      sp="NA"
    else:
      sp = tab["TN"]/(tab["TN"]+tab["FP"])
    
    mcc_nom = tab["TP"]*tab["TN"] - tab["FP"]*tab["FN"]
    mcc_den =  (tab["TP"]+tab["FP"])*(tab["TP"]+tab["FN"])*(tab["TN"]+tab["FP"])*(tab["TN"]+tab["FN"])
    
    if mcc_den == 0:
        mcc = 0
    else:
        mcc = mcc_nom /  math.sqrt(mcc_den)
    
    return [ss,sp,mcc]

def get_method(P):
    if  re.search("_X_OXY$", P):
        res = "PCOC"
    elif re.search("_X_OX$", P):
        res = "OC"
    elif re.search("_X_XY$", P):
        res = "PC"
    elif re.search("_X_CX$", P):
        res = "Topological"
    else:
        res = "PCOC_var"
    return res

########################################################################
##                       PCOC model                                   ##
########################################################################

def dico_typechg_new(C1,C2, N, repest, name_AC, tree="", NbCat_Est=10, n_sites=1000, set_e1e2=[], lseuil=[0.7,0.80,0.85,0.90,0.95,0.99], ID="", dist_C1_C2=None):
### dico des proba a posteriori qd nb clades avec changement = nbe clades convergents

    ## 11 : sequences generees par C1-C1
    ## 12 : sequences generees par C1-C2  
    simu_i = name_AC.split("_")[1]

    lsites = range(n_sites)
    
    bilan = {11:{"pos":lsites}, 12:{"pos":lsites}}
    P_l = []
    test1_l = ["X_XY","OX_OXY","XY_OXY","X_OXY","X_OX"]
    test2_l = ["p_max", "p_mean", "p_opt"]
    
    for test1 in test1_l:
        for test2 in test2_l:
            P = test2 + "_" + test1
            P_l.append(P)
            for k in [11,12]:
                bilan[k][P]=[]

    if not DEV:
        test1_l = ["X_XY","X_OXY","X_OX"]
        test2_l = ["p_mean"]
        P_l = []
        for test1 in test1_l:
            for test2 in test2_l:
                P = test2 + "_" + test1
                P_l.append(P)
        
    for a in ["X", "OX", "XY", "OXY"]:
        for b in ["opt_", "mean_", "max_"]:
            for k in [11,12]:
                bilan[k][b+a] = []

    dli= {"withOneChange" : {11:{}, 12:{}}, "noOneChange" :{11:{}, 12:{}}}
    for O in ["withOneChange","noOneChange"]: 
        infosC1C1C1C1OPT_file = "%s/Scenario_%s_A%s_C%s_%s_%s_opt_%s.infos"%(repest,simu_i,C1,C1,C1,C1,O)
        infosC1C1C1C2OPT_file = "%s/Scenario_%s_A%s_C%s_%s_%s_opt_%s.infos"%(repest,simu_i,C1,C1,C1,C2,O)
        infosC1C2C1C1OPT_file = "%s/Scenario_%s_A%s_C%s_%s_%s_opt_%s.infos"%(repest,simu_i,C1,C2,C1,C1,O)
        infosC1C2C1C2OPT_file = "%s/Scenario_%s_A%s_C%s_%s_%s_opt_%s.infos"%(repest,simu_i,C1,C2,C1,C2,O)

        ### read infos
        dli[O][11][12] = read_info_exception(infosC1C1C1C2OPT_file, n_sites, O, 11, 12)
        dli[O][11][11] = read_info_exception(infosC1C1C1C1OPT_file, n_sites, O, 11, 11)
        dli[O][12][12] = read_info_exception(infosC1C2C1C2OPT_file, n_sites, O, 12, 12)
        dli[O][12][11] = read_info_exception(infosC1C2C1C1OPT_file, n_sites, O, 12, 11)

    for s in lsites:
        for k in [11,12]:
            # get all likelihood:
            opt_X   = dli["noOneChange"][k][11][s]
            bilan[k]["opt_X"].append(opt_X)
            opt_XY  = dli["noOneChange"][k][12][s]
            bilan[k]["opt_XY"].append(opt_XY)
            opt_OX  = dli["withOneChange"][k][11][s]
            bilan[k]["opt_OX"].append(opt_OX)
            opt_OXY = dli["withOneChange"][k][12][s]
            bilan[k]["opt_OXY"].append(opt_OXY)
            
            bilan[k]["p_opt_X_XY"].append(   prob_ap(opt_XY,opt_X)   )
            bilan[k]["p_opt_OX_OXY"].append( prob_ap(opt_OXY,opt_OX) )
            bilan[k]["p_opt_X_OXY"].append(  prob_ap(opt_OXY,opt_X)  )
            bilan[k]["p_opt_X_OX"].append(   prob_ap(opt_OX,opt_X)   )
            bilan[k]["p_opt_XY_OXY"].append( prob_ap(opt_OXY,opt_XY) )

    if not set_e1e2:
        for e1 in range(1, (NbCat_Est+1)):
            for e2 in range(1, (NbCat_Est+1)):
                if e2 != e1:
                   set_e1e2.append((e1,e2))
                   
    X_l   = []
    XY_l  = []
    OX_l  = []
    OXY_l = []

    for s in lsites:
        X_l.append({11:[], 12:[]})
        XY_l.append({11:[], 12:[]})
        OX_l.append({11:[], 12:[]})
        OXY_l.append({11:[], 12:[]})
            
    for (e1,e2) in set_e1e2:
        if e2 == e1:
                continue
        dli = {"withOneChange" : {11:{}, 12:{}}, "noOneChange" :{11:{}, 12:{}}}
        for O in ["withOneChange","noOneChange"]: 
            #logger.debug ("Count e1: %s e2: %s O: %s", e1, e2,O)
            infosC1C1E1E1_file = "%s/Scenario_%s_A%s_C%s_%s_%s_%s.infos"%(repest,simu_i,C1,C1,e1,e1,O)
            infosC1C1E1E2_file = "%s/Scenario_%s_A%s_C%s_%s_%s_%s.infos"%(repest,simu_i,C1,C1,e1,e2,O)
            infosC1C2E1E1_file = "%s/Scenario_%s_A%s_C%s_%s_%s_%s.infos"%(repest,simu_i,C1,C2,e1,e1,O)
            infosC1C2E1E2_file = "%s/Scenario_%s_A%s_C%s_%s_%s_%s.infos"%(repest,simu_i,C1,C2,e1,e2,O)

            ### read infos
            dli[O][11][12] = read_info_exception(infosC1C1E1E2_file, n_sites, O, 11, 12)
            dli[O][11][11] = read_info_exception(infosC1C1E1E1_file, n_sites, O, 11, 11)
            dli[O][12][12] = read_info_exception(infosC1C2E1E2_file, n_sites, O, 12, 12)
            dli[O][12][11] = read_info_exception(infosC1C2E1E1_file, n_sites, O, 12, 11)


        for s in lsites:
            for k in [11,12]:
               # get all likelihood:
               X   = dli["noOneChange"][k][11][s]
               XY  = dli["noOneChange"][k][12][s]
               OX  = dli["withOneChange"][k][11][s]
               OXY = dli["withOneChange"][k][12][s]

               X_l  [s][k].append(X)
               XY_l [s][k].append(XY)
               OX_l [s][k].append(OX)
               OXY_l[s][k].append(OXY)

    # Bilan
    for s in lsites:
        for k in [11,12]:
            # p_max:
            max_X   = likelihood_max(X_l  [s][k])
            bilan[k]["max_X"].append(max_X)
            max_XY  = likelihood_max(XY_l [s][k])
            bilan[k]["max_XY"].append(max_XY)
            max_OX  = likelihood_max(OX_l [s][k])
            bilan[k]["max_OX"].append(max_OX)
            max_OXY = likelihood_max(OXY_l[s][k])
            bilan[k]["max_OXY"].append(max_OXY)
            
            bilan[k]["p_max_X_XY"].append(   prob_ap(max_XY,max_X)   )
            bilan[k]["p_max_OX_OXY"].append( prob_ap(max_OXY,max_OX) )
            bilan[k]["p_max_X_OXY"].append(  prob_ap(max_OXY,max_X)  )
            bilan[k]["p_max_X_OX"].append(   prob_ap(max_OX,max_X)   )
            bilan[k]["p_max_XY_OXY"].append( prob_ap(max_OXY,max_XY) )
            
            # p_mean:
            mean_X   = likelihood_mean(X_l  [s][k])
            bilan[k]["mean_X"].append(mean_X)
            mean_XY  = likelihood_mean(XY_l [s][k])
            bilan[k]["mean_XY"].append(mean_XY)
            mean_OX  = likelihood_mean(OX_l [s][k])
            bilan[k]["mean_OX"].append(mean_OX)
            mean_OXY = likelihood_mean(OXY_l[s][k])
            bilan[k]["mean_OXY"].append(mean_OXY)
            
            bilan[k]["p_mean_X_XY"].append(   prob_ap_no_log(mean_XY, mean_X)  )
            bilan[k]["p_mean_OX_OXY"].append( prob_ap_no_log(mean_OXY,mean_OX) )
            bilan[k]["p_mean_X_OXY"].append(  prob_ap_no_log(mean_OXY,mean_X)  )
            bilan[k]["p_mean_X_OX"].append(   prob_ap_no_log(mean_OX, mean_X)  )
            bilan[k]["p_mean_XY_OXY"].append( prob_ap_no_log(mean_OXY,mean_XY) )
            
            #logger.debug("mean_X: %s max_X:%s", bilan[k]["mean_X"][s], bilan[k]["max_X"][s])
            #logger.debug("mean_XY: %s max_XY:%s", bilan[k]["mean_XY"][s], bilan[k]["max_XY"][s])
            #logger.debug("p_mean_X_XY: %s p_max_X_XY:%s", bilan[k]["p_mean_X_XY"][s], bilan[k]["p_max_X_XY"][s])
            #logger.debug("mean_X: %s", bilan[k]["mean_X"][s])
            #logger.debug("mean_XY: %s", bilan[k]["mean_XY"][s])
            #logger.debug("mean_OXY: %s", bilan[k]["mean_OXY"][s])
            #logger.debug("mean_OX: %s", bilan[k]["mean_OX"][s])
            #logger.debug("k: %s, mean_X: %s, mean_XY: %s, p_mean_X_XY: %s", k, bilan[k]["mean_X"][s],bilan[k]["mean_XY"][s], bilan[k]["p_mean_X_XY"][s])
            #logger.debug("k: %s, mean_X: %s, mean_OX: %s, p_mean_X_OX: %s", k, bilan[k]["mean_X"][s],bilan[k]["mean_OX"][s], bilan[k]["p_mean_X_OX"][s])
            #logger.debug("p_mean_X_XY: %s", bilan[k]["p_mean_X_XY"][s])
            #logger.debug("p_mean_X_OX: %s", bilan[k]["p_mean_X_OX"][s])
            #logger.debug("p_mean_X_OXY: %s", bilan[k]["p_mean_X_OXY"][s])
            #logger.debug("p_mean_OX_OXY: %s p_max_OX_OXY:%s", bilan[k]["p_mean_OX_OXY"][s], bilan[k]["p_max_OX_OXY"][s])
    # Calcul (VP, FP, FN, VN)
    res = []
    for P in P_l:
        for seuil in lseuil:
            tab={}
           
            tab["PosteriorProbabilityType"]= P
            tab["TP"]=float(len([s for s in bilan[12][P] if s >= seuil]))
            tab["FN"]=float(len([s for s in bilan[12][P] if s < seuil]))
            tab["FP"]=float(len([s for s in bilan[11][P] if s >= seuil]))
            tab["TN"]=float(len([s for s in bilan[11][P] if s < seuil]))
        
            tab["C1"] = C1
            tab["C2"] = C2
            tab["SimuCoupleID"] = "A%s_C%s" %(C1, C2)
            tab["DistanceSimuCouple"] = dist_C1_C2["C" + str(C1)]["C" + str(C2)]
            tab["NumberOfSites"] = n_sites
            tab["NumberOfConvergentEvents"] = N
            tab["Threshold"] = seuil
            tab["Method"] = get_method(P)
            tab["ScenarioID"] = "Scenario_%s" %simu_i
            if tree:
                tab["InputTree"] = tree
            if ID:
                tab["RunID"] = ID
            
            tab["Sensitivity"], tab["Specificity"], tab["MCC"] = sensspec(tab)
            res.append(tab)
    
    return res, bilan

def dico_typechg_het_det(N,repest,repseq,ali_filename, n_sites, tree="", NbCat_Est = 10, set_e1e2 = [], lseuil = [0.7,0.80,0.85,0.90,0.95,0.99], ID= ""):
### dico des proba a posteriori qd nb clades avec changement = nbe clades convergents

    lsites = range(n_sites)

    bilan = {11:{"pos":lsites}, 12:{"pos":lsites}}
    P_l = []
    test1_l = ["X_XY","OX_OXY","XY_OXY","X_OXY","X_OX"]
    test2_l = ["p_max", "p_mean"]

    for test1 in test1_l:
        for test2 in test2_l:
            P = test2 + "_" + test1
            P_l.append(P)
            for k in [11,12]:
                bilan[k][P]=[]
    for a in ["X", "OX", "XY", "OXY"]:
        for b in ["mean_", "max_"]:
            for k in [11,12]:
                bilan[k][b+a] = []
                
    if not DEV:
        test1_l = ["X_XY","X_OXY","X_OX"]
        test2_l = ["p_mean"]
        P_l = []
        for test1 in test1_l:
            for test2 in test2_l:
                P = test2 + "_" + test1
                P_l.append(P)
    
    if not set_e1e2:
        for e1 in range(1, (NbCat_Est+1)):
            for e2 in range(1, (NbCat_Est+1)):
                if e2 != e1:
                    set_e1e2.append((e1, e2))
    X_l   = []
    XY_l  = []
    OX_l  = []
    OXY_l = []
    
    for s in lsites:
        X_l.append({11:[], 12:[]})
        XY_l.append({11:[], 12:[]})
        OX_l.append({11:[], 12:[]})
        OXY_l.append({11:[], 12:[]})
    for (e1,e2) in set_e1e2:
        if e2 == e1:
            continue
        dli = {"withOneChange" : {11:{}, 12:{}}, "noOneChange" :{11:{}, 12:{}}}
        for O in ["withOneChange","noOneChange"]:
            logger.debug ("Count e1: %s e2: %s O: %s", e1, e2,O)
            infosC1C2E1E1_file = "%s/%s_%s_%s_%s.infos"%(repest,ali_filename,e1,e1,O)
            infosC1C2E1E2_file = "%s/%s_%s_%s_%s.infos"%(repest,ali_filename,e1,e2,O)

            ### read infos
            dli[O][12][12]=read_info_exception(infosC1C2E1E2_file, n_sites, O, 12, 12)
            dli[O][12][11]=read_info_exception(infosC1C2E1E1_file, n_sites, O, 12, 11)


        for s in lsites:
            for k in [12]:
               # get all likelihood:
               X   = dli["noOneChange"][k][11][s]
               XY  = dli["noOneChange"][k][12][s]
               OX  = dli["withOneChange"][k][11][s]
               OXY = dli["withOneChange"][k][12][s]

               X_l  [s][k].append(X)
               XY_l [s][k].append(XY)
               OX_l [s][k].append(OX)
               OXY_l[s][k].append(OXY)
    # Bilan
    for s in lsites:
        for k in [12]:
            # p_max:
            max_X   = likelihood_max(X_l  [s][k])
            bilan[k]["max_X"].append(max_X)
            max_XY  = likelihood_max(XY_l [s][k])
            bilan[k]["max_XY"].append(max_XY)
            max_OX  = likelihood_max(OX_l [s][k])
            bilan[k]["max_OX"].append(max_OX)
            max_OXY = likelihood_max(OXY_l[s][k])
            bilan[k]["max_OXY"].append(max_OXY)

            bilan[k]["p_max_X_XY"].append(   prob_ap(max_XY,max_X)   )
            bilan[k]["p_max_OX_OXY"].append( prob_ap(max_OXY,max_OX) )
            bilan[k]["p_max_X_OXY"].append(  prob_ap(max_OXY,max_X)  )
            bilan[k]["p_max_X_OX"].append(   prob_ap(max_OX,max_X)   )
            bilan[k]["p_max_XY_OXY"].append( prob_ap(max_OXY,max_XY) )

            # p_mean:
            mean_X   = likelihood_mean(X_l  [s][k])
            bilan[k]["mean_X"].append(mean_X)
            mean_XY  = likelihood_mean(XY_l [s][k])
            bilan[k]["mean_XY"].append(mean_XY)
            mean_OX  = likelihood_mean(OX_l [s][k])
            bilan[k]["mean_OX"].append(mean_OX)
            mean_OXY = likelihood_mean(OXY_l[s][k])
            bilan[k]["mean_OXY"].append(mean_OXY)

            bilan[k]["p_mean_X_XY"].append(   prob_ap_no_log(mean_XY, mean_X)  )
            logger.debug("s: %s, XY %s, X %s, p_mean_X_XY %s", s,  mean_XY, mean_X, bilan[k]["p_mean_X_XY"][-1])
            bilan[k]["p_mean_OX_OXY"].append( prob_ap_no_log(mean_OXY,mean_OX) )
            bilan[k]["p_mean_X_OXY"].append(  prob_ap_no_log(mean_OXY,mean_X)  )
            bilan[k]["p_mean_X_OX"].append(   prob_ap_no_log(mean_OX, mean_X)  )
            bilan[k]["p_mean_XY_OXY"].append( prob_ap_no_log(mean_OXY,mean_XY) )

    # Calcul (VP, FP, FN, VN)
    res = []
    for P in P_l:
        for seuil in lseuil:
            tab={}

            tab["PosteriorProbabilityType"] = P
            tab["P"] = float(len([s for s in bilan[12][P] if s >= seuil]))
            tab["N"] = float(len([s for s in bilan[12][P] if s < seuil]))

            tab["NumberOfSites"] = n_sites
            tab["NumberOfConvergentEvents"] = N
            tab["Threshold"] = seuil
            tab["Method"] = get_method(P)
            if tree:
                tab["InputTree"] = tree
            if ID:
                tab["RunID"] = ID

            res.append(tab)

    return res, bilan

########################################################################
##           Ancestral reconstruction method                          ##
########################################################################


def outdiff_new(tr,name,repseq,reptree,repest,c1,c2):
  laa=list("ARNDCQEGHILKMFPSTWYV")
  ##### Sequences
  fasta_file = "%s/%s.fa" %(repseq,name)
  ali = AlignIO.read(fasta_file, "fasta")

  ##### Tree
  tr = Tree("%s/annotated_tree.nhx" %(reptree))
  l_branch_with_transition = tr.search_nodes(T="True")


  #### output per site nb of clades with ...

  ln = len(ali[0].seq)

  diffone=[0]*ln
  diffall=[0]*ln
  diffdeb=[0]*ln
  pdiff_mdiff=[]
  for i in range(ln):
    pdiff_mdiff.append({x:0 for x in laa})

  for chg in l_branch_with_transition:
    fat=chg.up
    seqanc=ali[int(fat.ND)].seq
    seqici=ali[int(chg.ND)].seq
    #logger.debug(seqanc == seqici)
    #logger.debug(seqanc)
    #logger.debug(seqici)
    lconv=chg.get_leaves()

    pdiff_diffone=[]
    pdiff_diffall=[0]*ln

    for ch in lconv:
      seqconv=ali[int(ch.ND)].seq
      
      #logger.debug(seqconv)

      for p in range(ln):
        if seqanc[p]!=seqconv[p]:
            pdiff_diffall[p]+=1
            pdiff_diffone.append(p)

    sdiff_diffone=set(pdiff_diffone)
    sdiff_diffdeb = []
    sdiff_diffall = []

    for p in range(ln):
        if seqanc[p]!=seqici[p]:
            pdiff_mdiff[p][seqici[p]]+=1
            sdiff_diffdeb.append(p)

        if pdiff_diffall[p]==len(lconv):
            sdiff_diffall.append(p)


    for x in sdiff_diffone:
      diffone[x]+=1

    for x in sdiff_diffall:
      diffall[x]+=1

    for x in sdiff_diffdeb:
      diffdeb[x]+=1

  mdiff=[max(x.values()) for x in pdiff_mdiff]


  #~ #### output per site nb of clades with at least a change
  #~ print "diffone"
  #~ print diffone

  #~ #### output per site nb of clades with changes at all leaves
  #~ print "diffall"
  #~ print diffall

  #~ #### output per site nb of clades with a change at the start
  #~ print "diffdeb"
  #~ print diffdeb

  #~ #### output per site max nb of clades with similar aa at the start,
  #~ #### different from the ancestral one
  #~ print "mdiff"
  #~ print mdiff

  #logger.debug(mdiff)
  with open(repest+"/"+name+".conv","w") as f:
    f.write("nb_cl_1diff\tnb_cl_alldiff\tnb_cl_diffstart\tmax_same_aa_diffstart\n")
    for i in range(ln):
        f.write("%d\t%d\t%d\t%d\n"%(diffone[i],diffall[i],diffdeb[i],mdiff[i]))


def dico_typechg_obs_sub(C1,C2,N,repest,name_AC,tree="", n_sites = 1000,  ID= "", dist_C1_C2 = ""):
### dico des proba a posteriori qd nb clades avec changement = nbe clades convergents

## 11 : sequences generees par C1-C1
## 12 : sequences generees par C1-C2
    simu_i = name_AC.split("_")[1]
    lsites = range(n_sites)
    bilan = {11:{"pos":lsites}, 12:{"pos":lsites}}
    lseuil = [1]
    P_l = ["p_ident"]
    for P in P_l:
        for k in [11,12]:
            bilan[k][P]=[]

    conv_file_C1C2 = "%s/Scenario_%s_A%s_C%s.conv"%(repest,simu_i,C1,C2)
    conv_file_C1C1 = "%s/Scenario_%s_A%s_C%s.conv"%(repest,simu_i,C1,C1)
    
    dch={x:[] for x in [11,12]}

    if not os.path.exists(conv_file_C1C2):
      logger.error("%s n'existe pas"%(conv_file_C1C2))
    if not os.path.exists(conv_file_C1C1):
      logger.error("%s n'existe pas"%(conv_file_C1C1))

    ### filtre sites et compte same_aa
    #lsites={}
    for k in [11,12]:
      if k==11:
        fconv=open(conv_file_C1C1,"r")
      else:
        fconv=open(conv_file_C1C2,"r")

      dchk=dch[k]
      l=fconv.readline()
      l=fconv.readline()
      i=0
      while l:
        ls=l.split()
        if len(ls)<=3:
          break
        nsaa=int(ls[3])
        dchk.append(nsaa)

        l=fconv.readline()

      fconv.close()

    res=[]
    for s in lsites:
        for k in [11, 12]:
            bilan[k]["p_ident"].append(1 if (dch[k][s] == N) else 0)
    
    for P in P_l:
        for seuil in lseuil:
            tab={}
            tab["TP"]=float(len([s for s in bilan[12][P] if s >= seuil]))
            tab["FN"]=float(len([s for s in bilan[12][P] if s < seuil]))
            tab["FP"]=float(len([s for s in bilan[11][P] if s >= seuil]))
            tab["TN"]=float(len([s for s in bilan[11][P] if s < seuil]))
            #tab["TP"]=float(len([s for s in dch[12] if s == N]))
            #tab["FN"]=float(len([s for s in dch[12] if s < N]))
            #tab["FP"]=float(len([s for s in dch[11] if s == N]))
            #tab["TN"]=float(len([s for s in dch[11] if s < N]))

            tab["PosteriorProbabilityType"] = "p_ident"
            tab["C1"] = C1
            tab["C2"] = C2
            tab["SimuCoupleID"] = "A%s_C%s" %(C1, C2)
            tab["Threshold"] = "NA"
            tab["NumberOfSites"] = n_sites
            tab["NumberOfConvergentEvents"] = N
            tab["DistanceSimuCouple"] = dist_C1_C2["C" + str(C1)]["C" + str(C2)]
            tab["Method"] = "Identical"
            tab["ScenarioID"] = "Scenario_%s" %simu_i
            if tree:
                tab["InputTree"] = tree
            if ID:
                tab["RunID"] = ID
            
            tab["Sensitivity"], tab["Specificity"], tab["MCC"] = sensspec(tab)
            res.append(tab)
            
    return res, bilan


########################################################################
##                         Topological method                         ##
########################################################################



def dico_typechg_topo(C1, C2, N, repest, name_AC, tree = "", set_t1 = [], n_sites = 0, ID = "", lseuil = [0.7,0.80,0.85,0.90,0.95,0.99], NbCat_Est = 10, dist_C1_C2 = ""):
    ### 11 : sequences generees par C1-C1
    ### 12 : sequences generees par C1-C2

    simu_i = name_AC.split("_")[1]

    lsites = range(n_sites)
    
    bilan = {11:{"pos":lsites}, 12:{"pos":lsites}}
    P_l = []
    test1_l = ["X_CX"]
    test2_l = ["p_max", "p_mean", "p_opt"]

    for test1 in test1_l:
        for test2 in test2_l:
            P = test2 + "_" + test1
            P_l.append(P)
            for k in [11,12]:
                bilan[k][P]=[]
    for a in ["X", "CX"]:
        for b in ["opt_", "mean_", "max_"]:
            for k in [11,12]:
                bilan[k][b+a] = []

    if not DEV:
        test1_l = ["X_CX"]
        test2_l = ["p_mean"]
        P_l = []
        for test1 in test1_l:
            for test2 in test2_l:
                P = test2 + "_" + test1
                P_l.append(P)

    if not set_t1:
        for t1 in range(1, (NbCat_Est+1)):
            set_t1.append(t1)
    
    infosC1C1C1C1OPT_file = "%s/Scenario_%s_A%s_C%s_%s_%s_opt_noOneChange.infos"%(repest,simu_i,C1,C1,C1,C1)
    infosC1C1T1OPT_file = "%s/Scenario_%s_A%s_C%s_topo_t%s_opt.infos"%(repest,simu_i,C1,C1,C1)
    
    infosC1C2C1C1OPT_file = "%s/Scenario_%s_A%s_C%s_%s_%s_opt_noOneChange.infos"%(repest,simu_i,C1,C2,C1,C1)
    infosC1C2T1OPT_file = "%s/Scenario_%s_A%s_C%s_topo_t%s_opt.infos"%(repest,simu_i,C1,C2,C1)
  
    dli={}
    dli[11]={}  # spec
    dli[12]={}  # sens

    ### read infos
    dli[11][12] = read_info_exception(infosC1C1T1OPT_file,n_sites, "noOneChange_opt", 11, 12)
    dli[11][11] = read_info_exception(infosC1C1C1C1OPT_file,n_sites, "noOneChange_opt", 11, 11)
    dli[12][12] = read_info_exception(infosC1C2T1OPT_file,n_sites, "noOneChange_opt", 12, 12)
    dli[12][11] = read_info_exception(infosC1C2C1C1OPT_file,n_sites, "noOneChange_opt", 12, 11)

    for s in lsites:
        for k in [11,12]:
            # get all likelihood:
            opt_X   = dli[k][11][s]
            bilan[k]["opt_X"].append(opt_X)
            opt_CX  = dli[k][12][s]
            bilan[k]["opt_CX"].append(opt_CX)
            bilan[k]["p_opt_X_CX"].append(prob_ap(opt_CX,opt_X))
    
    X_l   = []
    CX_l  = []
    
    for s in lsites:
        X_l.append({11:[], 12:[]})
        CX_l.append({11:[], 12:[]})
    
    for t1 in set_t1:
        dli={}
        dli[11]={}  # spec
        dli[12]={}  # sens
        logger.debug ("Count T1: %s", t1)
        infosC1C1T1_file = "%s/Scenario_%s_A%s_C%s_topo_t%s.infos"%(repest,simu_i,C1,C1,t1)
        infosC1C1E1E1_file = "%s/Scenario_%s_A%s_C%s_%s_%s_noOneChange.infos"%(repest,simu_i,C1,C1,t1,t1)

        infosC1C2E1E1_file = "%s/Scenario_%s_A%s_C%s_%s_%s_noOneChange.infos"%(repest,simu_i,C1,C2,t1,t1)
        infosC1C2T1OPT_file = "%s/Scenario_%s_A%s_C%s_topo_t%s.infos"%(repest,simu_i,C1,C2,t1)
        

        ### read infos
        dli[11][12]=read_info_exception(infosC1C1T1_file,n_sites, "noOneChange", 11, 12)
        dli[11][11]=read_info_exception(infosC1C1E1E1_file,n_sites, "noOneChange", 11, 11)
        dli[12][12]=read_info_exception(infosC1C2T1OPT_file,n_sites, "noOneChange", 12, 12)
        dli[12][11]=read_info_exception(infosC1C2E1E1_file,n_sites, "noOneChange", 12, 11)

        for s in lsites:
            for k in [11,12]:
               # get all likelihood:
               X   = dli[k][11][s]
               CX  = dli[k][12][s]

               X_l  [s][k].append(X)
               CX_l [s][k].append(CX)


    # Bilan
    for s in lsites:
        for k in [11,12]:
            # p_max:
            max_X   = likelihood_max(X_l  [s][k])
            bilan[k]["max_X"].append(max_X)
            max_CX  = likelihood_max(CX_l [s][k])
            bilan[k]["max_CX"].append(max_CX)

            bilan[k]["p_max_X_CX"].append(   prob_ap(max_CX,max_X)   )

            # p_mean:
            mean_X   = likelihood_mean(X_l  [s][k])
            bilan[k]["mean_X"].append(mean_X)
            mean_CX  = likelihood_mean(CX_l [s][k])
            bilan[k]["mean_CX"].append(mean_CX)
            
            bilan[k]["p_mean_X_CX"].append(   prob_ap_no_log(mean_CX, mean_X)  )

            #logger.debug("mean_X: %s max_X:%s", bilan[k]["mean_X"][s], bilan[k]["max_X"][s])
            #logger.debug("mean_CX: %s max_CX:%s", bilan[k]["mean_CX"][s], bilan[k]["max_CX"][s])
            #logger.debug("p_mean_X_CX: %s p_max_X_CX:%s", bilan[k]["p_mean_X_CX"][s], bilan[k]["p_max_X_CX"][s])

    # Calcul (VP, FP, FN, VN)
    res = []
    for P in P_l:
        for seuil in lseuil:
            tab={}
           
            tab["PosteriorProbabilityType"]= P
            tab["TP"]=float(len([s for s in bilan[12][P] if s >= seuil]))
            tab["FN"]=float(len([s for s in bilan[12][P] if s < seuil]))
            tab["FP"]=float(len([s for s in bilan[11][P] if s >= seuil]))
            tab["TN"]=float(len([s for s in bilan[11][P] if s < seuil]))
            tab["C1"] = C1
            tab["C2"] = C2
            tab["SimuCoupleID"] = "A%s_C%s" %(C1, C2)
            tab["DistanceSimuCouple"] = dist_C1_C2["C" + str(C1)]["C" + str(C2)]
            tab["NumberOfSites"] = n_sites
            tab["NumberOfConvergentEvents"] = N
            tab["Threshold"] = seuil
            tab["Method"] = get_method(P)
            tab["ScenarioID"] = "Scenario_%s" %simu_i
            if tree:
                tab["InputTree"] = tree
            if ID:
                tab["RunID"] = ID
            
            tab["Sensitivity"], tab["Specificity"], tab["MCC"] = sensspec(tab)
            res.append(tab)
    
    return res, bilan

