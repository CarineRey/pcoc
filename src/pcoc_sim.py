#!/usr/bin/python
#  pcoc_sim.py
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

import datetime
import glob
import random
import sys
import argparse
import os
import re
import csv
import time
import logging
import subprocess

import bpp_lib
import events_placing
import estim_data
import plot_data

import multiprocessing

import pandas as pd
import numpy as np
from ete3 import Tree

import shutil

DEV = False

##########
# inputs #
##########
start_time = time.time()

### Option defining
parser = argparse.ArgumentParser(prog="pcoc_sim.py",
                                 description='')
parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')
parser._optionals.title = "Miscellaneous options"
parser.add_argument('-cpu', type=int,
                    help="Number of cpu to use. (default: 1)",
                    default=1)

##############
requiredOptions = parser.add_argument_group('Required options')
requiredOptions.add_argument('-td', "--tree_dir", type=str,
                             help='Directory name containing input trees.', required=True)
requiredOptions.add_argument('-o', '--output_dir', type=str,
                   help="Output directory name.", required=True)
##############


##############
Options_trees = parser.add_argument_group('Convergent scenarios simulation options')
Options_trees.add_argument('-n_sc', type=int,  metavar="INT",
                    help="Number of convergent scenarios picked in each input tree. (default: 1)",
                    default=1)
Options_trees.add_argument('-m', '--manual_mode', metavar="\"x/y,z/...\"", type=str,
                    help="User defined convergent transition/branches. Transition node must be the first number and independent events must be separed by a \"/\". ex: \"1,2,3/67/55,56\" (default: None)",
                    default="")
Options_trees.add_argument('-c_max', type=int, metavar="INT",
                    help="Maximum number of potential convergent transitions (=convergent events). (default: 7)", # First, this number of events are drawn up in the tree, then the wanted number of event is randomly chosen among them. It is to avoid any bias when you want to compared 2 and 7 events. (default: 7)",
                    default=7)
Options_trees.add_argument('-c_min', type=int, metavar="INT",
                    help="Minimum number of transition (=convergent events). (default: 2)",
                    default=2)
Options_trees.add_argument('-c', type=int, metavar="INT",
                    help="Number of convergent transitions (=convergent events). If not defined, random between min events and max events",
                    default=0)
Options_trees.add_argument('-cr', type=float, metavar="FLOAT",
                    help="Maximum ratio between the number of convergent/non convergent leaves. (default: No limits)",
                    default=1)
Options_trees.add_argument('-flg', type=float, metavar="FLOAT",
                    help="For each input tree, branch length multiplicator. (default: no modification)",
                    default=1)
Options_trees.add_argument('-bl_new', metavar="FLOAT", type=float,
                    help="For each input tree, replace all branch lengths by [FLOAT]. (default: no modification)",
                    default=-1)
Options_trees.add_argument('--ali_noise', action="store_true",
                    help="Add noisy events in the convergent scenario.",
                    default=False)
Options_trees.add_argument('--bl_noise', action="store_true",
                    help="Add noise in the branch lengths of of tree for the detection process.",
                    default=False)
##############


##############
Options_ali = parser.add_argument_group('Alignment simulation options')
Options_ali.add_argument('-nb_sampled_couple', type=int,  metavar="INT",
                    help="For each convergent scenario, number of simulated alignment with different sampled couple of profiles (Ancestral/Convergent). (default: 1)",
                    default=1)
Options_ali.add_argument('-n_sites', type=int,  metavar="INT",
                    help="Number of simulated sites per alignment. (default: 100)",
                    default=100)
Options_ali.add_argument('-CATX_sim', type=int, choices = [10,60],
                    help="Profile categories to simulate data (10->C10 or 60->C60). (default: 60)",
                    default=60)
Options_ali.add_argument('-min_dist_CAT', type=float, metavar="FLOAT",
                    help="Minimum distance between Ancestral and Convergent profiles to simulate the alignment (default: no limits)",
                    default=0)
Options_ali.add_argument('--plot_ali', action="store_true",
                    help="For each couple of profiles, plot a summary of the convergent scenario containing the tree and the alignment.",
                    default=False)
Options_ali.add_argument('--get_likelihood_summaries', action="store_true",
                    help="For each couple of profiles, write a summary of the likelihoods per site.",
                    default=False)
Options_ali.add_argument('--no_clean_seqs', action="store_true",
                    help="Do not cleanup the sequences after the run.",
                    default=False)
##############


##############
Options_ben = parser.add_argument_group('Detection options')
Options_ben.add_argument('--pcoc', action="store_true",
                    help="Use the PCOC model approach to detect sites under convergent evolution.",
                    default=False)
Options_ben.add_argument('--ident', action="store_true",
                    help="Use the ancestral reconstruction approach to detect sites under convergent evolution.",
                    default=False)
Options_ben.add_argument('--topo', action="store_true",
                    help="Use the topological approach to detect sites under convergent evolution.",
                    default=False)
Options_ben.add_argument('-CATX_est', type=int, choices = [10,60],
                    help="Profile categories to estimate data (10->C10 or 60->C60). (default: 10)",
                    default=10)
Options_ben.add_argument('--plot_event_repartition', action="store_true",
                    help="Plot chosen random convergent events repartition for each input tree.",
                    default=False)
#Options_ben.add_argument('--plot', action="store_true",
#                    help="Plot sensitivity and specificity of each appraoch for each input tree.",
#                    default=False)
##############


##############
Options_other = parser.add_argument_group('Other Options')
#Options_other.add_argument('-log', type=str, default="",
#                  help="a log filename to report avancement (default: no)")
Options_other.add_argument('--no_cleanup', action="store_true",
                    help="Do not cleanup the working directory after the run.",
                    default=False)
Options_other.add_argument("-LD_LIB", metavar='LD_LIBRARY_PATH', type=str, default="",
                   help="Redefine the LD_LIBRARY_PATH env variable, bppsuite library must be present in the $PATH and in the $LD_LIBRARY_PATH")
Options_other.add_argument('--debug', action="store_true",
                    help="debug mode", default=False)

### Option parsing
args = parser.parse_args()

metadata_run_dico = {}
date = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
OutDirName = "%s/RUN_%s/" %(args.output_dir, date)
OutDirName = OutDirName.replace("//","/")

metadata_run_dico["RunID"] = date

### Set up the output directory
if os.path.isdir(OutDirName):
    pass
    #logger.info("The output directory %s exists", OutDirName)
elif OutDirName: # if OutDirName is not a empty string we create the directory
    #logger.info("The output directory %s does not exist, it will be created", OutDirName)
    os.makedirs(OutDirName)

### Set up the log file
LogFile = OutDirName + "/pcoc_sim.log"

### Set up the logger
# create logger
logger = logging.getLogger("pcoc")
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler(LogFile)
# create console handler with a higher log level
ch = logging.StreamHandler()
if args.debug:
    ch.setLevel(logging.DEBUG)
    fh.setLevel(logging.DEBUG)
else:
    ch.setLevel(logging.INFO)
    fh.setLevel(logging.INFO)
# create formatter and add it to the handlers
formatter_fh = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
formatter_ch = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter_fh)
ch.setFormatter(formatter_ch)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

logger.debug(sys.argv)

cpu = args.cpu
try:
    cpus = multiprocessing.cpu_count()
except NotImplementedError:
    cpus = 1   # arbitrary default
logger.info("%s on %s cpus", cpu, cpus)

if args.LD_LIB:
    logger.info("$LD_LIBRARY_PATH will be change from %s to %s", os.environ.get("LD_LIBRARY_PATH", ""), args.LD_LIB)
    os.environ["LD_LIBRARY_PATH"]=args.LD_LIB
else:
    logger.debug("$LD_LIBRARY_PATH is %s", os.environ.get("LD_LIBRARY_PATH", ""))


repbppconfig = OutDirName +  "bpp_config"
if not os.path.exists(repbppconfig):
    os.mkdir(repbppconfig)


bpp_lib.write_config(repbppconfig, estim=True, NbCat = args.CATX_sim)
dist_C1_C2 =  pd.read_csv(repbppconfig+'/CATC'+str(args.CATX_sim)+'Distances.csv', index_col=0)

lnf=glob.glob(args.tree_dir+"/*")

flg=args.flg
if flg == 1:
    logger.info("Branch length multiplicator:\tno")
    metadata_run_dico["Branch length multiplicator"] = "no"
else:
    logger.info("Branch length multiplicator:\t%s", flg)
    metadata_run_dico["Branch length multiplicator"] = flg

bl_new=args.bl_new
if bl_new >0:
    logger.info("Branch length remplacement:\t%s", bl_new)
    metadata_run_dico["Branch length remplacement"] = bl_new
else:
    bl_new = -1
    metadata_run_dico["Branch length remplacement"] = "no"

#http://stackoverflow.com/questions/23172293/use-python-to-extract-branch-lengths-from-newick-format
pattern = re.compile(r"\b[0-9]+(?:\.[0-9]+)?\b")


logger.debug("trees:\n  * %s", "\n  * ".join(lnf))
nb_input_tree_before = len(lnf)
logger.debug("%s trees in %s", nb_input_tree_before, args.tree_dir)


for treefilename in lnf:
    # test if a tree
    try:
        t=Tree(treefilename)
    except:
        logger.warning("%s is not a newick tree, this tree will not be used",treefilename)
        lnf.remove(treefilename)
        t=""
    if t:
        treefile=open(treefilename,"r")
        tree=treefile.read().strip()
        treefile.close()
        #test if branch length
        branch_lengths = pattern.findall(tree)
        if branch_lengths == []:
            logger.warning("No branch length in %s, this tree will not be used",treefilename)
            lnf.remove(treefilename)

nb_input_tree_after = len(lnf)
if nb_input_tree_after!= nb_input_tree_before:
    logger.warning("%s trees in %s after checking (%s before)", nb_input_tree_after, args.tree_dir, nb_input_tree_before)

logger.debug("trees:\n  * %s", "\n  * ".join(lnf))

logger.info("Number of input trees:\t%s", nb_input_tree_after)
metadata_run_dico["Number of input trees"] = nb_input_tree_after

if len(lnf) == 0:
    logger.error("No tree. Bye.")
    sys.exit(1)

if len(lnf) > 1 and args.manual_mode:
    logger.error("Only 1 tree if manual mode.")
    sys.exit(1)

metadata_run_dico["Input trees"] = ",".join([os.path.basename(t) for t in lnf])

Nbsimul=args.n_sc
maxTrans=args.c_max
minTrans=args.c_min
maxConvRate=args.cr
Nsites=args.n_sites


manual_mode_nodes = {}
if args.manual_mode:
    manual_mode_nodes = {"T":[],"C":[]}
    p_events = args.manual_mode.strip().split("/")
    for e in p_events:
        l_e = map(int,e.split(","))
        manual_mode_nodes["T"].append(l_e[0])
        manual_mode_nodes["C"].extend(l_e[1:])
    if args.c_min > len(manual_mode_nodes["T"]):
        minTrans=len(manual_mode_nodes["T"])
    if args.c_max != len(manual_mode_nodes["T"]):
        maxTrans = len(manual_mode_nodes["T"])



metadata_run_dico["Number of scenarios per input tree"] = Nbsimul
metadata_run_dico["Maximum number of convergent events"] = maxTrans
metadata_run_dico["Minimum number of convergent events"] = minTrans
metadata_run_dico["Maximum rate of the number of Convergent/Non-convergent leafs"] = maxConvRate
metadata_run_dico["Number of simulated sites"] = Nsites

logger.info("Number of scenarios per input tree (= 1 tree and 1 set of convergent events):\t%s", Nbsimul)
logger.info("Number of simulated sites:\t%s", Nsites)
logger.info("Maximum number of convergent events:\t%s", maxTrans)
logger.info("Minimum number of convergent events:\t%s", minTrans)
logger.info("Maximum rate of the number of Convergent/Non-convergent leafs:\t%s", maxConvRate)


NbCat_Sim = args.CATX_sim
NbCat_Est = args.CATX_est

MinDistCAT = args.min_dist_CAT
Nb_sampled_couple = args.nb_sampled_couple

metadata_run_dico["Profile categories use during simulation"] = NbCat_Sim
metadata_run_dico["Profile categories use during estimation"] = NbCat_Est
metadata_run_dico["Number of sampled profile couples per scenario"] = Nb_sampled_couple
metadata_run_dico["Minimum distance between 2 profiles of a couple to be use for simulations"] = MinDistCAT

pcoc_str = "No"
topo_str = "No"
ident_str = "No"
if args.pcoc:
    pcoc_str = "Yes"
if args.topo:
    topo_str = "Yes"
if args.ident:
    ident_str = "Yes"
logger.info("Run PCOC method:\t%s", pcoc_str)
logger.info("Run topological method:\t%s", topo_str)
logger.info("Run identical method:\t%s", ident_str)
metadata_run_dico["Run PCOC method"] = pcoc_str
metadata_run_dico["Run topological method"] = topo_str
metadata_run_dico["Run identical method"] = ident_str

pd.Series(metadata_run_dico).to_csv(OutDirName + "/run_metadata.tsv", sep='\t')

### Choose CAT profiles

CATcouples = []
if MinDistCAT > 0:
    with open(repbppconfig+'/CATC'+str(NbCat_Sim)+'Distances.csv', 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        first_line = True
        for row in reader:
            #print row
            if first_line:
                first_line = False
                colnames = row
                colnames.pop(0)
            else:
                C1 = int(row.pop(0).replace("C",""))
                for C2 in colnames:
                    C2 = int(C2.replace("C",""))
                    d=float(row.pop(0))
                    if d >= MinDistCAT:
                        CATcouples.append((C1,C2))
                    else:
                        #print("rejecte %s %s" %(C1,C2))
                        pass
else:
    for C1 in range(1,NbCat_Sim+1):
        for C2 in range(1,NbCat_Sim+1):
            if C1 != C2:
                CATcouples.append((C1,C2))

if not Nb_sampled_couple:
    Nb_sampled_couple = len(CATcouples)
elif Nb_sampled_couple > len(CATcouples):
    Nb_sampled_couple = len(CATcouples)


logger.info("Profile category uses during simulation:\t%s", NbCat_Sim)
logger.info("Profile category uses during estimation:\t%s", NbCat_Est)
logger.info("Minimum distance between 2 profiles of a couple to be use for simulations:\t%s", MinDistCAT)
logger.info("Number of sampled profile couples per scenario:\t%s", Nb_sampled_couple)


def remove_folder(path):
    # check if folder exists
    if os.path.exists(path):
        # remove if exists
        shutil.rmtree(path)

def mk_simu((i, tree_filename, OutDirNamePrefixTree), n_try = 0) :
    start_simu = time.time()
    name0="Scenario_%d"%(i)
    name0_info="Scenario %d/%d"%(i, Nbsimul)
    logger.info("START: %s", name0_info)
    metadata_simu_dico = {}
    logger.debug("Tree: %s", os.path.basename(tree_filename))
    metadata_simu_dico["InputTree"] = os.path.basename(tree_filename)
    metadata_simu_dico["ScenarioID"] = name0
    metadata_simu_dico["RunID"] = date

    metadata_simu_dico["BranchLengthMultiplicator"] = metadata_run_dico["Branch length multiplicator"]
    metadata_simu_dico["BranchLengthRemplacement"] = metadata_run_dico["Branch length remplacement"]
    metadata_simu_dico["MaxNumberOfConvergentEvents"] = maxTrans
    metadata_simu_dico["MinNumberOfConvergentEvents"] = minTrans

    tree=events_placing.init_tree(tree_filename)

    if args.ali_noise:
        logger.info("%s - Addition of noise in the alignment: Yes", name0_info)
        metadata_simu_dico["AliNoise"] = "Yes"
        tree, cz_nodes, r_cz = events_placing.noise_tree(tree, NbCat = args.CATX_sim)
    else:
        logger.info("%s - Addition of noise in the alignment: No", name0_info)
        metadata_simu_dico["AliNoise"] = "No"
        cz_nodes = {}
        r_cz = 0

    if args.c == 0:
        n_events=random.randrange(minTrans,(maxTrans+1))
    else:
        n_events = args.c

    logger.info("%s - Number of convergent events:\t%s", name0_info, n_events)
    metadata_simu_dico["NumberOfConvergentEvents"] = n_events

    annotated_tree, nodesWithTransitions, nodesWithConvergentModel, nodesWithAncestralModel, numberOfLeafsWithTransitions, numberOfLeafs = events_placing.placeNTransitionsInTree_new(n_events,maxTrans,maxConvRate,tree, manual_mode_nodes = manual_mode_nodes,nf=tree_filename)

    numberOfNodes = len(annotated_tree.get_descendants())
    if len(nodesWithConvergentModel)==0:
        if n_try >=5:
            logger.warning("%s - Impossible to put the convergent events", name0_info)
            return {},[],[],[],[]
        else:
            logger.info("%s - Start again (%s/5)", name0, n_try)
            R = mk_simu((i, tree_filename, OutDirNamePrefixTree), n_try = n_try + 1)
            return R
    else:
        logger.info("%s - Assignment of the convergent events: OK", name0_info)

    nodesWithAncestralModel_sim = [x for x in nodesWithAncestralModel]
    if cz_nodes:
        for (cz, nodes) in cz_nodes.items():
            t_node = nodes[0]
            nodesWithAncestralModel_sim = list(set(nodesWithAncestralModel_sim) - set(nodes))
            if not t_node in set(nodesWithTransitions + nodesWithConvergentModel):
                nodes = nodes[1:]
                cz_nodes[cz] = [t_node] + list(set(nodes) - set(nodesWithTransitions + nodesWithConvergentModel))
            else:
                cz_nodes[cz] = list(set(nodes) - set(nodesWithTransitions + nodesWithConvergentModel))
          #logger.debug("nodesWithNoisyModel: %s\n", "\n\t-".join(["%s: %s" %(cz,",".join(nodes)) for (cz,nodes) in cz_nodes.items()]))
        logger.debug("nodesWithNoisyModel: \n%s", "\n\t-".join(["%s: %s" %(cz, nodes) for (cz,nodes) in cz_nodes.items()]))
        metadata_simu_dico["RatioNoisyNonNoisyLeaves"] = float(len(events_placing.unlist(cz_nodes.values())))/float(numberOfNodes)

    logger.debug("nodesWithAncestralModel: %s", nodesWithAncestralModel)
    logger.debug("nodesWithAncestralModel_sim: %s", nodesWithAncestralModel_sim)
    logger.debug("nodesWithConvergentModel: %s", nodesWithConvergentModel)
    logger.debug("nodesWithTransitions: %s", nodesWithTransitions)
    logger.debug("numberOfNodes: %s", numberOfNodes)
    logger.debug("len: %s tot : %s ", len(nodesWithTransitions+nodesWithConvergentModel+nodesWithAncestralModel_sim+events_placing.unlist(cz_nodes.values())),numberOfNodes)

    if len(nodesWithTransitions+nodesWithConvergentModel+nodesWithAncestralModel_sim+events_placing.unlist(cz_nodes.values())) != numberOfNodes:
        logger.debug("annotated_tree: %s", annotated_tree.get_ascii(attributes=["Cz","C"]))

        sys.exit(1)

    metadata_simu_dico["numberOfLeafs"] = numberOfLeafs
    metadata_simu_dico["numberOfLeafsWithTransitions"] = numberOfLeafsWithTransitions

    repseq="%s/%s"%(repseq0,name0)
    repseq=repseq.replace("//","/")
    repest="%s/%s"%(repest0,name0)
    repest=repest.replace("//","/")
    reptree="%s/%s"%(reptree0,name0)
    reptree=reptree.replace("//","/")
    repplottreeali="%s/%s"%(repplottreeali0,name0)
    repplottreeali=repplottreeali.replace("//","/")
    replikelihoodsummary="%s/%s"%(replikelihoodsummary0,name0)
    replikelihoodsummary=replikelihoodsummary.replace("//","/")

    if not os.path.exists(repseq) and not os.path.exists(repest) \
       and not os.path.exists(reptree) and not os.path.exists(repplottreeali):
        os.mkdir(repseq)
        os.mkdir(repest)
        os.mkdir(reptree)
        if args.plot_ali:
            os.mkdir(repplottreeali)
        if args.get_likelihood_summaries:
            os.mkdir(replikelihoodsummary)
    else:
        logger.error("%s/%s must not exist", repseq0,name0)
        sys.exit(1)

    logger.debug("repseq: %s", repseq)
    logger.debug("repest: %s", repest)
    logger.debug("reptree: %s", repest)

    if args.bl_noise:
        logger.info("%s - Addition of noise in the branch lengths: Yes", name0_info)
        metadata_simu_dico["BLNoise"] = "Yes"
        (node_ID, bl_before, bl_err, bl_after) = events_placing.noise_bl(annotated_tree, reptree, vnodes=nodesWithConvergentModel+nodesWithTransitions, topo_met=args.topo)
        #logger.debug("node_ID: %s\nbl_before: %s\nbl_err: %s\nbl_after: %s",node_ID, bl_before, bl_err, bl_after)
        median_err = np.median(bl_err)
        sd_err = np.std(bl_err)
        mean_err = np.mean(bl_err)
        metadata_simu_dico["Mean_Bl_err"] = mean_err
        metadata_simu_dico["Sd_Bl_err"] = sd_err
        metadata_simu_dico["Median_Bl_err"] = median_err
        logger.info("Median bl err: %s", median_err)
        logger.info("Sd bl err: %s", median_err)
        logger.info("Mean bl err: %s", mean_err)
    else:
        logger.info("%s - Addition of noise in the branch lengths: No", name0_info)
        metadata_simu_dico["BLNoise"] = "No"


    ### construit les arbres d'etude :
    (tree_annotated, tree_conv_annotated,  allbranchlength, convbranchlength) = events_placing.mk_tree_4_simu_new(annotated_tree,reptree,nodesWithConvergentModel,nodesWithTransitions, flg = flg, bl_new = bl_new, topo_met=True, plot=args.plot_ali, cz_nodes=cz_nodes)

    metadata_simu_dico["AllBranchLength"] = allbranchlength
    metadata_simu_dico["ConvBranchLength"] = convbranchlength

    SampledCATcouples = random.sample(CATcouples, Nb_sampled_couple)

    l_TPFPFNTN_mod_het = []
    l_TPFPFNTN_topo = []
    l_TPFPFNTN_obs_sub = []
    k_couple = 0
    for (c1,c2) in SampledCATcouples:
        k_couple+=1
        logger.info("Scenario %s/%s - Couple %s/%s, C1: %s, C2: %s",i, Nbsimul, k_couple, Nb_sampled_couple, c1, c2)
        dict_benchmark = {11:{}, 12:{}}

        if args.ident:
            outputInternalSequences = "yes"
        else:
            outputInternalSequences = "no"

        tree_simu = reptree + "/tree.nhx"
        if args.bl_noise:
            tree_fn_estim = reptree + "/noisy_tree.nhx"
            treeconv_fn_estim = reptree + "/noisy_tree_conv.nhx"
        else:
            tree_fn_estim = reptree + "/tree.nhx"
            treeconv_fn_estim = reptree + "/tree_conv.nhx"

        logger.debug("Tree estim : %s", tree_fn_estim)
        logger.debug("Tree_conv estim: %s", treeconv_fn_estim)

        ### on simule des sequences
        if c1!=c2:
            # Positif
            nameAC="%s_A%d_C%d"%(name0,c1,c2)
            bpp_lib.make_simul(nameAC,nodesWithAncestralModel_sim,nodesWithTransitions,nodesWithConvergentModel,c1,c2,repseq,tree_simu,repbppconfig,outputInternalSequences=outputInternalSequences, number_of_sites=Nsites, nbCAT=NbCat_Sim,cz_nodes=cz_nodes)
            AC_fasta_file = "%s/%s.fa" %(repseq, nameAC)
            if not os.path.isfile(AC_fasta_file):
                logger.error("%s does not exist", AC_fasta_file)
                sys.exit(1)
            else:
                logger.debug("%s exists", AC_fasta_file)

            # Negatif
            nameA="%s_A%d_C%d"%(name0,c1,c1)
            bpp_lib.make_simul(nameA,nodesWithAncestralModel_sim,nodesWithTransitions,nodesWithConvergentModel,c1,c1,repseq,tree_simu,repbppconfig,outputInternalSequences=outputInternalSequences, number_of_sites=Nsites, nbCAT=NbCat_Sim,cz_nodes=cz_nodes)
            A_fasta_file = "%s/%s.fa" %(repseq, nameA)
            if not os.path.isfile(A_fasta_file):
                logger.error("%s does not exist", A_fasta_file)
                sys.exit(1)
            else:
                logger.debug("%s exists", A_fasta_file)

            if args.pcoc:
                # Test optimal
                # Positif
                bpp_lib.make_estim(nameAC,nodesWithAncestralModel,nodesWithTransitions,nodesWithConvergentModel,c1,c1,repseq,tree_fn_estim,repest,repbppconfig, NBCATest=NbCat_Sim, suffix="_opt_withOneChange", OneChange = True)
                bpp_lib.make_estim(nameAC,nodesWithAncestralModel,nodesWithTransitions,nodesWithConvergentModel,c1,c2,repseq,tree_fn_estim,repest,repbppconfig, NBCATest=NbCat_Sim, suffix="_opt_withOneChange", OneChange = True)
                # Negatif
                bpp_lib.make_estim(nameA,nodesWithAncestralModel,nodesWithTransitions,nodesWithConvergentModel,c1,c1,repseq,tree_fn_estim,repest,repbppconfig, NBCATest=NbCat_Sim, suffix="_opt_withOneChange", OneChange = True)
                bpp_lib.make_estim(nameA,nodesWithAncestralModel,nodesWithTransitions,nodesWithConvergentModel,c1,c2,repseq,tree_fn_estim,repest,repbppconfig, NBCATest=NbCat_Sim, suffix="_opt_withOneChange", OneChange = True)
                # Test optimal No One Change
                # Positif
                bpp_lib.make_estim(nameAC,nodesWithAncestralModel,nodesWithTransitions,nodesWithConvergentModel,c1,c1,repseq,tree_fn_estim,repest,repbppconfig, NBCATest=NbCat_Sim, suffix="_opt_noOneChange", OneChange = False)
                bpp_lib.make_estim(nameAC,nodesWithAncestralModel,nodesWithTransitions,nodesWithConvergentModel,c1,c2,repseq,tree_fn_estim,repest,repbppconfig, NBCATest=NbCat_Sim, suffix="_opt_noOneChange", OneChange = False)
                # Negatif
                bpp_lib.make_estim(nameA,nodesWithAncestralModel,nodesWithTransitions,nodesWithConvergentModel,c1,c1,repseq,tree_fn_estim,repest,repbppconfig, NBCATest=NbCat_Sim, suffix="_opt_noOneChange", OneChange = False)
                bpp_lib.make_estim(nameA,nodesWithAncestralModel,nodesWithTransitions,nodesWithConvergentModel,c1,c2,repseq,tree_fn_estim,repest,repbppconfig, NBCATest=NbCat_Sim, suffix="_opt_noOneChange", OneChange = False)

            if args.topo:
                # Test optimal

                # Positif
                bpp_lib.make_estim(nameAC,nodesWithAncestralModel,nodesWithTransitions,nodesWithConvergentModel,c1,c1,repseq,tree_fn_estim,repest,repbppconfig, NBCATest=NbCat_Sim, suffix="_opt_noOneChange", OneChange = False)
                # Negatif
                bpp_lib.make_estim(nameA,nodesWithAncestralModel,nodesWithTransitions,nodesWithConvergentModel,c1,c1,repseq,tree_fn_estim,repest,repbppconfig, NBCATest=NbCat_Sim, suffix="_opt_noOneChange",  OneChange = False)

            if args.topo or args.pcoc:
                ### estimations C1C2 et C1C1 ave E1: 1->10 et E2: 1-> 10

                #CATcouples_est = []
                #for e1 in range(1, (NbCat_Est+1)):
                #    for e2 in range(1, (NbCat_Est+1)):
                #            CATcouples_est.append((e1,e2))

                #SampledCATcouples_est = random.sample(CATcouples_est, 5)
                #for (e1,e2) in SampledCATcouples_est:
                #    if (e1,e1) not in SampledCATcouples_est:
                #        SampledCATcouples_est.append((e1,e1))

                #for (e1,e2) in SampledCATcouples_est:
                set_e1e2 = []
                for e1 in range(1, (NbCat_Est+1)):
                    for e2 in range(1, (NbCat_Est+1)):
                        set_e1e2.append((e1,e2))
                        #e2=1
                        if (e1 == e2) and (args.topo or args.pcoc):
                            logger.debug ("Estime e1: %s e2: %s", e1, e2)
                            # Positif
                            bpp_lib.make_estim(nameAC, nodesWithAncestralModel, nodesWithTransitions, nodesWithConvergentModel, e1, e2, repseq, tree_fn_estim, repest, repbppconfig, NBCATest=NbCat_Est, suffix="_noOneChange",  OneChange = False)
                            if args.pcoc:
                                bpp_lib.make_estim(nameAC, nodesWithAncestralModel, nodesWithTransitions, nodesWithConvergentModel, e1, e2, repseq, tree_fn_estim, repest, repbppconfig, NBCATest=NbCat_Est, suffix="_withOneChange",  OneChange = True)
                            # Negatif
                            bpp_lib.make_estim(nameA, nodesWithAncestralModel, nodesWithTransitions, nodesWithConvergentModel,e1, e2, repseq, tree_fn_estim, repest, repbppconfig, NBCATest=NbCat_Est, suffix="_noOneChange",  OneChange = False)
                            if args.pcoc:
                                bpp_lib.make_estim(nameA, nodesWithAncestralModel, nodesWithTransitions, nodesWithConvergentModel,e1, e2, repseq, tree_fn_estim, repest, repbppconfig, NBCATest=NbCat_Est, suffix="_withOneChange",  OneChange = True)

                        if (e1 != e2) and args.pcoc:
                            logger.debug ("Estime e1: %s e2: %s", e1, e2)
                            # Positif
                            bpp_lib.make_estim(nameAC, nodesWithAncestralModel, nodesWithTransitions, nodesWithConvergentModel, e1 , e2, repseq, tree_fn_estim, repest, repbppconfig, NBCATest=NbCat_Est, suffix="_noOneChange",  OneChange = False)
                            bpp_lib.make_estim(nameAC, nodesWithAncestralModel, nodesWithTransitions, nodesWithConvergentModel, e1 , e2, repseq, tree_fn_estim, repest, repbppconfig, NBCATest=NbCat_Est, suffix="_withOneChange",  OneChange = True)
                            # Negatif
                            bpp_lib.make_estim(nameA, nodesWithAncestralModel, nodesWithTransitions, nodesWithConvergentModel, e1, e2, repseq, tree_fn_estim, repest, repbppconfig, NBCATest=NbCat_Est, suffix="_noOneChange",  OneChange = False)
                            bpp_lib.make_estim(nameA, nodesWithAncestralModel, nodesWithTransitions, nodesWithConvergentModel, e1, e2, repseq, tree_fn_estim, repest, repbppconfig, NBCATest=NbCat_Est, suffix="_withOneChange",  OneChange = True)

                if args.pcoc:
                    ### Calcul VP FP FN VN Model het
                    res, bilan = estim_data.dico_typechg_new(c1,c2,n_events,repest,nameAC,tree = os.path.basename(tree_filename), set_e1e2 = set_e1e2 , NbCat_Est = NbCat_Est, n_sites = Nsites, ID = date, dist_C1_C2 = dist_C1_C2)
                    l_TPFPFNTN_mod_het.extend(res)

                    for k in [11,12]:
                        dict_benchmark[k]["PCOC"] = bilan[k]["p_mean_X_OXY"]
                        dict_benchmark[k]["PC"] = bilan[k]["p_mean_X_XY"]
                        dict_benchmark[k]["OC"] = bilan[k]["p_mean_X_OX"]

                    if args.get_likelihood_summaries:
                        for p in ["p_max_OX_OXY","p_max_XY_OXY","p_mean_OX_OXY","p_mean_XY_OXY"]:
                            del bilan[12][p]

                        for k in [11,12]:
                            if k == 11:
                                c2_k = c1
                            else:
                                c2_k = c2
                            df_bilan = pd.DataFrame.from_dict(bilan[k], orient='columns', dtype=None)
                            df_bilan["pos"] = df_bilan["pos"] + 1
                            df_bilan.to_csv(replikelihoodsummary + '/likelihood_summary_A%s_C%s.pcoc.tsv' %(c1,c2_k), index=False, sep='\t')

            if args.ident:
                # Estim
                estim_data.outdiff_new(tree_annotated, nameAC, repseq, reptree, repest, c1, c2)
                estim_data.outdiff_new(tree_annotated, nameA, repseq, reptree, repest, c1, c1)

                ### Calcul VP FP FN VN obs sub
                res_sub, bilan_sub = estim_data.dico_typechg_obs_sub(c1,c2,n_events,repest,nameAC,tree = os.path.basename(tree_filename), n_sites = Nsites, ID = date, dist_C1_C2 = dist_C1_C2)
                l_TPFPFNTN_obs_sub.extend(res_sub)

                for k in [11,12]:
                    dict_benchmark[k]["Identical"] = bilan_sub[k]["p_ident"]

            if args.topo:

                nodes = [n.ND for n in tree_conv_annotated.traverse() if not n.is_root()]
                logger.debug(nodes)

                #OPT
                # Positif
                bpp_lib.make_estim_conv(nameAC,nodes,c1,repseq,treeconv_fn_estim,repest,repbppconfig,suffix="_t"+str(c1)+"_opt", NBCATest=NbCat_Sim)
                # Negatif
                bpp_lib.make_estim_conv(nameA,nodes,c1,repseq,treeconv_fn_estim,repest,repbppconfig,suffix="_t"+str(c1)+"_opt", NBCATest=NbCat_Sim)

                set_t1 = []
                for t1 in range(1, (NbCat_Est+1)):
                    set_t1.append(t1)
                    logger.debug ("Estime t1: %s ", t1)
                    # Positif
                    bpp_lib.make_estim_conv(nameAC,nodes,t1,repseq,treeconv_fn_estim,repest,repbppconfig,suffix="_t"+str(t1), NBCATest=NbCat_Est)
                    # Negatif
                    bpp_lib.make_estim_conv(nameA,nodes,t1,repseq,treeconv_fn_estim,repest,repbppconfig,suffix="_t"+str(t1), NBCATest=NbCat_Est)


                res_topo, bilan_topo = estim_data.dico_typechg_topo(c1,c2,n_events,repest,nameAC, set_t1=set_t1, tree=os.path.basename(tree_filename), n_sites=Nsites, ID=date, NbCat_Est=NbCat_Est, dist_C1_C2=dist_C1_C2)
                l_TPFPFNTN_topo.extend(res_topo)

                for k in [11,12]:
                    dict_benchmark[k]["Topological"] = bilan_topo[k]["p_mean_X_CX"]

                if args.get_likelihood_summaries:
                    for k in [11,12]:
                        if k == 11:
                            c2_k = c1
                        else:
                            c2_k = c2
                        df_bilan_topo = pd.DataFrame.from_dict(bilan_topo[k], orient='columns', dtype=None)
                        df_bilan_topo["pos"] = df_bilan_topo["pos"] + 1
                        df_bilan_topo.to_csv(replikelihoodsummary + '/likelihood_summary_A%s_C%s.topo.tsv' %(c1,c2_k), index=False, sep='\t')

        if args.plot_ali and dict_benchmark != {11:{}, 12:{}}:
            Out_11 = "%s/tree_ali_%s_%s_negative_sites.pdf"%(repplottreeali, c1, c1)
            Out_12 = "%s/tree_ali_%s_%s_positive_sites.pdf"%(repplottreeali, c1, c2)
            plot_data.make_tree_ali_detect_combi(reptree, AC_fasta_file, Out_12, cz_nodes = cz_nodes, dict_values_pcoc = dict_benchmark[12])
            plot_data.make_tree_ali_detect_combi(reptree, A_fasta_file, Out_11, cz_nodes = cz_nodes, dict_values_pcoc = dict_benchmark[11])



    if not args.no_cleanup:
        remove_folder(repest)
        if not args.no_clean_seqs:
            remove_folder(repseq)

    metadata_simu_dico["Execution_time"] = str(time.time() - start_simu)

    logger.info("END: %s", name0_info)
    return metadata_simu_dico, l_TPFPFNTN_mod_het, l_TPFPFNTN_topo, l_TPFPFNTN_obs_sub, nodesWithTransitions


num_tree = 1
metadata_tree_dico = {}
for tree_filename in lnf:
    start_tree_time = time.time()
    logger.info("START: %s", os.path.basename(tree_filename))
    OutDirNamePrefixTree = "%s/Tree_%s/" %(OutDirName, num_tree)
    repseq0 = OutDirNamePrefixTree + "/sequences"
    repest0 = OutDirNamePrefixTree + "/estimations"
    reptree0 = OutDirNamePrefixTree + "/nw_trees"
    repplottreeali0 = OutDirNamePrefixTree + "/plot_tree_ali"
    replikelihoodsummary0 = OutDirNamePrefixTree + "/likelihood_summaries"

    if not os.path.exists(OutDirNamePrefixTree):
        os.mkdir(OutDirNamePrefixTree)
    if not os.path.exists(repseq0):
        os.mkdir(repseq0)
    if not os.path.exists(repest0):
        os.mkdir(repest0)
    if not os.path.exists(reptree0):
        os.mkdir(reptree0)
    if not os.path.exists(repplottreeali0) and args.plot_ali:
        os.mkdir(repplottreeali0)
    if not os.path.exists(replikelihoodsummary0) and args.get_likelihood_summaries:
        os.mkdir(replikelihoodsummary0)

    list_to_map = [(i+1, tree_filename, OutDirNamePrefixTree,) for i in range(Nbsimul)]
    if cpu == 1:
        r=[]
        for x in list_to_map:
            r.append(mk_simu(x))
    else:
        p = multiprocessing.Pool(processes=cpus)
        pool_results = p.map_async(mk_simu, list_to_map)
        pool_results.wait()
        r = pool_results.get()

    metada_simu_global = []
    metada_simu_per_couple_het = []
    metada_simu_per_couple_topo = []
    metada_simu_per_couple_sub = []

    bilan_nodesWithTransitions = []
    for z in r:
        if z != 0:
            x, y_het, y_topo, y_sub, nodesWithTransitions = z
            metada_simu_global.append(x)
            metada_simu_per_couple_het.extend(y_het)
            metada_simu_per_couple_topo.extend(y_topo)
            metada_simu_per_couple_sub.extend(y_sub)
            bilan_nodesWithTransitions.extend(nodesWithTransitions)

    ### Write metada on the tree
    metada_simu = pd.DataFrame(metada_simu_global)
    sorted_colums = ["ScenarioID"] + [c for c in metada_simu.columns if c != "ScenarioID"]
    metada_simu = metada_simu.reindex_axis(sorted_colums, axis=1)
    metada_simu.to_csv(OutDirNamePrefixTree + "/MetadataScenarios.tsv", sep='\t', index=False)


    metada_simu_het = pd.DataFrame()
    metada_simu_topo = pd.DataFrame()
    metada_simu_sub = pd.DataFrame()
    if args.pcoc:
        metada_simu_het = pd.DataFrame(metada_simu_per_couple_het)
        #metada_simu.to_csv(OutDirNamePrefixTree + "/Scenarios_metadata_per_couple_mod_het.tsv", sep='\t', index=False)
    if args.topo:
        metada_simu_topo = pd.DataFrame(metada_simu_per_couple_topo)
        #metada_simu.to_csv(OutDirNamePrefixTree + "/Scenarios_metadata_per_couple_topo.tsv", sep='\t', index=False)
    if args.ident:
        metada_simu_sub = pd.DataFrame(metada_simu_per_couple_sub)
        #metada_simu.to_csv(OutDirNamePrefixTree + "/Scenarios_metadata_per_couple_obs_sub.tsv", sep='\t', index=False)

    df_concat = [df for df in [metada_simu_het, metada_simu_topo, metada_simu_sub] if not df.empty]
    if df_concat:
        df_cat = pd.concat(df_concat)
        #df_cat = df_cat["RunID","InputTree","ScenarioID","SimuCoupleID","C1","C2","DistanceSimuCouple","Method","Threshold","FN","FP","TN","TP","Sensitivity","Specificity","MCC","NumberOfConvergentEvents","NumberOfSites","PosteriorProbabilityType"]]
        df_cat = df_cat[["RunID","InputTree","ScenarioID","SimuCoupleID","C1","C2","DistanceSimuCouple","Method","Threshold","FN","FP","TN","TP","Sensitivity","Specificity","MCC","NumberOfConvergentEvents","NumberOfSites"]]
        df_cat.to_csv(OutDirNamePrefixTree + "/BenchmarkResults.tsv", sep='\t', index=False)

    if args.plot_event_repartition:
        plot_data.mk_bilan_tree(events_placing.init_tree(tree_filename), bilan_nodesWithTransitions, OutDirNamePrefixTree + "/Tree_"+str(num_tree)+".pdf")

        #script_dirname = os.path.dirname(os.path.abspath(__file__))
        #R_command="Rscript %s/rscripts/mk_sens_spe_MCC_plot.R %s %s" %(script_dirname, os.environ['PWD'] + "/" + OutDirNamePrefixTree, os.environ['PWD'] + "/" +OutDirNamePrefixTree)
        #logger.info(R_command)
        #p = subprocess.Popen(R_command.split(" "))
        #(out_p, err_p) = p.communicate()
        #logger.debug("out: %s \n err: %s", out_p, err_p)

    metadata_tree_dico["Tree_%s" %(num_tree)] = [os.path.basename(tree_filename), str(time.time() - start_tree_time)]

    logger.info("END: %s (%s seconds)", os.path.basename(tree_filename), str(time.time() - start_tree_time))
    logger.info("%s/%s scenario succeded", len(r), Nbsimul)
    num_tree +=1

### Write metada on the run
## Reorder metadata dict
restructured_metadata_tree_dico = {}
restructured_metadata_tree_dico["Tree_#"] = metadata_tree_dico.keys()
restructured_metadata_tree_dico["Tree_filename"] = [metadata_tree_dico[k][0] for k in restructured_metadata_tree_dico["Tree_#"]]
restructured_metadata_tree_dico["Execution time"] = [metadata_tree_dico[k][1] for k in restructured_metadata_tree_dico["Tree_#"]]
df_tree = pd.DataFrame.from_dict(restructured_metadata_tree_dico)
## Reorder columns
df_tree = df_tree[["Tree_#", "Tree_filename", "Execution time"]]

df_tree.to_csv(OutDirName + "/tree_metadata.tsv", sep='\t', index=False)

if not args.no_cleanup:
    remove_folder(repbppconfig)
    remove_folder(repest0)
    if not args.no_clean_seqs:
        remove_folder(repseq0)


end_time = str(time.time() - start_time)
logger.info("--- %s seconds ---", end_time)
logger.info("Output dir: %s", OutDirName)

