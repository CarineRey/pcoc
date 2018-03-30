#!/usr/bin/python
#  pcoc_det.py
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
import sys
import argparse
import os
import re
import time
import logging
import math

import bpp_lib
import events_placing
import plot_data
import estim_data

#import multiprocessing

import pandas as pd
from ete3 import Tree

from Bio import AlignIO, SeqIO, Seq, SeqRecord

import shutil

##########
# inputs #
##########
start_time = time.time()

### Option defining
parser = argparse.ArgumentParser(prog="pcoc_det.py",
                                 description='')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

parser._optionals.title = "MISCELLANEOUS"
#parser.add_argument('-cpu', type=int,
#                    help="Number of cpu to use. (default: 1)",
#                    default=1)

##############
requiredOptions = parser.add_argument_group('REQUIRED OPTIONS')
requiredOptions.add_argument('-t', "--tree", type=str,
                             help='Input tree filename', required=True)
requiredOptions.add_argument('-aa', "--ali", type=str,
                             help='Input amino-acid aligment filename', required=True)
requiredOptions.add_argument('-m', '--manual_mode', type=str, metavar="\"x/y,z/...\"",
                    help="User defined convergent transition/branches. Transition node must be the first number and independent events must be separed by a \"/\". ex: \"1,2,3/67/55,56\" (default: None)",
                    required=True)
requiredOptions.add_argument('-o', '--output_dir', type=str,
                   help="Output directory name", required=True)
##############


##############
BasicOptions = parser.add_argument_group('OPTIONS FOR BASIC USAGE')

BasicOptions.add_argument('-f', '--filter_t',type=float,
                    help="ALL model: Posterior probability threshold to put result in \"filtered\" results. (default: 0.99)",
                    default=0.99)
BasicOptions.add_argument('-f_pcoc', '--filter_t_pcoc',type=float,
                    help="PCOC model: Posterior probability threshold to put result in \"filtered\" results. If = -1, take the value of -f, if > 1, discard this model. (default: -1)",
                    default=-1)
BasicOptions.add_argument('-f_pc', '--filter_t_pc',type=float,
                    help="PC model: Posterior probability threshold to put result in \"filtered\" results. If = -1, take the value of -f, if > 1, discard this model.(default: -1)",
                    default=-1)
BasicOptions.add_argument('-f_oc', '--filter_t_oc',type=float,
                    help="OC model: Posterior probability threshold to put result in \"filtered\" results. If = -1, take the value of -f, if > 1, discard this model.(default: -1)",
                    default=-1)
BasicOptions.add_argument('-ph', type=str,
                    help="Add these positions in the filtered position and highlight them with a star in the plot",
                    default=False)
BasicOptions.add_argument('--plot', action="store_true",
                    help="Plot the tree and the filtered sites of the alignment with their corresponding score.",
                    default=False)
BasicOptions.add_argument('--plot_complete_ali', action="store_true",
                    help="Plot the tree and each site of alignment with its corresponding score. (Can take time to be openned)",
                    default=False)
BasicOptions.add_argument('--reorder', action="store_true",
                    help="reorder the filtered plot by score.categories (>= 0.99, >=0.9, >= 0.8, < 0.8)",
                    default=False)
BasicOptions.add_argument('--svg', action="store_true",
                    help="additional svg output plots.",
                    default=False)
BasicOptions.add_argument('--no_cleanup_fasta', action="store_true",
                    help="Do not cleanup the fasta directory after the run.",
                    default=False)
##############


##############
AdvancedOptions = parser.add_argument_group('OPTIONS FOR ADVANCED USAGE')
AdvancedOptions.add_argument('-CATX_est', type=int, choices = [10,60],
                    help="Profile categorie to estimate data (10->C10 or 60->C60). (default: 10)",
                    default=10)
AdvancedOptions.add_argument('--gamma', action="store_true",
                    help="Use rate_distribution=Gamma(n=4) instead of Constant()",
                    default=False)
AdvancedOptions.add_argument('--inv_gamma', action="store_true",
                    help="Use rate_distribution=Gamma(n=4) instead of Constant()",
                    default=False)
AdvancedOptions.add_argument('--max_gap_allowed', type=int,
                    help="max gap allowed to take into account a site (in %%), must be between 0 and 100",
                    default=30)
AdvancedOptions.add_argument('--no_cleanup', action="store_true",
                    help="Do not cleanup the working directory after the run.",
                    default=False)
AdvancedOptions.add_argument("-LD_LIB", metavar='LD_LIBRARY_PATH', type=str, default="",
                   help="Redefine the LD_LIBRARY_PATH env variable, bppsuite library must be present in the $PATH and in the LD_LIBRARY_PATH")
AdvancedOptions.add_argument('--debug', action="store_true",
                    help="debug mode",
                    default=False)

### Option parsing
args = parser.parse_args()

positions_to_highlight = []
if args.ph:
    positions_to_highlight = args.ph.split(",")
    if all([x.isdigit() for x in positions_to_highlight]):
        positions_to_highlight = map(int, positions_to_highlight)
    else:
        print ("ERROR: %s can not be convert in numeric, are positions spaced by comma ?" %(positions_to_highlight))
        sys.exit(1)

metadata_run_dico = {}
date = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
OutDirName = "%s/RUN_%s/" %(args.output_dir, date)
OutDirName = OutDirName.replace("//","/")


metadata_run_dico["date"] = date

### Set up the output directory
if os.path.isdir(OutDirName):
    pass
    #logger.info("The output directory %s exists", OutDirName)
elif OutDirName: # if OutDirName is not a empty string we create the directory
    #logger.info("The output directory %s does not exist, it will be created", OutDirName)
    os.makedirs(OutDirName)

### Set up the log file
LogFile = OutDirName + "/pcoc_det.log"

### Set up the logger
# create logger
logger = logging.getLogger("pcoc")
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
#fh = logging.FileHandler(LogFile)
# create console handler with a higher log level
ch = logging.StreamHandler()
if args.debug:
    ch.setLevel(logging.DEBUG)
#    fh.setLevel(logging.DEBUG)
else:
    ch.setLevel(logging.INFO)
#    fh.setLevel(logging.INFO)
# create formatter and add it to the handlers
formatter_fh = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
formatter_ch = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
#fh.setFormatter(formatter_fh)
ch.setFormatter(formatter_ch)
# add the handlers to the logger
#logger.addHandler(fh)
logger.addHandler(ch)

logger.debug(sys.argv)

#cpu = args.cpu
#try:
#    cpus = multiprocessing.cpu_count()
#except NotImplementedError:
#    cpus = 1   # arbitrary default
#logger.info("%s on %s cpus", cpu, cpus)

if args.LD_LIB:
    logger.info("$LD_LIBRARY_PATH will be change from %s to %s", os.environ.get("LD_LIBRARY_PATH", ""), args.LD_LIB)
    os.environ["LD_LIBRARY_PATH"]=args.LD_LIB
else:
    logger.debug("$LD_LIBRARY_PATH is %s", os.environ.get("LD_LIBRARY_PATH", ""))


repbppconfig = OutDirName +  "bpp_config"
if not os.path.exists(repbppconfig):
    os.mkdir(repbppconfig)

bpp_lib.write_config(repbppconfig, estim=True)

#http://stackoverflow.com/questions/23172293/use-python-to-extract-branch-lengths-from-newick-format
pattern = re.compile(r"\b[0-9]+(?:\.[0-9]+)?\b")


logger.info("alignment: %s", args.ali)
ali_filename = args.ali

l_n_sites = []

if os.path.isfile(ali_filename):
    ali = AlignIO.read(ali_filename, "fasta")
else:
    logger.error("%s does not exist", ali_filename)
    sys.exit(1)

logger.info("alignment ok after checking")


n_sites = ali.get_alignment_length()
nb_seq = len(ali)


logger.info("tree: %s", args.tree)
tree_filename = args.tree
# test if a tree
try:
    t=Tree(tree_filename)
except:

    logger.error("%s is not a newick tree, this tree can not be used",tree_filename)
    t=""
if t:
    treefile=open(tree_filename,"r")
    tree=treefile.read().strip()
    treefile.close()
    #test if branch length
    branch_lengths = pattern.findall(tree)
    if branch_lengths == []:
        logger.error("No branch length in %s, this tree can not be used",tree_filename)

if not (t and branch_lengths):
    sys.exit(1)

logger.info("tree ok after checking")

manual_mode_nodes = {}
if args.manual_mode:
    manual_mode_nodes = {"T":[],"C":[]}
    p_events = args.manual_mode.strip().split("/")
    for e in p_events:
        l_e = map(int,e.split(","))
        manual_mode_nodes["T"].append(l_e[0])
        manual_mode_nodes["C"].extend(l_e[1:])

NbCat_Est = args.CATX_est
metadata_run_dico["NbCat_Est"] = NbCat_Est
metadata_run_dico["Tree"] = os.path.basename(tree_filename)
metadata_run_dico["Alignment"] = os.path.basename(ali_filename)

metadata_run_dico["Position_to_highlight"] = positions_to_highlight
metadata_run_dico["Convergent_branches"] = p_events

metadata_run_dico["gamma"] = args.gamma
metadata_run_dico["inv_gamma"] = args.inv_gamma
metadata_run_dico["max_gap_allowed"] = args.max_gap_allowed

if not (0 <= args.max_gap_allowed <= 100):
    logger.error("max_gap_allowed (%s) must be between 0 and 100", max_gap_allowed)
    sys.error(1)

p_filter_threshold = args.filter_t

dict_p_filter_threshold = {}
dict_p_filter_threshold["PCOC"] = p_filter_threshold
dict_p_filter_threshold["PC"] = p_filter_threshold
dict_p_filter_threshold["OC"] = p_filter_threshold

if args.filter_t_pcoc >= 0:
    dict_p_filter_threshold["PCOC"] = args.filter_t_pcoc
if args.filter_t_pc >= 0:
    dict_p_filter_threshold["PC"] = args.filter_t_pc
if args.filter_t_oc >= 0:
    dict_p_filter_threshold["OC"] = args.filter_t_oc


metadata_run_dico["pp_threshold_PCOC"] = dict_p_filter_threshold["PCOC"]
metadata_run_dico["pp_threshold_PC"] = dict_p_filter_threshold["PC"]
metadata_run_dico["pp_threshold_OC"] = dict_p_filter_threshold["OC"]



prefix_out = OutDirName + "/" + os.path.splitext(os.path.basename(ali_filename))[0]


pd.Series(metadata_run_dico).to_csv(OutDirName + "/run_metadata.tsv", sep='\t')

def remove_folder(path):
    # check if folder exists
    if os.path.exists(path):
        # remove if exists
        shutil.rmtree(path)

def filter_l(l, pos):
    new_l = [None]*len(pos)
    for i in range(len(pos)):
        p = pos[i]
        new_l[i] = l[p-1]
    return new_l

def reorder_l(l, order):
    new_l = [None]*len(order)
    for i in range(len(order)):
        p = order[i]
        new_l[p] = l[i]
    return new_l

def mk_detect(tree_filename, ali_basename, OutDirName):
    start_detec = time.time()
    metadata_simu_dico = {}
    logger.debug("Tree: %s", os.path.basename(tree_filename))
    metadata_simu_dico["tree"] = os.path.basename(tree_filename)

    g_tree = events_placing.gene_tree(tree_filename, manual_mode_nodes)
    g_tree.init_tree_det(n_sites)

    metadata_simu_dico["numberOfLeaves"] = g_tree.numberOfLeafs

    g_tree.init_inter_dir_det(repest0, reptree0, repfasta0, repbppconfig, repseq)

    logger.debug("repfasta: %s", g_tree.repfasta)
    logger.debug("repest: %s", g_tree.repest)
    logger.debug("reptree: %s", g_tree.reptree)

    ### construit les arbres d'etude :
    (allbranchlength, convbranchlength) = g_tree.mk_tree_for_simu(plot=args.plot)

    metadata_simu_dico["allbranchlength"] = allbranchlength
    metadata_simu_dico["convbranchlength"] = convbranchlength


    l_TPFPFNTN_mod_het = []
    l_TPFPFNTN_topo = []
    l_TPFPFNTN_obs_sub = []

    if not os.path.isfile(g_tree.repseq + "/" + ali_basename):
        logger.error("%s does not exist", g_tree.repseq + "/" + ali_basename)
        sys.exit(1)

    c1 = 1  # useless but compatibility
    c2 = 2  # useless but compatibility

    set_e1e2 = []
    for e1 in range(1, (NbCat_Est+1)):
        for e2 in range(1, (NbCat_Est+1)):
            set_e1e2.append((e1,e2))
    for (e1,e2) in set_e1e2:
        logger.debug ("Estime e1: %s e2: %s", e1, e2)
        # Positif
#$      bpp_lib.make_estim(ali_basename, nodesWithAncestralModel, nodesWithTransitions, nodesWithConvergentModel, e1, e2, repseq, tree_fn, repest, repbppconfig, NBCATest=NbCat_Est, suffix="_noOneChange",  OneChange = False, ext="", max_gap_allowed=args.max_gap_allowed, gamma=args.gamma, inv_gamma=args.inv_gamma)
        bpp_lib.make_estim(ali_basename, e1, e2, g_tree, NBCATest=NbCat_Est, suffix="_noOneChange",  OneChange = False, ext="", max_gap_allowed=args.max_gap_allowed, gamma=args.gamma, inv_gamma=args.inv_gamma)
#$      bpp_lib.make_estim(ali_basename, nodesWithAncestralModel, nodesWithTransitions, nodesWithConvergentModel, e1, e2, repseq, tree_fn, repest, repbppconfig, NBCATest=NbCat_Est, suffix="_withOneChange",  OneChange = True, ext="", max_gap_allowed=args.max_gap_allowed, gamma=args.gamma, inv_gamma=args.inv_gamma)
        bpp_lib.make_estim(ali_basename, e1, e2, g_tree, NBCATest=NbCat_Est, suffix="_withOneChange",  OneChange = True, ext="", max_gap_allowed=args.max_gap_allowed, gamma=args.gamma, inv_gamma=args.inv_gamma)


    ### post proba
    res, bilan = estim_data.dico_typechg_het_det(ali_basename,  g_tree, set_e1e2 = set_e1e2 , NbCat_Est = NbCat_Est, ID = date)
    l_TPFPFNTN_mod_het.extend(res)

    for p in ["p_max_OX_OXY","p_max_XY_OXY","p_mean_OX_OXY","p_mean_XY_OXY"]:
        if bilan[12].has_key(p):
            del bilan[12][p]

    dict_values_pcoc = {}
    dict_values_pcoc["PCOC"] = bilan[12]["p_mean_X_OXY"]
    dict_values_pcoc["PC"] = bilan[12]["p_mean_X_XY"]
    dict_values_pcoc["OC"] = bilan[12]["p_mean_X_OX"]

    # filter position:

    bilan_f = {}
    all_pos = range(1, n_sites +1)

    dict_pos_filtered = {}
    for model in ["PCOC", "PC", "OC"]:
        dict_pos_filtered[model] = [p for p in all_pos if dict_values_pcoc[model][p-1] >= dict_p_filter_threshold[model] ]
        if positions_to_highlight:
            dict_pos_filtered[model].extend(positions_to_highlight)
            dict_pos_filtered[model] = list(set(dict_pos_filtered[model]))
            dict_pos_filtered[model].sort()

    # filter dict_values_pcoc
    dict_values_pcoc_filtered = {}
    all_filtered_position = list(set(events_placing.unlist(dict_pos_filtered.values())))
    all_filtered_position.sort()
    dict_pos_filtered["union"] = all_filtered_position

    if args.reorder:
        for model in ["PCOC", "PC", "OC", "union"]:
            m_list = [model]
            if model == "union":
                m_list = ["PCOC", "PC", "OC"]
            nb_filtered_pos = len(dict_pos_filtered[model])
            new_order = [0]*nb_filtered_pos
            j = 0
            # 0.99
            for i in range(nb_filtered_pos):
                p = dict_pos_filtered[model][i]
                if any([dict_values_pcoc[m] [p-1] >= 0.99 for m in m_list]):
                    new_order[i] = j
                    j+=1
            # 0.9
            for i in range(nb_filtered_pos):
                p = dict_pos_filtered[model][i]
                if any([0.99 > dict_values_pcoc[m] [p-1] >= 0.9 for m in m_list]) and \
                   all([0.99 > dict_values_pcoc[m] [p-1] for m in m_list]):
                    new_order[i] = j
                    j+=1
            # 0.8
            for i in range(nb_filtered_pos):
                p = dict_pos_filtered[model][i]
                if any([ 0.9 > dict_values_pcoc[m] [p-1] >= 0.8 for m in m_list]) and \
                   all([ 0.9 > dict_values_pcoc[m] [p-1] for m in m_list]):
                    new_order[i] = j
                    j+=1
            # other
            for i in range(nb_filtered_pos):
                p = dict_pos_filtered[model][i]
                if all([ dict_values_pcoc[m] [p-1] < 0.8 for m in m_list]):
                    new_order[i] = j
                    j+=1
            dict_pos_filtered[model] = reorder_l(dict_pos_filtered[model], new_order)


    # filtered ali:
    ## Per model
    for model in ["PCOC", "PC", "OC", "union"]:
        filtered_ali = []
        for seq in ali:
            new_seq = SeqRecord.SeqRecord(Seq.Seq("".join(filter_l(list(seq.seq),dict_pos_filtered[model]))), seq.id, "", "")
            filtered_ali.append(new_seq)
        SeqIO.write(filtered_ali, g_tree.repfasta+"/filtered_ali."+model+".faa", "fasta")
        if model == "union":
            modelstr = "union"
        else:
            modelstr = model
        logger.info("%s model: # filtered position: %s/%s", modelstr.upper(), len(dict_pos_filtered[model]), n_sites)

    ## Output
    ### Table
    #### complete:
    df_bilan = pd.DataFrame.from_dict(dict_values_pcoc, orient='columns', dtype=None)
    df_bilan["pos"] = all_pos
    #### filtered:
    df_bilan_f = df_bilan[df_bilan.pos.isin(all_filtered_position)]
    df_bilan_f = df_bilan_f.copy()

    df_bilan.to_csv(prefix_out +  ".results.tsv", index=False, sep='\t')
    df_bilan_f.to_csv(prefix_out + ".filtered_results.tsv", index=False, sep='\t')


    ### Plot
    if args.plot:
        if args.plot_complete_ali:
            plot_data.make_tree_ali_detect_combi(g_tree, g_tree.repseq + "/" + ali_basename, prefix_out+"_plot_complete.pdf", dict_benchmark=dict_values_pcoc, hp=positions_to_highlight)
            if args.svg:
                plot_data.make_tree_ali_detect_combi(g_tree, g_tree.repseq + "/" + ali_basename, prefix_out+"_plot_complete.svg", dict_benchmark=dict_values_pcoc, hp=positions_to_highlight)


        for model in ["PCOC", "PC", "OC"]:
            if dict_pos_filtered[model] and dict_p_filter_threshold[model] <=1:
                dict_values_pcoc_filtered_model = {}
                for (key, val) in dict_values_pcoc.items():
                    dict_values_pcoc_filtered_model[key] = filter_l(val, dict_pos_filtered[model])
                plot_data.make_tree_ali_detect_combi(g_tree, g_tree.repfasta+"/filtered_ali."+model+".faa", prefix_out+"_plot_filtered_"+model+".pdf", hist_up = model, dict_benchmark = dict_values_pcoc_filtered_model, x_values= dict_pos_filtered[model], hp=positions_to_highlight, reorder = args.reorder)
                if args.svg:
                     plot_data.make_tree_ali_detect_combi(g_tree, g_tree.repfasta+"/filtered_ali."+model+".faa", prefix_out+"_plot_filtered_"+model+".svg", hist_up = model, dict_benchmark = dict_values_pcoc_filtered_model, x_values= dict_pos_filtered[model], hp=positions_to_highlight, reorder = args.reorder)

        # all model
        if dict_pos_filtered["union"]:
            model = "union"
            dict_values_pcoc_filtered_model = {}
            for (key, val) in dict_values_pcoc.items():
                dict_values_pcoc_filtered_model[key] = filter_l(val, dict_pos_filtered[model])
            plot_data.make_tree_ali_detect_combi(g_tree, g_tree.repfasta+"/filtered_ali."+model+".faa", prefix_out+"_plot_filtered_"+model+".pdf", dict_benchmark = dict_values_pcoc_filtered_model, x_values= dict_pos_filtered[model], hp=positions_to_highlight, reorder = False)
            if args.svg:
                plot_data.make_tree_ali_detect_combi(g_tree, g_tree.repfasta+"/filtered_ali."+model+".faa", prefix_out+"_plot_filtered_"+model+".svg", dict_benchmark = dict_values_pcoc_filtered_model, x_values= dict_pos_filtered[model], hp=positions_to_highlight, reorder = False)


    if not args.no_cleanup:
        remove_folder(g_tree.repest)
        remove_folder(g_tree.repbppconfig)
        remove_folder(g_tree.reptree)
        if not args.no_cleanup_fasta:
            remove_folder(g_tree.repfasta)

    metadata_simu_dico["time"] = str(time.time() - start_detec)

    return metadata_simu_dico, l_TPFPFNTN_mod_het, l_TPFPFNTN_topo, l_TPFPFNTN_obs_sub, g_tree.conv_events.nodesWithTransitions


if __name__ == "__main__":
    num_tree = 1
    metadata_tree_dico = {}
    start_tree_time = time.time()
    #logger.info("START: %s", os.path.basename(tree_filename))
    repseq = os.path.dirname(ali_filename)
    if not repseq:
        repseq = "."

    #logger.info("alignment directory: %s", repseq)

    ali_basename = os.path.basename(ali_filename)
    repest0 = OutDirName + "/Estimations"
    reptree0 = OutDirName + "/Trees"
    repfasta0 = OutDirName + "/fasta"

    if not os.path.exists(repest0):
        os.mkdir(repest0)
    if not os.path.exists(reptree0):
        os.mkdir(reptree0)
    if not os.path.exists(repfasta0):
        os.mkdir(repfasta0)

    mk_detect(tree_filename, ali_basename, OutDirName)

    logger.info("--- %s seconds ---", str(time.time() - start_time))
    logger.info("Output dir: %s", OutDirName)


