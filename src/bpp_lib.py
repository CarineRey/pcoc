#  bpp_lib.py
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

#import subprocess
import commands
import os
import sys,re

import logging

from profile_tools import CATC10Distances, CATC60Distances
logger = logging.getLogger("pcoc.bpp_lib")

import pandas as pd

debug_mode_bpp = False

########################################################################
##                    Functions                                       ##
########################################################################

def write_global_config(d, estim=True, sim_profiles_name=""):
    bpp_config_files_estim = [
                   ("CATseq_estim.bpp", CATseq_estim),
                   ("CATseq_conv.bpp", CATseq_conv_bpp)]
    
    bpp_config_files_sim = [("CATseq_sim.bpp", CATseq_sim_bpp)]


    files_list = bpp_config_files_sim
    if estim:
        files_list.extend(bpp_config_files_estim)

    if sim_profiles_name == "C10":
        files_list.append(("CATC10Distances.csv", CATC10Distances))
    elif sim_profiles_name == "C60":
        files_list.append(("CATC60Distances.csv", CATC60Distances))

    for (f, s) in files_list:

        with open(d+"/"+f, "w") as F:
            F.write(s)


########################################################################
##                    BPP simulations                                 ##
########################################################################

def make_simul(name, c1, c2, g_tree, sim_profiles,
               number_of_sites=1000,
               outputInternalSequences="yes",
               gamma = True,
               cz_nodes={}, CzOneChange=True):

    nodesWithAncestralModel  = g_tree.conv_events.nodesWithAncestralModel_sim
    nodesWithTransitions     = g_tree.conv_events.nodesWithTransitions_sim
    nodesWithConvergentModel = g_tree.conv_events.nodesWithConvergentModel_sim
    repseq                   = g_tree.repseq
    repbppconfig             = g_tree.repbppconfig
    tree_fn                  = g_tree.tree_fn_sim
    cz_nodes                 = g_tree.cz_nodes

    if not os.path.isfile(tree_fn):
        logger.error("%s is not a file", tree_fn)

    fasta_outfile = "%s/%s%s" %(repseq.replace("//","/"), name, ".fa")

    if outputInternalSequences != "yes":
        outputInternalSequences = "no"

    command_bppseqgen = "bppseqgen FASTA_OUT=%s TREE=%s NUMBER_OF_SITES=%s OUTPUT_INTERNAL_SEQUENCES=%s "\
            %(fasta_outfile, tree_fn, number_of_sites, outputInternalSequences)
    
    command_bppseqgen+=" param=%s Ne1=%d Ne2=%d " %(repbppconfig+"/CATseq_sim.bpp",c1,c2)

    if gamma:
        command_bppseqgen+=" RATE_DISTRIBUTION=\'Gamma(n=4)\' "
    else:
        command_bppseqgen+=" RATE_DISTRIBUTION=\'Constant()\' "
    
    number_of_models = 0
    nonhomogeneous_models=""

    if sim_profiles.name in ["C10","C60"]:
        command_bppseqgen += " NBCAT=%s " %(sim_profiles.nb_cat)
    else:
        if not os.path.isfile(sim_profiles.formatted_frequencies_filename):
            logger.error("%s is not a file", sim_profiles.formatted_frequencies_filename)
        command_bppseqgen += " PROFILE_F=%s " %(sim_profiles.formatted_frequencies_filename)
    
    if sim_profiles.name in ["C10","C60"]:
        command_bppseqgen+=" MODEL_A=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(c1)
        command_bppseqgen+=" MODEL_C=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(c2)
        command_bppseqgen+=" MODEL_OC=\'OneChange(model=$(MODEL_C))\' "
    else:
        command_bppseqgen+=" MODEL_A=\'LG08+F(frequencies=Empirical(file=$(PROFILE_F), col=$(Ne1)))\' "
        command_bppseqgen+=" MODEL_C=\'LG08+F(frequencies=Empirical(file=$(PROFILE_F), col=$(Ne2)))\' "
        command_bppseqgen+=" MODEL_OC=\'OneChange(model=$(MODEL_C))\' "


    command_bppseqgen += " \'ROOT_FREQ=FromModel(model=$(MODEL_A))\' "

    if c1!=c2:
        if nodesWithAncestralModel:
            number_of_models +=1
            command_bppseqgen+=" model%s=\'$(MODEL_A)\' " %(number_of_models)
            nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
            nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str, nodesWithAncestralModel)))

        if nodesWithTransitions:
            number_of_models +=1
            command_bppseqgen+=" model%s=\'$(MODEL_OC)\' " %(number_of_models)
            nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
            nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str, nodesWithTransitions)))

        if nodesWithConvergentModel:
            number_of_models +=1
            command_bppseqgen+=" model%s=\'$(MODEL_C)\' " %(number_of_models)
            nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
            nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str, nodesWithConvergentModel)))

    else:
        allNodes = nodesWithConvergentModel+nodesWithTransitions+nodesWithAncestralModel
        number_of_models +=1
        command_bppseqgen+=" model%s=\'$(MODEL_A)\' " %(number_of_models)
        nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
        nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str, allNodes)))
        

    # If noisy profiles
    sup_command_bppseqgen = ""
    if cz_nodes:
        for (cz, nodes) in cz_nodes.items():
            if nodes:
                if CzOneChange:
                    number_of_models +=1

                    if sim_profiles.name in ["C10","C60"]:
                        sup_command_bppseqgen+=" model%s=\'OneChange(model=LGL08_CAT_C%s(nbCat=$(NBCAT)))\' " %(number_of_models, cz)
                    else:
                        sup_command_bppseqgen+=" model%s=\'OneChange(model=LG08+F(frequencies=Empirical(file=$(PROFILE_F), col=%s)))\' " %(number_of_models, cz)
                    t_node = nodes[0]
                    nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
                    nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,str(t_node))

                    if len(nodes) > 1:
                        number_of_models +=1

                        if sim_profiles.name in ["C10","C60"]:
                            sup_command_bppseqgen+=" model%s=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(number_of_models, cz)
                        else:
                            sup_command_bppseqgen+=" model%s=\'LG08+F(frequencies=Empirical(file=$(PROFILE_F), col=%s))\' " %(number_of_models, cz)
                        nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
                        nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str,nodes[1:])))
                else:
                    number_of_models +=1
                    if sim_profiles.name in ["C10","C60"]:
                        sup_command_bppseqgen+=" model%s=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(number_of_models, cz)
                    else:
                        sup_command_bppseqgen+=" model%s=\'LG08+F(frequencies=Empirical(file=$(PROFILE_F), col=%s))\' " %(number_of_models, cz)
                    nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
                    nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str, nodes)))

    
    command_bppseqgen+=" %s NONHOMOGENEOUS_MODELS=\"%s\" " %(sup_command_bppseqgen, nonhomogeneous_models)

    out = commands.getoutput(command_bppseqgen)

    if debug_mode_bpp:
        logger.debug("%s\n%s\n%s", command_bppseqgen, out, command_bppseqgen)

    if not os.path.exists(fasta_outfile):
        logger.error("%s does not exist", fasta_outfile)
        logger.error("command_bppseqgen: %s\nout:\n%s\ncommand_bppseqgen: %s\n", command_bppseqgen, out,command_bppseqgen)
        sys.exit(42)


########################################################################
##                    BPP estimations                                 ##
########################################################################

def make_estim(name, c1, c2, g_tree, est_profiles, suffix="",
               OneChange=True, ext=".fa", gamma=False,
               max_gap_allowed=90, inv_gamma=False):

    nodesWithAncestralModel  = g_tree.conv_events.nodesWithAncestralModel_est
    nodesWithTransitions     = g_tree.conv_events.nodesWithTransitions_est
    nodesWithConvergentModel = g_tree.conv_events.nodesWithConvergentModel_est
    tree_fn                  = g_tree.tree_fn_est
    repseq                   = g_tree.repseq
    repest                   = g_tree.repest
    repbppconfig             = g_tree.repbppconfig

    fasta_file = "%s/%s%s" %(repseq, name, ext)

    if not os.path.isfile(tree_fn):
        logger.error("%s is not a file", tree_fn)
    if not os.path.isfile(fasta_file):
        logger.error("%s is not a file", fasta_file)

    output_infos     =  "%s/%s_%s_%s%s.infos"  %(repest, name , c1, c2, suffix)
    output_params    = "%s/%s_%s_%s%s.params" %(repest, name , c1, c2, suffix)
    output_tree      = "%s/%s_%s.dnd" %(repest, name , suffix)

    command_bppml = "bppml param=%s output.infos=%s output.estimates=%s output.tree.file=%s  TREE=%s \'FILESEQ=%s\' " \
        %(repbppconfig+"/CATseq_estim.bpp", output_infos, output_params, output_tree, tree_fn, fasta_file)

    number_of_models = 0
    nonhomogeneous_models=""

    if est_profiles.name in ["C10","C60"]:
        command_bppml += " NBCAT=%s " %(est_profiles.nb_cat)
    else:
        if not os.path.isfile(est_profiles.formatted_frequencies_filename):
            logger.error("%s is not a file", est_profiles.formatted_frequencies_filename)
        command_bppml += " PROFILE_F=%s " %(est_profiles.formatted_frequencies_filename)


    command_bppml += "Ne1=%d Ne2=%d" %(c1, c2)

    if est_profiles.name in ["C10","C60"]:
        command_bppml+=" MODEL_A=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(c1)
        command_bppml+=" MODEL_C=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(c2)
        command_bppml+=" MODEL_OC=\'OneChange(model=$(MODEL_C))\' "
    else:
        command_bppml+=" MODEL_A=\'LG08+F(frequencies=Empirical(file=$(PROFILE_F), col=$(Ne1)))\' "
        command_bppml+=" MODEL_C=\'LG08+F(frequencies=Empirical(file=$(PROFILE_F), col=$(Ne2)))\' "
        command_bppml+=" MODEL_OC=\'OneChange(model=$(MODEL_C))\' "

    command_bppml += " \'ROOT_FREQ=FromModel(model=$(MODEL_A))\' "

    if nodesWithAncestralModel:
        number_of_models +=1
        command_bppml+=" model%s=\'$(MODEL_A)\' " %(number_of_models)
        nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
        nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str, nodesWithAncestralModel)))


    if OneChange:
        if nodesWithConvergentModel:
            number_of_models +=1
            command_bppml+=" \'model%s=$(MODEL_C)\' " %(number_of_models)
            nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
            nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str, nodesWithConvergentModel)))

        if nodesWithTransitions:
            number_of_models +=1
            command_bppml+=" \'model%s=$(MODEL_OC)\' " %(number_of_models)
            nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
            nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str, nodesWithTransitions)))


    else:
        nodesWithTransitionsAndWithConvergentModel = nodesWithTransitions+nodesWithConvergentModel
        if nodesWithTransitionsAndWithConvergentModel:
            number_of_models +=1
            command_bppml+=" model%s=\'$(MODEL_C)\' " %(number_of_models)
            nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
            nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str, nodesWithTransitionsAndWithConvergentModel)))

    if gamma:
        command_bppml+=" RATE_DISTRIBUTION=\'Gamma(n=4)\' "
    elif inv_gamma:
        command_bppml+=" RATE_DISTRIBUTION=\'Invariant(dist=Gamma(n=4))\' "
    else:
        command_bppml+=" RATE_DISTRIBUTION=\'Constant()\' "

    if 0 <= max_gap_allowed <=100:
        command_bppml+=" MAX_GAP_ALLOWED=%s " %(max_gap_allowed)
    else:
        logger.error("max_gap_allowed (%s) must be between 0 and 100", max_gap_allowed)
        sys.error(1)

    command_bppml += "Ne1=%d Ne2=%d NONHOMOGENEOUS_MODELS=\"%s\"" %(c1, c2, nonhomogeneous_models)

    if debug_mode_bpp:
        logger.debug("%s", command_bppml)

    out = commands.getoutput(command_bppml)

    output_infos = "%s_1" %(output_infos)

    if re.search("^Number of sites retained.*: 0$",out,re.MULTILINE) or \
       re.search("^Number of sites.*: 0$",out,re.MULTILINE):
        logger.error("No site retained for %s (too many gaps), you can use the \"--max_gap_allowed\" option.", name)
        f_infos = open(output_infos,"w")
        f_infos.close()

    if not os.path.exists(output_infos):
        logger.error("%s does not exist", output_infos)
        logger.error("command_bppml: %s\nout:\n%s", command_bppml, out)
        sys.exit(42)
        
    ### Read outputs ###
    logger.debug("Read and save likelihoods (%s, %s)", c1, c2)
    ## bppml ##
    logger.debug("read output_infos from %s", output_infos)
    col_list=["Sites", "lnL"]
    df_bppml = pd.read_csv(output_infos, sep = '\s+',  usecols = col_list, header = 0)
    #logger.debug("bppml: %s", df_bppml.head().to_string())
    #logger.debug("bppml: %s", list(df_bppml))
    
    df_bppml = df_bppml.rename(columns={"lnL": "lnl"})
    #logger.debug("bppml: %s", df_bppml.head().to_string())
    df_bppml["C1"] = c1
    df_bppml["C2"] = c2
    df_bppml["OneChange"] = OneChange
    df_bppml["T_C1"] = None
    df_bppml["T_C2"] = None
    #logger.debug("bppml: %s", df_bppml.head().to_string())
    df_bppml["Sites"] = df_bppml["Sites"].str.replace("[","").str.replace("]","")
    df_bppml["Sites"] = pd.to_numeric(df_bppml["Sites"]) + 1
    
    return(df_bppml)

def make_estim_conv_topo(name, c1, g_tree, est_profiles, suffix="", gamma = False, max_gap_allowed=90):

    repseq        = g_tree.repseq
    repest        = g_tree.repest
    repbppconfig  = g_tree.repbppconfig

    tree_fn       = g_tree.treeconv_fn_est
    
    allNodes = [n.ND for n in g_tree.tree_conv_annotated.traverse() if not n.is_root()]
    logger.debug(allNodes)

    fasta_file = "%s/%s%s" %(repseq, name, ".fa")
    if not os.path.isfile(tree_fn):
        logger.error("%s is not a file", tree_fn)
    if not os.path.isfile(fasta_file):
        logger.error("%s is not a file", fasta_file)

    output_infos  = "%s/%s_topo%s.infos"   %(repest, name, suffix)
    output_params = "%s/%s_topo%s.params" %(repest, name, suffix)
    output_tree      = "%s/%s_%s.dnd" %(repest, name , suffix)

    command_bppml="bppml param=%s output.infos=%s output.estimates=%s output.tree.file=%s TREE=%s \'FILESEQ=%s\' " \
        %(repbppconfig + "/CATseq_conv.bpp", output_infos, output_params, output_tree, tree_fn, fasta_file)
    
    number_of_models = 0

    if est_profiles.name in ["C10","C60"]:
        command_bppml += " NBCAT=%s " %(est_profiles.nb_cat)
    else:
        if not os.path.isfile(est_profiles.formatted_frequencies_filename):
            logger.error("%s is not a file", est_profiles.formatted_frequencies_filename)
        command_bppml += " PROFILE_F=%s " %(est_profiles.formatted_frequencies_filename)

    command_bppml += " Ne1=%d " %(c1)

    if est_profiles.name in ["C10","C60"]:
        command_bppml+=" MODEL_A=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(c1)
    else:
        command_bppml+=" MODEL_A=\'LG08+F(frequencies=Empirical(file=$(PROFILE_F), col=$(Ne1)))\' "

    command_bppml += " \'ROOT_FREQ=FromModel(model=$(MODEL_A))\' "

    if allNodes:
        number_of_models +=1
        command_bppml+=" model%s=\'$(MODEL_A)\' " %(number_of_models)

    if gamma:
        command_bppml+=" RATE_DISTRIBUTION=\'Gamma(n=4)\' "
    else:
        command_bppml+=" RATE_DISTRIBUTION=\'Constant()\' "

    if 0 <= max_gap_allowed <=100:
        command_bppml+=" MAX_GAP_ALLOWED=%s " %(max_gap_allowed)
    else:
        logger.error("max_gap_allowed (%s) must be between 0 and 100", max_gap_allowed)
        sys.error(1)

    out = commands.getoutput(command_bppml)

    if debug_mode_bpp:
        logger.debug("%s\n%s\n%s", command_bppml, out, command_bppml)
    
    output_infos = "%s_1" %(output_infos)

    if not os.path.exists(output_infos):
        logger.error("%s does not exist", output_infos)
        logger.error("command_bppml: %s\nout:\n%s", command_bppml, out)
        sys.exit(42)

    if re.search("^Number of sites retained.*: 0$",out,re.MULTILINE) or \
       re.search("^Number of sites.*: 0$",out,re.MULTILINE):
        logger.warning("No site retained for %s", output_infos)
        f_infos = open(output_infos,"w")
        f_infos.close()


########################################################################
##                     Configuration files                            ##
########################################################################

#==> CATseq_estim.bpp <==
CATseq_estim= """alphabet=Protein

input.tree1=user(file=$(TREE), format=NHX)
input.data1=alignment(file=$(FILESEQ),\
                      sites_to_use=all,\
                      max_gap_allowed=$(MAX_GAP_ALLOWED)\
                      )
input.sequence.remove_saturated_sites = yes

root_freq1=$(ROOT_FREQ)

rate_distribution1=$(RATE_DISTRIBUTION)

optimization.ignore_parameters=BrLen*

process1=NonHomogeneous($(NONHOMOGENEOUS_MODELS), tree=1, root_freq=1, rate=1)
phylo1=Single(process=1, data=1)

result=phylo1

"""

#==> CATseq_conv.bpp <==
CATseq_conv_bpp = """alphabet=Protein

input.tree1=user(file=$(TREE), format=NHX)

input.data1=alignment(file=$(FILESEQ),\
                          sites_to_use=all,\
                          max_gap_allowed=$(MAX_GAP_ALLOWED)\
                          )
input.sequence.remove_saturated_sites = yes

root_freq1=$(ROOT_FREQ)
rate_distribution1=$(RATE_DISTRIBUTION)

process1=Homogeneous(model=1, tree=1, root_freq=1, rate=1)
#process{number} = {Homogeneous|OnePerBranch}(tree={number}, model={number}, rate={number} [, root_freq={number}])

phylo1=Single(process=1, data=1)

result=phylo1

### estimation
optimization.ignore_parameters=BrLen*
"""

#==> CATseq_sim.bpp <==
CATseq_sim_bpp = """alphabet=Protein

input.tree1=user(file=$(TREE), format=NHX)

root_freq1=$(ROOT_FREQ)
rate_distribution1=$(RATE_DISTRIBUTION)

process1=NonHomogeneous($(NONHOMOGENEOUS_MODELS), tree=1, root_freq=1, rate=1)

### simulation
#simul{int}={Simulation type}(process={number}, output.sequence.file={file path}, number_of_sites = {int>0}[,output.sequence.format={alignement format}, output.internal.sequences = true])
simul1=simul(process=1, output.sequence.file=$(FASTA_OUT), number_of_sites = $(NUMBER_OF_SITES),\
             output.sequence.format=Fasta, output.internal.sequences = $(OUTPUT_INTERNAL_SEQUENCES))

"""
