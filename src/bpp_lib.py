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
logger = logging.getLogger("pcoc.bpp_lib")


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

    number_of_models = 0

    command="bppseqgen FASTA_OUT=%s TREE=%s number_of_sites=%s output.internal.sequences=%s " %(fasta_outfile, tree_fn, number_of_sites, outputInternalSequences)

    if sim_profiles.name in ["C10","C60"]:
        command += " NBCAT=%s " %(sim_profiles.nb_cat)
    else:
        if not os.path.isfile(sim_profiles.formatted_frequencies_filename):
            logger.error("%s is not a file", sim_profiles.formatted_frequencies_filename)
        command += " PROFILE_F=%s " %(sim_profiles.formatted_frequencies_filename)
    
    command+=" param=%s.bpp Ne1=%d Ne2=%d " %(repbppconfig+"/CATseq_sim",c1,c2)

    if sim_profiles.name in ["C10","C60"]:
        command+=" modelA=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(c1)
        command+=" modelC=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(c2)
        command+=" modelOC=\'OneChange(model=LGL08_CAT_C%s(nbCat=$(NBCAT)))\' " %(c2)
    #else:
        #command+=" modelA=\'LG08+F(name = Fixed ,fitness=Empirical(file=$(PROFILE_F), col=$(Ne1)))\' " %(c1)
        #command+=" modelC=\'LG08+F(name = Fixed ,fitness=Empirical(file=$(PROFILE_F), col=$(Ne1)))\' " %(c2)

    command += " \'nonhomogeneous.root_freq=FromModel(model=$(modelA))\' "

    if c1!=c2:
        if nodesWithAncestralModel:
            number_of_models +=1
            command+=" model%s=\'$(modelA)\' " %(number_of_models)
            command+=" model%s.nodes_id=\'%s\' " %(number_of_models,",".join(map(str, nodesWithAncestralModel)))

        if nodesWithTransitions:
            number_of_models +=1
            command+=" model%s=\'$(modelOC)\' " %(number_of_models)
            command+=" model%s.nodes_id=\'%s\' " %(number_of_models,",".join(map(str, nodesWithTransitions)))

        if nodesWithConvergentModel:
            number_of_models +=1
            command+=" model%s=\'$(modelC)\' " %(number_of_models)
            command+=" model%s.nodes_id=\'%s\' " %(number_of_models,",".join(map(str, nodesWithConvergentModel)))

    else:
        allNodes = nodesWithConvergentModel+nodesWithTransitions+nodesWithAncestralModel
        number_of_models +=1
        command+=" model%s=\'$(modelA)\' " %(number_of_models)
        command+=" model%s.nodes_id=\'%s\' " %(number_of_models,",".join(map(str, allNodes)))

    # If noisy profiles
    sup_command = ""
    if cz_nodes:
        for (cz, nodes) in cz_nodes.items():
            if nodes:
                if CzOneChange:
                    number_of_models +=1
                    if sim_profiles.name in ["C10","C60"]:
                        sup_command+=" model%s=\'OneChange(model=LGL08_CAT_C%s(nbCat=$(NBCAT)))\' " %(number_of_models, cz)
                    #else:
                        #
                    t_node = nodes[0]
                    sup_command+=" model%s.nodes_id=\'%s\' " %(number_of_models,str(t_node))
                    if len(nodes) > 1:
                        number_of_models +=1
                        if sim_profiles.name in ["C10","C60"]:
                            sup_command+=" model%s=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(number_of_models, cz)
                        #else:
                            #
                        sup_command+=" model%s.nodes_id=\'%s\' " %(number_of_models,",".join(map(str,nodes[1:])))
                else:
                    number_of_models +=1
                    if sim_profiles.name in ["C10","C60"]:
                        sup_command+=" model%s=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(number_of_models, cz)
                    #else:
                        #
                    sup_command+=" model%s.nodes_id=\'%s\' " %(number_of_models,",".join(map(str, nodes)))

    command+=sup_command + " nonhomogeneous.number_of_models=%s " %(number_of_models)

    out = commands.getoutput(command)

    if debug_mode_bpp:
        logger.debug("%s\n%s\n%s", command, out, command)


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
    #logger.debug("fasta_file: %s",fasta_file )

    if not os.path.isfile(tree_fn):
        logger.error("%s is not a file", tree_fn)
    if not os.path.isfile(fasta_file):
        logger.error("%s is not a file", fasta_file)

    command = "bppml param=%s NAME=%s SUFFIX=%s REP_SEQ=%s TREE=%s REP_EST=%s \'input.sequence.file=%s\' " %(repbppconfig+"/CATseq_estim.bpp", name, suffix, repseq, tree_fn, repest, fasta_file)
    number_of_models = 0

    if est_profiles.name in ["C10","C60"]:
        command += " NBCAT=%s " %(est_profiles.nb_cat)
    else:
        if not os.path.isfile(est_profiles.formatted_frequencies_filename):
            logger.error("%s is not a file", est_profiles.formatted_frequencies_filename)
        command += " PROFILE_F=%s " %(est_profiles.formatted_frequencies_filename)


    command += "Ne1=%d Ne2=%d" %(c1, c2)
    if est_profiles.name in ["C10","C60"]:
        command+=" modelA=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(c1)
        command+=" modelC=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(c2)
        command+=" modelOC=\'OneChange(model=LGL08_CAT_C%s(nbCat=$(NBCAT)))\' " %(c2)
    #else:
        #
        #
        #

    command += " \'nonhomogeneous.root_freq=FromModel(model=$(modelA))\' "

    if nodesWithAncestralModel:
        number_of_models +=1
        command+=" model%s=\'$(modelA)\' " %(number_of_models)
        command+=" model%s.nodes_id=\'%s\' " %(number_of_models,",".join(map(str, nodesWithAncestralModel)))

    if OneChange:
        if nodesWithConvergentModel:
            number_of_models +=1
            command+=" \'model%s=$(modelC)\' " %(number_of_models)
            command+=" model%s.nodes_id=\'%s\' " %(number_of_models,",".join(map(str, nodesWithConvergentModel)))

        if nodesWithTransitions:
            number_of_models +=1
            # Mixture
            relproba1 = 0.9
            command+=" \'model%s=Mixture(model1=$(modelOC),model2=$(modelC),relproba1=%d)\' " %(number_of_models, relproba1)
            command+=" model%s.nodes_id=\'%s\' " %(number_of_models,",".join(map(str, nodesWithTransitions)))

    else:
        nodesWithTransitionsAndWithConvergentModel = nodesWithTransitions+nodesWithConvergentModel
        if nodesWithTransitionsAndWithConvergentModel:
            number_of_models +=1
            command+=" model%s=\'$(modelC)\' " %(number_of_models)
            command+=" model%s.nodes_id=\'%s\' " %(number_of_models,",".join(map(str, nodesWithTransitionsAndWithConvergentModel)))

    if gamma:
        command+=" rate_distribution=\'Gamma(n=4)\' "
    elif inv_gamma:
        command+=" rate_distribution=\'Invariant(dist=Gamma(n=4))\' "
    else:
        command+=" rate_distribution=\'Constant()\' "

    if 0 <= max_gap_allowed <=100:
        command+=" MAX_GAP_ALLOWED=%s " %(max_gap_allowed)
    else:
        logger.error("max_gap_allowed (%s) must be between 0 and 100", max_gap_allowed)
        sys.error(1)

    command += "nonhomogeneous.number_of_models=%d " %(number_of_models)

    if debug_mode_bpp:
        logger.debug("%s", command)

    out = commands.getoutput(command)

    if debug_mode_bpp:
        logger.debug("%s\n%s", out, command)

    info_filename = "%s/%s_%s_%s%s.infos" %(repest, name , c1, c2, suffix)

    if re.search("^Number of sites retained.*: 0$",out,re.MULTILINE) or \
       re.search("^Number of sites.*: 0$",out,re.MULTILINE):
        logger.warning("No site retained for %s (too much gaps), you can use the \"--max_gap_allowed\" option.", name)
        f_infos = open(info_filename,"w")
        f_infos.close()

    if not os.path.exists(info_filename):
        logger.error("%s does not exist", info_filename)
        logger.error("command: %s\nout:\n%s", command, out)
        sys.exit(42)

def make_estim_conv(name, c1, g_tree, est_profiles, suffix="", gamma = False, max_gap_allowed=90):

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

    command="bppml param=%s NAME=%s SUFFIX=%s REP_SEQ=%s TREE=%s REP_EST=%s \'input.sequence.file=%s\' "%(repbppconfig + "/CATseq_conv.bpp", name, suffix, repseq, tree_fn, repest, fasta_file)
    number_of_models = 0

    if est_profiles.name in ["C10","C60"]:
        command += " NBCAT=%s " %(est_profiles.nb_cat)
    else:
        if not os.path.isfile(est_profiles.formatted_frequencies_filename):
            logger.error("%s is not a file", est_profiles.formatted_frequencies_filename)
        command += " PROFILE_F=%s " %(est_profiles.formatted_frequencies_filename)

    command += " Ne1=%d " %(c1)

    if est_profiles.name in ["C10","C60"]:
        command+=" modelA=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(c1)
    #else:
        #

    command += " \'nonhomogeneous.root_freq=FromModel(model=$(modelA))\' "

    if allNodes:
        number_of_models +=1
        if est_profiles.name in ["C10","C60"]:
            command+=" model%s=\'$(modelA)\' " %(number_of_models)
        command+=" model%s.nodes_id=\'%s\' " %(number_of_models,"\'"+ ",".join(map(str, allNodes))+"\'")

    if gamma:
        command+=" rate_distribution=\'Gamma(n=4)\' "
    else:
        command+=" rate_distribution=\'Constant()\' "

    if 0 <= max_gap_allowed <=100:
        command+=" MAX_GAP_ALLOWED=%s " %(max_gap_allowed)
    else:
        logger.error("max_gap_allowed (%s) must be between 0 and 100", max_gap_allowed)
        sys.error(1)

    command += " nonhomogeneous.number_of_models=%d " %(number_of_models)

    out = commands.getoutput(command)

    if debug_mode_bpp:
        logger.debug("%s\n%s\n%s", command, out, command)

    if re.search("^Number of sites retained.*: 0$",out,re.MULTILINE) or \
       re.search("^Number of sites.*: 0$",out,re.MULTILINE):
        info_filename = "%s/%s_topo%s.infos" %(repest, name , suffix)
        logger.warning("No site retained for %s", info_filename)
        f_infos = open(info_filename,"w")
        f_infos.close()


########################################################################
##                     Configuration files                            ##
########################################################################

#==> CATseq_estim.bpp <==
CATseq_estim= """alphabet=Protein
input.tree.file=$(TREE)
input.tree.format=Nhx

input.sequence.sites_to_use=all
input.sequence.max_gap_allowed=$(MAX_GAP_ALLOWED)%

nonhomogeneous = general

output.tree.file=$(REP_EST)/$(NAME)$(SUFFIX).dnd

### estimation
input.sequence.remove_saturated_sites=yes
optimization.ignore_parameters=*

output.infos=$(REP_EST)/$(NAME)_$(Ne1)_$(Ne2)$(SUFFIX).infos
output.estimates=$(REP_EST)/$(NAME)_$(Ne1)_$(Ne2)$(SUFFIX).params
"""

#==> CATseq_conv.bpp <==
CATseq_conv_bpp = """alphabet=Protein
input.tree.file=$(TREE)
input.tree.format=Nhx

input.sequence.sites_to_use=all
input.sequence.max_gap_allowed=$(MAX_GAP_ALLOWED)%

nonhomogeneous = general

### estimation
input.sequence.remove_saturated_sites=yes
optimization.ignore_parameters=*

output.infos=$(REP_EST)/$(NAME)_topo$(SUFFIX).infos
output.estimates=$(REP_EST)/$(NAME)_topo$(SUFFIX).params
"""

#==> CATseq_sim.bpp <==
CATseq_sim_bpp = """alphabet=Protein

input.tree.file=$(TREE)
input.tree.format=Nhx

nonhomogeneous = general

#rate_distribution=Constant()
rate_distribution=Gamma(n=4)

### simulation
output.sequence.file=$(FASTA_OUT)
"""
