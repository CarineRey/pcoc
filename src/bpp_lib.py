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
            %(fasta_outfile, tree_fn, nbCAT,number_of_sites, outputInternalSequences)
    
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
        command_bppseqgen+=" modelA=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(c1)
        command_bppseqgen+=" modelC=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(c2)
        command_bppseqgen+=" modelOC=\'OneChange(model=LGL08_CAT_C%s(nbCat=$(NBCAT)))\' " %(c2)
    #else:
        #command+=" modelA=\'LG08+F(name = Fixed ,fitness=Empirical(file=$(PROFILE_F), col=$(Ne1)))\' " %(c1)
        #command+=" modelC=\'LG08+F(name = Fixed ,fitness=Empirical(file=$(PROFILE_F), col=$(Ne1)))\' " %(c2)

    command_bppseqgen += " \'ROOT_FREQ=FromModel(model=$(modelA))\' "

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
    sup_command_bppseqgen = ""
    if cz_nodes:
        for (cz, nodes) in cz_nodes.items():
            if nodes:
                if CzOneChange:
                    number_of_models +=1

                    if sim_profiles.name in ["C10","C60"]:
                        sup_command_bppseqgen+=" model%s=\'OneChange(model=LGL08_CAT_C%s(nbCat=$(NBCAT)))\' " %(number_of_models, cz)
                    #else:
                            #

                    t_node = nodes[0]
                    nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
                    nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,str(t_node))

                    if len(nodes) > 1:
                        number_of_models +=1

                        if sim_profiles.name in ["C10","C60"]:
                            sup_command_bppseqgen+=" model%s=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(number_of_models, cz)
                        #else:
                            #
                        
                        nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
                        nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str, nodes[1:])))
                else:
                    number_of_models +=1
                    if sim_profiles.name in ["C10","C60"]:
                        sup_command_bppseqgen+=" model%s=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(number_of_models, cz)
                    #else:
                        #
                    
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
    #logger.debug("fasta_file: %s",fasta_file )

    if not os.path.isfile(tree_fn):
        logger.error("%s is not a file", tree_fn)
    if not os.path.isfile(fasta_file):
        logger.error("%s is not a file", fasta_file)

    command_bppml = "bppml param=%s NAME=%s SUFFIX=%s REP_SEQ=%s TREE=%s REP_EST=%s FILESEQ=%s  " \
        %(repbppconfig+"/CATseq_estim.bpp", name, suffix, repseq, tree_fn, repest, fasta_file)
    
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
        command_bppml+=" modelA=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(c1)
        command_bppml+=" modelC=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(c2)
        command_bppml+=" modelOC=\'OneChange(model=LGL08_CAT_C%s(nbCat=$(NBCAT)))\' " %(c2)
    #else:
        #
        #
        #

    command_bppml += " \'ROOT_FREQ=FromModel(model=$(modelA))\' "

    if nodesWithAncestralModel:
        number_of_models +=1
<<<<<<< HEAD
        command_bppml+=" model%s=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(number_of_models, c1)
        nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
        nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str, nodesWithAncestralModel)))
=======
        command+=" model%s=\'$(modelA)\' " %(number_of_models)
        command+=" model%s.nodes_id=\'%s\' " %(number_of_models,",".join(map(str, nodesWithAncestralModel)))
>>>>>>> 79d8f9b... big refactoring

    if OneChange:
        if nodesWithConvergentModel:
            number_of_models +=1
            command+=" \'model%s=$(modelC)\' " %(number_of_models)
            command+=" model%s.nodes_id=\'%s\' " %(number_of_models,",".join(map(str, nodesWithConvergentModel)))

        if nodesWithTransitions:
            number_of_models +=1
<<<<<<< HEAD
            command_bppml+=" model%s=\'OneChange(model=LGL08_CAT_C%s(nbCat=$(NBCAT)))\' " %(number_of_models, c2)
            nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
            nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str, nodesWithTransitions)))

        if nodesWithConvergentModel:
            number_of_models +=1
            command_bppml+=" model%s=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(number_of_models, c2)
            nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
            nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str, nodesWithConvergentModel)))
=======
            # Mixture
            relproba1 = 0.9
            command+=" \'model%s=Mixture(model1=$(modelOC),model2=$(modelC),relproba1=%d)\' " %(number_of_models, relproba1)
            command+=" model%s.nodes_id=\'%s\' " %(number_of_models,",".join(map(str, nodesWithTransitions)))

>>>>>>> 79d8f9b... big refactoring
    else:
        nodesWithTransitionsAndWithConvergentModel = nodesWithTransitions+nodesWithConvergentModel
        if nodesWithTransitionsAndWithConvergentModel:
            number_of_models +=1
<<<<<<< HEAD
            command_bppml+=" model%s=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(number_of_models, c2)
            nonhomogeneous_models+="model%s=%s," %(number_of_models, number_of_models)
            nonhomogeneous_models+="model%s.nodes_id=(%s)," %(number_of_models,",".join(map(str, nodesWithTransitionsAndWithConvergentModel)))
=======
            command+=" model%s=\'$(modelC)\' " %(number_of_models)
            command+=" model%s.nodes_id=\'%s\' " %(number_of_models,",".join(map(str, nodesWithTransitionsAndWithConvergentModel)))
>>>>>>> 79d8f9b... big refactoring

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

<<<<<<< HEAD
    command_bppml += "Ne1=%d Ne2=%d NONHOMOGENEOUS_MODELS=\"%s\"" %(c1, c2, nonhomogeneous_models)

    if debug_mode_bpp:
        logger.debug("%s", command_bppml)
    out = commands.getoutput(command_bppml)
=======
    command += "nonhomogeneous.number_of_models=%d " %(number_of_models)

    if debug_mode_bpp:
        logger.debug("%s", command)

    out = commands.getoutput(command)
>>>>>>> 79d8f9b... big refactoring

    if debug_mode_bpp:
        logger.debug("%s\n%s", out, command_bppml)

    info_filename = "%s/%s_%s_%s%s.infos_1" %(repest, name , c1, c2, suffix)

    if re.search("^Number of sites retained.*: 0$",out,re.MULTILINE) or \
       re.search("^Number of sites.*: 0$",out,re.MULTILINE):
        logger.error("No site retained for %s (too many gaps), you can use the \"--max_gap_allowed\" option.", name)
        f_infos = open(info_filename,"w")
        f_infos.close()

    if not os.path.exists(info_filename):
        logger.error("%s does not exist", info_filename)
        logger.error("command_bppml: %s\nout:\n%s", command_bppml, out)
        sys.exit(42)

<<<<<<< HEAD
def make_estim_conv_topo(name, c1, g_tree, suffix="", NBCATest=10, gamma = False, max_gap_allowed=90):

=======
def make_estim_conv(name, c1, g_tree, est_profiles, suffix="", gamma = False, max_gap_allowed=90):
>>>>>>> 79d8f9b... big refactoring

    repseq        = g_tree.repseq
    repest        = g_tree.repest
    repbppconfig  = g_tree.repbppconfig

    tree_fn       = g_tree.treeconv_fn_est

    fasta_file = "%s/%s%s" %(repseq, name, ".fa")
    if not os.path.isfile(tree_fn):
        logger.error("%s is not a file", tree_fn)
    if not os.path.isfile(fasta_file):
        logger.error("%s is not a file", fasta_file)

<<<<<<< HEAD
    command="bppml param=%s NAME=%s SUFFIX=%s REP_SEQ=%s TREE=%s REP_EST=%s Ne1=%d NBCAT=%s FILESEQ=%s " \
        %(repbppconfig + "/CATseq_conv.bpp", name, suffix, repseq, tree_fn, repest, c1, NBCATest, fasta_file)
=======
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
>>>>>>> 79d8f9b... big refactoring

    if gamma:
        command+=" RATE_DISTRIBUTION=\'Gamma(n=4)\' "
    else:
        command+=" RATE_DISTRIBUTION=\'Constant()\' "

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

input.tree1=user(file=$(TREE), format=NHX)
input.data1=alignment(file=$(FILESEQ),\
                      sites_to_use=all,\
                      max_gap_allowed=$(MAX_GAP_ALLOWED)\
                      )
input.sequence.remove_saturated_sites = yes

root_freq1=$(ROOT_FREQ)

rate_distribution1=$(RATE_DISTRIBUTION)

process1=NonHomogeneous($(NONHOMOGENEOUS_MODELS), tree=1, root_freq=1, rate=1)
phylo1=Single(process=1, data=1)

result=phylo1

### estimation
optimization.ignore_parameters=BrLen*

output.tree.file=$(REP_EST)/$(NAME)$(SUFFIX).dnd
output.infos=$(REP_EST)/$(NAME)_$(Ne1)_$(Ne2)$(SUFFIX).infos
output.estimates=$(REP_EST)/$(NAME)_$(Ne1)_$(Ne2)$(SUFFIX).params
"""

#==> CATseq_conv.bpp <==
CATseq_conv_bpp = """alphabet=Protein
<<<<<<< HEAD

input.tree1=user(file=$(TREE), format=NHX)

input.data1=alignment(file=$(FILESEQ),\
                          sites_to_use=all,\
                          max_gap_allowed=$(MAX_GAP_ALLOWED)\
                          )
input.sequence.remove_saturated_sites = yes

root_freq1=FromModel(model=LGL08_CAT_C$(Ne1)(nbCat=$(NBCAT)))
rate_distribution1=$(RATE_DISTRIBUTION)

model1=LGL08_CAT_C$(Ne1)(nbCat=$(NBCAT))

process1=Homogeneous(model=1, tree=1, root_freq=1, rate=1)
#process{number} = {Homogeneous|OnePerBranch}(tree={number}, model={number}, rate={number} [, root_freq={number}])

phylo1=Single(process=1, data=1)

result=phylo1
=======
input.tree.file=$(TREE)
input.tree.format=Nhx

input.sequence.sites_to_use=all
input.sequence.max_gap_allowed=$(MAX_GAP_ALLOWED)%

nonhomogeneous = general
>>>>>>> 79d8f9b... big refactoring

### estimation
optimization.ignore_parameters=BrLen*

output.infos=$(REP_EST)/$(NAME)_topo$(SUFFIX).infos
output.estimates=$(REP_EST)/$(NAME)_topo$(SUFFIX).params
"""

#==> CATseq_sim.bpp <==
CATseq_sim_bpp = """alphabet=Protein

<<<<<<< HEAD
input.tree1=user(file=$(TREE), format=NHX)

root_freq1=FromModel(model=LGL08_CAT_C$(Ne1)(nbCat=$(NBCAT)))
rate_distribution1=$(RATE_DISTRIBUTION)

process1=NonHomogeneous($(NONHOMOGENEOUS_MODELS), tree=1, root_freq=1, rate=1)
=======
input.tree.file=$(TREE)
input.tree.format=Nhx

nonhomogeneous = general

#rate_distribution=Constant()
rate_distribution=Gamma(n=4)
>>>>>>> 79d8f9b... big refactoring

### simulation
#simul{int}={Simulation type}(process={number}, output.sequence.file={file path}, number_of_sites = {int>0}[,output.sequence.format={alignement format}, output.internal.sequences = true])
simul1=simul(process=1, output.sequence.file=$(FASTA_OUT), number_of_sites = $(NUMBER_OF_SITES),\
             output.sequence.format=Fasta, output.internal.sequences = $(OUTPUT_INTERNAL_SEQUENCES))

"""
<<<<<<< HEAD

#==> CATseq_sim_noconv.bpp <==
CATseq_sim_noconv_bpp = """alphabet=Protein

input.tree1=user(file=$(TREE), format=NHX)

root_freq1=FromModel(model=LGL08_CAT_C$(Ne1)(nbCat=$(NBCAT)))
rate_distribution1=$(RATE_DISTRIBUTION)

process1=NonHomogeneous($(NONHOMOGENEOUS_MODELS), tree=1, root_freq=1, rate=1)

### simulation
simul1=simul(process=1, output.sequence.file=$(FASTA_OUT), number_of_sites = $(NUMBER_OF_SITES),\
             output.sequence.format=Fasta, output.internal.sequences = $(OUTPUT_INTERNAL_SEQUENCES))
"""


########################################################################
##                     CAT profiles distance                          ##
########################################################################

CATC10Distances = '''"","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10"
"C1",0,0.394936382844135,0.647911483590955,0.504741739857853,0.649936426685023,0.462502041970525,0.455330902048981,0.680835730034952,0.500267443565405,0.354142373558561
"C2",0.394936382844135,0,0.366929020923609,0.399520579337121,0.388681433281979,0.347808726198259,0.26742130880874,0.5467207096324,0.328990117070772,0.308902043578295
"C3",0.647911483590955,0.366929020923609,0,0.646797515574788,0.484219935483386,0.624318750987495,0.552790608645897,0.7100404194913,0.612602719459113,0.616142326500867
"C4",0.504741739857853,0.399520579337121,0.646797515574788,0,0.61498560351171,0.288643152604794,0.334011322172709,0.624833201460655,0.393308051004086,0.357862515199004
"C5",0.649936426685023,0.388681433281979,0.484219935483386,0.61498560351171,0,0.596077269209793,0.482582338428004,0.621079788639518,0.572826917164944,0.592224093993954
"C6",0.462502041970525,0.347808726198259,0.624318750987495,0.288643152604794,0.596077269209793,0,0.336112756102295,0.620605426684893,0.33971233033462,0.365907677902604
"C7",0.455330902048981,0.26742130880874,0.552790608645897,0.334011322172709,0.482582338428004,0.336112756102295,0,0.435041188980131,0.294206074959596,0.330589548366017
"C8",0.680835730034952,0.5467207096324,0.7100404194913,0.624833201460655,0.621079788639518,0.620605426684893,0.435041188980131,0,0.609005758193588,0.615419950246057
"C9",0.500267443565405,0.328990117070772,0.612602719459113,0.393308051004086,0.572826917164944,0.33971233033462,0.294206074959596,0.609005758193588,0,0.374899594976228
"C10",0.354142373558561,0.308902043578295,0.616142326500867,0.357862515199004,0.592224093993954,0.365907677902604,0.330589548366017,0.615419950246057,0.374899594976228,0'''

CATC60Distances = '''"","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25","C26","C27","C28","C29","C30","C31","C32","C33","C34","C35","C36","C37","C38","C39","C40","C41","C42","C43","C44","C45","C46","C47","C48","C49","C50","C51","C52","C53","C54","C55","C56","C57","C58","C59","C60"
"C1",0,0.468943978251089,0.619511922258355,0.32347011811786,0.445737874237941,0.31335445949712,0.606726332690357,0.444551051309446,0.372059461519443,0.268617844518868,0.495078736444999,0.35315693411213,0.254395263725802,0.454699047384902,0.165111214335744,0.275429947834589,0.378270406240024,0.607527928544704,0.272778696198201,0.464520962229137,0.922904076761031,0.355915670018763,0.497210289333075,0.519613666536877,0.395831778071236,0.262267634048403,0.182394187978795,0.293360968911815,0.451497305131059,0.336982989891302,0.258333776269069,0.63482529361811,0.395228491377792,0.672316921435624,0.447465543711839,0.336930328032101,0.73795235532303,0.191688240475808,0.525529840194857,0.642760642989431,0.406278615153105,0.191890213701116,0.525490195468797,0.26045452596257,0.314338398215097,0.422469938234662,0.272260510094909,0.49916708408605,0.356829269016759,0.1611209933261,0.707703642001276,0.162426809915323,0.325281292459037,0.315977926132789,0.244399925474771,0.461349065525982,0.329864217978574,0.272855751053774,0.401024644052398,0.606190439625681
"C2",0.468943978251089,0,0.703605662603301,0.567139429192973,0.574328670142162,0.493689489094861,0.701164798323324,0.671674038736365,0.523681319136652,0.446481958367255,0.653834729807725,0.361943428905081,0.591813216820769,0.645676011255589,0.38296217638241,0.477235256717359,0.535973764688221,0.692534365773873,0.254622944112007,0.645138302283394,0.985455653195309,0.440088797891112,0.612978653642644,0.624536553117282,0.682789077003093,0.432895061459954,0.411602547913554,0.211689481572857,0.588906295898231,0.497133135773248,0.582757282686316,0.719001661544569,0.420206918087509,0.755908253953561,0.615088635249367,0.689368904605218,0.928666808278951,0.438364577199164,0.615622628465977,0.724735221146926,0.622745599385895,0.583802724212561,0.130241905555658,0.533413102275447,0.593335553821491,0.438368181572418,0.462269717991225,0.688719383412591,0.344127673408234,0.456812522342311,0.791995081913527,0.516810469026954,0.427349119927189,0.538597676778906,0.375466548975976,0.585879273467451,0.535713504998622,0.513971652196677,0.695483967143157,0.692518216381593
"C3",0.619511922258355,0.703605662603301,0,0.494062077168436,0.41687771356073,0.611044644516638,0.692178363905846,0.718191111633515,0.561871729485542,0.439036466101009,0.747742589316601,0.512220228558914,0.632119860234518,0.455764058053165,0.608382243572234,0.641535956453588,0.480220443607987,0.298072784762931,0.629071044587284,0.732786418868623,1.01781986064006,0.345969417118892,0.661134308362314,0.151550352639521,0.761091614722705,0.553145074029308,0.573146973516457,0.576501616463269,0.225878923100233,0.40861959659323,0.675855640235749,0.568692357456592,0.39312475789478,0.264895811762571,0.703994216037174,0.754619078596702,0.978139487223776,0.612247756514232,0.708136363070036,0.53718454294183,0.716781395473256,0.679429445577233,0.750273374475198,0.496862935656964,0.699248806340843,0.617326615358582,0.631491229152073,0.787922284912928,0.626669796656057,0.603430069636162,0.825316236083171,0.57710938231357,0.633469741649031,0.427878453350798,0.542649284934117,0.711063013597631,0.680569021917312,0.629262663250619,0.740619508725804,0.388212048625072
"C4",0.32347011811786,0.567139429192973,0.494062077168436,0,0.369451965233533,0.435981050211927,0.593339749718081,0.327169712971302,0.393071607850036,0.274341741645925,0.600076675029399,0.394325910841101,0.259722323860086,0.255779804837107,0.341726107943767,0.44822397111334,0.354727641089797,0.510236439142285,0.424958374275872,0.591982575039002,0.931937304639618,0.302963467150492,0.51513223045837,0.393811440770416,0.508920388212989,0.340499296373118,0.322113139244766,0.386274363416391,0.298458542570305,0.291829693341455,0.444168445894821,0.579067730641874,0.367617399097779,0.536028271551357,0.510328462301058,0.391417021535622,0.846474965044117,0.396423039797892,0.578882910828008,0.577030661380971,0.558866268555434,0.406939733640189,0.619294526315412,0.18165384724313,0.47691200510977,0.464102254096015,0.400199441717554,0.649681492845493,0.466990403225362,0.29880418791937,0.711316520668075,0.232747643897451,0.457405422046103,0.262589216927737,0.329284512863754,0.564313078730144,0.510533503235342,0.407365341837343,0.446071160720811,0.524005443404456
"C5",0.445737874237941,0.574328670142162,0.41687771356073,0.369451965233533,0,0.437618428938742,0.298178582442575,0.570155110585478,0.283423226209989,0.230554475204126,0.620572471559088,0.368704587299198,0.487434134429667,0.447254934157088,0.431576596091154,0.483996194599205,0.374293596465489,0.377920430173037,0.467825627857516,0.6027823796088,0.886230088785716,0.32324653241015,0.321229632560589,0.304442694746683,0.624343843046928,0.28451817453965,0.415063641627758,0.421331760981034,0.312480597357034,0.290798386617972,0.520167700099171,0.335277951101073,0.267172505556907,0.524436701686654,0.566512724227427,0.608457004114464,0.879486204128319,0.44408685304625,0.522926114107663,0.462953442970312,0.580334659603419,0.521087648608552,0.621873208450086,0.384168037619774,0.550203847895309,0.416256622946525,0.473079095555507,0.665666310367608,0.475313078437224,0.4355344347563,0.466757551081269,0.354073250324663,0.458070986272189,0.350548152673448,0.389118182357529,0.580869805987346,0.531525046522946,0.470065489847547,0.605558623940919,0.427187049839298
"C6",0.31335445949712,0.493689489094861,0.611044644516638,0.435981050211927,0.437618428938742,0,0.57942854779518,0.559900523782482,0.392636425975046,0.275311449133405,0.268113149832369,0.352717498708046,0.463486437872754,0.532625040201485,0.249567739088519,0.11649767203036,0.418412122704464,0.596626069423967,0.27831171482003,0.269506527070446,0.920663650289505,0.358941597509656,0.459552192529263,0.517958428142256,0.575319343542149,0.25588707879816,0.289645608543058,0.352352554750139,0.475400006223926,0.378664813629489,0.406196939649671,0.616642468649653,0.383841615376691,0.671068236101062,0.513122044523329,0.568996777568023,0.848782269134433,0.216785690655304,0.451721716057755,0.633962712953464,0.328945167632753,0.414333851129122,0.557164610950304,0.389177597396146,0.461962803371258,0.409270430654839,0.182230774779725,0.563877698778723,0.296726725935504,0.287748737882351,0.675223016606337,0.358978585557201,0.30066415930272,0.406531062287049,0.230491271011104,0.249696218504311,0.370810246684122,0.211695821233226,0.583169149197816,0.599679097828673
"C7",0.606726332690357,0.701164798323324,0.692178363905846,0.593339749718081,0.298178582442575,0.57942854779518,0,0.711772007427947,0.40013068741843,0.457328825547921,0.736431636324324,0.56468930617703,0.662682818936594,0.684039794707485,0.59330522430943,0.626605528939434,0.583837964510244,0.623304256127225,0.616945035723713,0.718802526154398,0.925567110689886,0.581372417228134,0.267489387301424,0.5927581592906,0.753345334326865,0.423705811024579,0.592879047010023,0.589999798771308,0.601776823210756,0.540991270314464,0.656474065323847,0.432590083628541,0.503981993645517,0.775143856880362,0.700145085283037,0.749401576311975,0.967189202561296,0.597048615736821,0.609992363060696,0.641090674085538,0.700828117562171,0.663979660180279,0.735629760129771,0.605027623484982,0.683422141758384,0.535168330614542,0.625983356973929,0.773777387182285,0.621533096936597,0.606337118445915,0.224547245139235,0.515825023689743,0.594529524613339,0.58667325647008,0.575941777387981,0.704926244501409,0.661060313134029,0.621713067808978,0.746415395919635,0.659016754996741
"C8",0.444551051309446,0.671674038736365,0.718191111633515,0.327169712971302,0.570155110585478,0.559900523782482,0.711772007427947,0,0.543069177867481,0.486347349686748,0.691756002595733,0.565251824681309,0.447033083629495,0.549949431121168,0.480350601104602,0.554980813714162,0.551435542007668,0.713506734307437,0.546876019267369,0.676459711657326,0.996808274226744,0.543935881376424,0.628843260779125,0.634351973270816,0.581185662295279,0.492302279847897,0.480353633705456,0.531246049584541,0.573751303033905,0.520596907469917,0.500333726456977,0.738161996438313,0.569259139618578,0.756861299654531,0.617900486719104,0.491417369443434,0.88955636810667,0.507760139965615,0.678047566337205,0.747218446886628,0.635436680205864,0.482042340028739,0.708291136043944,0.459131215371924,0.561023543481146,0.58070021921764,0.535045918217404,0.71107756717598,0.590579000529638,0.462571983392586,0.800509464110627,0.381457802216891,0.564887107608785,0.501766038720405,0.507242494850124,0.667301222586974,0.594594673777067,0.530974806998285,0.564331074580251,0.716534024521658
"C9",0.372059461519443,0.523681319136652,0.561871729485542,0.393071607850036,0.283423226209989,0.392636425975046,0.40013068741843,0.543069177867481,0,0.255845299921277,0.577593039287877,0.35868569505834,0.450244499876,0.506733101522743,0.358560676442525,0.422615289745944,0.382741793296419,0.512046978878238,0.403546573152664,0.557245893281019,0.656797567536004,0.380800012704793,0.342082424559482,0.451586150841352,0.574422998984556,0.249821701590595,0.366183153541822,0.372602479546771,0.429928103062708,0.332706491473197,0.463210157097297,0.479366436323747,0.323639647376321,0.649853293578737,0.512305484892613,0.552212992579668,0.834066890237607,0.38047416276833,0.493864828028631,0.541279785654817,0.531076850319492,0.452578084618671,0.56694938663169,0.382876943936348,0.484937458959249,0.242233360023974,0.417178219756003,0.612139027823236,0.428167845709379,0.371499072450497,0.532173755721412,0.312599606560049,0.402396396326535,0.386203857688164,0.352728101867672,0.54069836852979,0.469360476718262,0.407344457699293,0.557319213217732,0.523617923689854
"C10",0.268617844518868,0.446481958367255,0.439036466101009,0.274341741645925,0.230554475204126,0.275311449133405,0.457328825547921,0.486347349686748,0.255845299921277,0,0.495699562272933,0.218396994396507,0.353320443626171,0.390050639082499,0.238190024812502,0.315712417298246,0.281856077795867,0.407256814667629,0.289505875415068,0.480992144926937,0.879610882401652,0.203357850004495,0.361956083895468,0.317176682784106,0.503240493272718,0.140181240729087,0.224872902670048,0.255057053631077,0.27982476752755,0.16685416439635,0.379354819872822,0.450953254667842,0.174573310040803,0.537092666883157,0.466291926514768,0.49628377406297,0.822242742021967,0.263838363073264,0.427944447381569,0.468702884014876,0.453380218627235,0.373014881359612,0.505312920687681,0.253634245623358,0.435566177562718,0.3155806084574,0.29938197376346,0.575875290002891,0.308386775249859,0.261477687009535,0.58901055493649,0.229290237979817,0.295181328043234,0.235187107186444,0.196444834553932,0.444154038041896,0.402966707315303,0.311417211453606,0.514005417126351,0.421301193522285
"C11",0.495078736444999,0.653834729807725,0.747742589316601,0.600076675029399,0.620572471559088,0.268113149832369,0.736431636324324,0.691756002595733,0.577593039287877,0.495699562272933,0,0.546918334663554,0.615073139418015,0.671281461590765,0.447207560961383,0.267342323994985,0.577750014434041,0.73988221109944,0.456319569983299,0.301765781169973,1.0087312521244,0.53414693804818,0.654716969894024,0.673967436558625,0.712557074736071,0.49611323330574,0.480672263953929,0.552945286474951,0.63530197846252,0.56359161313221,0.605970727867706,0.761046176016852,0.581508125318891,0.793768881423387,0.632743778744556,0.690108322885644,0.927372975804954,0.436187767533964,0.678318738891652,0.771979798173737,0.474895637964203,0.571503698070367,0.713848092592111,0.549991791366543,0.591550423954319,0.614889673008844,0.284257663514504,0.67068259889406,0.475704059067776,0.43786172799246,0.820250341424467,0.550411313704897,0.554259800334421,0.587126275984016,0.410144852447543,0.168213079770385,0.526824389053471,0.315459410850532,0.696781639649939,0.740699788718268
"C12",0.35315693411213,0.361943428905081,0.512220228558914,0.394325910841101,0.368704587299198,0.352717498708046,0.56468930617703,0.565251824681309,0.35868569505834,0.218396994396507,0.546918334663554,0,0.465071143451012,0.511525216722596,0.283436843968238,0.369113749678537,0.358656230861484,0.422964446390257,0.261758122330903,0.542714520440843,0.918199828793461,0.29491209873364,0.494301581157657,0.389706021935912,0.588954013054174,0.279533209307178,0.303464446275943,0.228135147480816,0.388526836344666,0.187069819714376,0.471333139028564,0.474703015330613,0.168538647825081,0.640670854849212,0.522093771818986,0.583862855910219,0.858215802349474,0.334590437459933,0.50840457977527,0.416742440719539,0.525320167277005,0.471423650424705,0.43448625759003,0.380987971982154,0.494319212440531,0.346990633819313,0.346209149975707,0.613436859841688,0.252213063163846,0.33682843153268,0.703670543629005,0.358718025532318,0.353765234411601,0.369412505409737,0.227236142884967,0.471814310429044,0.446790583528168,0.389357118078659,0.587478709874173,0.390773755797755
"C13",0.254395263725802,0.591813216820769,0.632119860234518,0.259722323860086,0.487434134429667,0.463486437872754,0.662682818936594,0.447033083629495,0.450244499876,0.353320443626171,0.615073139418015,0.465071143451012,0,0.368127495084612,0.337179529175986,0.455258064217555,0.424703441432407,0.639646243842381,0.431946930850561,0.610326533432289,0.963529807676529,0.395352647466874,0.568917150267139,0.535397321153356,0.338913949207648,0.382996598478367,0.2643532968931,0.413087696368735,0.439748843653356,0.385498400857674,0.400507654163012,0.682809214925338,0.476938347233878,0.661560517595125,0.509494169281271,0.285321687655063,0.858376278190219,0.38395413830716,0.617868288935916,0.6849068075278,0.560417403869192,0.280300948381639,0.645841135017446,0.213155985082615,0.473950935974115,0.52235456986836,0.398981239884322,0.670394467941121,0.497245813455184,0.277812003366818,0.759796755655022,0.245095950788344,0.479737617567608,0.271758909320733,0.357806832096919,0.584876891377853,0.531113780163513,0.404704119133087,0.433044014438861,0.64088130688729
"C14",0.454699047384902,0.645676011255589,0.455764058053165,0.255779804837107,0.447254934157088,0.532625040201485,0.684039794707485,0.549949431121168,0.506733101522743,0.390050639082499,0.671281461590765,0.511525216722596,0.368127495084612,0,0.455194582307774,0.546589685069592,0.419004796712696,0.566427623397469,0.524938930547575,0.668934286198461,0.986216949894065,0.320661235802792,0.609487070305006,0.401650817611115,0.630875509552028,0.46350960016732,0.440831636073038,0.488499313387306,0.26453360873442,0.407815706950498,0.570426009182863,0.682161089209834,0.464829721834953,0.410069842107403,0.584840822642755,0.47983316126792,0.897862611400354,0.511758442349993,0.66190673203791,0.68179514763532,0.647639698415138,0.541084756739226,0.697085084006434,0.225181000898801,0.566611461898011,0.571692802878778,0.493375996796884,0.720995325235635,0.559363549700433,0.404544352061691,0.786728504634223,0.392584930208902,0.559397381253563,0.328859665211022,0.425841858542097,0.637037755423762,0.599599141745524,0.503635053196221,0.501673471166873,0.601552697495693
"C15",0.165111214335744,0.38296217638241,0.608382243572234,0.341726107943767,0.431576596091154,0.249567739088519,0.59330522430943,0.480350601104602,0.358560676442525,0.238190024812502,0.447207560961383,0.283436843968238,0.337179529175986,0.455194582307774,0,0.232056200813537,0.373797160760779,0.594409575319843,0.178447966633,0.443454211036772,0.912663824550714,0.323346908008492,0.47346355771296,0.508507957415185,0.513915222522437,0.216429699168664,0.210123191116234,0.209139039777318,0.444878438313788,0.320950886650078,0.357348495260703,0.621595520676703,0.340706872090782,0.664687942759284,0.459594352131801,0.410557723733412,0.753036721877741,0.206526160581451,0.467374153799532,0.629910925070826,0.425842589530226,0.339636851389865,0.440412379890125,0.276839723380629,0.337457380939558,0.347609157680571,0.200192626626384,0.503604281634953,0.2384319775176,0.133088028611607,0.694958459240215,0.218260498886747,0.272119465447813,0.356973449268605,0.147241519568328,0.390001044204822,0.31749300244279,0.257880227913497,0.421155428386801,0.594433033790373
"C16",0.275429947834589,0.477235256717359,0.641535956453588,0.44822397111334,0.483996194599205,0.11649767203036,0.626605528939434,0.554980813714162,0.422615289745944,0.315712417298246,0.267342323994985,0.369113749678537,0.455258064217555,0.546589685069592,0.232056200813537,0,0.436788262504781,0.631295643996048,0.254605180485584,0.24491255746018,0.93260573649485,0.380260622623572,0.521200259503461,0.5519066131163,0.560435651373258,0.304620565675738,0.281176922714582,0.349886003146793,0.503542063509295,0.402473883505863,0.384337406966974,0.655707851345909,0.42016107695246,0.695191753165725,0.510622165691906,0.546472263579869,0.792013845294214,0.181491772133776,0.529579904558293,0.666920170862358,0.297190021908688,0.380394017497455,0.544194166414571,0.393583528731511,0.404450490780091,0.449812712186142,0.169893472559265,0.502295990173124,0.31224346221312,0.264158691166392,0.723254333664237,0.356656718726676,0.33056083405453,0.420981037445664,0.240000234412063,0.260507542304688,0.312324752719078,0.196945580855326,0.54293722802066,0.631578089245295
"C17",0.378270406240024,0.535973764688221,0.480220443607987,0.354727641089797,0.374293596465489,0.418412122704464,0.583837964510244,0.551435542007668,0.382741793296419,0.281856077795867,0.577750014434041,0.358656230861484,0.424703441432407,0.419004796712696,0.373797160760779,0.436788262504781,0,0.471036065536498,0.409539279528644,0.575602482244944,0.927911049224207,0.315824644953957,0.513367766320939,0.387291930133632,0.582558372732104,0.340323740936286,0.361265843368956,0.379522700448739,0.34192397526982,0.304849845111876,0.492116105698663,0.544864398653198,0.329366400265435,0.559143805462344,0.243197025247424,0.541245300580819,0.863403653246417,0.403101077736933,0.558394657011458,0.539124568776952,0.557266736786312,0.460200744721362,0.588875813414587,0.284136659275024,0.508526839832287,0.450850082027128,0.409370356917924,0.646701249011766,0.434977683082074,0.330362683342582,0.707903957299348,0.355126278796164,0.444153873333882,0.34002831863995,0.332012444638915,0.539081316342633,0.505116468896852,0.307232896532445,0.552644051404323,0.487356135135748
"C18",0.607527928544704,0.692534365773873,0.298072784762931,0.510236439142285,0.377920430173037,0.596626069423967,0.623304256127225,0.713506734307437,0.512046978878238,0.407256814667629,0.73988221109944,0.422964446390257,0.639646243842381,0.566427623397469,0.594409575319843,0.631295643996048,0.471036065536498,0,0.617801671991061,0.72505675231368,1.00551859423844,0.429587262171959,0.629352929445174,0.21738171778581,0.754011312945154,0.523795964724138,0.573586123801221,0.562565164579122,0.340399268493857,0.316360199836968,0.66738005903009,0.352488036293827,0.304193662684218,0.547364361502396,0.693487738131092,0.748760096561289,0.971808527311716,0.601166330712022,0.690376097523133,0.269376329858476,0.709126642141422,0.669790446913477,0.738295800703081,0.534218396238306,0.688122650038705,0.582453533581569,0.62236968521049,0.780136454362113,0.609994374079796,0.595092479026951,0.791319299553853,0.556452936880184,0.619794337086452,0.477807664480688,0.532755566558482,0.704161190377643,0.670277684097521,0.620052782170961,0.737593550402074,0.177501040623811
"C19",0.272778696198201,0.254622944112007,0.629071044587284,0.424958374275872,0.467825627857516,0.27831171482003,0.616945035723713,0.546876019267369,0.403546573152664,0.289505875415068,0.456319569983299,0.261758122330903,0.431946930850561,0.524938930547575,0.178447966633,0.254605180485584,0.409539279528644,0.617801671991061,0,0.46831440979517,0.927025769861893,0.322506930188222,0.509445123338381,0.535559814509209,0.55303156961854,0.272036284270613,0.227879833178296,0.156220759174676,0.483076674723587,0.368247018580634,0.418144550244407,0.644361509103396,0.344443838981682,0.683477941682317,0.492480482284909,0.534233757345129,0.834999614754497,0.238445811132037,0.511188924355839,0.653189164245983,0.458065529798934,0.413147894321575,0.34670795767969,0.370244449410866,0.444890346592307,0.376501490755417,0.226179101498779,0.569687312667964,0.180846635384636,0.249939604569539,0.716057235487017,0.34418499694003,0.294828541347715,0.400847666530179,0.161633208201771,0.379995843659285,0.382665697817542,0.300760722449478,0.553367164763465,0.617758769320473
"C20",0.464520962229137,0.645138302283394,0.732786418868623,0.591982575039002,0.6027823796088,0.269506527070446,0.718802526154398,0.676459711657326,0.557245893281019,0.480992144926937,0.301765781169973,0.542714520440843,0.610326533432289,0.668934286198461,0.443454211036772,0.24491255746018,0.575602482244944,0.72505675231368,0.46831440979517,0,0.997026509597963,0.539148818514298,0.629068329740728,0.658612692270487,0.688622887859886,0.474143265531691,0.474437535105141,0.540773669823445,0.624038717065606,0.552828183160291,0.5005771898907,0.74571215092136,0.569637087811258,0.780641552438602,0.635545428939213,0.69052286627019,0.883988684691822,0.341055794669743,0.636749468779813,0.757767148946119,0.230480876973337,0.519352845224457,0.696041718874709,0.553321116978876,0.547869598643733,0.592034659601226,0.387136460467716,0.589431991434023,0.519452482417006,0.465413928350513,0.80323519507791,0.527212965711451,0.47045044558441,0.563823802403422,0.459212221223372,0.39527515859644,0.432682998872874,0.351631159061793,0.681246442279234,0.726048799308373
"C21",0.922904076761031,0.985455653195309,1.01781986064006,0.931937304639618,0.886230088785716,0.920663650289505,0.925567110689886,0.996808274226744,0.656797567536004,0.879610882401652,1.0087312521244,0.918199828793461,0.963529807676529,0.986216949894065,0.912663824550714,0.93260573649485,0.927911049224207,1.00551859423844,0.927025769861893,0.997026509597963,0,0.92026080449921,0.905155781516104,0.964005697452757,1.02504534293252,0.873131539478423,0.915861257928854,0.917114268917167,0.949338478291611,0.915412856343254,0.960062735188031,0.99374877962651,0.904480397616091,1.05647715518105,0.98672729990775,1.02204857563364,1.18565750682898,0.918373901373996,0.97461514267106,1.02623604663995,0.985479337250967,0.958929514714455,1.00893333209683,0.929032546241457,0.972137995180423,0.691205491883363,0.933812048339074,1.03534811982955,0.93478867594912,0.922113557780756,1.01064765433977,0.901095442378223,0.927249144401921,0.927660851113743,0.90963239261685,0.987295543138433,0.957190139171425,0.930642793879994,1.01967554398889,1.01183746599561
"C22",0.355915670018763,0.440088797891112,0.345969417118892,0.302963467150492,0.32324653241015,0.358941597509656,0.581372417228134,0.543935881376424,0.380800012704793,0.203357850004495,0.53414693804818,0.29491209873364,0.395352647466874,0.320661235802792,0.323346908008492,0.380260622623572,0.315824644953957,0.429587262171959,0.322506930188222,0.539148818514298,0.92026080449921,0,0.496759742458928,0.267911283490708,0.539575148462019,0.301390149431214,0.264331647949984,0.28630256984152,0.186405219984738,0.256391541104002,0.466323639268975,0.561049991750436,0.245877171553343,0.387783120262003,0.51657434398263,0.56790775830591,0.864349772556204,0.346327857384401,0.531670480217084,0.5570031838795,0.524414681365123,0.444856938969505,0.507145709222375,0.266317947943463,0.50483113704041,0.407614582030843,0.350475408691788,0.627088381427999,0.343255836340188,0.33403772840506,0.701794034515378,0.35244608494717,0.389311983509621,0.200968849072137,0.241452973520867,0.472816865980692,0.470790064678127,0.377717688835849,0.573446088746122,0.461988656692588
"C23",0.497210289333075,0.612978653642644,0.661134308362314,0.51513223045837,0.321229632560589,0.459552192529263,0.267489387301424,0.628843260779125,0.342082424559482,0.361956083895468,0.654716969894024,0.494301581157657,0.568917150267139,0.609487070305006,0.47346355771296,0.521200259503461,0.513367766320939,0.629352929445174,0.509445123338381,0.629068329740728,0.905155781516104,0.496759742458928,0,0.563749206975622,0.668016439072078,0.278551768145731,0.484118051202934,0.485461624633941,0.546620497161902,0.484014726195577,0.5430283005977,0.550677423502039,0.443721984534518,0.72289023590682,0.618546748513802,0.661497253557029,0.910977905270097,0.479196753342908,0.408001840703331,0.660358046245435,0.600900462629687,0.56410915804371,0.646535261046312,0.512551894432011,0.593716589070892,0.417112724602884,0.524712149586222,0.695556413436217,0.514263982707578,0.50222864088989,0.280783123157254,0.417752400991542,0.433743348893553,0.502079633163896,0.474372084983094,0.619060459915402,0.559006986674598,0.52156765782347,0.667551524220771,0.649804892741607
"C24",0.519613666536877,0.624536553117282,0.151550352639521,0.393811440770416,0.304442694746683,0.517958428142256,0.5927581592906,0.634351973270816,0.451586150841352,0.317176682784106,0.673967436558625,0.389706021935912,0.535397321153356,0.401650817611115,0.508507957415185,0.5519066131163,0.387291930133632,0.21738171778581,0.535559814509209,0.658612692270487,0.964005697452757,0.267911283490708,0.563749206975622,0,0.673411648609782,0.441296657128217,0.471910061796978,0.479325030847051,0.148336849712515,0.268763206656793,0.588397602104912,0.455267679902818,0.27012181998341,0.349918043439633,0.623342924891687,0.669504390581545,0.922205914384708,0.516374205988384,0.625871047604382,0.416864167188504,0.639325258673129,0.586202637781788,0.675552676274757,0.404615002187007,0.61537698494644,0.521285554665185,0.538463569906722,0.719013884329151,0.535869736392912,0.504248679991568,0.740346885283587,0.468328751165605,0.541692354209239,0.335546173539632,0.44159094598852,0.634017083665719,0.597953717461128,0.537513939013243,0.659976813011697,0.283915186693157
"C25",0.395831778071236,0.682789077003093,0.761091614722705,0.508920388212989,0.624343843046928,0.575319343542149,0.753345334326865,0.581185662295279,0.574422998984556,0.503240493272718,0.712557074736071,0.588954013054174,0.338913949207648,0.630875509552028,0.513915222522437,0.560435651373258,0.582558372732104,0.754011312945154,0.55303156961854,0.688622887859886,1.02504534293252,0.539575148462019,0.668016439072078,0.673411648609782,0,0.516171136039328,0.348944435005501,0.537645518638282,0.610289584481125,0.543475678674711,0.42221327133689,0.777549279187319,0.606760096268271,0.79798488217306,0.637325942062039,0.541117732788221,0.950441245945755,0.46851076860913,0.707598285950763,0.782978040917477,0.613489881488598,0.245238033720581,0.72956919213213,0.492586005151705,0.610179210068826,0.628614005917004,0.556009051766219,0.756184167890742,0.610948840553355,0.504045993520841,0.837364915684415,0.456835961245835,0.573951345697087,0.381837392817694,0.524548734673293,0.689615533871171,0.63196575975533,0.541077183372523,0.679414070003811,0.749363720684664
"C26",0.262267634048403,0.432895061459954,0.553145074029308,0.340499296373118,0.28451817453965,0.25588707879816,0.423705811024579,0.492302279847897,0.249821701590595,0.140181240729087,0.49611323330574,0.279533209307178,0.382996598478367,0.46350960016732,0.216429699168664,0.304620565675738,0.340323740936286,0.523795964724138,0.272036284270613,0.474143265531691,0.873131539478423,0.301390149431214,0.278551768145731,0.441296657128217,0.516171136039328,0,0.244423060065947,0.255377961370487,0.399380496662667,0.287541773100799,0.362673824207936,0.512277298622182,0.270135041716659,0.625249424475355,0.46603033583887,0.498421465679387,0.815905701698716,0.247032379383085,0.347513536527177,0.56091691912169,0.440758683393344,0.373092760732272,0.484471891762578,0.312199098631419,0.426125225095416,0.282322460833754,0.300060657029581,0.56407230513406,0.29395221875772,0.266238009698841,0.524487483635394,0.225660499940735,0.232179544707664,0.321663817221695,0.223892047759916,0.444841509030003,0.382723507534019,0.310577321836364,0.519414551704725,0.53348141328341
"C27",0.182394187978795,0.411602547913554,0.573146973516457,0.322113139244766,0.415063641627758,0.289645608543058,0.592879047010023,0.480353633705456,0.366183153541822,0.224872902670048,0.480672263953929,0.303464446275943,0.2643532968931,0.440831636073038,0.210123191116234,0.281176922714582,0.361265843368956,0.573586123801221,0.227879833178296,0.474437535105141,0.915861257928854,0.264331647949984,0.484118051202934,0.471910061796978,0.348944435005501,0.244423060065947,0,0.237288877381455,0.400364805906871,0.30794547510742,0.319105523636971,0.615063339116454,0.341763396467537,0.621265286373063,0.469754722926182,0.459582122920977,0.835916488206502,0.20453807232228,0.506620790748898,0.621027537607354,0.432096382118773,0.244504680082339,0.482187609909647,0.271063924867126,0.433555487384763,0.392103036204565,0.258340026554264,0.584144981954395,0.29761664444089,0.231952917899641,0.69782666056414,0.259879706210098,0.31307766487275,0.221399459214981,0.190087891102868,0.426158535549681,0.405756013959942,0.28729478195908,0.524907110332623,0.575010975473217
"C28",0.293360968911815,0.211689481572857,0.576501616463269,0.386274363416391,0.421331760981034,0.352352554750139,0.589999798771308,0.531246049584541,0.372602479546771,0.255057053631077,0.552945286474951,0.228135147480816,0.413087696368735,0.488499313387306,0.209139039777318,0.349886003146793,0.379522700448739,0.562565164579122,0.156220759174676,0.540773669823445,0.917114268917167,0.28630256984152,0.485461624633941,0.479325030847051,0.537645518638282,0.255377961370487,0.237288877381455,0,0.430060095591185,0.322089187413107,0.431054898933899,0.600307113472454,0.273826377533175,0.640196557320922,0.493551514814792,0.532943153643131,0.842378564752852,0.289799561642991,0.493042380368135,0.60422810819762,0.508369168843089,0.419249130060024,0.270631499424548,0.348960096502402,0.460792461262308,0.311662050254735,0.325248562129919,0.595470513436138,0.241908414920233,0.286377844323989,0.697705462884124,0.326972405864186,0.296507040085987,0.358347609592515,0.207521512948674,0.486225336710659,0.421169353728456,0.368234019191426,0.553646142548348,0.564653579278115
"C29",0.451497305131059,0.588906295898231,0.225878923100233,0.298458542570305,0.312480597357034,0.475400006223926,0.601776823210756,0.573751303033905,0.429928103062708,0.27982476752755,0.63530197846252,0.388526836344666,0.439748843653356,0.26453360873442,0.444878438313788,0.503542063509295,0.34192397526982,0.340399268493857,0.483076674723587,0.624038717065606,0.949338478291611,0.186405219984738,0.546620497161902,0.148336849712515,0.610289584481125,0.399380496662667,0.400364805906871,0.430060095591185,0,0.268211595069952,0.53863556588014,0.524987525692414,0.291709235526796,0.304502567441694,0.571309640142899,0.589617945565534,0.89199963804343,0.464398042100051,0.601761894138972,0.505368952420986,0.602298383236675,0.52277358730872,0.643102687998476,0.29525008506275,0.564506473481567,0.495200329897359,0.478151614580435,0.686922102934629,0.495316350176727,0.428796560622695,0.734761426577364,0.40242534595499,0.501482099936212,0.249388396619155,0.381154622685559,0.593585996059346,0.558058565130016,0.479487691574107,0.589234012243131,0.393981708947593
"C30",0.336982989891302,0.497133135773248,0.40861959659323,0.291829693341455,0.290798386617972,0.378664813629489,0.540991270314464,0.520596907469917,0.332706491473197,0.16685416439635,0.56359161313221,0.187069819714376,0.385498400857674,0.407815706950498,0.320950886650078,0.402473883505863,0.304849845111876,0.316360199836968,0.368247018580634,0.552828183160291,0.915412856343254,0.256391541104002,0.484014726195577,0.268763206656793,0.543475678674711,0.287541773100799,0.30794547510742,0.322089187413107,0.268211595069952,0,0.448947672310388,0.414907158201902,0.160319118684277,0.555179891742119,0.509048007378812,0.515206134189077,0.840481704028436,0.356650269426897,0.529666135319935,0.34528532925621,0.526818762197503,0.431324955266686,0.557092949006847,0.292541870917605,0.476922560157763,0.395089203641271,0.375772928292223,0.617095640294103,0.37836382496765,0.324111755595511,0.692797993769956,0.297948084705988,0.401774464330188,0.275806865751016,0.272113648664954,0.513240110929121,0.466859934444576,0.389240080616537,0.52802668452317,0.291317428871639
"C31",0.258333776269069,0.582757282686316,0.675855640235749,0.444168445894821,0.520167700099171,0.406196939649671,0.656474065323847,0.500333726456977,0.463210157097297,0.379354819872822,0.605970727867706,0.471333139028564,0.400507654163012,0.570426009182863,0.357348495260703,0.384337406966974,0.492116105698663,0.66738005903009,0.418144550244407,0.5005771898907,0.960062735188031,0.466323639268975,0.5430283005977,0.588397602104912,0.42221327133689,0.362673824207936,0.319105523636971,0.431054898933899,0.53863556588014,0.448947672310388,0,0.689595021727618,0.493078939910316,0.724721809336184,0.554069272858707,0.47508346048068,0.8443778961398,0.2314599794901,0.540446669053057,0.70024513301759,0.340122400730196,0.239104328420478,0.627301916109117,0.423158675221858,0.460155935456361,0.503346212646366,0.439599363819062,0.583429421158168,0.487217936956372,0.386340507101658,0.746917609374612,0.320415366969492,0.318742406181261,0.40357333636689,0.417022953808906,0.583355095977002,0.411495456408508,0.418267075314188,0.575461808163806,0.66669561696544
"C32",0.63482529361811,0.719001661544569,0.568692357456592,0.579067730641874,0.335277951101073,0.616642468649653,0.432590083628541,0.738161996438313,0.479366436323747,0.450953254667842,0.761046176016852,0.474703015330613,0.682809214925338,0.682161089209834,0.621595520676703,0.655707851345909,0.544864398653198,0.352488036293827,0.644361509103396,0.74571215092136,0.99374877962651,0.561049991750436,0.550677423502039,0.455267679902818,0.777549279187319,0.512277298622182,0.615063339116454,0.600307113472454,0.524987525692414,0.414907158201902,0.689595021727618,0,0.394399435854976,0.757480969874195,0.719559741573102,0.775114239580656,0.986988358127568,0.627755779875079,0.692024992371477,0.319263500922143,0.730052273633461,0.693040565264197,0.758708238534587,0.611009612513933,0.709765851845594,0.589835126989428,0.651334567305948,0.7985713910129,0.638399192541888,0.628823936455558,0.640124959984822,0.562396553200945,0.638450999953018,0.578151028509855,0.579574574679549,0.728789978467713,0.690768348696266,0.64825118097583,0.767165569047348,0.399258944577706
"C33",0.395228491377792,0.420206918087509,0.39312475789478,0.367617399097779,0.267172505556907,0.383841615376691,0.503981993645517,0.569259139618578,0.323639647376321,0.174573310040803,0.581508125318891,0.168538647825081,0.476938347233878,0.464829721834953,0.340706872090782,0.42016107695246,0.329366400265435,0.304193662684218,0.344443838981682,0.569637087811258,0.904480397616091,0.245877171553343,0.443721984534518,0.27012181998341,0.606760096268271,0.270135041716659,0.341763396467537,0.273826377533175,0.291709235526796,0.160319118684277,0.493078939910316,0.394399435854976,0,0.547461708628235,0.535991356637101,0.605212251533544,0.87030418834125,0.377943889302608,0.490340850585528,0.370834959225565,0.548995591568073,0.497444057296324,0.470805898927992,0.368804620044216,0.520209162552343,0.33493216501991,0.404046190522197,0.635793195795689,0.345064230733882,0.378942627007027,0.652784808818616,0.362498270528592,0.374293057951088,0.332827100733322,0.282923862597511,0.52241131265194,0.481050861415013,0.425050649605606,0.602089133994835,0.322552749384722
"C34",0.672316921435624,0.755908253953561,0.264895811762571,0.536028271551357,0.524436701686654,0.671068236101062,0.775143856880362,0.756861299654531,0.649853293578737,0.537092666883157,0.793768881423387,0.640670854849212,0.661560517595125,0.410069842107403,0.664687942759284,0.695191753165725,0.559143805462344,0.547364361502396,0.683477941682317,0.780641552438602,1.05647715518105,0.387783120262003,0.72289023590682,0.349918043439633,0.79798488217306,0.625249424475355,0.621265286373063,0.640196557320922,0.304502567441694,0.555179891742119,0.724721809336184,0.757480969874195,0.547461708628235,0,0.751650507541166,0.788966111373702,1.0138071959356,0.668436923364168,0.764787689905048,0.758649392915683,0.765063240471482,0.725313679880568,0.800194209310676,0.521591150223483,0.748399572185362,0.693305394502899,0.682746591537805,0.833249581013668,0.688490090055084,0.653952739382453,0.87303153451203,0.637495273938257,0.691222370526757,0.467655465761648,0.606215652141994,0.758573258611221,0.732862226244628,0.681271042753241,0.775812819922794,0.616299752861611
"C35",0.447465543711839,0.615088635249367,0.703994216037174,0.510328462301058,0.566512724227427,0.513122044523329,0.700145085283037,0.617900486719104,0.512305484892613,0.466291926514768,0.632743778744556,0.522093771818986,0.509494169281271,0.584840822642755,0.459594352131801,0.510622165691906,0.243197025247424,0.693487738131092,0.492480482284909,0.635545428939213,0.98672729990775,0.51657434398263,0.618546748513802,0.623342924891687,0.637325942062039,0.46603033583887,0.469754722926182,0.493551514814792,0.571309640142899,0.509048007378812,0.554069272858707,0.719559741573102,0.535991356637101,0.751650507541166,0,0.569566489586758,0.898475106315924,0.487164723975162,0.646360609632523,0.727967442777778,0.618603202031804,0.51019898272532,0.658555484162286,0.421983909927373,0.56263961103669,0.566541095510533,0.487629410276036,0.699336581171078,0.534356248918418,0.402451617575494,0.790802379459916,0.4563079795335,0.534602769718138,0.515353322625515,0.458195053012401,0.608237537224522,0.573202876723561,0.339167105577533,0.600474103402251,0.697057840973735
"C36",0.336930328032101,0.689368904605218,0.754619078596702,0.391417021535622,0.608457004114464,0.568996777568023,0.749401576311975,0.491417369443434,0.552212992579668,0.49628377406297,0.690108322885644,0.583862855910219,0.285321687655063,0.47983316126792,0.410557723733412,0.546472263579869,0.541245300580819,0.748760096561289,0.534233757345129,0.69052286627019,1.02204857563364,0.56790775830591,0.661497253557029,0.669504390581545,0.541117732788221,0.498421465679387,0.459582122920977,0.532943153643131,0.589617945565534,0.515206134189077,0.47508346048068,0.775114239580656,0.605212251533544,0.788966111373702,0.569566489586758,0,0.852492514775502,0.495262407127831,0.708488301291122,0.780149699959622,0.644117765463845,0.409451133495558,0.730228138271506,0.342587201811713,0.486473528669282,0.626908657741179,0.486040827531019,0.705704161050135,0.605207210009703,0.334855321914748,0.834259896929464,0.302607223480069,0.58154560913728,0.506771227874937,0.480339840383051,0.672437326442707,0.589019729461297,0.489737415394397,0.367450300692441,0.748945084417619
"C37",0.73795235532303,0.928666808278951,0.978139487223776,0.846474965044117,0.879486204128319,0.848782269134433,0.967189202561296,0.88955636810667,0.834066890237607,0.822242742021967,0.927372975804954,0.858215802349474,0.858376278190219,0.897862611400354,0.753036721877741,0.792013845294214,0.863403653246417,0.971808527311716,0.834999614754497,0.883988684691822,1.18565750682898,0.864349772556204,0.910977905270097,0.922205914384708,0.950441245945755,0.815905701698716,0.835916488206502,0.842378564752852,0.89199963804343,0.840481704028436,0.8443778961398,0.986988358127568,0.87030418834125,1.0138071959356,0.898475106315924,0.852492514775502,0,0.79363433953854,0.937459061400891,0.996032498391887,0.870512267011418,0.797317795059667,0.945220327149697,0.822006478932225,0.429874964817721,0.877319922172855,0.835453010304886,0.322892750105378,0.871406704556867,0.758594020283384,1.03298663337379,0.748041520328626,0.840369316659988,0.859686725523428,0.828903425054883,0.918696581087334,0.52202012748556,0.826425400624533,0.53498536883079,0.971806976149503
"C38",0.191688240475808,0.438364577199164,0.612247756514232,0.396423039797892,0.44408685304625,0.216785690655304,0.597048615736821,0.507760139965615,0.38047416276833,0.263838363073264,0.436187767533964,0.334590437459933,0.38395413830716,0.511758442349993,0.206526160581451,0.181491772133776,0.403101077736933,0.601166330712022,0.238445811132037,0.341055794669743,0.918373901373996,0.346327857384401,0.479196753342908,0.516374205988384,0.46851076860913,0.247032379383085,0.20453807232228,0.289799561642991,0.464398042100051,0.356650269426897,0.2314599794901,0.627755779875079,0.377943889302608,0.668436923364168,0.487164723975162,0.495262407127831,0.79363433953854,0,0.479245297195504,0.637939810334335,0.254842960718193,0.280468058905178,0.501171673227938,0.349805290216482,0.387507373097155,0.405861479547131,0.265784586268636,0.500923210313823,0.325180987243366,0.260910599191852,0.696630804744268,0.283992665560328,0.224846238261243,0.345509450395001,0.250730719864976,0.412740120523042,0.294766036635817,0.263128638241721,0.525005357822291,0.600703376924904
"C39",0.525529840194857,0.615622628465977,0.708136363070036,0.578882910828008,0.522926114107663,0.451721716057755,0.609992363060696,0.678047566337205,0.493864828028631,0.427944447381569,0.678318738891652,0.50840457977527,0.617868288935916,0.66190673203791,0.467374153799532,0.529579904558293,0.558394657011458,0.690376097523133,0.511188924355839,0.636749468779813,0.97461514267106,0.531670480217084,0.408001840703331,0.625871047604382,0.707598285950763,0.347513536527177,0.506620790748898,0.493042380368135,0.601761894138972,0.529666135319935,0.540446669053057,0.692024992371477,0.490340850585528,0.764787689905048,0.646360609632523,0.708488301291122,0.937459061400891,0.479245297195504,0,0.718764598312063,0.595291094514426,0.599663186645921,0.643257932254072,0.561833050144338,0.62201593141854,0.420247411394198,0.544928677269052,0.707670228118814,0.476879752990754,0.534501162112757,0.663507810914775,0.49832482778655,0.302622141569481,0.552941155602683,0.490756969055794,0.630476203212427,0.552419110019252,0.549817475827844,0.711783088037604,0.695053070844941
"C40",0.642760642989431,0.724735221146926,0.53718454294183,0.577030661380971,0.462953442970312,0.633962712953464,0.641090674085538,0.747218446886628,0.541279785654817,0.468702884014876,0.771979798173737,0.416742440719539,0.6849068075278,0.68179514763532,0.629910925070826,0.666920170862358,0.539124568776952,0.269376329858476,0.653189164245983,0.757767148946119,1.02623604663995,0.5570031838795,0.660358046245435,0.416864167188504,0.782978040917477,0.56091691912169,0.621027537607354,0.60422810819762,0.505368952420986,0.34528532925621,0.70024513301759,0.319263500922143,0.370834959225565,0.758649392915683,0.727967442777778,0.780149699959622,0.996032498391887,0.637939810334335,0.718764598312063,0,0.742153411767058,0.702021933402975,0.767896445201577,0.611663688643266,0.719737918046952,0.609482065796407,0.659682451926414,0.809835103660505,0.642430466729857,0.635403507675994,0.814613840063102,0.592854546732904,0.654087150771035,0.574448081862852,0.580794777231691,0.738005188706619,0.703312171240533,0.65849420627348,0.772572826900325,0.186202426601806
"C41",0.406278615153105,0.622745599385895,0.716781395473256,0.558866268555434,0.580334659603419,0.328945167632753,0.700828117562171,0.635436680205864,0.531076850319492,0.453380218627235,0.474895637964203,0.525320167277005,0.560417403869192,0.647639698415138,0.425842589530226,0.297190021908688,0.557266736786312,0.709126642141422,0.458065529798934,0.230480876973337,0.985479337250967,0.524414681365123,0.600900462629687,0.639325258673129,0.613489881488598,0.440758683393344,0.432096382118773,0.508369168843089,0.602298383236675,0.526818762197503,0.340122400730196,0.730052273633461,0.548995591568073,0.765063240471482,0.618603202031804,0.644117765463845,0.870512267011418,0.254842960718193,0.595291094514426,0.742153411767058,0,0.434827001015213,0.66980170631986,0.52671586847864,0.518994022664169,0.565411288276532,0.441714159766972,0.564115799432015,0.528648262843566,0.463378244887938,0.786682001048261,0.47489249345468,0.376687609879756,0.514592233133282,0.469669779852515,0.518193681327891,0.389538327100213,0.397927313012186,0.659046929855938,0.709901360778531
"C42",0.191890213701116,0.583802724212561,0.679429445577233,0.406939733640189,0.521087648608552,0.414333851129122,0.663979660180279,0.482042340028739,0.452578084618671,0.373014881359612,0.571503698070367,0.471423650424705,0.280300948381639,0.541084756739226,0.339636851389865,0.380394017497455,0.460200744721362,0.669790446913477,0.413147894321575,0.519352845224457,0.958929514714455,0.444856938969505,0.56410915804371,0.586202637781788,0.245238033720581,0.373092760732272,0.244504680082339,0.419249130060024,0.52277358730872,0.431324955266686,0.239104328420478,0.693040565264197,0.497444057296324,0.725313679880568,0.51019898272532,0.409451133495558,0.797317795059667,0.280468058905178,0.599663186645921,0.702021933402975,0.434827001015213,0,0.631726528362551,0.369295041129946,0.414983808510213,0.516869978080127,0.405191444572083,0.572566869572041,0.488519362320021,0.331652916918369,0.757071812390957,0.284889742423444,0.421435142830143,0.334104616508906,0.389437922578989,0.558119829743011,0.43077594101046,0.371505687039148,0.51707664011488,0.667484037590324
"C43",0.525490195468797,0.130241905555658,0.750273374475198,0.619294526315412,0.621873208450086,0.557164610950304,0.735629760129771,0.708291136043944,0.56694938663169,0.505312920687681,0.713848092592111,0.43448625759003,0.645841135017446,0.697085084006434,0.440412379890125,0.544194166414571,0.588875813414587,0.738295800703081,0.34670795767969,0.696041718874709,1.00893333209683,0.507145709222375,0.646535261046312,0.675552676274757,0.72956919213213,0.484471891762578,0.482187609909647,0.270631499424548,0.643102687998476,0.557092949006847,0.627301916109117,0.758708238534587,0.470805898927992,0.800194209310676,0.658555484162286,0.730228138271506,0.945220327149697,0.501171673227938,0.643257932254072,0.767896445201577,0.66980170631986,0.631726528362551,0,0.591131358326235,0.628657417786601,0.454224949236742,0.537468955921391,0.718677542767682,0.42904013558015,0.5187234898479,0.819741309378551,0.564618242820323,0.469443300682332,0.598169463735777,0.45413409384883,0.65735569006303,0.578893693620671,0.57733197758032,0.731099880467167,0.739116676057197
"C44",0.26045452596257,0.533413102275447,0.496862935656964,0.18165384724313,0.384168037619774,0.389177597396146,0.605027623484982,0.459131215371924,0.382876943936348,0.253634245623358,0.549991791366543,0.380987971982154,0.213155985082615,0.225181000898801,0.276839723380629,0.393583528731511,0.284136659275024,0.534218396238306,0.370244449410866,0.553321116978876,0.929032546241457,0.266317947943463,0.512551894432011,0.404615002187007,0.492586005151705,0.312199098631419,0.271063924867126,0.348960096502402,0.29525008506275,0.292541870917605,0.423158675221858,0.611009612513933,0.368804620044216,0.521591150223483,0.421983909927373,0.342587201811713,0.822006478932225,0.349805290216482,0.561833050144338,0.611663688643266,0.52671586847864,0.369295041129946,0.591131358326235,0,0.433241966220583,0.456624777799817,0.331785203464942,0.617911152893947,0.422245742666795,0.203656082703586,0.714751873743664,0.215174804910972,0.424326507329213,0.240474020651235,0.267076921797713,0.512040670236659,0.470163486423695,0.322848870519691,0.400810301921442,0.549579891262311
"C45",0.314338398215097,0.593335553821491,0.699248806340843,0.47691200510977,0.550203847895309,0.461962803371258,0.683422141758384,0.561023543481146,0.484937458959249,0.435566177562718,0.591550423954319,0.494319212440531,0.473950935974115,0.566611461898011,0.337457380939558,0.404450490780091,0.508526839832287,0.688122650038705,0.444890346592307,0.547869598643733,0.972137995180423,0.50483113704041,0.593716589070892,0.61537698494644,0.610179210068826,0.426125225095416,0.433555487384763,0.460792461262308,0.564506473481567,0.476922560157763,0.460155935456361,0.709765851845594,0.520209162552343,0.748399572185362,0.56263961103669,0.486473528669282,0.429874964817721,0.387507373097155,0.62201593141854,0.719737918046952,0.518994022664169,0.414983808510213,0.628657417786601,0.433241966220583,0,0.534783399620873,0.439045029809679,0.227756357889304,0.503668071689989,0.341831027175029,0.773930351814164,0.340959917296155,0.466726158712327,0.497370059085422,0.429510893449376,0.573518162334401,0.19598523820724,0.432349758430596,0.270660509381645,0.687635302910806
"C46",0.422469938234662,0.438368181572418,0.617326615358582,0.464102254096015,0.416256622946525,0.409270430654839,0.535168330614542,0.58070021921764,0.242233360023974,0.3155806084574,0.614889673008844,0.346990633819313,0.52235456986836,0.571692802878778,0.347609157680571,0.449812712186142,0.450850082027128,0.582453533581569,0.376501490755417,0.592034659601226,0.691205491883363,0.407614582030843,0.417112724602884,0.521285554665185,0.628614005917004,0.282322460833754,0.392103036204565,0.311662050254735,0.495200329897359,0.395089203641271,0.503346212646366,0.589835126989428,0.33493216501991,0.693305394502899,0.566541095510533,0.626908657741179,0.877319922172855,0.405861479547131,0.420247411394198,0.609482065796407,0.565411288276532,0.516869978080127,0.454224949236742,0.456624777799817,0.534783399620873,0,0.445398475358396,0.647083346285784,0.361213464235771,0.420355737345394,0.639156443978873,0.398382971005873,0.344299981474179,0.45182333159002,0.357432065006837,0.557149807348132,0.491065522863742,0.465115825896079,0.628862930157556,0.587129293130081
"C47",0.272260510094909,0.462269717991225,0.631491229152073,0.400199441717554,0.473079095555507,0.182230774779725,0.625983356973929,0.535045918217404,0.417178219756003,0.29938197376346,0.284257663514504,0.346209149975707,0.398981239884322,0.493375996796884,0.200192626626384,0.169893472559265,0.409370356917924,0.62236968521049,0.226179101498779,0.387136460467716,0.933812048339074,0.350475408691788,0.524712149586222,0.538463569906722,0.556009051766219,0.300060657029581,0.258340026554264,0.325248562129919,0.478151614580435,0.375772928292223,0.439599363819062,0.651334567305948,0.404046190522197,0.682746591537805,0.487629410276036,0.486040827531019,0.835453010304886,0.265784586268636,0.544928677269052,0.659682451926414,0.441714159766972,0.405191444572083,0.537468955921391,0.331785203464942,0.439045029809679,0.445398475358396,0,0.580804141958494,0.262587083276251,0.191915591400073,0.72464900542827,0.329744190349168,0.376675305756723,0.402148838387991,0.153849403570434,0.212053673039551,0.402980208642839,0.192090813522108,0.516309993500715,0.623105041443264
"C48",0.49916708408605,0.688719383412591,0.787922284912928,0.649681492845493,0.665666310367608,0.563877698778723,0.773777387182285,0.71107756717598,0.612139027823236,0.575875290002891,0.67068259889406,0.613436859841688,0.670394467941121,0.720995325235635,0.503604281634953,0.502295990173124,0.646701249011766,0.780136454362113,0.569687312667964,0.589431991434023,1.03534811982955,0.627088381427999,0.695556413436217,0.719013884329151,0.756184167890742,0.56407230513406,0.584144981954395,0.595470513436138,0.686922102934629,0.617095640294103,0.583429421158168,0.7985713910129,0.635793195795689,0.833249581013668,0.699336581171078,0.705704161050135,0.322892750105378,0.500923210313823,0.707670228118814,0.809835103660505,0.564115799432015,0.572566869572041,0.718677542767682,0.617911152893947,0.227756357889304,0.647083346285784,0.580804141958494,0,0.616723486893548,0.532221925382019,0.853707762649427,0.536988455538342,0.560175971717097,0.637780039151776,0.578541020160058,0.664563224176448,0.215201217571734,0.569279377190744,0.450918915635847,0.780414335400981
"C49",0.356829269016759,0.344127673408234,0.626669796656057,0.466990403225362,0.475313078437224,0.296726725935504,0.621533096936597,0.590579000529638,0.428167845709379,0.308386775249859,0.475704059067776,0.252213063163846,0.497245813455184,0.559363549700433,0.2384319775176,0.31224346221312,0.434977683082074,0.609994374079796,0.180846635384636,0.519452482417006,0.93478867594912,0.343255836340188,0.514263982707578,0.535869736392912,0.610948840553355,0.29395221875772,0.29761664444089,0.241908414920233,0.495316350176727,0.37836382496765,0.487217936956372,0.638399192541888,0.345064230733882,0.688490090055084,0.534356248918418,0.605207210009703,0.871406704556867,0.325180987243366,0.476879752990754,0.642430466729857,0.528648262843566,0.488519362320021,0.42904013558015,0.422245742666795,0.503668071689989,0.361213464235771,0.262587083276251,0.616723486893548,0,0.321681713618348,0.721717526124649,0.402558111566604,0.327070191027084,0.437796309276182,0.180470170956563,0.35953033666062,0.435833873019832,0.357076621370934,0.613188414532196,0.609947090314072
"C50",0.1611209933261,0.456812522342311,0.603430069636162,0.29880418791937,0.4355344347563,0.287748737882351,0.606337118445915,0.462571983392586,0.371499072450497,0.261477687009535,0.43786172799246,0.33682843153268,0.277812003366818,0.404544352061691,0.133088028611607,0.264158691166392,0.330362683342582,0.595092479026951,0.249939604569539,0.465413928350513,0.922113557780756,0.33403772840506,0.50222864088989,0.504248679991568,0.504045993520841,0.266238009698841,0.231952917899641,0.286377844323989,0.428796560622695,0.324111755595511,0.386340507101658,0.628823936455558,0.378942627007027,0.653952739382453,0.402451617575494,0.334855321914748,0.758594020283384,0.260910599191852,0.534501162112757,0.635403507675994,0.463378244887938,0.331652916918369,0.5187234898479,0.203656082703586,0.341831027175029,0.420355737345394,0.191915591400073,0.532221925382019,0.321681713618348,0,0.709577150785425,0.191266698984074,0.362769571514659,0.345812558284819,0.178175537340372,0.394392100717976,0.367742013120968,0.212540197776207,0.368751176291128,0.596819477504742
"C51",0.707703642001276,0.791995081913527,0.825316236083171,0.711316520668075,0.466757551081269,0.675223016606337,0.224547245139235,0.800509464110627,0.532173755721412,0.58901055493649,0.820250341424467,0.703670543629005,0.759796755655022,0.786728504634223,0.694958459240215,0.723254333664237,0.707903957299348,0.791319299553853,0.716057235487017,0.80323519507791,1.01064765433977,0.701794034515378,0.280783123157254,0.740346885283587,0.837364915684415,0.524487483635394,0.69782666056414,0.697705462884124,0.734761426577364,0.692797993769956,0.746917609374612,0.640124959984822,0.652784808818616,0.87303153451203,0.790802379459916,0.834259896929464,1.03298663337379,0.696630804744268,0.663507810914775,0.814613840063102,0.786682001048261,0.757071812390957,0.819741309378551,0.714751873743664,0.773930351814164,0.639156443978873,0.72464900542827,0.853707762649427,0.721717526124649,0.709577150785425,0,0.62295627967359,0.682694473578846,0.70269199681198,0.688696662966652,0.7927956024176,0.751917565767482,0.720346102853114,0.832701412813695,0.821149847699333
"C52",0.162426809915323,0.516810469026954,0.57710938231357,0.232747643897451,0.354073250324663,0.358978585557201,0.515825023689743,0.381457802216891,0.312599606560049,0.229290237979817,0.550411313704897,0.358718025532318,0.245095950788344,0.392584930208902,0.218260498886747,0.356656718726676,0.355126278796164,0.556452936880184,0.34418499694003,0.527212965711451,0.901095442378223,0.35244608494717,0.417752400991542,0.468328751165605,0.456835961245835,0.225660499940735,0.259879706210098,0.326972405864186,0.40242534595499,0.297948084705988,0.320415366969492,0.562396553200945,0.362498270528592,0.637495273938257,0.4563079795335,0.302607223480069,0.748041520328626,0.283992665560328,0.49832482778655,0.592854546732904,0.47489249345468,0.284889742423444,0.564618242820323,0.215174804910972,0.340959917296155,0.398382971005873,0.329744190349168,0.536988455538342,0.402558111566604,0.191266698984074,0.62295627967359,0,0.355766877982088,0.309686808808016,0.278511648875085,0.515856681845954,0.385189309767734,0.328484387506723,0.36342193572136,0.56121754085793
"C53",0.325281292459037,0.427349119927189,0.633469741649031,0.457405422046103,0.458070986272189,0.30066415930272,0.594529524613339,0.564887107608785,0.402396396326535,0.295181328043234,0.554259800334421,0.353765234411601,0.479737617567608,0.559397381253563,0.272119465447813,0.33056083405453,0.444153873333882,0.619794337086452,0.294828541347715,0.47045044558441,0.927249144401921,0.389311983509621,0.433743348893553,0.541692354209239,0.573951345697087,0.232179544707664,0.31307766487275,0.296507040085987,0.501482099936212,0.401774464330188,0.318742406181261,0.638450999953018,0.374293057951088,0.691222370526757,0.534602769718138,0.58154560913728,0.840369316659988,0.224846238261243,0.302622141569481,0.654087150771035,0.376687609879756,0.421435142830143,0.469443300682332,0.424326507329213,0.466726158712327,0.344299981474179,0.376675305756723,0.560175971717097,0.327070191027084,0.362769571514659,0.682694473578846,0.355766877982088,0,0.417506699302395,0.318327589981013,0.50527308731264,0.361555239500733,0.389272559968361,0.591119826838224,0.621624490227438
"C54",0.315977926132789,0.538597676778906,0.427878453350798,0.262589216927737,0.350548152673448,0.406531062287049,0.58667325647008,0.501766038720405,0.386203857688164,0.235187107186444,0.587126275984016,0.369412505409737,0.271758909320733,0.328859665211022,0.356973449268605,0.420981037445664,0.34002831863995,0.477807664480688,0.400847666530179,0.563823802403422,0.927660851113743,0.200968849072137,0.502079633163896,0.335546173539632,0.381837392817694,0.321663817221695,0.221399459214981,0.358347609592515,0.249388396619155,0.275806865751016,0.40357333636689,0.578151028509855,0.332827100733322,0.467655465761648,0.515353322625515,0.506771227874937,0.859686725523428,0.345509450395001,0.552941155602683,0.574448081862852,0.514592233133282,0.334104616508906,0.598169463735777,0.240474020651235,0.497370059085422,0.45182333159002,0.402148838387991,0.637780039151776,0.437796309276182,0.345812558284819,0.70269199681198,0.309686808808016,0.417506699302395,0,0.314561542125568,0.545972535885533,0.491689687486932,0.399961816474001,0.54699520830946,0.493294791871022
"C55",0.244399925474771,0.375466548975976,0.542649284934117,0.329284512863754,0.389118182357529,0.230491271011104,0.575941777387981,0.507242494850124,0.352728101867672,0.196444834553932,0.410144852447543,0.227236142884967,0.357806832096919,0.425841858542097,0.147241519568328,0.240000234412063,0.332012444638915,0.532755566558482,0.161633208201771,0.459212221223372,0.90963239261685,0.241452973520867,0.474372084983094,0.44159094598852,0.524548734673293,0.223892047759916,0.190087891102868,0.207521512948674,0.381154622685559,0.272113648664954,0.417022953808906,0.579574574679549,0.282923862597511,0.606215652141994,0.458195053012401,0.480339840383051,0.828903425054883,0.250730719864976,0.490756969055794,0.580794777231691,0.469669779852515,0.389437922578989,0.45413409384883,0.267076921797713,0.429510893449376,0.357432065006837,0.153849403570434,0.578541020160058,0.180470170956563,0.178175537340372,0.688696662966652,0.278511648875085,0.318327589981013,0.314561542125568,0,0.324236831330463,0.395986692333113,0.245529119653936,0.504681307780223,0.536993086904029
"C56",0.461349065525982,0.585879273467451,0.711063013597631,0.564313078730144,0.580869805987346,0.249696218504311,0.704926244501409,0.667301222586974,0.54069836852979,0.444154038041896,0.168213079770385,0.471814310429044,0.584876891377853,0.637037755423762,0.390001044204822,0.260507542304688,0.539081316342633,0.704161190377643,0.379995843659285,0.39527515859644,0.987295543138433,0.472816865980692,0.619060459915402,0.634017083665719,0.689615533871171,0.444841509030003,0.426158535549681,0.486225336710659,0.593585996059346,0.513240110929121,0.583355095977002,0.728789978467713,0.52241131265194,0.758573258611221,0.608237537224522,0.672437326442707,0.918696581087334,0.412740120523042,0.630476203212427,0.738005188706619,0.518193681327891,0.558119829743011,0.65735569006303,0.512040670236659,0.573518162334401,0.557149807348132,0.212053673039551,0.664563224176448,0.35953033666062,0.394392100717976,0.7927956024176,0.515856681845954,0.50527308731264,0.545972535885533,0.324236831330463,0,0.509153573317376,0.309381840073154,0.679124913288166,0.705323078062751
"C57",0.329864217978574,0.535713504998622,0.680569021917312,0.510533503235342,0.531525046522946,0.370810246684122,0.661060313134029,0.594594673777067,0.469360476718262,0.402966707315303,0.526824389053471,0.446790583528168,0.531113780163513,0.599599141745524,0.31749300244279,0.312324752719078,0.505116468896852,0.670277684097521,0.382665697817542,0.432682998872874,0.957190139171425,0.470790064678127,0.559006986674598,0.597953717461128,0.63196575975533,0.382723507534019,0.405756013959942,0.421169353728456,0.558058565130016,0.466859934444576,0.411495456408508,0.690768348696266,0.481050861415013,0.732862226244628,0.573202876723561,0.589019729461297,0.52202012748556,0.294766036635817,0.552419110019252,0.703312171240533,0.389538327100213,0.43077594101046,0.578893693620671,0.470163486423695,0.19598523820724,0.491065522863742,0.402980208642839,0.215201217571734,0.435833873019832,0.367742013120968,0.751917565767482,0.385189309767734,0.361555239500733,0.491689687486932,0.395986692333113,0.509153573317376,0,0.399275906462517,0.425260591852547,0.670692454507753
"C58",0.272855751053774,0.513971652196677,0.629262663250619,0.407365341837343,0.470065489847547,0.211695821233226,0.621713067808978,0.530974806998285,0.407344457699293,0.311417211453606,0.315459410850532,0.389357118078659,0.404704119133087,0.503635053196221,0.257880227913497,0.196945580855326,0.307232896532445,0.620052782170961,0.300760722449478,0.351631159061793,0.930642793879994,0.377717688835849,0.52156765782347,0.537513939013243,0.541077183372523,0.310577321836364,0.28729478195908,0.368234019191426,0.479487691574107,0.389240080616537,0.418267075314188,0.64825118097583,0.425050649605606,0.681271042753241,0.339167105577533,0.489737415394397,0.826425400624533,0.263128638241721,0.549817475827844,0.65849420627348,0.397927313012186,0.371505687039148,0.57733197758032,0.322848870519691,0.432349758430596,0.465115825896079,0.192090813522108,0.569279377190744,0.357076621370934,0.212540197776207,0.720346102853114,0.328484387506723,0.389272559968361,0.399961816474001,0.245529119653936,0.309381840073154,0.399275906462517,0,0.5207984917144,0.621905308047654
"C59",0.401024644052398,0.695483967143157,0.740619508725804,0.446071160720811,0.605558623940919,0.583169149197816,0.746415395919635,0.564331074580251,0.557319213217732,0.514005417126351,0.696781639649939,0.587478709874173,0.433044014438861,0.501673471166873,0.421155428386801,0.54293722802066,0.552644051404323,0.737593550402074,0.553367164763465,0.681246442279234,1.01967554398889,0.573446088746122,0.667551524220771,0.659976813011697,0.679414070003811,0.519414551704725,0.524907110332623,0.553646142548348,0.589234012243131,0.52802668452317,0.575461808163806,0.767165569047348,0.602089133994835,0.775812819922794,0.600474103402251,0.367450300692441,0.53498536883079,0.525005357822291,0.711783088037604,0.772572826900325,0.659046929855938,0.51707664011488,0.731099880467167,0.400810301921442,0.270660509381645,0.628862930157556,0.516309993500715,0.450918915635847,0.613188414532196,0.368751176291128,0.832701412813695,0.36342193572136,0.591119826838224,0.54699520830946,0.504681307780223,0.679124913288166,0.425260591852547,0.5207984917144,0,0.738908010975953
"C60",0.606190439625681,0.692518216381593,0.388212048625072,0.524005443404456,0.427187049839298,0.599679097828673,0.659016754996741,0.716534024521658,0.523617923689854,0.421301193522285,0.740699788718268,0.390773755797755,0.64088130688729,0.601552697495693,0.594433033790373,0.631578089245295,0.487356135135748,0.177501040623811,0.617758769320473,0.726048799308373,1.01183746599561,0.461988656692588,0.649804892741607,0.283915186693157,0.749363720684664,0.53348141328341,0.575010975473217,0.564653579278115,0.393981708947593,0.291317428871639,0.66669561696544,0.399258944577706,0.322552749384722,0.616299752861611,0.697057840973735,0.748945084417619,0.971806976149503,0.600703376924904,0.695053070844941,0.186202426601806,0.709901360778531,0.667484037590324,0.739116676057197,0.549579891262311,0.687635302910806,0.587129293130081,0.623105041443264,0.780414335400981,0.609947090314072,0.596819477504742,0.821149847699333,0.56121754085793,0.621624490227438,0.493294791871022,0.536993086904029,0.705323078062751,0.670692454507753,0.621905308047654,0.738908010975953,0'''
=======
>>>>>>> 79d8f9b... big refactoring
