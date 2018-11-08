#  events_placing.py
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

import random
import sys
import os
import numpy as np
from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import itertools


import logging
logger = logging.getLogger("pcoc.events_placing")



class conv_events(object):
    def __init__(self,nodesWithTransitions, nodesWithConvergentModel, nodesWithAncestralModel, all_possibilities_of_transitions, numberOfLeafsWithTransitions, numberOfLeafs):
        self.nodesWithTransitions = nodesWithTransitions
        self.nodesWithTransitions_sim = [x for x in nodesWithTransitions]
        self.nodesWithTransitions_est = [x for x in nodesWithTransitions]

        self.all_possibilities_of_transitions = all_possibilities_of_transitions

        self.nodesWithConvergentModel = nodesWithConvergentModel
        self.nodesWithConvergentModel_sim = [x for x in nodesWithConvergentModel]
        self.nodesWithConvergentModel_est = [x for x in nodesWithConvergentModel]

        self.nodesWithAncestralModel = nodesWithAncestralModel
        self.nodesWithAncestralModel_sim = [x for x in nodesWithAncestralModel]
        self.nodesWithAncestralModel_est = [x for x in nodesWithAncestralModel]

        self.numberOfLeafsWithTransitions = numberOfLeafsWithTransitions
        self.numberOfLeafsWithTransitions_sim = numberOfLeafsWithTransitions
        self.numberOfLeafsWithTransitions_est = numberOfLeafsWithTransitions

        self.numberOfLeafs = numberOfLeafs


    def get_nodes_sim_TorC(self):
        return(self.nodesWithTransitions_sim + self.nodesWithConvergentModel_sim)
    def get_nodes_est_TorC(self):
        return(self.nodesWithTransitions_est + self.nodesWithConvergentModel_est)

    def get_number_nodes_sim(self,cz_nodes):
        logger.info(unlist(cz_nodes.values()))

        return(len(self.nodesWithTransitions_sim+
                   self.nodesWithConvergentModel_sim+
                   self.nodesWithAncestralModel_sim+
                   unlist(cz_nodes.values()))
               )
    def get_number_nodes_est(self):
        return(len(self.nodesWithTransitions_est+
                   self.nodesWithConvergentModel_est+
                   self.nodesWithAncestralModel_est)
               )


class gene_tree(object):
    def __init__(self,ini_tree_fn, manual_mode_nodes):
        self.init_tree_fn = ini_tree_fn
        self.init_tree = None

        self.annotated_tree = None

        self.conv_events = None

        self.manual_mode_nodes = manual_mode_nodes


        self.cz_nodes = None

        self.numberOfLeafs = None

        self.auto_trim_tree = []
        self.reptree = None

    def init_tree_sim(self,flg, bl_new):
        self.flg = flg
        self.bl_new = bl_new
        self.init_tree = init_tree(self.init_tree_fn)
        self.numberOfLeafs = len(self.init_tree.get_tree_root().get_leaves()) + 1

    def init_tree_det(self,n_sites):
        self.flg = 1
        self.bl_new = None
        self.n_sites = n_sites

        if self.auto_trim_tree:
            self.init_tree_fn, self.manual_mode_nodes = trim_tree(self.init_tree_fn, self.manual_mode_nodes, self.auto_trim_tree, self.reptree)



        self.init_tree = init_tree(self.init_tree_fn)
        self.numberOfLeafs = len(self.init_tree.get_tree_root().get_leaves()) + 1

        self.annotated_tree, nodesWithTransitions, observedNumTransitions = manualTransitions_new(self.manual_mode_nodes, self.init_tree)


        nodesWithTransitions     = [k.ND for k in self.annotated_tree.search_nodes(T=1)]
        nodesWithConvergentModel = [k.ND for k in self.annotated_tree.search_nodes(C=1, T=0)]
        nodesWithAncestralModel  = [k.ND for k in self.annotated_tree.search_nodes(T=0,C=0)]

        nodesWithAncestralModel.sort()
        nodesWithConvergentModel.sort()
        nodesWithTransitions.sort()

        self.n_events = len(nodesWithTransitions)

        if nodesWithAncestralModel:
        #### remove root
            root_ND = self.annotated_tree.get_tree_root().ND
            if root_ND in nodesWithTransitions:
                nodesWithTransitions.remove(root_ND)
            else:
                nodesWithAncestralModel.remove(root_ND)

        numberOfNodes = len(self.annotated_tree.get_descendants())
        numberOfLeafs = len(self.annotated_tree.get_leaves())


        numberOfLeafsWithTransitions = 0
        for n_ND in nodesWithConvergentModel:
                node = self.annotated_tree.search_nodes(ND=n_ND)[0]
                if node.is_leaf():
                    numberOfLeafsWithTransitions += 1

        numberOfConvLeafs = 0
        for n_ND in nodesWithConvergentModel + nodesWithTransitions:
                node = self.annotated_tree.search_nodes(ND=n_ND)[0]
                if node.is_leaf():
                    numberOfConvLeafs += 1

        self.numberOfConvLeafs = numberOfConvLeafs

        self.conv_events = conv_events(nodesWithTransitions, nodesWithConvergentModel, nodesWithAncestralModel, nodesWithTransitions, numberOfLeafsWithTransitions, numberOfLeafs)

        self.outgroup_ND = self.annotated_tree.get_tree_root().ND

        logger.debug("tree: %s", self.annotated_tree.get_ascii(attributes=["ND","name","C","T"]))

        if self.conv_events.get_number_nodes_est() != numberOfNodes:
            logger.error("annotated_tree: %s", annotated_tree.get_ascii(attributes=["T","C"]))
            sys.exit(1)



    def init_inter_dir_det(self, repest0, reptree0, repfasta0, repbppconfig, repseq):
        repest=repest0.replace("//","/")
        reptree=reptree0.replace("//","/")
        repfasta=repfasta0.replace("//","/")

        if not os.path.exists(repest) and not os.path.exists(reptree):
            os.mkdir(repest)
            os.mkdir(reptree)

        self.repfasta              = repfasta
        self.repest                = repest
        self.repseq                = repseq
        self.reptree               = reptree
        self.repbppconfig          = repbppconfig

    def init_inter_dir(self,name0, repseq0, repest0, reptree0, repplottreeali0, replikelihoodsummary0, repbppconfig, plot_ali, get_likelihood_summaries):
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
            if plot_ali:
                os.mkdir(repplottreeali)
            if get_likelihood_summaries:
                os.mkdir(replikelihoodsummary)
        else:
            logger.error("%s/%s must not exist", repseq0,name0)
            sys.exit(1)


        self.repseq                = repseq
        self.repest                = repest
        self.reptree               = reptree
        self.repplottreeali        = repplottreeali
        self.replikelihoodsummary  = replikelihoodsummary
        self.repbppconfig          = repbppconfig


    def add_noisy_profils(self,NbCat):
       self.init_tree, self.cz_nodes, r_cz = noise_tree(self.init_tree, NbCat = NbCat)

    def placeNTransitionsInTree(self,n_events, minTrans, maxTrans, maxConvRate):
        self.n_events = n_events
        self.minTrans = minTrans
        self.maxTrans = maxTrans
        self.maxConvRate = maxConvRate

        self.annotated_tree, self.conv_events = placeNTransitionsInTree(n_events, maxTrans, maxConvRate, self.init_tree, manual_mode_nodes = self.manual_mode_nodes, nf=self.init_tree_fn)
        self.numberOfNodes = len(self.annotated_tree.get_descendants())

        self.outgroup_ND = self.annotated_tree.get_tree_root().ND


    def conv_events_is_ok(self):
        return (len(self.conv_events.nodesWithTransitions) >= 1)

    def resolve_conflicts_between_noisy_and_conv_profils(self):
        if self.cz_nodes:
            for (cz, nodes) in self.cz_nodes.items():
                t_node = nodes[0]
                self.conv_events.nodesWithAncestralModel_sim = list(set(self.conv_events.nodesWithAncestralModel_sim) - set(nodes))
                if not t_node in set(self.conv_events.get_nodes_sim_TorC()):
                    nodes = nodes[1:]
                    self.cz_nodes[cz] = [t_node] + list(set(nodes) - set(self.conv_events.get_nodes_sim_TorC()))
                else:
                    self.cz_nodes[cz] = list(set(nodes) - set(self.conv_events.get_nodes_sim_TorC()))
            logger.debug("nodesWithNoisyModel: \n\t-%s", "\n\t-".join(["C%s: %s" %(cz, nodes) for (cz,nodes) in self.cz_nodes.items()]))

        logger.debug("nodesWithAncestralModel: %s", self.conv_events.nodesWithAncestralModel)
        logger.debug("nodesWithAncestralModel_sim: %s", self.conv_events.nodesWithAncestralModel_sim)
        logger.debug("nodesWithConvergentModel: %s", self.conv_events.nodesWithConvergentModel)
        logger.debug("nodesWithConvergentModel_sim: %s", self.conv_events.nodesWithConvergentModel_sim)
        logger.debug("nodesWithTransitions: %s", self.conv_events.nodesWithTransitions)
        logger.debug("numberOfNodes: %s", self.numberOfNodes)
        logger.debug("len: %s tot_sim : %s ", self.conv_events.get_number_nodes_sim(self.cz_nodes), self.numberOfNodes)
        logger.debug("len: %s tot_est : %s ", self.conv_events.get_number_nodes_est(), self.numberOfNodes)

        if self.conv_events.get_number_nodes_sim(self.cz_nodes) != self.numberOfNodes:
            logger.debug("annotated_tree: %s", self.annotated_tree.get_ascii(attributes=["Cz","C"]))
            sys.exit(1)

    def mk_tree_for_simu(self,topo_met=True, plot=False):
        (self.tree_annotated, self.tree_conv_annotated, allbranchlength, convbranchlength) = mk_tree_4_simu(self, topo_met=topo_met, plot=plot)

        self.tree_fn_sim     = self.reptree + "/tree.nhx"
        self.tree_fn_est     = self.reptree + "/tree.nhx"
        self.treeconv_fn_est = self.reptree + "/tree_conv.nhx"

        self.annotated_tree_fn_sim = "%s/annotated_tree.nhx" %(self.reptree)
        self.annotated_tree_fn_est = "%s/annotated_tree.nhx" %(self.reptree)

        return (allbranchlength, convbranchlength)

    def add_noise_in_events_def(self, ev_noise):
        logger.debug("annotated_tree: %s", self.annotated_tree.get_ascii(attributes=["ND","T","C"]))
        logger.debug(self.conv_events.get_number_nodes_est())

        logger.debug("before T: %s", self.conv_events.nodesWithTransitions_est)
        logger.debug("before C: %s", self.conv_events.nodesWithConvergentModel_est)
        logger.debug("before A: %s", self.conv_events.nodesWithAncestralModel_est)
        logger.debug(self.conv_events.all_possibilities_of_transitions)

        new_possible_transitions = list(set(self.conv_events.all_possibilities_of_transitions) - set(self.conv_events.nodesWithTransitions_est))
        new_transitions = []
        discarded_transistions = []

        if ev_noise == "+1":
            new_transitions = random.sample(new_possible_transitions, 1)
        elif ev_noise == "-1":
            discarded_transistions = random.sample(self.conv_events.nodesWithTransitions_est, 1)
        elif ev_noise == "=1":
            discarded_transistions = random.sample(self.conv_events.nodesWithTransitions_est, 1)
            new_transitions = random.sample(new_possible_transitions, 1)

        new_C = []
        discarded_C = []

        for ND_nt in new_transitions:
            nt = self.annotated_tree.search_nodes(ND=ND_nt)
            logger.debug("annotated_tree: %s", nt[0].get_ascii(attributes=["ND","T","C"]))
            for nodes in nt[0].traverse():
                if nodes.ND != ND_nt:
                    new_C.append(nodes.ND)

        for ND_dt in discarded_transistions:
                nt = self.annotated_tree.search_nodes(ND=ND_dt)
                logger.debug("annotated_tree: %s", nt[0].get_ascii(attributes=["ND","T","C"]))
                for nodes in nt[0].traverse():
                    if nodes.ND != ND_dt:
                        discarded_C.append(nodes.ND)

        self.conv_events.nodesWithAncestralModel_est = list(set( self.conv_events.nodesWithAncestralModel_est) - set(new_C) - set(new_transitions)) + discarded_C + discarded_transistions
        self.conv_events.nodesWithTransitions_est = list(set(self.conv_events.nodesWithTransitions_est)  - set(discarded_transistions)) + new_transitions
        self.conv_events.nodesWithConvergentModel_est = list(set(self.conv_events.nodesWithConvergentModel_est) - set(discarded_C)) + new_C

        logger.debug("after T: %s", self.conv_events.nodesWithTransitions_est)
        logger.debug("after C: %s", self.conv_events.nodesWithConvergentModel_est)#self.nodesWithConvergentModel_est
        logger.debug("after A: %s", self.conv_events.nodesWithAncestralModel_est)#self.nodesWithAncestralModel_est
        logger.debug(self.conv_events.get_number_nodes_est())#self.nodesWithTransitions_est


    def add_noise_in_bl_tree_est(self,topo_met):
        (bl_before, bl_err, bl_after) = noise_bl(self.annotated_tree, self.reptree, vnodes=self.conv_events.get_nodes_est_TorC(), topo_met=topo_met)
        median_err = np.median(bl_err)
        sd_err = np.std(bl_err)
        mean_err = np.mean(bl_err)
        logger.info("Median bl err: %s", median_err)
        logger.info("Sd bl err: %s", median_err)
        logger.info("Mean bl err: %s", mean_err)

        self.tree_fn_est     = self.reptree + "/noisy_tree.nhx"
        self.treeconv_fn_est = self.reptree + "/noisy_tree_conv.nhx"

        return (mean_err, sd_err, median_err)

    def add_noise_in_root_tree_est(self, root_noise, topo_met):
        vnodes=self.conv_events.get_nodes_est_TorC()

        noisy_tree = self.annotated_tree.copy(method="deepcopy")
        logger.debug("sim root %s",self.outgroup_ND)


        root_children = noisy_tree.get_tree_root().children

        if root_noise == "ll":
            new_outgroup = root_children[0].children[0]
            self.outgroup_ND = root_children[1].ND
        elif root_noise == "lr":
            new_outgroup = root_children[0].children[1]
            self.outgroup_ND = root_children[1].ND
        elif root_noise == "rr":
            new_outgroup = root_children[1].children[1]
            self.outgroup_ND = root_children[0].ND
        elif root_noise == "rl":
            new_outgroup = root_children[1].children[0]
            self.outgroup_ND = root_children[0].ND

        noisy_tree.set_outgroup(new_outgroup)



        logger.debug("det root %s",self.outgroup_ND)

        new_nodesWithTransitions_est = []
        new_nodesWithConvergentModel_est = []
        new_nodesWithAncestralModel_est = []

        if not topo_met:
            tconv = None
        else:
            noisy_tconv = build_conv_topo(noisy_tree, vnodes)
            noisy_tconv.write(format=1, features=["ND"],outfile="%s/noisy_root_tree_conv.nhx"%(self.reptree), format_root_node=True)

        noisy_tree.write(format=1, features=["ND"],outfile="%s/noisy_root_tree.nhx"%(self.reptree), format_root_node=True)

        if self.cz_nodes:
            noisy_tree.write(format=1, features=["ND","T","C","Cz"],outfile="%s/annotated_tree_est.nhx"%(self.reptree),format_root_node=True)
        else:
            noisy_tree.write(format=1, features=["ND","T","C"],outfile="%s/annotated_tree_est.nhx"%(self.reptree),format_root_node=True)


        self.tree_fn_est     = self.reptree + "/noisy_root_tree.nhx"
        self.treeconv_fn_est = self.reptree + "/noisy_root_tree_conv.nhx"
        self.annotated_tree_fn_est = "%s/annotated_tree_est.nhx" %(self.reptree)

        return

# Basic tree style
tree_style = TreeStyle()
tree_style.show_leaf_name = True
tree_style.show_branch_length = True
tree_style.min_leaf_separation  = 4
tree_style.branch_vertical_margin   = 2

nstyle_T = NodeStyle()
nstyle_T["fgcolor"] = "orange"
#nstyle_T["shape"] = "square"
nstyle_T["size"] = 5
nstyle_T["vt_line_color"] = "orange"
nstyle_T["hz_line_color"] = "orange"
nstyle_T["vt_line_width"] = 2
nstyle_T["hz_line_width"] = 2


nstyle_C = NodeStyle()
nstyle_C["fgcolor"] = "orange"
nstyle_C["size"] = 5
nstyle_C["vt_line_color"] = "orange"
nstyle_C["hz_line_color"] = "orange"
nstyle_C["vt_line_width"] = 2
nstyle_C["hz_line_width"] = 2

nstyle = NodeStyle()
nstyle["fgcolor"] = "blue"
nstyle["size"] = 5
nstyle["hz_line_width"] = 2
nstyle["vt_line_width"] = 2
nstyle["vt_line_color"] = "#0052FF"
nstyle["hz_line_color"] = "#0052FF"



def add_t(node):
    nd = TextFace("-")
    nd.fsize = 4
    nd.background.color = "black"
    nd.margin_right = 0
    nd.margin_top = 0
    nd.margin_left = 0
    nd.margin_bottom = 0
    nd.border.width = 1
    nd2 = TextFace(" ")
    nd2.fsize = 4

    node.add_face(nd, column=0, position = "float")
    node.add_face(nd2, column=1, position = "float")

def build_conv_topo(annotated_tree, vnodes):

      tconv = annotated_tree.copy(method="deepcopy")
      for n in tconv.iter_leaves():
        n.add_features(L=1)
      for n in tconv.traverse():
        n.add_features(COPY=0)
      # get the most recent ancestral node of all the convergent clades
      l_convergent_clades = tconv.search_nodes(T=True)
      common_anc_conv=tconv.get_common_ancestor(l_convergent_clades)

      # duplicate it at its same location (branch lenght = 0). we get
      # a duplicated subtree with subtrees A and B (A == B)

      dist_dup = common_anc_conv.dist
      if not common_anc_conv.is_root():
        dup_point = common_anc_conv.add_sister(name="dup_point",dist=0.000001)
        dup_point_root = False
      else:
        dup_point = Tree()
        dup_point_root = True
        dup_point.dist=0.000001

      dup_point.add_features(ND=0,T=False, C=False, Cz=False)

      common_anc_conv.detach()
      common_anc_conv_copy = common_anc_conv.copy(method="deepcopy")

      # tag duplicated nodes:

      for n in common_anc_conv_copy.traverse():
        n.COPY=1
        if n.ND not in vnodes and not n.is_root():
            n.dist=0.000001

      # pruned A from all branches not leading to any convergent clade
      l_leaves_to_keep_A = common_anc_conv.search_nodes(COPY=0, C=False, L=1)
      #logger.debug("A: %s",l_leaves_to_keep_A)
      common_anc_conv.prune(l_leaves_to_keep_A, preserve_branch_length=True)

      # pruned B from all branches not leading to any non-convergent clade
      l_leaves_to_keep_B = common_anc_conv_copy.search_nodes(COPY=1, C=True, L=1)
      #logger.debug("B : %s", l_leaves_to_keep_B)
      common_anc_conv_copy.prune(l_leaves_to_keep_B, preserve_branch_length=True)


      dup_point.add_child(common_anc_conv_copy)
      dup_point.add_child(common_anc_conv)

      tconv = dup_point.get_tree_root()

      nodeId = 0
      for node in tconv.traverse("postorder"):
          node.ND = nodeId
          nodeId += 1

      return tconv

def mk_tree_4_simu(g_tree, topo_met = False, plot = False):

  conv_combi_events = g_tree.conv_events
  cz_nodes = g_tree.cz_nodes
  annotated_tree = g_tree.annotated_tree
  reptree = g_tree.reptree
  flg = g_tree.flg
  bl_new = g_tree.bl_new

  vnodes = conv_combi_events.get_nodes_sim_TorC()

  ## multiply branch lengths
  if flg != 1:
    logger.warning("multiply branch lengths by: %s", flg)
    for n in annotated_tree.traverse("postorder"):
        n.dist = n.dist * flg
  ## multiply branch lengths
  if bl_new >= 0:
    logger.warning("Replace all branch lengths by: %s", bl_new )
    for n in annotated_tree.traverse("postorder"):
        n.dist = bl_new

  allbranchlength = 0
  convbranchlength = 0
  lnodes=map(str,vnodes)


  ## numerote nodes
  lab_esp={}
  for n in annotated_tree.traverse("postorder"):
    if n.is_leaf():
      lab_esp[n.ND] = n.name
    allbranchlength += n.dist
    if n.ND in vnodes:
        convbranchlength += n.dist

  if cz_nodes:
    annotated_tree.write(format=1, features=["ND","T","C","Cz"],outfile="%s/annotated_tree.nhx"%(reptree),format_root_node=True)
  else:
    annotated_tree.write(format=1, features=["ND","T","C"],outfile="%s/annotated_tree.nhx"%(reptree),format_root_node=True)
  annotated_tree.write(format=1, features=["ND"],outfile="%s/tree.nhx"%(reptree),format_root_node=True)


  if plot:
    cz_nodes_s = {}
    if cz_nodes:
        cols = ["#008000","#800080","#007D80","#9CA1A2","#A52A2A","#ED8585","#FF8EAD","#8EB1FF","#FFE4A1","#ADA1FF"]
        col_i = 0

        for Cz in cz_nodes.keys():
            cz_nodes_s[Cz] = NodeStyle()
            cz_nodes_s[Cz]["fgcolor"] = cols[col_i]
            cz_nodes_s[Cz]["size"] = 5
            cz_nodes_s[Cz]["hz_line_width"] = 2
            cz_nodes_s[Cz]["vt_line_width"] = 2
            cz_nodes_s[Cz]["vt_line_color"] = cols[col_i]
            cz_nodes_s[Cz]["hz_line_color"] = cols[col_i]
            col_i +=1
    for n in annotated_tree.traverse():
        if n.T:
            n.set_style(nstyle_T)
            add_t(n)

        elif n.C:
            n.set_style(nstyle_C)
        elif cz_nodes_s and n.Cz:
            n.set_style(cz_nodes_s[n.Cz])
            if n.ND == cz_nodes[n.Cz][0]:
                add_t(n)
        else:
            n.set_style(nstyle)

    annotated_tree.render("%s/tree.pdf"%(reptree), tree_style=tree_style)

  if not topo_met:
      tconv = None
  else:
      tconv = build_conv_topo(annotated_tree, conv_combi_events.get_nodes_sim_TorC())

      if cz_nodes:
          tconv.write(format=1, features=["ND","T","C","Cz"],outfile="%s/annotated_tree_conv.nhx"%(reptree),format_root_node=True)
      else:
          tconv.write(format=1, features=["ND","T","C"],outfile="%s/annotated_tree_conv.nhx"%(reptree),format_root_node=True)


      tconv.write(format=1, features=["ND"],outfile="%s/tree_conv.nhx"%(reptree),format_root_node=True)


      if plot:
          for n in tconv.traverse():
            #print n.get_ascii(attributes=["ND","name","C","T"])
            if n.T:
                n.set_style(nstyle_T)
                add_t(n)
            elif n.C:
                n.set_style(nstyle_C)
            elif cz_nodes_s and n.Cz:
                n.set_style(cz_nodes_s[n.Cz])
                if n.ND == cz_nodes[int(n.Cz)][0]:
                    add_t(n)
            else:
                n.set_style(nstyle)

          tconv.render("%s/tree_conv.pdf"%(reptree), tree_style=tree_style)

  return (annotated_tree, tconv, allbranchlength, convbranchlength)



def unlist(list2d):
    return [item for sublist in list2d for item in sublist]


def init_tree(nf):
  t = Tree(nf)

  #Alternatively we could read a tree from a file into a string "line", and then use:
  # t =  Tree( line )

  nodeId = 0
  for n in t.traverse("postorder"):
    n.add_features(ND=nodeId)
    nodeId = nodeId + 1

  return t


## recursively set all the values of the nodes under this one
## (included) to this value
def setsubtree_new(node, value):
  node.C=value
  for son in node.traverse("levelorder"):
    son.C=value

def setsubtree_Cz(node, value):
  node.C=value
  for son in node.traverse("levelorder"):
    son.C=value


# We could use some dynamic programming to be able to generate paths that yield n transitions exactly.
# INstead we randomly generate transitions on the tree until we get the desired number.
# We have two states: ancestral (0) and convergent (1).
# We count the numbers of transitions

def randomTransitions_new(numTransitions, tree):
  numberOfNodes = len(tree.get_tree_root().get_descendants()) + 1
  rate = float(numTransitions)/float(numberOfNodes)
  totalNumberOfTransitions = 0
  nodesWithTransitions = list()

  for n in tree.traverse():
      n.add_features(T=False)
      n.add_features(C=False)
  for node in tree.traverse("levelorder"):
    if node.is_root() :
        pass
    elif (node.up.C == True):
      node.C=True
    else :
      sisterHasAlreadyTransitioned=node.get_sisters()[0].T #Here we assume binary trees!

      #randomly draw whether we do a transition or not
      transitionBool = random.uniform(0,1)<rate
      if (transitionBool and not sisterHasAlreadyTransitioned):
        setsubtree_new(node,True)
        nodesWithTransitions.append(node.ND)
        node.T=True
        node.C=True
        totalNumberOfTransitions = totalNumberOfTransitions + 1

  return tree, nodesWithTransitions, totalNumberOfTransitions

def manualTransitions_new(manual_mode_nodes, tree):
  numberOfNodes = len(tree.get_tree_root().get_descendants()) + 1


  totalNumberOfTransitions = 0
  nodesWithTransitions = list()

  for n in tree.traverse():
      n.add_features(T=False)
      n.add_features(C=False)
  for node in tree.traverse("levelorder"):
    if node.is_root() :
        pass
    elif node.ND in manual_mode_nodes["T"]:
      totalNumberOfTransitions = totalNumberOfTransitions + 1
      nodesWithTransitions.append(node.ND)
      node.T=True
      node.C=True
    elif node.ND in manual_mode_nodes["C"]:
      node.C=True

  return tree, nodesWithTransitions, totalNumberOfTransitions


def randomTransitions_cz(Cz, tree):

  numberOfNodes = len(tree.get_tree_root().get_descendants()) + 1
  rate = float(len(Cz))/float(numberOfNodes)
  totalNumberOfTransitions = 0
  nodesWithNoisyModel = {}
  nodesWithTransition = {}

  for n in tree.traverse():
      n.add_features(Cz=False)
  for node in tree.traverse("levelorder"):
    if node.is_root() :
        pass
    elif (node.up.Cz):
      node.Cz=node.up.Cz
    else :
      sisterHasAlreadyTransitioned=node.get_sisters()[0].Cz #Here we assume binary trees!

      #randomly draw whether we do a transition or not
      transitionBool = random.uniform(0,1)<rate
      if (transitionBool and (not sisterHasAlreadyTransitioned) and totalNumberOfTransitions < len(Cz)):
        setsubtree_Cz(node,Cz[totalNumberOfTransitions])
        node.Cz=Cz[totalNumberOfTransitions]
        totalNumberOfTransitions = totalNumberOfTransitions + 1
        nodesWithTransition[node.Cz] = node.ND

  for Cz, t_node in nodesWithTransition.items():
      nodesWithNoisyModel[Cz] = [t_node]

  for node in tree.traverse():
      if node.Cz:
          if not nodesWithTransition[node.Cz] == node.ND:
              nodesWithNoisyModel[node.Cz].append(node.ND)

  return tree, nodesWithNoisyModel, totalNumberOfTransitions

def noise_tree(tree_ini, NbCat = 10):

    numberOfLeafs = len(tree_ini.get_tree_root().get_leaves()) + 1
    logger.debug( "(numberOfLeafs start: %s)", numberOfLeafs)

    numberOfNodes = len(tree_ini.get_tree_root().get_descendants()) + 1
    len_cz_nodes = numberOfNodes
    r_cz = float(len_cz_nodes)/float(numberOfNodes)
    ntry = 0

    while  (r_cz < 0.1 or r_cz > 0.30) and ntry < 500:
        N = random.randrange(0,10)
        Cz = random.sample(range(1,NbCat+1), N)
        tree = tree_ini.copy(method="deepcopy")
        tree, cz_nodes,  ncz = randomTransitions_cz(Cz, tree)
        r_cz = float(len(unlist(cz_nodes.values())))/float(numberOfNodes)
        ntry +=1
        logger.debug( "(ntry: %s)", ntry)
        logger.debug( "(cz_nodes: %s)", cz_nodes)
        logger.debug( "(ncz: %s)", ncz)
        logger.debug( "(Cz: %s)", Cz)
        logger.debug( "(len(unlist(cz_nodes.values())): %s)", len(unlist(cz_nodes.values())))
        logger.debug( "(cz_nodes.values(): %s)", cz_nodes.values())
        logger.debug( "(r_cz: %s)", r_cz)
        logger.debug( "(r_cz > 0.1: %s)", r_cz > 0.1)
        logger.debug( "(not r_cz == 0: %s)", not r_cz == 0)


    if ntry == 500:
        tree = tree_ini.copy(method="deepcopy")
        cz_nodes = {}
        r_cz = 0

    numberOfLeafs = len(tree.get_tree_root().get_leaves()) + 1
    logger.debug( "(numberOfLeafs end: %s)", numberOfLeafs)

    return tree, cz_nodes, r_cz


def noise_bl(tree, reptree, vnodes=None, topo_met=False, cz_nodes= {} ):
    noisy_tree =tree.copy(method="deepcopy")

    bl_before = []
    bl_err = []
    bl_after = []

    for n in noisy_tree.traverse("postorder"):
        bl_before.append(n.dist)
        random_err = np.random.gamma(10,0.1)
        bl_err.append(random_err)
        n.dist = n.dist * random_err
        bl_after.append(n.dist)

    if not topo_met:
        tconv = None
    else:
        noisy_tconv = build_conv_topo(noisy_tree, vnodes)
        noisy_tconv.write(format=1, features=["ND"],outfile="%s/noisy_tree_conv.nhx"%(reptree), format_root_node=True)


    noisy_tree.write(format=1, features=["ND"],outfile="%s/noisy_tree.nhx"%(reptree), format_root_node=True)
    return (bl_before, bl_err, bl_after)


def placeNTransitionsInTree(numTransitions, maxTransitions, maxConvRate, tree_ini, manual_mode_nodes = {}, nf=""):

  logger.debug( "Wanted Transitions: %s (Max: %s)", numTransitions, maxTransitions)

  observedNumTransitions = -1
  nodesWithTransitions_filtered = []

  numTries = 0
  numTries2 = 0
  numTries3 = 1000
  maxConvRateObs = 1000

  while (maxConvRate < maxConvRateObs or observedNumTransitions != maxTransitions) and numTries2 < 200:
    numTries2 += 1
    numTries = 0
    observedNumTransitions = -1

    #logger.debug("numTries2: %s", numTries2)
    while observedNumTransitions != maxTransitions and numTries < 500:
        numTries = numTries + 1
        tree_all_transitions = tree_ini.copy(method="deepcopy")

        if manual_mode_nodes:
            tree_all_transitions, nodesWithTransitions, observedNumTransitions = manualTransitions_new(manual_mode_nodes, tree_all_transitions)
            #logger.debug ("Observed Number of Transitions (MANUEL): "+ str(observedNumTransitions ) + " compared to "+ str(maxTransitions) + " wanted in:" + nf )
        else:
            tree_all_transitions, nodesWithTransitions, observedNumTransitions = randomTransitions_new(maxTransitions, tree_all_transitions)
            #logger.debug ("Observed Number of Transitions: "+ str(observedNumTransitions ) + " compared to "+ str(maxTransitions) + " wanted in:" + nf )

    #logger.debug("numTries: %s", numTries)
    if numTries < 500:
        numTries3 = 0
        nodesWithTransitions_allCombinations = list(itertools.combinations(nodesWithTransitions,maxTransitions-numTransitions))
        nodesWithTransitions_list = random.sample(nodesWithTransitions_allCombinations, len(nodesWithTransitions_allCombinations))

        while maxConvRate < maxConvRateObs and nodesWithTransitions_list and numTries3 < 1000:
            #logger.debug("%s < %s and %s and %s < 1000", maxConvRate , maxConvRateObs ,nodesWithTransitions_list ,numTries3)
            numTries3 = numTries3 + 1
            #logger.debug("numTries3: %s", numTries3)
            #delete maxTransitions-numTransitions random convergence events
            tree_final = tree_all_transitions.copy(method="deepcopy")
            lr=nodesWithTransitions_list.pop()
            for n_ND in lr:
                #print n_ND
                node = tree_final.search_nodes(ND=n_ND)[0]
                #print node.get_ascii(attributes=["ND","name","C","T"])
                setsubtree_new(node,False)
                node.T=False
                #print node.get_ascii(attributes=["ND","name","C","T"])

            nodesWithTransitions_filtered=[k for k in nodesWithTransitions if not k in lr]

            numberOfLeafsWithTransitions = 0
            for n_ND in nodesWithTransitions_filtered:
                node = tree_final.search_nodes(ND=n_ND)[0]
                #logger.debug("nodesWithTransitions: %s", n.get_ascii())
                numberOfLeafsWithTransitions += len(node.get_leaves())

            numberOfLeafs = len(tree_final.get_tree_root().get_leaves()) + 1
            numberOfLeafsWithoutTransitions = numberOfLeafs - numberOfLeafsWithTransitions

            maxConvRateObs = float(numberOfLeafsWithTransitions)/ float(numberOfLeafs)

            #logger.debug("OK numberOfLeafs: %s", numberOfLeafs)
            #logger.debug("numberOfLeafsWithTransitions: %s", numberOfLeafsWithTransitions)
            #logger.debug("numberOfLeafsWithoutTransitions: %s", numberOfLeafsWithoutTransitions)
            #logger.debug("rate: %s", maxConvRateObs)


  #logger.debug("numTries2: %s", numTries2)
  if numTries3 >=  1000 and maxConvRateObs > maxConvRate:
    logger.warning("It seems like it is too difficult to place "+ str(maxTransitions) + " events in this tree:" + nf)
    numberOfLeafsWithTransitions = 0
    numberOfLeafs = 0
    all_possibilities_of_transitions = []
    nodesWithTransitions = []
    nodesWithConvergentModel = []
    nodesWithAncestralModel = []
    tree_final = Tree()
  else:
    all_possibilities_of_transitions = nodesWithTransitions
    nodesWithTransitions = [k.ND for k in tree_final.search_nodes(T=1)]
    nodesWithConvergentModel = [k.ND for k in tree_final.search_nodes(C=1, T=0)]
    nodesWithAncestralModel =  [k.ND for k in tree_final.search_nodes(T=0,C=0)]
    for n in nodesWithTransitions_filtered:
        #logger.debug("nodesWithTransitions: %s", tree_final.search_nodes(ND=n)[0].get_ascii(attributes=["ND","name"]))
        pass

  #logger.debug("tree: %s", tree_final.get_ascii(attributes=["ND","name","C","T"]))


  nodesWithAncestralModel.sort()
  nodesWithConvergentModel.sort()
  nodesWithTransitions.sort()

  if nodesWithAncestralModel:
  #### remove root
    root_ND = tree_final.get_tree_root().ND
    if root_ND in nodesWithTransitions:
        nodesWithTransitions.remove(root_ND)
    else:
        nodesWithAncestralModel.remove(root_ND)

  conv_combi_events = conv_events(nodesWithTransitions, nodesWithConvergentModel, nodesWithAncestralModel, all_possibilities_of_transitions, numberOfLeafsWithTransitions, numberOfLeafs)

  return tree_final, conv_combi_events

def trim_tree(tree_fn, manual_mode_nodes, sp_present, reptree):
    trim_t = init_tree(tree_fn)
    trim_t.prune(sp_present)
    trim_nodeId =0

    for n in trim_t.traverse("postorder"):
        n.add_features(TND=trim_nodeId)
        trim_nodeId = trim_nodeId + 1

    trim_manual_mode_nodes = {"T":[],"C":[]}
    for k in ["C","T"]:
        for ND in manual_mode_nodes[k]:
            TND_l = trim_t.search_nodes(ND=ND)
            if TND_l:
                TND = TND_l[0].TND
                trim_manual_mode_nodes[k].append(TND)

    trim_t.write(format=1,outfile="%s/trim_tree.nw"%(reptree), format_root_node=True)
    return("%s/trim_tree.nw"%(reptree),trim_manual_mode_nodes)

