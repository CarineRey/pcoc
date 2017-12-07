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
import numpy as np
from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import itertools


import logging
logger = logging.getLogger("simul_CAT_rand.placing")


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

def mk_tree_4_simu_new(annotated_tree, reptree, nodesWithConvergentModel,nodesWithTransitions, flg = 1, bl_new = -1, topo_met = False, plot = False, cz_nodes = {}):

  vnodes = nodesWithConvergentModel+nodesWithTransitions

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
      tconv = build_conv_topo(annotated_tree, nodesWithConvergentModel+nodesWithTransitions)
      
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
    logger.debug( "(numberOfLeafs debut: %s)", numberOfLeafs)

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
    logger.debug( "(numberOfLeafs fin: %s)", numberOfLeafs)

    return tree, cz_nodes, r_cz


def noise_bl(tree, reptree, vnodes=None, topo_met=False):
    noisy_tree =tree.copy(method="deepcopy")
    
    for n in noisy_tree.traverse("postorder"):
        random_err = np.random.exponential(1)
        n.dist = n.dist * random_err
    
    if not topo_met:
        tconv = None
    else:
        noisy_tconv = build_conv_topo(noisy_tree, vnodes)
        noisy_tconv.write(format=1, features=["ND"],outfile="%s/noisy_tree_conv.nhx"%(reptree), format_root_node=True)
    
    noisy_tree.write(format=1, features=["ND"],outfile="%s/noisy_tree.nhx"%(reptree), format_root_node=True)

def placeNTransitionsInTree_new(numTransitions, maxTransitions, maxConvRate, tree_ini, manual_mode_nodes = {}, nf=""):

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
        else:
            tree_all_transitions, nodesWithTransitions, observedNumTransitions = randomTransitions_new(maxTransitions, tree_all_transitions)
            #logger.debug ("Observed Number of Transitions: "+ str(observedNumTransitions ) + " compared to "+ str(maxTransitions) + " wanted in:" + nf )

    #logger.debug("numTries: %s", numTries)
    if numTries < 500:
        numTries3 = 0
        nodesWithTransitions_allCombinations = list(itertools.combinations(nodesWithTransitions,maxTransitions-numTransitions))
        nodesWithTransitions_list = random.sample(nodesWithTransitions_allCombinations, len(nodesWithTransitions_allCombinations))

        while maxConvRate < maxConvRateObs and nodesWithTransitions_list and numTries3 < 1000:
            numTries3 = numTries3 + 1
            #logger.debug("numTries3: %s", numTries3)
            #delete maxTransitions-numTransitions random convergence events
            tree_final = tree_all_transitions.copy(method="deepcopy")
            lr=nodesWithTransitions_list.pop()
            if not nodesWithTransitions_list:
                numTries3 = 1000
            for n_ND in lr:
                #print n_ND
                node = tree_final.search_nodes(ND=n_ND)[0]
                #print node.get_ascii(attributes=["ND","name","C","T"])
                setsubtree_new(node,False)
                node.T=False
                #print node.get_ascii(attributes=["ND","name","C","T"])

            #print nodesWithTransitions
            #print lr
            nodesWithTransitions_filtered=[k for k in nodesWithTransitions if not k in lr]
            #print nodesWithTransitions_filtered
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
    nodesWithTransitions = []
    nodesWithConvergentModel = []
    nodesWithAncestralModel = []
    tree_final = Tree()
  else:
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

  return tree_final, nodesWithTransitions, nodesWithConvergentModel, nodesWithAncestralModel, numberOfLeafsWithTransitions, numberOfLeafs
