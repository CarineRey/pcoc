from ete3 import Tree
import commands
import random
import itertools
import comp_seq
import copy


import logging
logger = logging.getLogger("simul_CAT_rand.placing")


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

def make_simul(name,nodesWithAncestralModel,nodesWithTransitions,nodesWithConvergentModel,c1,c2,repseq,reptree,repbppconfig,number_of_sites=1000,outputInternalSequences="yes",nbCAT=10,cz_nodes={},CzOneChange=True):
  if outputInternalSequences != "yes":
      outputInternalSequences = "no"

  if c1 != c2:
      number_of_models = 3
  else:
      number_of_models = 1
  sup_command = ""


  if cz_nodes:
      for (cz, nodes) in cz_nodes.items():
          if nodes:
            if CzOneChange:
                number_of_models +=1
                sup_command+=" model%s=\'OneChange(model=LGL08_CAT_C%s(nbCat=$(NBCAT)))\' " %(number_of_models,cz)
                t_node = nodes[0]
                sup_command+=" model%s.nodes_id=\"%s\" " %(number_of_models,str(t_node))
                if len(nodes) > 1:
                    number_of_models +=1
                    sup_command+=" model%s=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(number_of_models,cz)
                    sup_command+=" model%s.nodes_id=\"%s\" " %(number_of_models,",".join(map(str,nodes[1:])))
            else:
                number_of_models +=1
                sup_command+=" model%s=\'LGL08_CAT_C%s(nbCat=$(NBCAT))\' " %(number_of_models,cz)
                sup_command+=" model%s.nodes_id=\"%s\" " %(number_of_models,",".join(map(str,nodes)))
            


  command="bppseqgen NAME=%s REP_SEQ=%s REP_TREE=%s NBCAT=%s number_of_sites=%s output.internal.sequences=%s " %(name,repseq.replace("//","/"), reptree, nbCAT,number_of_sites,outputInternalSequences)
  if c1!=c2:

    n1="\""+ ",".join(map(str,nodesWithAncestralModel))+"\""
    n2="\""+ ",".join(map(str,nodesWithTransitions))+"\""
    n3="\""+ ",".join(map(str,nodesWithConvergentModel))+"\""

    command+="param=%s.bpp mod1Nodes=%s mod2Nodes=%s mod3Nodes=%s Ne1=%d Ne2=%d" %(repbppconfig+"/CATseq_sim_conv",n1,n2,n3,c1,c2)

  else:
    allNodes = nodesWithConvergentModel+nodesWithTransitions+nodesWithAncestralModel
    n="\""+ ",".join(map(str,allNodes))+"\""
    command+="param=%s.bpp mod1Nodes=%s Ne1=%d" %(repbppconfig+"/CATseq_sim_noconv", n,c1)

  command+=sup_command + " nonhomogeneous.number_of_models=%s " %(number_of_models)
  #logger.debug(command)
  out = commands.getoutput(command)
  #logger.debug(out)

def make_estim(name,nodesWithAncestralModel,nodesWithTransitions,nodesWithConvergentModel,c1,c2,repseq,reptree,repest,repbppconfig, test_c1c1 = False,NBCATest=10,suffix="", OneChange = True, ext = ".fa"):

  n1="\""+ ",".join(map(str,nodesWithAncestralModel))+"\""
  fasta_file = "%s/%s%s" %(repseq,name,ext)

  #logger.debug("fasta_file: %s",fasta_file )
  if OneChange:
    n2="\""+ ",".join(map(str,nodesWithConvergentModel))+"\""
    n3="\""+ ",".join(map(str,nodesWithTransitions))+"\""
    raw_command = "bppml param=%s NAME=%s SUFFIX=%s REP_SEQ=%s REP_TREE=%s REP_EST=%s mod1Nodes=%s mod2Nodes=%s  mod3Nodes=%s \"input.sequence.file=%s\" NBCAT=%s "%(repbppconfig+"/CATseq_OneChange.bpp",name,suffix,repseq,reptree,repest,n1,n2,n3,fasta_file,NBCATest)


  else:
    n2="\""+ ",".join(map(str,nodesWithConvergentModel+nodesWithTransitions))+"\""
    n3 = ""
    raw_command = "bppml param=%s NAME=%s SUFFIX=%s REP_SEQ=%s REP_TREE=%s REP_EST=%s mod1Nodes=%s mod2Nodes=%s \"input.sequence.file=%s\" NBCAT=%s "%(repbppconfig+"/CATseq.bpp",name,suffix,repseq,reptree,repest,n1,n2,fasta_file,NBCATest)

  #command= raw_command + "Ne1=%d Ne2=%d" %(c1,c2)
  command= raw_command + "Ne1=%d Ne2=%d" %(c1,c2)
  #logger.debug(command)
  out = commands.getoutput(command)
  #logger.debug(out)

  if c1!=c2 and test_c1c1:
    command = raw_command + "Ne1=%d Ne2=%d" %(c1,c1)
    logger.debug(command)
    out = commands.getoutput(command)
    logger.debug(out)


def make_estim_conv(name,nodes,c1,repseq,reptree,repest,repbppconfig,suffix="",NBCATest=10):

  allNodes=nodes

  n1="\""+ ",".join(map(str,allNodes))+"\""

  command="bppml param=%s NAME=%s  SUFFIX=%s REP_SEQ=%s REP_TREE=%s REP_EST=%s mod1Nodes=%s Ne1=%d  NBCAT=%s  \"input.sequence.file=%s/%s.fa\" "%(repbppconfig + "/CATseq_conv.bpp",name, suffix, repseq,reptree,repest,n1,c1,NBCATest,repseq,name)

  #logger.debug(command)
  out = commands.getoutput(command)
  #logger.debug(out)

