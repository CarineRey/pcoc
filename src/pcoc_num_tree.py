#!/usr/bin/python
#  pcoc_num_tree.py
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

import argparse
import time
from ete3 import Tree, NodeStyle, TreeStyle, TextFace

##########
# inputs #
##########
start_time = time.time()

### Option defining
parser = argparse.ArgumentParser(prog="pcoc_num_tree.py",
                                 description='')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

##############
requiredOptions = parser.add_argument_group('Required arguments')
requiredOptions.add_argument('-t', "--tree", type=str,
                             help='input tree name', required=True)
requiredOptions.add_argument('-o', '--output', type=str,
                   help="Output file (pdf)", required=True)
##############
Options = parser.add_argument_group('Options')
Options.add_argument('-m', '--manual_mode_test', type=str, metavar="\"x/y,z/...\"",
                    help="Highlight, in the output plot, user defined convergent transition/branches. Transition node must be the first number and independent events must be separed by a \"/\". ex: \"1,2,3/67/55,56\" (default: None)",
                    default="")
##############

### Option parsing
args = parser.parse_args()


def init_tree(nf):
    t = Tree(nf)

    #Alternatively we could read a tree from a file into a string "line", and then use:
    # t =  Tree( line )

    nodeId = 0
    for n in t.traverse("postorder"):
        n.add_features(ND=nodeId)
        nodeId = nodeId + 1

    return t


# Basic tree style
tree_style = TreeStyle()
tree_style.show_leaf_name = False
tree_style.show_branch_length = False
tree_style.draw_guiding_lines = True
tree_style.complete_branch_lines_when_necessary = True

nstyle = NodeStyle()
nstyle["fgcolor"] = "black"
nstyle["size"] = 0

nstyle_L = NodeStyle()
nstyle_L["fgcolor"] = "black"
nstyle_L["size"] = 0

tree = init_tree(args.tree)

manual_mode_nodes = {}
if args.manual_mode_test:
    manual_mode_nodes = {"T":[], "C":[]}
    p_events = args.manual_mode_test.strip().split("/")
    for e in p_events:
        l_e = map(int, e.split(","))
        manual_mode_nodes["T"].append(l_e[0])
        manual_mode_nodes["C"].extend(l_e[1:])

for n in tree.traverse():
    if n.is_leaf():
        n.set_style(nstyle_L)
        n.add_face(TextFace(str(n.name)), column=0, position="aligned")
    else:
        n.set_style(nstyle)
    nd = TextFace(str(n.ND))

    if manual_mode_nodes:
        if n.ND in manual_mode_nodes["T"]:
            nd.background.color = "red"
        elif n.ND in manual_mode_nodes["C"]:
            nd.background.color = "orange"
        else:
            nd.background.color = "white"
    else:
        nd.background.color = "white"
    nd.margin_right = 2
    nd.margin_top = 1
    nd.margin_left = 2
    nd.margin_bottom = 1
    nd.border.width = 1
    n.add_face(nd, column=0, position="float")
    n.add_face(TextFace("       "), column=0, position="branch-bottom")

tree.render(args.output, tree_style=tree_style)
print args.output
