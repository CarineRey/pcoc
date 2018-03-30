#!/usr/bin/python
#  plot_data.py
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

from ete3 import NodeStyle, TreeStyle, faces, CircleFace, PhyloTree, \
                 TextFace, add_face_to_node, SeqMotifFace, TreeFace, \
                 Tree, StaticItemFace

from PyQt4.QtGui import (QGraphicsRectItem, QGraphicsLineItem,
                         QGraphicsPolygonItem, QGraphicsEllipseItem,
                         QPen, QColor, QBrush, QPolygonF, QFont,
                         QPixmap, QFontMetrics, QPainter,
                         QRadialGradient, QGraphicsSimpleTextItem, \
                         QGraphicsTextItem, QGraphicsItem)

import types
import sys

import logging
logger = logging.getLogger("pcoc.plot_data")


#### SequencePlotFace methods modification
class SequencePlotFace_mod(faces.SequencePlotFace):
    def draw_x_axis(self):
            lineItem = QGraphicsLineItem(self.col_w/2,
                                               self.coordY(self.ylim[0])+2,
                                               self.width-self.col_w/2,
                                               self.coordY(self.ylim[0])+2,
                                               parent=self.item)
            lineItem.setPen(QPen(QColor('black')))
            lineItem.setZValue(10)
            all_vals = list(range(0, len(self.values), self.x_inter_values))
            if (len(self.values)-1)%self.x_inter_values:
                all_vals += [len(self.values)-1]

            hp_x = []
            if self.hp:
                for x in list(range(0, len(self.values))):
                    if self.x_values[x] in self.hp:
                        hp_x.append(x)
                        if not x in all_vals:
                            all_vals += [x]
            all_vals.sort()

            for x in all_vals:
                lineItem = QGraphicsLineItem(0, self.coordY(self.ylim[0])+2,
                                                   0, self.coordY(self.ylim[0])+6,
                                                   parent=self.item)
                lineItem.setX(x*self.col_w + self.col_w/2)
                lineItem.setPen(QPen(QColor('black')))
                lineItem.setZValue(10)
                if x in hp_x:
                    text = QGraphicsSimpleTextItem("*" + str(self.x_values[x]))
                    qfont = QFont("Arial", self.fsize -1)
                    #qfont.setBold(True)
                    text.setFont(qfont)
                else:
                    text = QGraphicsSimpleTextItem(" " + str(self.x_values[x]))
                    text.setFont(QFont("Arial", self.fsize - 1))
                text.rotate(-90)
                text.setParentItem(self.item)
                text.setZValue(10)
                tw = text.boundingRect().width()
                th = text.boundingRect().height()
                # Center text according to masterItem size
                text.setPos( x*self.col_w - th/2  + self.col_w/2, tw + self.coordY(self.ylim[0])+ 7)


    def draw_y_axis(self):
            lineItem = QGraphicsLineItem(0, self.coordY(self.ylim[0]),
                                         0, self.coordY(self.ylim[1]),
                                         parent=self.item)
            lineItem.setPen(QPen(QColor('black')))
            lineItem.setZValue(10)
            max_w = 0
            for y in set(self.hlines + list(self.ylim)):
                if y in list(self.ylim):
                    lineItem = QGraphicsLineItem(0, self.coordY(y),
                                                   -5, self.coordY(y),
                                                   parent=self.item)
                    lineItem.setPen(QPen(QColor('black')))
                    lineItem.setZValue(10)
                    text = QGraphicsSimpleTextItem(str(y))
                    text.setFont(QFont("Arial", self.fsize-2))
                    text.setParentItem(self.item)
                    tw = text.boundingRect().width()
                    max_w = tw if tw > max_w else max_w
                    th = text.boundingRect().height()
                    # Center text according to masterItem size
                    text.setPos(-tw - 5, self.coordY(y)-th/2)
                else:
                    text = QGraphicsSimpleTextItem(str(y))
                    text.setFont(QFont("Arial", self.fsize-4))
                    text.setParentItem(self.item)
                    tw = text.boundingRect().width()
                    max_w = tw if tw > max_w else max_w
                    th = text.boundingRect().height()
                    # Center text according to masterItem size
                    text.setPos(self.width + 5, self.coordY(y)-th/2)

            if self.ylabel:
                text = QGraphicsSimpleTextItem(self.ylabel)
                text.setFont(QFont("Arial", self.fsize-1))
                text.setParentItem(self.item)
                text.rotate(-90)
                tw = text.boundingRect().width()
                th = text.boundingRect().height()
                # Center text according to masterItem size
                text.setPos(-th -5-max_w, tw/2+self.coordY(sum(self.ylim)/2))

    def set_sticks_color(self):
        stick_colors = []
        y_values = self.values
        for y_val in y_values:
            if y_val >= 0.99:
                stick_colors.append("red")
            elif y_val >= 0.90 :
                stick_colors.append("orange")
            elif y_val >= 0.80 :
                stick_colors.append("#EFDB00")
            else:
                stick_colors.append("gray")
        self.colors = stick_colors

    def draw_colored_boxes(self, h):

            stick_colors = self.colors

            num_col_red = stick_colors.count("red")
            num_col_orange = stick_colors.count("orange")
            num_col_yellow = stick_colors.count("#EFDB00")
            num_col_gray = stick_colors.count("gray")

            colored_box_h = self.height + h - 5

            if num_col_red:
                rect_red = QGraphicsRectItem(self.col_w * 0 - 1 , self.coordY(self.ylim[1])- 5, self.col_w * num_col_red , colored_box_h)
                qpen = QPen(QColor('red'))
                qpen.setWidth(1)
                qpen.setStyle(5) # dash line : http://doc.qt.io/qt-4.8/qt.html#PenStyle-enum
                rect_red.setPen(qpen)
                rect_red.setBrush(QColor("#FFD1D1"))
                rect_red.setZValue(-1)
                rect_red.setParentItem(self.item)
                #rect_red.setOpacity(0.5)

            if num_col_orange:
                rect_orange = QGraphicsRectItem(self.col_w * num_col_red + 1 , self.coordY(self.ylim[1])- 5, self.col_w * num_col_orange - 2 , colored_box_h)
                qpen = QPen(QColor('orange'))
                qpen.setWidth(1)
                qpen.setStyle(5) # dash line : http://doc.qt.io/qt-4.8/qt.html#PenStyle-enum
                rect_orange.setPen(qpen)
                rect_orange.setBrush(QColor("#FFE3B0"))
                rect_orange.setZValue(-1)
                rect_orange.setParentItem(self.item)
                #rect_orange.setOpacity(0.5)

            if num_col_yellow:
                rect_yellow = QGraphicsRectItem(self.col_w * (num_col_orange + num_col_red) + 1 , self.coordY(self.ylim[1])- 5, self.col_w * num_col_yellow - 2 , colored_box_h)
                qpen = QPen(QColor('#EFDB00'))
                qpen.setWidth(1)
                qpen.setStyle(5) # dash line : http://doc.qt.io/qt-4.8/qt.html#PenStyle-enum
                rect_yellow.setPen(qpen)
                rect_yellow.setBrush(QColor("#FBFFA5"))
                rect_yellow.setZValue(-1)
                rect_yellow.setParentItem(self.item)
                #rect_yellow.setOpacity(0.5)

            if num_col_gray:
                rect_gray = QGraphicsRectItem(self.col_w * (num_col_orange + num_col_red + num_col_yellow)  + 1 , self.coordY(self.ylim[1])- 5, self.col_w * num_col_gray, colored_box_h)
                qpen = QPen(QColor('gray'))
                qpen.setWidth(1)
                qpen.setStyle(5) # dash line : http://doc.qt.io/qt-4.8/qt.html#PenStyle-enum
                rect_gray.setPen(qpen)
                rect_gray.setBrush(QColor("#E2E2E2"))
                rect_gray.setZValue(-1)
                rect_gray.setParentItem(self.item)
                #rect_gray.setOpacity(0.5)

    def update_items(self):
            self.item =  QGraphicsRectItem(-30, 0, self.width+40, self.height+70)
            #self.item.setPen(QPen(QColor('gray')))
            self.item.setPen(QPen(QColor('white')))

            try:
                put_colored_boxes = self.put_colored_boxes
            except AttributeError:
                put_colored_boxes = (False, 0)

            if put_colored_boxes[0]:
                (_, tree_h) = put_colored_boxes
                self.draw_colored_boxes(tree_h)

            # draw lines
            for line, col in zip(self.hlines, self.hlines_col):
                self.draw_hlines(line, col)
            # draw plot
            width = self.col_w
            for i, val in enumerate(self.values):
                self.draw_fun(width * i + self.col_w / 2 , val, i)
            # draw error bars
            if self.errors:
                for i in range(len(self.errors)):
                    self.draw_errors(width * i + self.col_w / 2 , i)
            # draw x axis
            self.draw_x_axis()
            # draw y axis
            self.draw_y_axis()
            # put header
            self.write_header()


### New StaticItemFace class

class SequenceScoreFace(StaticItemFace):
    def __init__(self, dict_values_pcoc, col_height=10, col_width=10, fsize=9):
        faces.Face.__init__(self)
        self.type = "item"
        self.item = None
        self.dict_values_pcoc = dict_values_pcoc
        self.nb_values = len(dict_values_pcoc.values()[1])
        self.col_w = float(col_width)
        self.col_h = float(col_height)
        self.fsize = fsize
        self.nb_models = len(dict_values_pcoc.keys())
        self.height = (self.nb_models + 1) * self.col_h
        self.width = self.nb_values * self.col_w
        self.x_axis = False
        self.hp = []

    def draw_fun(self, x, y, val, col_width = 10, col_height = 10, color = "gray"):
        color = get_corr_color(val)
        rect = QGraphicsRectItem(x, y, col_width, col_height, parent = self.item)
        rect.setPen(QPen(QColor('black')))
        rect.setBrush(QColor(color))

    def draw_legend(self):
        legend_h = self.height * ((self.nb_models - 1) / float(self.nb_models))
        if legend_h < 35:
            legend_h = 35
        legend_rect = QGraphicsRectItem(-20, 0, 10, legend_h, parent=self.item)
        x0 = -20
        n_cat = 6.
        for y, str_y in [(1,1),(1/n_cat*5,0.99), (1/n_cat*4, 0.9), (1/n_cat*3, 0.8), (1/n_cat*2, 0.7), (1/n_cat*1, 0.5) , (1/n_cat*0,0)]:
            y_stick = legend_h - y * legend_h
            lineItem = QGraphicsLineItem(x0 - 5, y_stick , x0, y_stick,
                                               parent=self.item)
            lineItem.setPen(QPen(QColor('black')))
            text = QGraphicsSimpleTextItem(str(str_y))
            text.setFont(QFont("Arial", self.fsize-4))
            text.setParentItem(self.item)
            tw = text.boundingRect().width()
            th = text.boundingRect().height()
            # Center text according to masterItem size
            text.setPos(x0 - tw - 7, y_stick - th/2)

        for (y1,y2,c) in [(1,1/n_cat*5, 0.99), (1/n_cat*5, 1/n_cat*4, 0.9), (1/n_cat*4, 1/n_cat*3, 0.8), (1/n_cat*3, 1/n_cat*2, 0.7), (1/n_cat*2, 1/n_cat*1, 0.5), (1/n_cat*1, 1/n_cat*0, 0)]:
            y1_stick = legend_h - y1 * legend_h
            y2_stick = legend_h - y2 * legend_h
            self.draw_fun(x0, y1_stick, c, col_width = 10, col_height = y2_stick-y1_stick)


    def draw_x_axis(self):
        y0 = self.nb_models * self.col_h + 5
        lineItem = QGraphicsLineItem(self.col_w/2,
                                           y0,
                                           self.width-self.col_w/2,
                                           y0,
                                           parent=self.item)
        lineItem.setPen(QPen(QColor('black')))
        lineItem.setZValue(10)
        all_vals = list(range(0, self.nb_values, self.x_inter_values))
        if (self.nb_values-1)%self.x_inter_values:
            all_vals += [self.nb_values-1]

        hp_x = []
        if self.hp:
            for x in list(range(0, self.nb_values)):
                if self.x_values[x] in self.hp:
                    hp_x.append(x)
                    if not x in all_vals:
                        all_vals += [x]
        all_vals.sort()

        for x in all_vals:
            lineItem = QGraphicsLineItem(0, y0,
                                               0, y0+6,
                                               parent=self.item)
            lineItem.setX(x*self.col_w + self.col_w/2)
            lineItem.setPen(QPen(QColor('black')))
            lineItem.setZValue(10)
            if x in hp_x:
                text = QGraphicsSimpleTextItem("*" + str(self.x_values[x]))
                qfont = QFont("Arial", self.fsize -1)
                #qfont.setBold(True)
                text.setFont(qfont)
            else:
                text = QGraphicsSimpleTextItem(" " + str(self.x_values[x]))
                text.setFont(QFont("Arial", self.fsize - 1))
            text.rotate(-90)
            text.setParentItem(self.item)
            text.setZValue(10)
            tw = text.boundingRect().width()
            th = text.boundingRect().height()
            # Center text according to masterItem size
            text.setPos( x*self.col_w - th/2 + self.col_w/2, tw + y0 + 7)


    def update_items(self):
        rect_h = self.height
        if self.x_axis:
            rect_h += 30
        self.item = QGraphicsRectItem(0, 0, self.width + 40 , rect_h)
        self.item.setPen(QPen(QColor('white')))

        #X axis

        if self.x_axis:
            self.draw_x_axis()

        # Legend

        self.draw_legend()

        # Y axes and colo rect
        yi = -1
        for model in ["PCOC", "PC", "OC", "Topological", "Identical"]:
            if self.dict_values_pcoc.has_key(model):
                yi += 1
                y = yi * self.col_w

                # Y axes
                ## Stick
                yaxis = (yi + 0.5) * self.col_w
                lineItem = QGraphicsLineItem(self.width, yaxis,
                                             self.width + 5, yaxis,
                                             parent=self.item)
                lineItem.setPen(QPen(QColor('black')))
                ## Text
                text = QGraphicsSimpleTextItem(model)
                text.setFont(QFont("Arial", self.fsize-2))
                text.setParentItem(self.item)
                tw = text.boundingRect().width()
                th = text.boundingRect().height()
                ## Center text according to masterItem size
                text.setPos(self.width + 5, yaxis - th/2)

                # Color rect for each model
                values = self.dict_values_pcoc[model]
                for i, val in enumerate(values):
                    self.draw_fun( i * self.col_w , y, val,
                                   col_width=self.col_w)


### Function

def get_corr_color(x):
    score = [1,0.99,0.9,0.8,0.7,0.5,0]
    color = ["red", "red",  "orange", "#EFDB00", "#6BAC00", "#7174D0",  "#A3A3A3"]

    for i in range(len(score)):
        if x >= score[i]:
            res = color[i]
            break
    return(res)



### Tree style

# Basic tree style
tree_style = TreeStyle()
tree_style.show_leaf_name = True
tree_style.show_branch_length = True
tree_style.min_leaf_separation  = 4
tree_style.branch_vertical_margin   = 2

nstyle_T_sim = NodeStyle()
nstyle_T_sim["fgcolor"] = "orange"
nstyle_T_sim["size"] = 5
nstyle_T_sim["vt_line_color"] = "orange"
nstyle_T_sim["hz_line_color"] = "orange"
nstyle_T_sim["vt_line_width"] = 2
nstyle_T_sim["hz_line_width"] = 2

nstyle_T_sim_est = NodeStyle()
nstyle_T_sim_est["fgcolor"] = "orange"
nstyle_T_sim_est["bgcolor"] = "#CECECE"
nstyle_T_sim_est["size"] = 5
nstyle_T_sim_est["vt_line_color"] = "orange"
nstyle_T_sim_est["hz_line_color"] = "orange"
nstyle_T_sim_est["vt_line_width"] = 2
nstyle_T_sim_est["hz_line_width"] = 2

nstyle_T_est = NodeStyle()
nstyle_T_est["fgcolor"] = "#0052FF"
nstyle_T_est["bgcolor"] = "#CECECE"
nstyle_T_est["size"] = 5
nstyle_T_est["vt_line_color"] = "#0052FF"
nstyle_T_est["hz_line_color"] = "#0052FF"
nstyle_T_est["vt_line_width"] = 2
nstyle_T_est["hz_line_width"] = 2


nstyle_C_sim = NodeStyle()
nstyle_C_sim["fgcolor"] = "orange"
nstyle_C_sim["size"] = 5
nstyle_C_sim["vt_line_color"] = "orange"
nstyle_C_sim["hz_line_color"] = "orange"
nstyle_C_sim["vt_line_width"] = 2
nstyle_C_sim["hz_line_width"] = 2

nstyle_C_sim_est = NodeStyle()
nstyle_C_sim_est["fgcolor"] = "orange"
nstyle_C_sim_est["bgcolor"] = "#FFE9AD"
nstyle_C_sim_est["size"] = 5
nstyle_C_sim_est["vt_line_color"] = "orange"
nstyle_C_sim_est["hz_line_color"] = "orange"
nstyle_C_sim_est["vt_line_width"] = 2
nstyle_C_sim_est["hz_line_width"] = 2


nstyle_C_est = NodeStyle()
nstyle_C_est["fgcolor"] = "#0052FF"
nstyle_C_est["bgcolor"] = "#FFE9AD"
nstyle_C_est["size"] = 5
nstyle_C_est["vt_line_color"] = "#0052FF"
nstyle_C_est["hz_line_color"] = "#0052FF"
nstyle_C_est["vt_line_width"] = 2
nstyle_C_est["hz_line_width"] = 2


nstyle = NodeStyle()
nstyle["fgcolor"] = "#0052FF"
nstyle["size"] = 5
nstyle["hz_line_width"] = 2
nstyle["vt_line_width"] = 2
nstyle["vt_line_color"] = "#0052FF"
nstyle["hz_line_color"] = "#0052FF"


nstyle_T2 = NodeStyle()
nstyle_T2["fgcolor"] = "orange"
nstyle_T2["size"] = 5
nstyle_T2["vt_line_color"] = "#0052FF"
nstyle_T2["hz_line_color"] = "orange"
nstyle_T2["vt_line_width"] = 2
nstyle_T2["hz_line_width"] = 2



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

def add_sim_root(node):
    nd = TextFace("R")
    nd.fsize = 4
    nd.background.color = "#ABABFF"
    nd.margin_right = 0
    nd.margin_top = 0
    nd.margin_left = 0
    nd.margin_bottom = 0
    nd.border.width = 1
    nd2 = TextFace(" ")
    nd2.fsize = 4
    

    node.add_face(nd, column=0, position = "float")
    node.add_face(nd2, column=1, position = "float")

def make_tree_ali_detect_combi(g_tree, ali_nf, Out,
                               hist_up = "",
                               x_values=[], hp = [],
                               dict_benchmark = {},
                               reorder = False,
                               det_tool=False):
    reptree = g_tree.reptree
    cz_nodes = g_tree.cz_nodes
    ### Tree

    ## Tree style
    phylotree_style = TreeStyle()
    phylotree_style.show_leaf_name = False
    phylotree_style.show_branch_length = False
    phylotree_style.draw_guiding_lines  = True
    phylotree_style.min_leaf_separation  = 1

    ## For noisy tree
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
    
    sim_root_ND = g_tree.outgroup_ND

    def my_layout(node):
        ## Sequence name
        F = TextFace(node.name, tight_text=True)
        add_face_to_node(F, node, column=0, position="aligned")

        ## Sequence motif
        if node.is_leaf():
            motifs_n = []
            box_color = "black"
            opacity = 1
            if node.T == "True" or node.C == "True":
                motifs_n.append([0, len(node.sequence), "[]", 10, 12, box_color, box_color, None])

            motifs_n.append([0, len(node.sequence), "seq", 10, 10, None, None, None])
            seq_face = SeqMotifFace(seq=node.sequence, seqtype='aa',
                                    seq_format='seq',
                                    fgcolor = box_color,
                                    motifs=motifs_n)
            seq_face.overlaping_motif_opacity = opacity
            add_face_to_node(seq_face, node, column=1, position='aligned')

        ## Nodes style
        if det_tool and node.T == "True":
            node.set_style(nstyle_T_sim)
            add_t(node)
        elif det_tool and node.C == "True":
            node.set_style(nstyle_C_sim)
        elif det_tool:
            node.set_style(nstyle)
        #if not det_tool no background
        elif node.T == "True" and not int(node.ND) in g_tree.conv_events.nodesWithTransitions_est:
            node.set_style(nstyle_T_sim)
            add_t(node)
        elif node.T == "True" and int(node.ND) in g_tree.conv_events.nodesWithTransitions_est:
            node.set_style(nstyle_T_sim_est)
            add_t(node)
        elif int(node.ND) in g_tree.conv_events.nodesWithTransitions_est:
            node.set_style(nstyle_T_est)
        elif node.C == "True" and not int(node.ND) in g_tree.conv_events.nodesWithConvergentModel_est:
            node.set_style(nstyle_C_sim)
        elif node.C == "True" and int(node.ND) in g_tree.conv_events.nodesWithConvergentModel_est:
            node.set_style(nstyle_C_sim_est)
        elif int(node.ND) in g_tree.conv_events.nodesWithConvergentModel_est:
            node.set_style(nstyle_C_est)
        elif cz_nodes_s and node.Cz != "False":
            node.set_style(cz_nodes_s[int(node.Cz)])
            if int(node.ND) == int(cz_nodes[int(node.Cz)][0]):
                add_t(node)
        else:
            node.set_style(nstyle)
        
        if int(node.ND) == sim_root_ND and not det_tool:
            add_sim_root(node)

    phylotree_style.layout_fn = my_layout

    # Get tree dimensions
    tree_nf = g_tree.annotated_tree_fn_est
    logger.debug("tree_nf: %s",tree_nf)
    t = PhyloTree(tree_nf)
    t.link_to_alignment (ali_nf)
    t.render(Out)
    tree_face = TreeFace(t , phylotree_style)
    tree_face.update_items()
    tree_h = tree_face._height()
    tree_w = tree_face._width()

    ### X axes:
    if not x_values: # Complete representation
        x_values_up = [x+1 for x in range(0,len(dict_benchmark.values()[1]))]
        inter = 5
    else: # Filtered representation
        x_values_up = x_values
        inter = 1


    ### Histogram up
    if hist_up in ["PCOC", "PC", "OC", "Topological", "Identical"]:
        header_hist_value_up = 'Posterior probability (' + hist_up.upper() +' model)'
        hist_value_up = dict_benchmark[hist_up]
        # Define emphased lines
        hlines=[0.8, 0.9,0.99]
        hlines_col=['#EFDB00', 'orange', 'red']

        # Type of representation
        kind= 'stick' # bar/curve/stick



        y_values_up = hist_value_up

        hist = SequencePlotFace_mod(y_values_up, hlines=hlines,
                                        hlines_col=hlines_col,  kind=kind,
                                        header=header_hist_value_up,
                                        height=40, col_width=10, ylim = [0,1])

        hist.hp = hp
        hist.x_values = x_values_up
        hist.x_inter_values = inter
        hist.set_sticks_color()


        if reorder:
            # draw colored boxes
            hist.put_colored_boxes = (True, tree_h + 40)

        phylotree_style.aligned_header.add_face(hist, 1)


    ### Rect all model:
    sequencescorebox = SequenceScoreFace(dict_benchmark, col_width=10)
    sequencescorebox.hp = hp
    sequencescorebox.x_values = x_values_up
    sequencescorebox.x_inter_values = inter
    if not hist_up in ["PCOC", "PC", "OC", "Topological", "Identical"]:
        sequencescorebox.x_axis=True

    phylotree_style.aligned_header.add_face(sequencescorebox, 1)


    tree_nf = reptree + "/annotated_tree.nhx"
    logger.debug("tree_nf: %s",tree_nf)

    res = t.render(Out, tree_style=phylotree_style)
    del t
    return(res)

def mk_bilan_tree(tree, bilan_nodesWithTransitions,output):
    d_nodesWithTransitions = {}
    for x in bilan_nodesWithTransitions:
        if x in d_nodesWithTransitions:
            d_nodesWithTransitions[x]+=1
        else:
            d_nodesWithTransitions[x]=1

    for n in tree.traverse():
        if n.ND in d_nodesWithTransitions:
            n.set_style(nstyle_T2)
            C = CircleFace(radius=d_nodesWithTransitions[n.ND] * 2 , color="RoyalBlue", style="sphere", label = str(d_nodesWithTransitions[n.ND]))
            C.opacity = 0.5
            n.add_face(C, 0, position="float")
        else:
            n.set_style(nstyle)

    tree.render(output, tree_style=tree_style)
