#! /bin/bash
set -e

mkdir -p test_det
cd test_det
ROOT_DIR=$PWD
OUTPUT_DIR=$ROOT_DIR/out_test_det/
rm -r $OUTPUT_DIR || ls $ROOT_DIR
mkdir -p  $OUTPUT_DIR

mkdir -p example
cd example

## Get the data
rm -r data || echo ""
mkdir -p data
### tree
wget https://raw.githubusercontent.com/gilles-didier/Convergence/master/data/phylogeny.nwk -P data/
tree=$PWD/data/phylogeny.nwk
tree_with_cond=$PWD/data/phylogeny_with_cond.nwk
echo "((((Microbat:0.152207[&&NHX:Condition=1:Transition=1],Megabat:0.118463[&&NHX:Condition=0]):0.0219285[&&NHX:Condition=0],((Dolphin:0.0710076[&&NHX:Condition=1],Cow:0.119764[&&NHX:Condition=1]):0.0182287[&&NHX:Condition=1:Transition=1],Alpaca:0.119871[&&NHX:Condition=0]):0.0352354[&&NHX:Condition=0]):0.0262836[&&NHX:Condition=0],((Marmoset:0.0667873[&&NHX:Condition=0],Human:0.0427597[&&NHX:Condition=0]):0.0844445[&&NHX:Condition=0],Mouse:0.344614[&&NHX:Condition=0]):0.00810167[&&NHX:Condition=0]):0.00560883[&&NHX:Condition=0],Elephant:0.168097[&&NHX:Condition=0]);" > $tree_with_cond

### alignment
wget https://raw.githubusercontent.com/gilles-didier/Convergence/master/data/alignments/SLC26A5.fasta -P data/
ali=$PWD/data/SLC26A5.fasta

DOCKER_CMD="docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc"

## First identify the nodes with convergent transitions with pcoc_num_tree.py
$DOCKER_CMD pcoc_num_tree.py -t $tree -o $PWD/data/num_tree.pdf

## In the num_tree.pdf, get the numbers of the branches with the convergent transitions.
## Here we are interested in the branches leading to the Microbat and the Dolphin, so nodes 0 and 3:
## You can check your combination using the -m parameter, which will produce a pdf of the tree where the convergent branches lead to red nodes.
$DOCKER_CMD pcoc_num_tree.py -t $tree -o $PWD/data/colored_num_tree.pdf -m 0/3

## Then run pcoc_det on your alignment uding the tree and the scenario
$DOCKER_CMD pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3

## Then run pcoc_det on your alignment using the annotated tree
$DOCKER_CMD pcoc_det.py -t $tree_with_cond -aa $ali -o $PWD/output_pcoc_det -m "-"

## Then run pcoc_det on your alignment
$DOCKER_CMD pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --gamma

## if you want to have a nice plot of the sites detected as convergent you can use the option --plot with or without the option --reorder
## With the --reorder option:
$DOCKER_CMD pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --plot --reorder
## Without the --reorder option:
$DOCKER_CMD pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --plot -plot_title "Test PCOC"
## if you want to have a nice plot of your complete alignment you can use --plot_complete_ali
$DOCKER_CMD pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --plot_complete_ali --plot --reorder -plot_title "Test PCOC"
$DOCKER_CMD pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --plot_complete_ali --plot --reorder --inv_gamma -plot_title "Test PCOC"
$DOCKER_CMD pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --plot_complete_ali --plot --reorder --gamma -plot_title "Test PCOC"


## another example with several convergent branches, with OC applied on branches 11, 6, and 15, and PC applied on all of 11,9,10,6,15:
$DOCKER_CMD pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 11,9,10/6/15 --plot -plot_title "Test PCOC"


## Test with lacking seqs:
ali=$PWD/data/SLC26A5.fasta
ali_trim=$PWD/data/SLC26A5_trim.fasta

tail -n 56  $ali > $ali_trim #remove elephant (15)
$DOCKER_CMD pcoc_det.py -t $tree -aa $ali_trim -o $PWD/output_pcoc_det -m 11,9,10/6/15 --plot -plot_title "Test PCOC" --auto_trim_tree

