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
### alignment
wget https://raw.githubusercontent.com/gilles-didier/Convergence/master/data/alignments/SLC26A5.fasta -P data/
ali=$PWD/data/SLC26A5.fasta

## First identify the nodes with convergent transitions with pcoc_num_tree.py
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_num_tree.py -t $tree -o $PWD/data/num_tree.pdf

## In the num_tree.pdf, get the numbers of the branches with the convergent transitions.
## Here we are interested in the branches leading to the Microbat and the Dolphin, so nodes 0 and 3:
## You can check your combination using the -m parameter, which will produce a pdf of the tree where the convergent branches lead to red nodes.
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_num_tree.py -t $tree -o $PWD/data/colored_num_tree.pdf -m 0/3

## Then run pcoc_det on your alignment
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3

## Then run pcoc_det on your alignment
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --gamma

## if you want to have a nice plot of the sites detected as convergent you can use the option --plot with or without the option --reorder
## With the --reorder option:
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --plot --reorder
## Without the --reorder option:
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --plot
## if you want to have a nice plot of your complete alignment you can use --plot_complete_ali
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --plot_complete_ali --plot --reorder
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --plot_complete_ali --plot --reorder --inv_gamma
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --plot_complete_ali --plot --reorder --gamma


## another example with several convergent branches, with OC applied on branches 11, 6, and 15, and PC applied on all of 11,9,10,6,15:
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 11,9,10/6/15 --plot
