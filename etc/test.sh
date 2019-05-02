#! /bin/bash

set -e

mkdir -p test
cd test
ROOT_DIR=$PWD
TREES_DIR=$ROOT_DIR/../etc/data/tree_ensembl
DATA_DIR=$ROOT_DIR/../example/data/
OUTPUT_DIR=$ROOT_DIR/out_test/
rm -rf $OUTPUT_DIR || ls $ROOT_DIR
mkdir -p  $OUTPUT_DIR


meth="--pcoc"
meth="--pcoc --ident --topo"
debug="--no_cleanup --debug"
debug=""


DOCKER_CMD="docker run --rm -i -e LOCAL_USER_ID=`id -u $USER` -v $ROOT_DIR:$ROOT_DIR -v $TREES_DIR:$TREES_DIR -v $DATA_DIR:$DATA_DIR  carinerey/pcoc:v1.1.0"

<<<<<<< HEAD
=======

# test: P CONV
$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test1 -n_sc 1 -nb_sampled_couple 1 -sim_profiles $DATA_DIR/aa_order_checking.csv -est_profiles C10  -n_sites 10 -c 5 -c_max 7 -cpu 1  $meth $debug  --no_clean_seqs -p_conv 0.8"
#$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test1 -n_sc 1 -nb_sampled_couple 1 -est_profiles $DATA_DIR/aa_per_properties.csv -n_sites 10 -c 5 -c_max 7 -cpu 1  $meth $debug --no_clean_seqs -p_conv 0.8"


>>>>>>> 8117b8e... reduce test time
# test7: PROFILE INPUT
#$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test1 -n_sc 1 -nb_sampled_couple 1 -sim_profiles $DATA_DIR/C10_aa_frequencies.csv -n_sites 1 -c 5 -c_max 7 -cpu 1  $meth $debug --plot_ali -min_dist_CAT 0 --no_clean_seqs"
#$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test1 -n_sc 1 -nb_sampled_couple 1 -sim_profiles $DATA_DIR/aa_per_properties.csv -n_sites 1 -c 5 -c_max 7 -cpu 1  $meth $debug --plot_ali -min_dist_CAT 0 --no_clean_seqs"

# test1: TOPO IDENT PCOC NO_NOISE
$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test1 -n_sc 1 -nb_sampled_couple 1 -sim_profiles C10 -n_sites 1 -c 5 -c_max 7 -cpu 1  $meth $debug --plot_ali -min_dist_CAT 1 --no_clean_seqs"

#test2: TOPO IDENT PCOC BL_NOISE
$DOCKER_CMD bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test2 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1 $meth $debug --bl_noise" 

#test3: TOPO IDENT PCOC ALI_NOISE
$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test3 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --ali_noise --plot_ali"

#test4: TOPO IDENT PCOC EV_NOISE
$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test4 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --ev_noise +1 --plot_ali -min_dist_CAT 0.8"
$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test4 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --ev_noise -1 --plot_ali -min_dist_CAT 0.8"
$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test4 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --ev_noise =1 --plot_ali -min_dist_CAT 0.8"

#test5: TOPO IDENT PCOC ROOT_NOISE
$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test5 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --root_noise ll --plot_ali -min_dist_CAT 0.8"
$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test5 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --root_noise rl --plot_ali -min_dist_CAT 0.8"
$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test5 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --root_noise rr --ev_noise +1 --plot_ali -min_dist_CAT 0.8"
$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test5 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --root_noise lr --bl_noise --ali_noise --ev_noise =1 --plot_ali -min_dist_CAT 0.8"

#test6: IDENT MANUAL
$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test6 -n_sc 1 -nb_sampled_couple 1 -n_sites 10 -c 2 -m 0/33 -cpu 1  --ident $debug "
$DOCKER_CMD  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test6 -n_sc 1 -nb_sampled_couple 1 -n_sites 10 -c 3 -m 0/33/47,45,46 -cpu 1  --ident $debug "

echo TOPO IDENT PCOC NO_NOISE
cat $OUTPUT_DIR/test1/RUN*/Tree_1/BenchmarkResults.tsv | cut -f 4,8,10-12,13,9 | grep -e Threshold -e "0.9" -e "NA"

echo TOPO IDENT PCOC BL_NOISE
cat $OUTPUT_DIR/test2/RUN*/Tree_1/BenchmarkResults.tsv | cut -f 4,8,10-12,13,9 | grep -e Threshold -e "0.9" -e "NA"

echo TOPO IDENT PCOC ALI_NOISE
cat $OUTPUT_DIR/test3/RUN*/Tree_1/BenchmarkResults.tsv | cut -f 4,8,10-12,13,9 | grep -e Threshold -e "0.9" -e "NA"

echo TOPO IDENT PCOC EV_NOISE
cat $OUTPUT_DIR/test4/RUN*/Tree_1/BenchmarkResults.tsv | cut -f 4,8,10-12,13,9 | grep -e Threshold -e "0.9" -e "NA"

echo TOPO IDENT PCOC ROOT_NOISE
cat $OUTPUT_DIR/test5/RUN*/Tree_1/BenchmarkResults.tsv | cut -f 4,8,10-12,13,9 | grep -e Threshold -e "0.9" -e "NA"

echo IDENT MANUAL
cat $OUTPUT_DIR/test6/RUN*/Tree_1/BenchmarkResults.tsv | cut -f 4,8,10-12,13,9 | grep -e Threshold -e "0.9" -e "NA"

