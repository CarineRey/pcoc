#! /bin/bash


set -e


mkdir -p test
cd test
ROOT_DIR=$PWD
TREES_DIR=$ROOT_DIR/../etc/data/tree_ensembl
OUTPUT_DIR=$ROOT_DIR/out_test/
rm -r $OUTPUT_DIR || ls $ROOT_DIR
mkdir -p  $OUTPUT_DIR


meth="--pcoc"
meth="--pcoc --ident --topo"
debug="--no_cleanup --debug"
debug=""
# test1: TOPO IDENT PCOC NO_NOISE
#docker run --rm -i -e LOCAL_USER_ID=`id -u $USER` -v $ROOT_DIR:$ROOT_DIR -v $TREES_DIR:$TREES_DIR  carinerey/pcoc  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test1 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --plot_ali"

#test2: TOPO IDENT PCOC BL_NOISE
#docker run --rm -i -e LOCAL_USER_ID=`id -u $USER` -v $ROOT_DIR:$ROOT_DIR -v $TREES_DIR:$TREES_DIR  carinerey/pcoc  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test2 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1 $meth $debug --bl_noise" 

#test3: TOPO IDENT PCOC ALI_NOISE
#docker run --rm -i -e LOCAL_USER_ID=`id -u $USER` -v $ROOT_DIR:$ROOT_DIR -v $TREES_DIR:$TREES_DIR  carinerey/pcoc  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test3 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --ali_noise --plot_ali"

#test4: TOPO IDENT PCOC EV_NOISE
#docker run --rm -i -e LOCAL_USER_ID=`id -u $USER` -v $ROOT_DIR:$ROOT_DIR -v $TREES_DIR:$TREES_DIR  carinerey/pcoc  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test4 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --ev_noise +1 --plot_ali -min_dist_CAT 0.8"
#docker run --rm -i -e LOCAL_USER_ID=`id -u $USER` -v $ROOT_DIR:$ROOT_DIR -v $TREES_DIR:$TREES_DIR  carinerey/pcoc  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test4 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --ev_noise -1 --plot_ali -min_dist_CAT 0.8"
#docker run --rm -i -e LOCAL_USER_ID=`id -u $USER` -v $ROOT_DIR:$ROOT_DIR -v $TREES_DIR:$TREES_DIR  carinerey/pcoc  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test4 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --ev_noise =1 --plot_ali -min_dist_CAT 0.8"

#test5: TOPO IDENT PCOC ROOT_NOISE
docker run --rm -i -e LOCAL_USER_ID=`id -u $USER` -v $ROOT_DIR:$ROOT_DIR -v $TREES_DIR:$TREES_DIR  carinerey/pcoc  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test5 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --root_noise ll --plot_ali -min_dist_CAT 0.8"
docker run --rm -i -e LOCAL_USER_ID=`id -u $USER` -v $ROOT_DIR:$ROOT_DIR -v $TREES_DIR:$TREES_DIR  carinerey/pcoc  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test5 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --root_noise rl --plot_ali -min_dist_CAT 0.8"
docker run --rm -i -e LOCAL_USER_ID=`id -u $USER` -v $ROOT_DIR:$ROOT_DIR -v $TREES_DIR:$TREES_DIR  carinerey/pcoc  bash -c "pcoc_sim.py -td $TREES_DIR -o $OUTPUT_DIR/test5 -n_sc 1 -nb_sampled_couple 1 -n_sites 100 -c 5 -c_max 7 -cpu 1  $meth $debug --root_noise rl --ev_noise +1 --plot_ali -min_dist_CAT 0.8"


echo TOPO IDENT PCOC NO_NOISE
#cat $OUTPUT_DIR/test1/RUN*/Tree_1/BenchmarkResults.tsv | cut -f 4,8,10-12,13,9 | grep -e Threshold -e "0.9" -e "NA"

echo TOPO IDENT PCOC BL_NOISE
#cat $OUTPUT_DIR/test2/RUN*/Tree_1/BenchmarkResults.tsv | cut -f 4,8,10-12,13,9 | grep -e Threshold -e "0.9" -e "NA"

echo TOPO IDENT PCOC ALI_NOISE
#cat $OUTPUT_DIR/test3/RUN*/Tree_1/BenchmarkResults.tsv | cut -f 4,8,10-12,13,9 | grep -e Threshold -e "0.9" -e "NA"

echo TOPO IDENT PCOC EV_NOISE
#cat $OUTPUT_DIR/test4/RUN*/Tree_1/BenchmarkResults.tsv | cut -f 4,8,10-12,13,9 | grep -e Threshold -e "0.9" -e "NA"

echo TOPO IDENT PCOC ROOT_NOISE
cat $OUTPUT_DIR/test5/RUN*/Tree_1/BenchmarkResults.tsv | cut -f 4,8,10-12,13,9 | grep -e Threshold -e "0.9" -e "NA"

