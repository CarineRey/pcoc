# PCOC-TOOLKIT

[[etc/PCOC_logo.pdf || logo]]

# Tools:

The pcoc-toolkit is composed of 3 tools:

* pcoc_sim.py: [(usage)](#pcoc_sim.py-usage)

simulates an alignment and a scenario of convergent evolution according to the PCOC model for each tree of a directory.


* pcoc_det.py: [(usage)](#pcoc_det.py-usage)

detects sites of an alignment under convergent evolution according to the PCOC model.
It requires an input alignment, an input tree and the number of nodes where the convergent transitions occurred.



* pcoc_num_tree.py: [(usage)](#pcoc_num_tree.py-usage)

allows a user to number the nodes of an input tree so as to set the nodes where convergent transitions occurred in pcoc_sim.py and pcoc_det.py



# Installation and Usage:

docker is the easiest way to use the pcoc-toolkit locally.

For example to use pcoc_sim.py, without installation you just have to type:

 * If you want to have user permission on output files:

```
# run as user (to have user permission on output files)
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc [tool].py [options]
```

 * If you want to have root permission on output files:

```
# run as root (to have root permission on output files)
docker run -v $PWD:$PWD --rm carinerey/pcoc [tool].py [options]
```

To get all options:

```docker run --rm carinerey/pcoc```


or only for a specific tool:

```docker run --rm carinerey/pcoc [tool].py -h```


Docker update after an update of pcoc-toolkit:

```docker pull carinerey/pcoc```


# Example:


## Simulation:

```{sh}
mkdir -p example
cd example
wget https://raw.githubusercontent.com/Ensembl/ensembl-compara/release/89/scripts/pipeline/species_tree.39mammals.branch_len.nw -P trees

docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_sim.py -td $PWD/trees -o $PWD/output_pcoc --pcoc --topo --ident
```

Results will be found in $PWD/output_pcoc/RUN_yyyymmdd_hhmmss/

You will find:

  * Run_metadata.tsv: a tabular file with metadata about this run
  * Tree_metadata.tsv: a tabular file with metadata about all input trees
  * Tree_#: for each input tree a directory containing results for this input tree:

    * Metadata_Scenarios.tsv: a tabular file with metadata about this input tree
    * BenchmarkResults.tsv: a tabular file which summarizes the number of TP/FP/TN/FN/... for each used methods (--pcoc/--ident/--topo), for different thresholds, and for different couples of profiles used to simulate the data ( see --nb_sampled_couple, min_dist_CAT)
    * nw_trees:  a directory containing for each scenario the newick formated tree
    * sequences:  a directory containing for each scenario, all simulated alignments (the option --no_clean_seqs must be used, otherwise the sequences are discarded)
    * likelihood_summaries:  a directory containing for each scenario, likelihood summaries (the option --get_likelihood_summaries must be used)
    * plot_tree_ali:  a directory containing for each scenario, a plot with the tree, the alignment and the posterior probabilities per site for each used method (the option --plot_ali must be used)

## Detection:

```{sh}
mkdir -p example
cd example

## Get the data
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

## if you want to have a nice plot of the sites detected as convergent you can use the option --plot with or without the option --reorder
## With the --reorder option:
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --plot --reorder
## Without the --reorder option:
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --plot
## if you want to have a nice plot of your complete alignment you can use --plot_complete_ali
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 0/3 --plot_complete_ali --plot --reorder


## another example with several convergent branches, with OC applied on branches 11, 6, and 15, and PC applied on all of 11,9,10,6,15:
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD carinerey/pcoc pcoc_det.py -t $tree -aa $ali -o $PWD/output_pcoc_det -m 11,9,10/6/15 --plot


```
Results will be found in the folder $PWD/output_pcoc_det/RUN_yyyymmdd_hhmmss/ where "yyyymmdd_hhmmss" corresponds to the time of the analysis.

You will find:

  * Run_metadata.tsv: a tabular file with metadata about this run
  * pcoc_det.log: a log file
  * ali.results.tsv: a tabular file with posterior probabilities for each position according the PCOC models
  * ali.filtered_results.tsv: the same file but filtered for positions with a posterior probability superior to the threshold of one of the PCOC models (see -f/-f_pcoc/-f_pc/-f_oc)
  * fasta: a directory containing filtered position fasta (only if --no_cleanup_fasta has been used)
  * ali_plot_complete_ali.pdf: a plot of the tree, the complete alignment and the posterior probability according to PCOC models (only if --plot_complete_ali has been used)
  * ali_plot_filtered_ali_PCOC.pdf: the same plot as ali_plot_complete_ali.pdf but filtered for the positions with a posterior probability according to the PCOC model above its threshold (only if --plot)
  * ali_plot_filtered_ali_PC.pdf: the same plot as ali_plot_filtered_ali_PCOC.pdf but for the PC model (only if --plot)
  * ali_plot_filtered_ali_OC.pdf: the same plot as ali_plot_filtered_ali_PCOC.pdf but for the OC model (only if --plot)
  * ali_plot_filtered_ali_union.pdf: the same plot as ali_plot_filtered_ali_PCOC.pdf but filtered for the positions with a posterior probability according one of the PCOC models above its threshold (only if --plot)

# Options by tool:

* [pcoc_sim.py usage](#pcoc_sim.py-usage)
* [pcoc_det.py usage](#pcoc_det.py-usage)
* [pcoc_num_tree.py usage](#pcoc_num_tree.py-usage)


## pcoc_sim.py usage:

```
usage: pcoc_sim.py [-h] [--version] [-cpu CPU] -td TREE_DIR -o OUTPUT_DIR
                   [-n_sc INT] [-m "x/y,z/..."] [-c_max INT] [-c INT]
                   [-c_min INT] [-cr FLOAT] [-flg FLOAT] [-bl_new FLOAT]
                   [--ali_noise] [--bl_noise]
                   [-nb_sampled_couple NB_SAMPLED_COUPLE] [-n_sites N_SITES]
                   [-CATX_sim {10,60}] [-min_dist_CAT FLOAT] [--plot_ali]
                   [--get_likelihood_summaries] [--no_clean_seqs]
                   [-CATX_est {10,60}] [--pcoc] [--ident] [--topo]
                   [--plot_event_repartition] [--no_cleanup]
                   [-LD_LIB LD_LIBRARY_PATH] [--debug]

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -cpu CPU              Number of cpu to use. (default: 1)

required arguments:
  -td TREE_DIR, --tree_dir TREE_DIR
                        Directory name containing input trees.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory name. (Default output)

Convergent scenarios simulation options:
  -n_sc INT             Number of convergent scenarios draw up from each input
                        tree. (default: 1)
  -m "x/y,z/...", --manual_mode "x/y,z/..."
                        User defined convergent transition/branches.
                        Transition node must be the first number and
                        independent events must be separed by a "/". ex:
                        "1,2,3/67/55,56" (default: None)
  -c_max INT            Maximum number of potential convergent transitions
                        (=convergent events). First, this number of events are
                        draw up in the tree, then the wanted number of event
                        is randomly chosen among them. It is to avoid any bias
                        when you want to compared 2 and 7 events. (default: 8)
  -c INT                Number of convergent transitions (=convergent events).
                        If not defined, random between min events and max
                        events
  -c_min INT            Minimum number of transition (=convergent events).
                        (default: 2)
  -cr FLOAT             Maximum ratio between the number of convergent/non
                        convergent leaves. (default: No limits)
  -flg FLOAT            For each input tree, branch length multiplicator.
                        (default: no modification)
  -bl_new FLOAT         For each input tree, replace all branch lengths by
                        [FLOAT]. (default: no modification)
  --ali_noise           Add noisy events in the convergent scenario.
  --bl_noise            Add noise in the branch lengths of of tree for the
                        detection process.

Alignment simulation options:
  -nb_sampled_couple NB_SAMPLED_COUPLE
                        For each convergent scenario, number of simulated
                        alignment with different sampled couple of profiles
                        (Ancestral/Convergent). (default: 1)
  -n_sites N_SITES      Number of simulated sites per alignment. (default:
                        100)
  -CATX_sim {10,60}     Profile categories to simulate data (10->C10 or
                        60->C60). (default: 60)
  -min_dist_CAT FLOAT   Minimum distance between Ancestral and Convergent
                        profiles to simulate the alignment (default: no
                        limits)
  --plot_ali            For each couple of profiles, plot a summary of the
                        convergent scenario containing the tree and the
                        alignment.
  --get_likelihood_summaries
                        For each couple of profiles, write a summary of the
                        likelihoods per site.
  --no_clean_seqs       Do not cleanup the sequences after the run.

Benchmark options:
  -CATX_est {10,60}     Profile categories to estimate data (10->C10 or
                        60->C60). (default: 10)
  --pcoc                Use the PCOC model approach to detect site under
                        convergent evolution.
  --ident               Use the ancestral reconstruction approach to detect
                        site under convergent evolution.
  --topo                Use the topological approach to detect site under
                        convergent evolution.
  --plot_event_repartition
                        Plot chosen random convergent events repartition for
                        each input tree.

Options:
  --no_cleanup          Do not cleanup the working directory after the run.
  -LD_LIB LD_LIBRARY_PATH
                        Redefine the LD_LIBRARY_PATH env variable, bppsuite
                        library must be present in the $PATH and in the
                        LD_LIBRARY_PATH
  --debug               debug mode
```

## pcoc_det.py usage:

```
usage: pcoc_det.py [-h] [--version] -t TREE -aa ALI -m "x/y,z/..." -o
                   OUTPUT_DIR [-CATX_est {10,60}] [-f FILTER_T]
                   [-f_pcoc FILTER_T_PCOC] [-f_pc FILTER_T_PC]
                   [-f_oc FILTER_T_OC] [--plot] [-ph PH] [--reorder] [--gamma]
                   [--max_gap_allowed MAX_GAP_ALLOWED] [--plot_complete_ali]
                   [--no_cleanup_fasta] [--no_cleanup]
                   [-LD_LIB LD_LIBRARY_PATH] [--debug]

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

Required arguments:
  -t TREE, --tree TREE  Input tree filename
  -aa ALI, --ali ALI    Input amino-acid aligment filename
  -m "x/y,z/...", --manual_mode "x/y,z/..."
                        User defined convergent transition/branches.
                        Transition node must be the first number and
                        independent events must be separed by a "/". ex:
                        "1,2,3/67/55,56" (default: None)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory name

Options:
  -CATX_est {10,60}     Profile categorie to estimate data (10->C10 or
                        60->C60). (default: 10)
  -f FILTER_T, --filter_t FILTER_T
                        ALL model: Posterior probability threshold to put
                        result in "filtered" results. (default: 0.99)
  -f_pcoc FILTER_T_PCOC, --filter_t_pcoc FILTER_T_PCOC
                        PCOC model: Posterior probability threshold to put
                        result in "filtered" results. If = -1, take the value
                        of -f, if > 1, discard this model. (default: -1)
  -f_pc FILTER_T_PC, --filter_t_pc FILTER_T_PC
                        PC model: Posterior probability threshold to put
                        result in "filtered" results. If = -1, take the value
                        of -f, if > 1, discard this model.(default: -1)
  -f_oc FILTER_T_OC, --filter_t_oc FILTER_T_OC
                        OC model: Posterior probability threshold to put
                        result in "filtered" results. If = -1, take the value
                        of -f, if > 1, discard this model.(default: -1)
  --plot                Plot the tree and the filtered sites of the alignment
                        with their corresponding score.
  -ph PH                Add these positions in the filtered position and
                        highlight them with a star in the plot
  --reorder             reorder the filtered plot by score.categories (>=
                        0.99, >=0.9, >= 0.8, < 0.8)
  --gamma               Use rate_distribution=Gamma(n=4) instead of Constant()
  --max_gap_allowed MAX_GAP_ALLOWED
                        max_gap_allowed, must be between 0 and 100
  --plot_complete_ali   Plot the tree and which site of alignment with its
                        corresponding score. (Can take time to be open)
  --no_cleanup_fasta    Do not cleanup the fasta directory after the run.
  --no_cleanup          Do not cleanup the working directory after the run.
  -LD_LIB LD_LIBRARY_PATH
                        Redefine the LD_LIBRARY_PATH env variable, bppsuite
                        library must be present in the $PATH and in the
                        LD_LIBRARY_PATH
  --debug               debug mode
```

## pcoc_num_tree.py usage:

```
usage: pcoc_num_tree.py [-h] [--version] -t TREE -o OUTPUT [-m "x/y,z/..."]

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

Required arguments:
  -t TREE, --tree TREE  input tree name
  -o OUTPUT, --output OUTPUT
                        Output file (pdf)

Options:
  -m "x/y,z/...", --manual_mode_test "x/y,z/..."
                        Highlight, in the output plot, user defined convergent
                        transition/branches. Transition node must be the first
                        number and independent events must be separed by a
                        "/". ex: "1,2,3/67/55,56" (default: None)
```
