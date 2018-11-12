# PCOC

[![Build Status](https://travis-ci.org/CarineRey/pcoc.svg?branch=master)](https://travis-ci.org/CarineRey/pcoc)

(For full explanation of the pipeline, see the [**PCOC** paper](https://doi.org/10.1101/247296)).

**PCOC** is a convergent substitution detection tool with two main facets:

* [I.  **PCOC** allows users to detect convergence in an empirical dataset (species tree and amino-acid alignment).](#i-using-pcoc-to-detect-convergence-in-an-empirical-dataset)
     * Conditions to use the **PCOC** *detection* pipeline:
         * you already have an empirical dataset, *i.e.* (at least) an amino-acid alignment and a corresponding gene tree


* [II.  **PCOC** allows users to simulate a dataset and compare the power of methods that detect convergence in it.](https://github.com/CarineRey/pcoc/wiki#ii-using-pcoc-to-simulate-a-dataset-and-compare-convergent-detection-methods-in-it)
     * Conditions to use the **PCOC** *simulation* pipeline:

        * you want to start a project on convergence and you wonder how many samples are necessary to have sufficient power to detect convergent sites. With the **PCOC** *simulation* pipeline, you can simulate if detection of convergence is possible with a dataset of a given size.

        * you want to test the statistical power of **PCOC** on your dataset.

        * you want to test the underlying model of **PCOC** and compare it to other algorithms used to detect convergence


***

**Table of contents:**

   * [I. Using <strong>PCOC</strong> to detect convergence in an empirical dataset](#i-using-pcoc-to-detect-convergence-in-an-empirical-dataset)
      * [1. Prerequisite on input data](#1-prerequisite-on-input-data)
      * [2. <strong>PCOC</strong> Installation](#2-pcoc-installation)
      * [3. Usage](#3-usage)
         * [A. Prepare tree for detection analysis: set putative convergent leaves/nodes](#a-prepare-tree-for-detection-analysis-set-putative-convergent-leavesnodes)
         * [B. Run detection analysis](#b-run-detection-analysis)
      * [4. Output files](#4-output-files)
      * [5. Using simulation to calculate the statistical power of PCOC on your dataset](#5-using-simulation-to-calculate-the-statistical-power-of-pcoc-on-your-dataset)
   * [II. Using <strong>PCOC</strong> to simulate a dataset and compare methods for detecting convergent evolution on it](#ii-using-pcoc-to-simulate-a-dataset-and-compare-methods-for-detecting-convergent-evolution-on-it)
      * [1. Usage](#1-usage)
         * [A. Define the convergent scenarios](#a-define-the-convergent-scenarios)
         * [B. Simulate the data according to the convergent scenario](#b-simulate-the-data-according-to-the-convergent-scenario)
         * [C. Test for convergence in the simulated data using different methods (including <strong>PCOC</strong>)](#c-test-for-convergence-in-the-simulated-data-using-different-methods-including-pcoc)
      * [2. Output files](#2-output-files)
      * [3. Some examples taken from the PCOC paper](#3-some-examples-taken-from-the-pcoc-paper)


***

# I. Using **PCOC** to detect convergence in an empirical dataset

## 1. Prerequisite on input data

The user needs:

* an amino-acid alignment
* a gene tree with values for branch lengths (these values are used by the program)

Once you have these data, you can install **PCOC**.

If you don't have your data yet, you can use the dataset from [(Besnard *et al.*, 2009)](https://doi.org/10.1093/molbev/msp103) used in the PCOC paper and courteously provided by the authors:

> Besnard, G., Muasya, A. M., Russier, F., Roalson, E. H., Salamin, N., and Christin, P.-A. 2009. Phylogenomics of C4 Photosynthesis in Sedges (Cyperaceae): Multiple Appearances and Genetic Convergence. Molecular Biology and Evolution, 26(8): 1909–1919, https://doi.org/10.1093/molbev/msp103

```{sh}
#get the PEPC protein alignment for sedges (plant species at C3/C4 transition), and rename this alignment with shorter name
wget  https://raw.githubusercontent.com/CarineRey/pcoc/master/data/det/cyp_coding.aa.coor_mays.fa
mv cyp_coding.aa.coor_mays.fa ali.fa

#get the corresponding tree, rename it and move it to a new folder
mkdir -p tree_dir
wget https://raw.githubusercontent.com/CarineRey/pcoc/master/data/det/cyp_coding.phy_phyml_tree.txt -P tree_dir
mv tree_dir/cyp_coding.phy_phyml_tree.txt tree_dir/tree.nw
```

**Of note**: the gene tree needs to have the **exact** same labels (same number of labels, same names) as the alignment.

## 2. **PCOC** Installation

### Docker

PCOC relies on the `Bpp` suite (https://github.com/BioPP/bppsuite) which is in constant development. In order to avoid installation and compilation of the latest version on your machine (this can take a long time), we suggest you to use `Docker`. Of course, you can also use **PCOC** without Docker but `Docker` is the easiest way to use the **PCOC** toolkit locally. `Docker` will create a local environment on your computer that will contain **PCOC** dependencies  (`Bpp` and `python` with some modules (`ete3` and `Biopython`)).They will all be packaged within the **PCOC** `Docker` image.

If you don't have Docker on your machine, get it [here](https://docs.docker.com/engine/installation/) first. (*Be aware that installation might differ if you're a Linux, a Mac or a Windows user.*)

To download or update the docker image you have to type the following command line:

```{sh}
docker pull carinerey/pcoc
```
*Note:* Mac users may first need to start Docker using the Docker application (if it's not running in background). If not, you could get the error message `Cannot connect to the Docker daemon at unix:///var/run/docker.sock. Is the docker daemon running?`.

Then, to use **PCOC** through a docker container (containing the **PCOC**  image), you have to type:

```{sh}
# run as user (to have user permission on output files) ##NOT WORKING
docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD -e CWD=$PWD carinerey/pcoc [SOME PCOC TOOL] [SOME PCOC OPTIONS]
```

What does this command line mean?

* `docker run [container options] [image] [cmd]`  is the command to launch a `[cmd]` in a docker container with the specified `[container options]` containing a given `[image]`.
* `` -e LOCAL_USER_ID=`id -u $USER` ``: allows you to get read and write permissions on the ouput files (if you remove this part of the command line, files will belong to the root user).
* `-v $PWD:$PWD`: allows sharing the current directory between your computer and the container.
* `--rm`: allows automaticlly removing the container when the job is finished (this limits disk space usage)
* `-e CWD=$PWD` : allows specifying the current working directory and use relative path from this directory (here `$PWD`)
*  `carinerey/pcoc ` is the name of the docker image hosted on DockerHub
*  `[SOME PCOC TOOL] [SOME PCOC OPTIONS] ` is the command run in the container (see below)

Below, we will shorten this command line by passing all the Docker options into an environment variable.

```{sh}
CMD_PCOC_DOCKER="docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD -e CWD=$PWD carinerey/pcoc"
```
Thus, to run **PCOC**, simply type:

```{sh}
$CMD_PCOC_DOCKER [SOME PCOC TOOL] [SOME PCOC OPTIONS]
```

### Singularity

**PCOC** can also be run under [Singularity](https://www.sylabs.io/), a container system designed specifically for use on shared compute clusters.  By default, Singularity containers run as the current user, not root, so the environment settings and current directory mappings are not necessary.

Unlike docker, Singularity "pull" will create an image in the current working directory that is referenced directly.

```{sh}
singularity pull docker://carinerey/pcoc
singularity run ./pcoc.simg [SOME PCOC TOOL] [SOME PCOC OPTIONS]
```

## 3. Usage


### A. Prepare tree for detection analysis: set putative convergent leaves/nodes

Based on biological relevance of the studied taxon, each user determines manually where the convergent transitions occurred, and which nodes are assumed to be in the convergent phenotypic state. This is done with the help of the `pcoc_num_tree.py` script.

 *Note that PCOC aims to detect sites with convergent substitutions in the case of a given set of transitions but not to identify the set of convergent transitions.*

#### i. The 1st step is to "number" the tree: each node in the input tree is labeled and given a number.


```{sh}
$CMD_PCOC_DOCKER pcoc_num_tree.py -t tree_dir/tree.nw -o num_tree.pdf ##takes a few seconds; replace file name of tree.nw with your input tree
```

#### ii. In the 2nd step, you use the numbered tree (output of the 1st step) to write a string corresponding to the convergent scenario you are interested in.


The string looks like `1,2,3/55,56/67`. The syntax is defined as such:

* a clade of nodes corresponding to a same convergent transition is written as a list of node IDs separated by a `,`. **The node ID corresponding to the initial convergent transition must be the first number** (here, `1`, `55` and `67`).

* each independent set of nodes corresponding to independent convergent clades are separated by a `/`.


In the example dataset used in the **PCOC** paper, species with the convergent phenotypes are (see the figure in num_tree.pdf for the node IDs):

| 1st convergent clade                                  | 2nd convergent clade          | 3rd convergent clade       | 4th convergent clade        | 5th convergent clade          |
|--------------------------------------------------------|-----------------------------------|----------------------------------|----------------------------------|----------------------------------|
| Ele.bald (node 0)                                      | Bulbostyl (38)                    | Fimb.li2 (45)                    | Cyp.capi (78)                    | Rhy.rubr (133)                   |
| Ele.bal2 (node 1)                                      | all internal nodes in this group  | Fimb.di2 (46)                    | Volkiell (79)                    | Rhy.glob (134)                   |
| Ele.bal4 (node 3)                                      |                                   | Fimb.fer (47)                    | Cyp.ust2 (80)                    | Rhy.glo2 (135)                   |
| Ele.vivi (node 5)                                      |                                   | all internal nodes in this group | Remirea  (81)                    | all internal nodes in this group |
| Ele.vivA (node 6)                                      |                                   |                                  | Cyp.iria (82)                    |                                  |
| all internal nodes in this group (nodes 2, 4, 7 and 8) |                                   |                                  | Killinga (83)                    |                                  |
|                                                        |                                   |                                  | Pycreus  (84)                    |                                  |
|                                                        |                                   |                                  | Cyp.long (86)                    |                                  |
|                                                        |                                   |                                  | Cyp.rotu (87)                    |                                  |
|                                                        |                                   |                                  | Cyp.papy (89)                    |                                  |
|                                                        |                                   |                                  | Cyp.ustu (90)                    |                                  |
|                                                        |                                   |                                  | all internal nodes in this group |                                  |

so the string corresponding to this scenario is:

`8,0,1,2,3,4,5,6,7/49,45,46,47,48/38/98,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97/137,133,134,135,136`.

(Remember that the common ancestor of each clade with a convergent transition has to be the first in the list.)

#### iii. (Optional) The 3rd step is to visualize the final tree with colored labeled branches (to ascertain that no branch has been forgotten during 2nd step)


```{sh}
# this takes a few seconds
scenario="8,0,1,2,3,4,5,6,7/49,45,46,47,48/38/98,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97/137,133,134,135,136"
$CMD_PCOC_DOCKER pcoc_num_tree.py -t tree_dir/tree.nw -o num_tree.pdf -m $scenario
```

### B. Run detection analysis

This is done by the `pcoc_det.py` script. Please see the man page for the full list of options: `$CMD_PCOC_DOCKER pcoc_det.py -h`

#### Minimal usage:

To run the detection tool of the PCOC toolkit, you need the amino-acid alignment of the studied gene, the "species" tree with values for branch lengths, and the string corresponding to the convergent scenario (obtained from the previous section).

With the option `-f`, you can fix the threshold at which a position is retained if its posterior probability is above this threshold for one of the 3 implemented models (PCOC, PC or OC).

```{sh}
# here we set the threshold to 0.8
# this takes one minute maximum
$CMD_PCOC_DOCKER pcoc_det.py -t tree_dir/tree.nw -aa ali.fa -o output_pcoc_det -m $scenario -f 0.8
```

## 4. Output files

Results are in the folder `output_pcoc_det/RUN_yyyymmdd_hhmmss/` where `yyyymmdd_hhmmss` corresponds to the date and time of the analysis.

You will find:

  * `ali.filtered_results.tsv`: a tabular file with posterior probabilities for each position with a posterior probability above the threshold of one of the PCOC models (PCOC, PC and OC; see -f/-f_pcoc/-f_pc/-f_oc)
  * `ali.results.tsv`: a tabular file with posterior probabilities for all positions according the **PCOC** models (PCOC, PC and OC)
  * `Run_metadata.tsv`: a tabular file with metadata about this run
  * `pcoc_det.log`: a log file

*Of note*: if you want to get the alignment of filtered positions only, add the `--no_cleanup_fasta` option to your command. This will create a directory named `fasta` with these particular positions for each model (PCOC, PC and OC) and the union of them. You will most probably only use the output from the PCOC model, but the PC and OC outputs are also provided.

#### Graphical output options:

* You can get a graphical output for all positions using  `--plot_complete_ali --plot ` (default is pdf output, add `--svg` for svg output)

    ```{sh}
    $CMD_PCOC_DOCKER pcoc_det.py -t tree_dir/tree.nw -aa ali.fa -o output_pcoc_det -m $scenario --plot_complete_ali --plot
    ```
    You will additionally find in the output directory `ali_plot_complete.pdf`. This file contains a plot of the tree, the complete alignment and the posterior probability according to PCOC models (only if `--plot_complete_ali` has been used)

* You can get a graphical output for the filtered positions only using `--plot` and reorder them in function of their posterior probability with `--reorder`

    ```{sh}
    $CMD_PCOC_DOCKER pcoc_det.py -t  tree_dir/tree.nw -aa ali.fa -o output_pcoc_det -m $scenario --plot --reorder
    ```

    You will additionnaly find in the output directory:

    * `ali_plot_filtered_PCOC.pdf`: the same plot as `ali_plot_complete.pdf` but filtered for the positions with a posterior probability according to the PCOC model above `-f` threshold (only if `--plot`)
    * `ali_plot_filtered_PC.pdf`: the same plot as ali_plot_filtered_PCOC.pdf but for the PC model
    * `ali_plot_filtered_OC.pdf`: the same plot as ali_plot_filtered_ali_PCOC.pdf but for the OC model
    * `ali_plot_filtered_union.pdf`: the same plot with the union of all positions wit the different models

 * You can also add user defined sites in the output and highlight them

    Based on previous work, you can have a list of candidate sites that are susceptible to be under convergent evolution. If you want to visualize how these sites go through the **PCOC** *detection* pipeline, you can choose to specifically display them in the output with the option `-ph`.

    In the example dataset, this is the case for 16 sites (based on Besnard et al. (2009)).

    ```{sh}
    candidates="505,540,572,573,623,665,731,733,749,751,761,770,780,810,839,906"
    ```
    You can visualize them by typing:

    ```{sh}
    $CMD_PCOC_DOCKER pcoc_det.py -t  tree_dir/tree.nw -aa ali.fa -o  output_pcoc_det -m $scenario  -f 0.8 --plot --reorder -ph $candidates
    ```
    The sites with a `*` in the image output correspond to the candidate sites (set with `-ph`).

## 5. Using simulation to calculate the statistical power of PCOC on your dataset

Once you got the candidate convergent sites for your gene, you may want
to known the sensitivity (True Positive Rate; TPR) and the specificity
or maybe the opposite, the False Positive Rate (FPR) of PCOC in your own
dataset.

For that, you can use the **PCOC** simulation and benchmark part using
the same tree and same parameters that were used for detection.
It will perform simulations of a large number of sites with convergent
evolution, and of sites without convergent evolution and provide the
amount of true positives and false negatives.

```{sh}
$CMD_PCOC_DOCKER pcoc_sim.py -t  tree_dir/ -o  output_pcoc_det_sim -m $scenario -c 5  -n_sc 1 --pcoc -nb_sampled_couple 10 -n_sites 100
```
*(For more details about pcoc_sim.py options see the next section.)*

In brief, for 10 random couple of profiles (`-nb_sampled_couple 10 `), pcoc_sim.py will simulate two alignments of 100 sites (`-n_sites 100`):

* one with convergent evolution to get the TPR of convergent site detection,
* another alignement without convergent evolution to get the FPR.

Then, pcoc_sim.py will detect convergent evolution in both alignments.

For each couple of ancestral and convergent profiles, the convergent evolution scenario will contain the 5 convergent transitions (`-c 5`) defined in `$scenario`.

*You can add `--ident` or `--topo` to test the Topological and Identical methods (see **PCOC** paper for more details).*

In the output directory, you fill find a tabular file `Tree_1/BenchmarkResults.tsv`, which contains the number of FP and TP for each method used, for each couple of profiles and for different thresholds.

You have just to average all the lines for a given method and a given threshold to get your TPR and FPR.

For example, we have:

* Sensitivity = TPR = mean(TP_simulations) / n_sites_simulations
* 1 - Specificity = FPR = mean(FP_simulations) / n_sites_simulations

In our example, if we take a threshold equal to 0.99 and as we used n_sites_simulations=100, you should found something like :

* TPR = 0.9976
* FPR = 0.0004

Then, to get the expected number of FP in your dataset (FP_ali), you have to calculate the expected number of positive and negative sites (P_ali and F_ali).
If you take, as in the PCOC paper, an expectd proportion of 2% of convergent sites in your data, you will have:

* P_ali = 0.02 x n_sites_dataset = 0.02 x 458 = 9
* N_ali = 0.98 x n_sites_dataset = 0.98 x 458 = 449

*This figure (2%) is extracted from the paper of Thomas and Hahn (2015) where they found 140 truly convergent sites among 6000 sites.*


And so, you will have:

* TP_ali = Sensitivity * P_ali = 9 * 0.9976 = 8.98
* FP_ali = FPR * N_ali = 449 * 0.0004 = 0.18

Finally, you should find a very low False Discovery Rate (FDR):

* FDR_ali = FP_ali / (TP_ali + FP_ali) = 0.02
___


# II. Using **PCOC** to simulate a dataset and compare methods for detecting convergent evolution on it

*If you don't have **PCOC** yet, see [the installation procedure in the first section](#2-pcoc-installation)*

**PCOC** allows users to simulate a dataset and compare **PCOC** with other convergence detection tools (using the simulated dataset). This **PCOC** feature can be used in 2 main situations:

* you want to start a project on convergence and you wonder how many samples are necessary to have enough power to detect convergent sites. With the **PCOC** *simulation* pipeline, you can simulate if detection of convergence is possible with a dataset of a given size

* you want to test the underlying model of **PCOC** and compare it to other algorithms used to detect convergence

In both cases, the options will be exactly the same. The only thing that differs between these approaches is the definition of the convergence scenario.


(i) In the first situation, you provide a candidate convergence scenario (relevant to the biology of the studied species). In practice, you have a species tree with given convergent events and you want to test the power of convergent detection tools on this particular topology.

(ii) In the second situation, the convergence scenario will be randomly generated. You have one (or more) species tree(s) and you want to compare the power of **PCOC** with other convergence detection tools under random convergent scenarios on this (these) topologies.


## 1. Usage

The **PCOC** simulation and benchmark part is done by the `pcoc_sim.py` script.

*(see man page for full list of options:* `$CMD_PCOC_DOCKER pcoc_sim.py -h`)

> *If you don't have a species tree yet, you can use the dataset from [Besnard *et al.*, 2009](https://doi.org/10.1093/molbev/msp103) used in the **PCOC** paper and courteously provided by the authors:*
>```{sh}
>#get the corresponding tree, rename it and move it to a new folder
>mkdir -p tree_dir
>wget https://raw.githubusercontent.com/CarineRey/pcoc/master/data/det/cyp_coding.phy_phyml_tree.txt -P tree_dir
>mv tree_dir/cyp_coding.phy_phyml_tree.txt tree_dir/tree.nw
>```
> The other species trees used in the **PCOC** paper are available here: https://github.com/CarineRey/pcoc/tree/master/data/sim/tree_dir

A typical **PCOC** simulation and detection pipeline is composed of three main steps:
* A. Define the convergent scenario(s)
* B. Simulate the data accordingly (this will be different under different scenarios)
* C. Test for convergence in the simulated data using different methods (including **PCOC**)

Thus, **PCOC** command-lines will be divided into "blocks" that correspond to the different steps in this pipeline. Prior to each analysis, you will have to define the directory containing your tree(s) (`-td` option) and the output directory (`-o` option).

 ```{sh}
$CMD_PCOC_DOCKER pcoc_sim.py -td tree_dir -o output_pcoc_sim [REST OF OPTIONS]
```
Below, we'll go through all the three "steps"... step by step.


### A. Define the convergent scenarios

In both situations (user-defined or random scenario), **PCOC** will *randomly* choose ONE number of convergent events between a minimum (`c_min`) and a maximum number (`c_max`). If you use `-c`, the number of convergent events will be fixed to `-c`. (`c` can still be used in combination with `c_max` to pick `-c` convergent events among `-c_max` events.)

To run **PCOC** multiple times and test different scenarios, use the `-n_sc` option. For instance, to repeat the test 10 times:
 ```{sh}
 pcoc_sim.py [...]  -n_sc 10  [...]
```

 * If you have a given convergent scenario, you have to set the scenario (the so-called `$scenario` variable [used above](#iii-optional-the-3rd-step-is-to-visualize-the-final-tree-with-colored-labeled-branches-to-ascertain-that-no-branch-has-been-forgotten-during-2nd-step)) using the same syntax defined in the [previous section](#--in-the-2nd-step-you-use-the-numbered-tree-output-of-the-1st-step-to-write-a-string-corresponding-to-the-convergent-scenario-you-are-interested-in).

    The `c_max` will be fixed as the total number of convergent events in the given scenario and the default `c_min` is 2 (except if you set `--c_min` to other value > 2). *Note that* `c_min` *has to be above 2 because convergence has to happen in at least 2 clades. Likewise,* `c_max` *has to be below the species number.*

    For instance, in the Besnard et al. (2009) dataset, `c_max` = 5. If you want to force **PCOC** to test these 5 events, use `-c 5` (If you set `c_max 5`, **PCOC** will randomly choose a number of events between 2 and 5 among the 5 possible).
     ```{sh}
     $scenario="8,0,1,2,3,4,5,6,7/49,45,46,47,48/38/98,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97/137,133,134,135,136"
     pcoc_sim.py [...]  -m $scenario  [...]
    ```

* If you want to test random convergent scenarios:


    You have to fix the `c_max` and the `c_min` of convergent event and **PCOC** will randomly place them in the tree. ( You can use the `--plot_event_repartition` option to see the distribution in the tree).

    *Note: To be fair when you compare 2 and 7 events in a topology, PCOC will first draw up the maximum number of convergent events and then randomly choose a number between those two (between the* `c_min` *and the* `c_max` *– except if you use* `-c`*) . In fact, it is more simple to place 2 events closest to the root of the tree than 7.*

    ```{sh}
    pcoc_sim.py [...]  -c_min 2 -c_max 7 [...] #this is overwritten by the -m option
    ```


**PCOC** can not modify the species tree topology but can modify the branch lengths. You can:

* multiply all the branch lengths by a factor X using `--flg X`
* change all the branch length by a same value X using `--bl_new X`


### B. Simulate the data according to the convergent scenario

Once **PCOC** has placed the convergent events on the species tree, it will use the tree to generate:

* an alignment corresponding to this convergent scenario

* another alignment corresponding to a non-convergent scenario (the same tree but without convergent evolution). *This allows for measuring the rate of false-positives.*

*Of note: the length of the alignments can be set with* `-n_sites`.

Sequence generation is based on the evolution of sequences from the root to the leaves of the tree according to amino-acid profiles (*i. e.* sets of amino-acid frequencies).
**PCOC** uses the C60 profiles (containing 60 amino-acid profiles) defined in [(Si Quang et al, 2009)](https://doi.org/10.1093/bioinformatics/btn445). Among those 60 profiles, 1 couple is picked and the first profile is qualified as "ancestor" and the second as "convergent".

In the convergent scenario, **PCOC** "traverses" the tree starting from the root with the ancestor profile. Once it reaches a node with a convergent transition, **PCOC** switches from the ancestor to the convergent profile. The other nodes keep the ancestral profile. *In the non-convergent scenario, the ancestor profile is applied to all nodes.*

*Alternatively, you can use the C10 profiles (containing 10 amino-acid profiles) instead of C60 using* `-CATX_sim 10`.

By default, only one couple of profiles is used for a scenario but you can use the `-nb_sampled_couple` to increase this number. This will generate multiple sets of independent alignments corresponding to the same convergent scenario.


To increase the realism of the simulation, you can add some noise in the simulated data. You can either:

* add noise in the alignment by using different profiles in the non-convergent part of the tree (see **PCOC** paper) using `--ali_noise`

* add noise in the scenario by modifying the branch lengths after simulating the data using `--bl_noise`. This will allow using a different tree for the simulation and the detection part by mimicking mis-estimated branch lengths.



### C. Test for convergence in the simulated data using different methods (including **PCOC**)

Finally, `pcoc_sim.py` can use and compare 3 methods to detect convergence in the simulated dataset (see **PCOC** paper for more details):
*  the **PCOC** models using `-pcoc`
*  a method detecting identical substitutions in convergent clades using `-ident` (first defined in [Zhang and Kumar (1997)](https://doi.org/10.1093/oxfordjournals.molbev.a025789))
*  a topological method using `-topo` (defined in [Parker et al. (2013)](https://www.nature.com/articles/nature12511))

For instance, to compare the **PCOC** method and the topological method only, use

```{sh}
pcoc_sim.py [...]  -pcoc -topo #do not include -ident
```

By default, C10 amino-acid profiles are used in the estimation part instead of the C60 amino-acid profiles in the simulation part. *Alternatively, you can use the C60 profiles instead of C10 using* `-CATX_est 60`. *Note that* `-CATX_est` *is different from* `-CATX_sim` *used in the simulation part.*

## 2. Output files

To have an example output, you can use the following command line adapted from the example (see below).

```{sh}
$CMD_PCOC_DOCKER pcoc_sim.py -td tree_dir -o output_pcoc_sim -n_sc 1 -nb_sampled_couple 1 -n_sites 50 -c_min 2 -c_max 7 -c 5 --pcoc --ident --topo
```

Results are in the folder `output_pcoc_sim/RUN_yyyymmdd_hhmmss/` where `yyyymmdd_hhmmss` corresponds to the date and time of the analysis.

You will find:

* 1. `run_metadata.tsv`: a tabular file containing a summary of the used options:

    >|                                              |                  |
    >|----------------------------------------------|------------------|
    >| RunID                                        | 20180111_174413  |
    >| Input trees                                  | tree.nw          |
    >| Run PCOC method                              | Yes              |
    >| Run identical method                         | Yes              |
    >| Run topological method                       | Yes              |
    >| Number of simulated sites                    | 50              |
    >|[...]                                         | [...]            |

* 2. `tree_metadata.tsv`: a correspondence table between the result directory name(s) `Tree_i` and the tree file name(s). (This contain also the execution time). *In this tutorial we have only used 1 tree but this table is useful if use more input trees (you would then have* `Tree_1`, `Tree_2`, *etc.).*

    >| Tree_# | Tree_filename | Execution time |
    >|--------|---------------|----------------|
    >| Tree_1 | tree.nw       | 32.9199180603  |
    >| [...]  |               |                |

* 3. `pcoc_sim.log`: a log file

* 4. a result directory `Tree_i` for each input tree containing:

    + a. `BenchmarkResults.tsv`: a tabular file containing metadata (columns "RunID", "InputTree", ..., "DistanceSimuCouple") and performance statistics ("FN", ... , "MCC") for each method ("Method") and for each threshold ("Threshold") used, for each simulated dataset. "NumberOfConvergentEvents" and "NumberOfSites" are other metadata related to the scenario and to the alignment.

        >| RunID           | InputTree | ScenarioID | SimuCoupleID | C1    | C2 | DistanceSimuCouple | Method      | Threshold | FN    | FP  | TN    | TP    | Sensitivity | Specificity | MCC    | NumberOfConvergentEvents | NumberOfSites|
        >|-----------------|-----------|------------|--------------|-------|----|--------------------|-------------|-----------|-------|-----|-------|-------|-------------|-------------|--------|--------------------------|---------------|
        >| 20180111_174413 | tree.nw   | Scenario_1 | A53_C44      | 53    | 44 | 0.424326507329     | PC          | 0.7       | 26.0  | 2.0 | 98.0  | 74.0  | 0.74        | 0.98        | 0.7416 | 5                        | 100           |
        >| 20180111_174413 | tree.nw   | Scenario_1 | A53_C44      | 53    | 44 | 0.424326507329     | PCOC        | 0.7       | 0.0   | 1.0 | 99.0  | 100.0 | 1.0         | 0.99        | 0.9900 | 5                        | 100           |
        >| 20180111_174413 | tree.nw   | Scenario_1 | A53_C44      | 53    | 44 | 0.424326507329     | OC          | 0.7       | 0.0   | 2.0 | 98.0  | 100.0 | 1.0         | 0.98        | 0.98019| 5                        | 100           |
        >| 20180111_174413 | tree.nw   | Scenario_1 | A53_C44      | 53    | 44 | 0.424326507329     | Topological | 0.7       | 16.0  | 3.0 | 97.0  | 84.0  | 0.84        | 0.97        | 0.81693| 5                        | 100           |
        >| [...]           | [...]     | [...]      | [...]        | [...] |[...]| [...]             | [...]       | [...]     | [...] |[...]| [...] | [...] | [...]       | [...]       | [...]  | [...]                    | [...]         |
        >| 20180111_174413 | tree.nw   | Scenario_1 | A53_C44      | 53    | 44 | 0.424326507329     | Identical   | NA        | 100.0 | 0.0 | 100.0 | 0.0   | 0.0         | 1.0         | 0.0    | 5                        | 100           |

    + b. `Metadata_Scenarios.tsv`: a tabular file containing (a lot of) metadata for each scenario (Number of convergent transitions,  Number of leaves, ... ). Column names should be explicit enough.

        >| ScenarioID | numberOfLeaves | Execution_time | InputTree | [...] |
        >|------------|----------------|----------------|-----------|-------|
        >| Scenario_1 | 41             | 32.8589060307  | tree.nw   | [...] |
        >| [...]        |                |                |           |       |

    + c. a directory `nw_trees` containing NHX formatted trees annotated with convergent transitions for each scenario. You will find 2 prefix:
        + `tree` corresponding to the input topology tree without modification
        + `tree_conv` corresponding to the topology tree with modifications link to the topological method

        Each node is labeled by a string starting with
        `&&NHX:ND=[NODE NUMBER]:T=True/False:C=True/False`. `T` refers to a node with convergent transition and `C` to nodes with a convergent profile. *In case you used* `--ali_noise` *option, a new field in node labels is added:*
        `Cz=False/[NUMBER OF NOISY PROFILE]` *(in the tree with prefix* `annotated_` *).*

## 3. Some examples taken from the PCOC paper

In this section, some command-lines used in the **PCOC** paper will be used as examples but option values have been changed to reduce the computational time.

The 2nd figure of the paper addresses the impact of:

1. the number of convergent events on the ability to detect convergence.

```{sh}
$CMD_PCOC_DOCKER pcoc_sim.py -td tree_dir -o output_pcoc_sim -n_sc 180 -nb_sampled_couple 1 -n_sites 50 -c_min 2 -c_max 7 -cpu 8 --pcoc --ident --topo -CATX_est 10 -CATX_sim 60
```

2. the impact of branch lengths on the ability to detect convergence.

```{sh}
bl_l="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1"
for bl in $bl_l
do
$CMD_PCOC_DOCKER pcoc_sim.py -td tree_dir -o output_pcoc_sim -n_sc 32 -nb_sampled_couple 1 -n_sites 50 -c_min 2 -c 5 -c_max 7 -cpu 8 --pcoc --ident --topo  -CATX_est 10 -CATX_sim 60 -bl_new $bl
done
```

Some figures of the supplementary data address the robustness of the **PCOC** models towards some noisy simulated dataset.

1. noise in the alignment (for example Supplementary figure S8)

```{sh}
$CMD_PCOC_DOCKER pcoc_sim.py -td tree_dir -o output_pcoc_sim -n_sc 180 -nb_sampled_couple 1 -n_sites 50 -c_min 2 -c_max 7 -cpu 8 --pcoc --ident --topo -CATX_est 10 -CATX_sim 60 --ali_noise
```
2. noise add in the branch length after simulate the alignment (for example Supplementary figure S9)

```{sh}
$CMD_PCOC_DOCKER pcoc_sim.py -td tree_dir -o output_pcoc_sim -n_sc 180 -nb_sampled_couple 1 -n_sites 50 -c_min 2 -c_max 7 -cpu 8 --pcoc --ident --topo -CATX_est 10 -CATX_sim 60 --bl_noise
```

___

Have fun with ![](etc/logo.png) !!
