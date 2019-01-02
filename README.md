[![Fold-U release](https://img.shields.io/badge/fold--u-v3.0-blue.svg)](https://github.com/meetU-MasterStudents/Fold_U/releases/tag/v3.0)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Python version](https://img.shields.io/badge/python-3-brightgreen.svg)
[![Documentation Status](https://readthedocs.org/projects/fold-u/badge/?version=latest)](https://fold-u.readthedocs.io/en/latest/?badge=latest)

<br>

# Fold U: A Protein Structure Prediction Program

<p align="center">
  <img width="400" src="img/logo_foldu.png" alt="logo_foldu"/>
</p>

Our program is the second step (downstream) of a **protein structure prediction project**. This step consists of threading a query sequence on different given templates.


This project is part of the **Meet-U 2018-2019** competition.
Meet-U is a collaborative pedagogical and research initiative between several Universities of Paris area. The course is intended to Master students (2nd year) in Bioinformatics. For more details, please refer to [http://www.meet-u.org/](http://www.meet-u.org/).

## Installation

### Clone the repository
```
git clone https://github.com/meetU-MasterStudents/Fold_U.git
cd Fold_U
```

### Requirements

1. A linux distribution.

2. Install the few required **python packages** :

```
pip3 install -r requirements.txt
# This command will install the following modules:
# docopt==0.6.2
# numpy==1.15.2
# biopython==1.72
# pandas==0.23.4
# schema==0.6.8
# tqdm==4.28.1
# matplotlib==2.2.2
# m2r # for Sphinx
```

3. **R** software environment. The full procedure is described [here](https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-18-04). R is used for the Machine Learning step. The required packages *boot*, *dplyr* and *readr* will be automatically installed if not already, you have nothing to do.
```
sudo apt-get install r-base
```

4. Download the **uniref90** database [here](https://www.uniprot.org/downloads), put it in `data/databases/uniref90` and format it with the following command line :
```
makeblastdb -in databases/uniref90.fasta -dbtype prot
```

5. Install **blast-legacy** with conda :
```
conda install -c bioconda blast-legacy
```

6. **MODELLER** is also required, and can be installed easily with Conda. You need to register to get a license key [here](https://salilab.org/modeller/registration.html), and follow instructions during installation to insert license key in the program.
```
conda install -c salilab modeller
```

7. If necessary, change the paths in the header of the following scripts : `bin/psipred.4.02/runpsipred` and `bin/psipred.4.02/runpsipred_single`

### Generate the PSSM files for each template

Run the following script before running `fold_u`:
```
./bin/salut_1.0/salut1.sh data/templates/ data/metafold.list
```

## Run the program

`fold_u` takes in input the protein sequence of the studied query (fasta format) and an uniref database.

It returns the mandatory outputs:
* `ranking.txt` file with 3 columns Rank | Template | Score
* The **top N pdb structures** (top 10 by default)

And additionnal outputs:
* `scores.csv`: this CSV formatted file contains the complete results with scores values.
* The MODELLER alignments.

### Toy example

```
./fold_u data/queries/Agglutinin/Agglutinin.fasta data/databases/uniref90/uniref90 -o results/Agglutinin
```

#### Get help

```
./fold_u -h

Usage:
    ./fold_u QUERY_FASTA UNIREF_DB [--nb_pdb NUM] [--nb_psiblast NUM] [--cpu NUM] [--dope FILE]
                                   [--metafold FILE] [--benchmark FILE] [--output PATH]

Arguments:
    QUERY_FASTA                           Path to the QUERY fasta sequence to predicted.
    UNIREF_DB                             Path to Uniref database.

Options:
    -h, --help                            Show this
    -p NUM, --nb_pdb NUM                  Number of pdb to create
                                          [default: 10]
    -t NUM, --nb_psiblast                 Round number for PSI-BLAST
                                          [default: 3]
    -c NUM, --cpu NUM                     Number of cpus to use for parallelisation. By default
                                          using all available (0).
                                          [default: 0]
    -d FILE, --dope FILE                  Path to the dope.par file
                                          [default: data/dope.par]
    -m FILE, --metafold FILE              Path to the metafold.list file
                                          [default: data/metafold.list]
    -b FILE, --benchmark FILE             Path to the benchmark.list file
                                          [default: data/benchmark.list]
    -o PATH, --output PATH                Path to the directory containing
                                          the result files (scores and pdb)
                                          [default: ./results]

```
## Run all the queries + Benchmarking

This program is also **benchmarked** using ROC style plots and **Top N** information to evaluate the power and the relevance of the different scores. The score results are generated for all queries. Each plot represents the cumulative sum of benchmarks encountered along the ranking (from rank 1 to rank 405) for each calculated scores. A top N results table is also generated showing the number of **"Family"**, **"Superfamily"** and **"Fold"** benchmarks found in the top N ranks.

We wrote a script that runs the `fold_u` program for each query if results are not yet generated. It returns a `results/plots` folder containing the generated plots and prints the **top N tables** in the terminal.

```
./scripts/benchmarking.py data/databases/uniref90/uniref90
```

#### Get help

```
./scripts/benchmarking.py -h

    Usage:
        ./script/benchmarking.py UNIREF_DB [--selected_score SCORE] [--cpu NUM] [--output PATH]

    Arguments:
        UNIREF_DB                             Path to Uniref database.        

    Options:
        -h, --help                            Show this
        -s SCORE, --selected_score SCORE      Score for which you wish to see the statistics:
                                              "alignment", "threading", "modeller",
                                              "secondary_structure", "solvent_access"
                                              or "sum_scores",
                                              or all of them at once: "all" [default: all]
        -c NUM, --cpu NUM                     Number of cpus to use for parallelisation. By default
                                              using all available (0).
                                              [default: 0]
        -o PATH, --output PATH                Path to the directory containing
                                              the result files (scores and plot)
                                              [default: ./results/plots]
```

### Machine Learning

The R script `script/machine_learning.R` uses logistic regression to find the best weights to apply to each type of score in order to optimize the benchmarking. That is to say it will learn the specificities of each scores according to the benchmarks (Fold, Superfamily and Family) in order to get the most information from each.

### Benchmark results

#### Top N tables

Results for the benchmark done with the merged program (upstream + downstream). Templates are ranked with the **combined (sum) score**.

```
Table summarizing the top N results.

         Family    Superfamily    Fold        Total

top 5    0/1       0/6            1/13        1/19
         0.0  %    0.0  %         7.7  %      5.3  %
----------------------------------------------------
top 10   0/1       2/6            1/13        3/19
         0.0  %    33.3 %         7.7  %      15.8 %
----------------------------------------------------
top 15   0/1       2/6            1/13        3/19
         0.0  %    33.3 %         7.7  %      15.8 %
----------------------------------------------------
top 20   0/1       2/6            1/13        3/19
         0.0  %    33.3 %         7.7  %      15.8 %
----------------------------------------------------
top 25   0/1       2/6            1/13        3/19
         0.0  %    33.3 %         7.7  %      15.8 %
----------------------------------------------------
top 50   0/1       3/6            3/13        6/19
         0.0  %    50.0 %         23.1 %      31.6 %
----------------------------------------------------
top 75   0/1       4/6            5/13        9/19
         0.0  %    66.7 %         38.5 %      47.4 %
----------------------------------------------------
top 100  0/1       5/6            6/13        11/19
         0.0  %    83.3 %         46.2 %      57.9 %
----------------------------------------------------
top 150  0/1       5/6            7/13        12/19
         0.0  %    83.3 %         53.8 %      63.2 %
----------------------------------------------------
top 200  0/1       5/6            7/13        12/19
         0.0  %    83.3 %         53.8 %      63.2 %
----------------------------------------------------
top 250  0/1       6/6            8/13        14/19
         0.0  %    100.0%         61.5 %      73.7 %
----------------------------------------------------
top 300  0/1       6/6            9/13        15/19
         0.0  %    100.0%         69.2 %      78.9 %
----------------------------------------------------
top 350  1/1       6/6            12/13       19/19
         100.0%    100.0%         92.3 %      100.0%
----------------------------------------------------
```


Templates ranked with the **weighted combined score**.
We can see that the program is less specific and more sensitive. The machine learning make the program find the benchmarks less fast than without, but all benchmarks are found faster than without machine learning ! That is important as it enables the program to work for more kinds of proteins.

```
Table summarizing the top N results.

         Family    Superfamily    Fold        Total

top 5    0/1       0/6            0/13        0/20
         0.0  %    0.0  %         0.0  %      0.0  %
----------------------------------------------------
top 10   0/1       0/6            0/13        0/20
         0.0  %    0.0  %         0.0  %      0.0  %
----------------------------------------------------
top 15   0/1       0/6            1/13        1/20
         0.0  %    0.0  %         7.7  %      5.0  %
----------------------------------------------------
top 20   0/1       1/6            3/13        4/20
         0.0  %    16.7 %         23.1 %      20.0 %
----------------------------------------------------
top 25   0/1       1/6            3/13        4/20
         0.0  %    16.7 %         23.1 %      20.0 %
----------------------------------------------------
top 50   0/1       3/6            3/13        6/20
         0.0  %    50.0 %         23.1 %      30.0 %
----------------------------------------------------
top 75   0/1       4/6            4/13        8/20
         0.0  %    66.7 %         30.8 %      40.0 %
----------------------------------------------------
top 100  0/1       4/6            6/13        10/20
         0.0  %    66.7 %         46.2 %      50.0 %
----------------------------------------------------
top 150  0/1       5/6            9/13        14/20
         0.0  %    83.3 %         69.2 %      70.0 %
----------------------------------------------------
top 200  0/1       6/6            12/13       18/20
         0.0  %    100.0%         92.3 %      90.0 %
----------------------------------------------------
top 250  1/1       6/6            13/13       20/20
         100.0%    100.0%         100.0%      100.0%
----------------------------------------------------
top 300  1/1       6/6            13/13       20/20
         100.0%    100.0%         100.0%      100.0%
----------------------------------------------------
top 350  1/1       6/6            13/13       20/20
         100.0%    100.0%         100.0%      100.0%
----------------------------------------------------
```


#### Generated plot
<p align="center">
  <img width="750" src="results/queries/plots/all_scores_plot.png" alt="Enrichment"/>
</p>

## Documentation

The documentation of our program is generated with Sphinx and and built on [Read The Docs](https://fold-u.readthedocs.io/en/latest/?badge=latest).


## Authors

We are master students in bioinformatics at Paris Diderot University.
- [Franz-Arnold Ake](https://github.com/franzx5)
- [Gabriel Cretin](https://github.com/gabrielctn)
- [Tom Gutman](https://github.com/tomgutman)
- [Hélène Kabbech](https://github.com/kabhel)
- [Flora Mikaeloff](https://github.com/FloraMika)


## Acknowledgment

Thanks to [Maïté Cretin](https://www.linkedin.com/in/maitewho/) for the nice logo.

## License

This project is licensed under the MIT License.
