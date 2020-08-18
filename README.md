# bbi-sci

## Intro
This pipeline is the under-construction pipe for processing 2-level and 3-level data. It uses the pipeline management system Nextflow to operate.

The pipeline is run in two parts, the first is [bbi-dmux](https://github.com/bbi-lab/bbi-dmux) which runs the demultiplexing, and the second is [bbi-sci](https://github.com/bbi-lab/bbi-sci/) which completes the preprocessing. The instructions below apply to both pipelines, and both pipelines can use the same configuration file.

## Prerequisites
1. This script requires Nextflow version >= 20.

2. As the Nextflow pipeline is run interactively, please use a terminal multiplexer such as tmux. tmux sessions are persistent which means that programs in tmux will continue to run even if you get disconnected. You can start a tmux session by using:
```
module load tmux/latest
tmux
```
If you get disconnected, you can return to the head node you were using (grid-head1 or grid-head2) and type:
```
tmux attach
```
which will return you to your session. See a handy tmux tutorial [here](https://www.hostinger.com/tutorials/tmux-beginners-guide-and-cheat-sheet/).

3. Always start with a qlogin session before you begin the pipeline. This can be done by
```
qlogin -l mfree=20G
```

## Installation

### modules
After starting a qlogin session:

First, you need to have python available. You should have version 3.6.4 in order to have nextflow work for you. Please make sure that this is the version you load in your ~/.bashrc file as this is the version that you will use to install the packages below. For example, in your ~/.bashrc file have:

```
module load python/3.6.4
```

You must also have a few modules other than python loaded:

```
module load drmaa/latest
module load git/latest
```

After loading the above modules, you must install the following python packages:

```
pip install --user drmaa
pip install --user biopython
pip install --user fmt
pip install --user pysam

git clone https://github.com/andrewhill157/barcodeutils.git
cd barcodeutils
python setup.py install --user
cd ..
```

Then, install monocle3 and garnett by running:

```
module load gcc/8.1.0
module load R/3.6.1
R
```
Then from within R, follow the installation instructions on the [monocle3 website](https://cole-trapnell-lab.github.io/monocle3/).
And the instructions for garnett on the [Garnett website](https://cole-trapnell-lab.github.io/garnett/docs_m3/#install-from-github).

You will also require scrublet, a tool used to detect doublets in single-cell RNA-seq data. You can install it from source by running:

```
git clone https://github.com/AllonKleinLab/scrublet.git
cd scrublet
pip install -r requirements.txt --user
python setup.py install --user
cd ..
```

Once monocle3 and scrublet are installed, install nextflow by typing:

```
curl -s https://get.nextflow.io | bash
```
You probably also want to add Nextflow to your path so you can access it from anywhere. Do this by adding the following to your .bashrc file (located in your home directory).
```
export PATH=/path/to/whereever/you/downloaded/nextflow:$PATH
```

Next, pull the pipeline to make sure you're on the latest version
```
nextflow pull bbi-lab/bbi-dmux
nextflow pull bbi-lab/bbi-sci
```

Check that it all worked by running:
```
nextflow run bbi-dmux --help
nextflow run bbi-sci --help
```
You should get some help info printed.

## Running the pipeline

To run the pipeline you need two files, a sample sheet and a configuration file.

#### Sample sheet:
The sample sheet should be a csv with 3 columns: RT Barcode, Sample ID, and Reference Genome. Here's an example:

```
RT Barcode,Sample ID,Reference Genome
2P9-A01,Sample1,Mouse
2P9-A02,Sample1,Mouse
2P9-A03,Sample2,Human
2P9-A04,Sample2,Human
```

#### Configuration file:

The second thing you need are configuration files that pass in arguments to the pipeline. These are the *experiment.config* and *nextflow.config* files.

##### *experiment.config* file

 The *experiment.config* file is helpful as it allows you to specify if your data is 2-level or 3-level, allocate memory requirements, process only a subset of your samples and use custom genomes to map your data. We highly recommend using this instead of passing arguments on the command line so that you have a record of the run you called.

Notes:
- an example experiment configuration file is included in the package and includes further information on usage

- for the Shendure lab cluster use `process.queue = "shendure-long.q"` in either the *experiment.config* or *nextflow.config* files

##### *nextflow.config* file

The *nextflow.config* file defines processing values such as the required modules, memory, and number of CPUs for each processing stage, which do not change typically from run-to-run. The file can be left in the bbi-\* installation directory where Nextflow searches for it automatically when the pipeline starts up. The supplied *nextflow.config* file has two profiles: the default profile, called *standard*, defines modules used by the pipeline on CentOS 7 systems in the UW Genome Sciences cluster, and the *centos_6* profile, which defines modules used by the pipeline on CentOS 6 systems in the UW Genome Sciences cluster. In order to run the pipelines with the *centos_6* profile, add the command line parameter `-profile centos_6` to the nextflow run command, for example

```
nextflow run bbi-dmux -profile centos_6 -c experiment.config
```

This *nextflow.config* file has comments that give additional information.

#### Run the pipeline:

After making your configuration file, run the pipeline by first running:

```
nextflow run bbi-dmux -c experiment.config
```

After running the dmuxing step, check out the dmuxing dashboard by downloading the dmux_dash folder and opening the html in a web browser. When you're satisfied that dmuxing was successful, run:
```
nextflow run bbi-sci -c experiment.config
```


For either piece of the pipeline, if there is an error, you can continue the pipeline where it left off with either
```
nextflow run bbi-dmux -c experiment.config -resume
```
or
```
nextflow run bbi-sci -c experiment.config -resume
```

#### The work folder:
Nextflow stores all of the intermediate files in its 'work' folder, which will be in the output directory you specified. This folder can get quite large, so after the pipeline is finished, you can delete it using:

```
rm -r work/
```

Warning: After you delete the work folder, -resume will no longer restart from the middle of the run, you'll have to start from the beginning if you need to regenerate any files. So please make sure that your pipeline has fully completed before deleting the work folder.

#### Questions and errors:
If you run into problems, please leave a detailed description in the issue tab above!

### Acknowledgements
Many members of the Shendure and Trapnell labs as well as the BBI team have worked on portions of this pipeline or its predecessors, especially Andrew Hill and Jonathan Packer. Many thanks to them!
