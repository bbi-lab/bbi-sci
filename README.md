# bbi-sci

## Intro
This pipeline is the under-construction pipe for processing 2-level and 3-level data. It uses the pipeline management system Nextflow to operate.

The pipeline is run in two parts, the first is [bbi-dmux](https://github.com/bbi-lab/bbi-dmux) which runs the demultiplexing, and the second is [bbi-sci](https://github.com/bbi-lab/bbi-sci/) which completes the preprocessing. The instructions below apply to both pipelines, and both pipelines can use the same configuration file.

## Prerequisites
1. As the Nextflow pipeline is run interactively, please use a terminal multiplexer such as tmux. tmux sessions are persistent which means that programs in tmux will continue to run even if you get disconnected. You can start a tmux session by using:
```
module load tmux/latest
tmux
```
If you get disconnected, you can return to the head node you were using (grid-head or shead) and type:
```
tmux attach
```
which will return you to your session. See a handy tmux tutorial [here](https://www.hostinger.com/tutorials/tmux-beginners-guide-and-cheat-sheet/).

2. Always start with a qlogin session before you begin the pipeline. This can be done by
```
qlogin -l mfree=10G
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

After loading the above modules, you must install the following python package:

```
pip install --user drmaa
pip install --user biopython
pip install --user fmt

git clone https://github.com/andrewhill157/barcodeutils.git
cd barcodeutils
python setup.py install --user
cd ..
```

Then, install monocle3 by running:

```
module load gcc/8.1.0
module load R/3.5.2
R
```
Then from within R, follow the installation instructions on the [monocle3 website](https://cole-trapnell-lab.github.io/monocle3/monocle3_docs/#installing-monocle-3).

Once monocle3 is installed, install nextflow by typing:

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
The second thing you need is a config file which passes in your arguments to the pipeline. This file is very helpful as it allows you to specify if your data is 2-level or 3-level, allocate memory requirements, process only a subset of your samples and use custom genomes to map your data. We highly recommend using this instead of passing arguments on the command line so that you have a record of the run you called.

An example configuration file is included in the package and includes further information on usage.

For Shendure lab cluster
```
process.queue = "ravana.q"
```

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
