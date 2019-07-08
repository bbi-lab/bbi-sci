# bbi-sci

## Intro
This pipeline is the under-construction pipe for processing 2-level and 3-level data. It uses the pipeline management system Nextflow to operate.

## Prerequisites
1. As the Nextflow pipeline is run interactively, please use a terminal multiplexer such as tmux. tmux sessions are persistent which means that programs in tmux will continue to run even if you get disconnected.
2. Always start with a qlogin session before you begin the pipeline. This can be done by
```
qlogin -l mfree=10G
```
3. Make sure you have R and Monocle 3 installed prior to running the pipeline.

## Installation
First, install nextflow by logging onto the cluster, starting a qlogin session and typing:

```
curl -s https://get.nextflow.io | bash
```
You probably also want to add Nextflow to your path so you can access it from anywhere. Do this by adding the following to your .bashrc file (located in your home directory).
```
export PATH=/path/to/whereever/you/downloaded/nextflow:$PATH
```

Next, pull the pipeline to make sure you're on the latest version
```
nextflow pull bbi-lab/bbi-sci -user [your github id]
```
Enter your github password when it prompts you.

Check that it all worked by running:
```
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

An example configuration file is included in the package.

For Shendure lab cluster
```
process.queue = "ravana.q"
```

#### Run the pipeline:

To run the pipeline, do

```
nextflow run bbi-sci -c experiment.config
```

If there is an error, you can continue the pipeline where it left off with

```
nextflow run bbi-sci -c experiment.config -resume
```

#### The work folder:
Nextflow stores all of the intermediate files in its 'work' folder, which will be in the output directory you specified. This folder can get quite large, so after the pipeline is finished, you can delete it using:

```
rm -r work/
```

Warning: After you delete the work folder, -resume will no longer restart from the middle of the run, you'll have to start from the beginning if you need to regenerate any files. So please make sure that your pipeline has fully completed before deleting the work folder.
