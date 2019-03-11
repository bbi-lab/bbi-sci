# bbi-sci

## Intro
This pipeline is the under-construction pipe for processing 2-level data. It uses the pipeline management system Nextflow to operate.

## Installation
First, install nextflow by logging onto the cluster, starting a qlogin session and typing:

```
curl -s https://get.nextflow.io | bash
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
The second thing you need is a config file which passes in your arguments to the pipeline. I highly recommend using this instead of passing arguments on the command line so that you have a record of the run you called.

An example configuration file is included in the package.

#### Run the pipeline:

To run the pipeline, do 

```
nextflow run bbi-sci -c experiment.config
```
