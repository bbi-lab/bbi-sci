# bbi-sci

## Intro
This pipeline is the under-construction pipe for processing 2-level and 3-level data. It uses the pipeline management system Nextflow to operate.

The pipeline is run in two parts, the first is [bbi-dmux](https://github.com/bbi-lab/bbi-dmux) which runs the demultiplexing, and the second is [bbi-sci](https://github.com/bbi-lab/bbi-sci/) which completes the preprocessing. The instructions below apply to both pipelines, and both pipelines can use the same configuration file.

## Prerequisites
1. This script requires Nextflow version >= 20.07.1 and <= 22.10.4.

2. As the Nextflow pipeline is run interactively, please use a terminal multiplexer such as tmux. tmux sessions are persistent which means that programs in tmux will continue to run even if you get disconnected. You can start a tmux session using the command:

```
tmux
```

If you get disconnected, you can return to the head node you were using (grid-head1 or grid-head2) and type:

```
tmux attach
```

which will return you to your session. See a handy tmux tutorial [here](https://www.hostinger.com/tutorials/tmux-beginners-guide-and-cheat-sheet/).

3. Always start with a qlogin session before you begin the pipeline. This can be done using

```
qlogin -l mfree=16G
```

## Installation


If you install the pipeline on a cluster with a mix of CPU architectures, when you qlogin to the cluster for the installation procedure,
request a node with the minimum CPU ID level on which you intend to run the pipeline.
For example, on the Shendure lab cluster use

```
qlogin -l mfree=20G -l cpuid_level=11
```

Omit `-l cpuid_level` when running the pipeline.

### Modules
After starting a qlogin session:

First, you need to have python available. You should have version 3.12.1 in order to have nextflow work for you. Please make sure that this is the version you load in your ~/.bashrc file as this is the version that you will use to install the packages below. For example, in your ~/.bashrc file have:

```
module load python/3.12.1
```

After loading the above modules, you must install the following python packages:

```
pip install --user biopython
pip install --user fmt
pip install --user pysam
pip install --user matplotlib

git clone https://github.com/andrewhill157/barcodeutils.git
pushd barcodeutils
python setup.py install --user
popd
```

Install monocle3, garnett, DropletUtils, and randomColoR by loading the R/4.3.2 module and running R:

```
module load R/4.3.2
R
```

Then from within R, follow the installation instructions for the following R packages:
- Monocle3: [monocle3 website](https://cole-trapnell-lab.github.io/monocle3/)
- Garnett: [Garnett website](https://cole-trapnell-lab.github.io/garnett/docs_m3/#install-from-github)
- DropletUtils: [DropletUtils website](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html)
- randomColoR: Run ```install.packages("randomcoloR")

You will also require scrublet, a tool used to detect doublets in single-cell RNA-seq data. You can install it from source by running:

```
git clone https://github.com/AllonKleinLab/scrublet.git
pushd scrublet
pip install -r requirements.txt --user
python setup.py install --user
popd
```

Please note: If you are doing a hashing experiment, you need scipy, which is installed already in python 3.12.1.


Once monocle3 and scrublet are installed, install nextflow by typing:

```
curl -s https://get.nextflow.io | bash
```
You probably also want to add Nextflow to your path so you can access it from anywhere. Do this by adding the following to your .bashrc file (located in your home directory).

```
export PATH=/path/to/where/you/installed/nextflow/:$PATH
```

Next, pull the pipeline to make sure you're on the latest version

```
nextflow pull bbi-lab/bbi-dmux
nextflow pull bbi-lab/bbi-sci
```

Build the pypy3 virtual environment required for the bbi-dmux pipeline. Do this in the directory ~/.nextflow/assets/bbi-lab/bbi-dmux using the commands

```
pushd ~/.nextflow/assets/bbi-lab/bbi-dmux
bash create_virtual_envs.sh
popd
```

The environment may need to be built on a node with the CPU architecture on which you will run the scripts. There are additional details in the create_virtual_envs.sh script.

Check that it works by running:

```
nextflow run bbi-dmux --help
nextflow run bbi-sci --help
```

You should get some help info printed.

Alternatively, you can install and run the pipeline from clones of the Github repositories, for example,

```
mkdir ~/git
cd ~/git
git clone https://github.com/bbi-lab/bbi-dmux
git clone https://github.com/bbi-lab/bbi-sci
cd bbi-dmux
bash create_virtual_envs.sh
```

In this case, you run the pipelines using the commands

```
nextflow run ~/git/bbi-dmux/main.nf -profile ubuntu_22_04 -c experiment.config
nextflow run ~/git/bbi-sci/main.nf -profile ubuntu_22_04 -c experiment.config
```

There are bash scripts at

```
~/git/bbi-dmux/scripts/run.scirna-demux.sh
~/git/bbi-dmux/scripts/run.scirna-analyze.sh
```

that you can edit and use to run the pipelines.

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

The second things you need are configuration files that pass in arguments to the pipeline. These are the *experiment.config* and *nextflow.config* files.

##### *experiment.config* file

 The *experiment.config* file is helpful as it allows you to specify if your data is 2-level or 3-level, allocate memory requirements, process only a subset of your samples and use custom genomes to map your data. We highly recommend using this instead of passing arguments on the command line so that you have a record of the run you called.

Notes:
- an example experiment configuration file is included in the package and includes further information on usage

- for the Shendure lab cluster use `process.queue = "shendure-long.q"` in either the *experiment.config* or *nextflow.config* files

##### *nextflow.config* file

The *nextflow.config* file defines processing values such as the required modules, memory, and number of CPUs for each processing stage, which do not change typically from run-to-run. The file can be left in the bbi-\* installation directory where Nextflow searches for it automatically when the pipeline starts up. The supplied *nextflow.config* file has two profiles: the default profile, called *standard*, defines modules used by the pipeline on CentOS 7 systems in the UW Genome Sciences cluster, and the *ubuntu_22_04* profile, which defines modules used by the pipeline on Ubuntu 22.04 systems in the UW Genome Sciences cluster. In order to run the pipelines with the *ubuntu_22_04* profile, add the command line parameter `-profile ubuntu_22_04` to the nextflow run command, for example

```
nextflow run bbi-dmux -profile ubuntu_22_04 -c experiment.config
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


For either piece of the pipelines, if there is an error, you can continue the pipeline where it left off with either

```
nextflow run bbi-dmux -c experiment.config -resume
```
or

```
nextflow run bbi-sci -c experiment.config -resume
```

#### Recovery mode:

After running the dmux part of the pipeline, you will find all unassigned reads in fastq files labeled "Undetermined...". You can run bbi-dmux in 'recovery mode' to generate a table with information about why each read wasn't assigned, and a summary file with percentages. There will be a table and summary file for each lane.

Run recovery mode (AFTER running bbi-dmux as above) like this:

```
nextflow run bbi-dmux --run_recovery true -c experiment.config
```
You should provide the same experiment.config file that you used above.

#### The work folder:
Nextflow stores all of the intermediate files in its 'work' folder, which will be in the output directory you specified. This folder can get quite large, so after the pipeline is finished, you can delete it using:

```
rm -r work/
```

Warning: After you delete the work folder, -resume will no longer restart from the middle of the run, you'll have to start from the beginning if you need to regenerate any files. So please make sure that your pipeline has fully completed before deleting the work folder.

#### Troubleshooting:

##### I got a "UTF encoding" error for my samplesheet
If you get this error, it means that the python script for reading the samplesheet encountered an issue with the UTC encoding
of your samplesheet. We haven't pinpointed why this sometimes happens (probably something to do with your Excel settings when making 
the sheet), but fixing it is pretty straight forward. You have two options:
- Option 1 (from the terminal): Open the samplesheet in the terminal (using your preferred editor, e.g. vim or nano) and 
copy the contents by highlighting and using Cmd + C or Ctrl + C. Open a new file in your editor and paste the contents. Save and 
use this new file as your samplesheet.

- Option 2 (from your local computer): Open the samplesheet in Excel and use Save As to save the sheet using the file format 
"Comma Separated Values (.csv)". Use this new file as your samplesheet. 

#### Questions and errors:
If you run into problems, please leave a detailed description in the issue tab above!

### Acknowledgements
Many members of the Shendure and Trapnell labs as well as the BBI team have worked on portions of this pipeline or its predecessors, especially Andrew Hill and Jonathan Packer. Many thanks to them!
