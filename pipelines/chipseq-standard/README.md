# ChIP-Seq Standard Pipeline

## 0. Create a new analysis

If a new analysis project has not already been created, do so with the following command:

```
~/pipeline-master/code/code.main/pipeline-new-analysis chipseq-standard /path/to/<project_name>
```

## 1. Set input files

Within the corresponding `<project_name>/inputs/fastq` or `<project_name>/inputs/bam` directory, subdirectories should be created with the name of each sample to be included in the analysis. The following naming scheme is preferable:

<Cell_line>-<ChIP>-<treatment>-<SampleID>

Each subdirectory should contain all fastq or bam files to be used for that sample through the analysis pipeline. Symlinks can be used if the files are not contained in the same location as the project analysis directory, and are preferable to save storage space. 

*QQ: Include sample directory file structures*
*QQ: Include sample code for creating symlinks and dirs?*
*QQ: treatment of special characters, etc., in file/dir names*

## 2. Create project sample sheet

A sample sheet must be created for the analysis project. Run the follow command to do so:

```
<project>/inputs/code/create-sample-sheet.tcsh <genome> <fragment-size>
```




INSTRUCTIONS

Go into inputs directory and read the README file.

Then, return to this directory and start the pipeline as follows:

./code.main/pipeline-execute PROJECT-NAME E-MAIL


