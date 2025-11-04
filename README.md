# MAGRec-TK
## Introduction
MAGRec-TK (**M**etagenome **A**ssembled **G**enome **Rec**overy **T**ool**K**it) is a pipeline developed to recover MAG from [MG-TK](https://github.com/hildebra/mg-tk) outputs.
MAGRec-TK is implemented by python, shell and nextflow framework.
Author: Daisuke Suzuki dsuzuki8forb@gmail.com

## Installation
1. Run command `git clone https://github.com/huminfo8/magrec-tk.git` (Might need to enter the node can download)
2. `cd magrec-tk`
4. Run ` bash installer.sh`
5. Check the environmental path (ex. /hpc-home/user/micromamba/envs/MAGRecTK)
## Prepare
### Fill out the **config.tsv** and **param.tsv**
Parameters and their explanations. These must be saved in tsv format.
#### config.tsv
| Name  | Requirement | Explanation | Example |
| ------------- | ------------- |------------- |------------- |
| genecat  | Mandatory| MG-TK genecat path | /ei/projects/c/c0293548-5d3d-418c-9de2-db0fcffd25e9/data/CRC.Selective/GeneCat |
| outdir  | Optional| Specify the path you want to put your results. If not specified, then the pipeline automatically put results in "genecat/BinSB/intra_phylo". |/ei/projects/c/c0293548-5d3d-418c-9de2-db0fcffd25e9/PipelineNxt/Out|
| tempdir  | Optional| Specify the path you want to put your intermediate files. If not specified, then the pipeline automatically put results in "genecat/BinSB/intra_phylo/MAGRecTK". |/ei/projects/c/c0293548-5d3d-418c-9de2-db0fcffd25e9/scratch|
| mgslist  | Optional| Specify the file such as "MGS.tsv" that has interested MGS. If you want to get all of the MGS MAG, then specify "ALL" |mgs.tsv|
| param  |Mandatory | Specify the file thresholds for bins quality in separared file "param.tsv" |param.tsv	|
| env  | Mandatory|Specify the path you got in the installation step |/hpc-home/user/micromamba/envs/MAGRecTK|

#### param.tsv
| Name   | Explanation | Example |
| ------------- | ------------- |------------- |
| comp | completeness |80 | 
| cont | contamination |10 |

#### mgs.tsv (Optional)
|MGS.10|
|---|
|MGS.8|

## Run
Just type `sbatch run.sh` and run! Log will be appeared in logs directory.
## Workflow
<img width="1290" height="1119" alt="Image" src="https://github.com/user-attachments/assets/615dc00e-ab6a-4f25-8955-9e79436f20bc" />
