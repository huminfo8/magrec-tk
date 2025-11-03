#!/bin/bash
#SBATCH -N 1                           # number of nodes
#SBATCH -n 1                           # number of tasks
#SBATCH -c 1                           # number of cores
#SBATCH -p "qib-long,ei-long"
#SBATCH --mem 24G                    # memory
#SBATCH --mail-type=ALL                # mail me everything
#SBATCH -o ./logs/%x_%A_%a_%j.otxt
#SBATCH -e ./logs/%x_%A_%a_%j.etxt
#SBATCH -J MAGRecTK

configFile="config.tsv"
py_1=SCRIPT/jsonConverter.py

if [ -f $configFile ];then
    echo -e ">>>LOG\tConfig File Exist! Continue"
    ENVPATH=`cat $configFile | grep -w "env"| cut -f 2`


    eval "$($(command -v micromamba) shell hook --shell bash)"
    micromamba activate ${ENVPATH}

    python $py_1 "`pwd`/"$configFile

    TMPDIR=`jq -r '.tempdir' config.json`
    OUTDIR=`jq -r '.outdir' config.json`

    if [ ! -d ${TMPDIR} ];then
        mkdir ${TMPDIR}
    fi
    if [ ! -d ${OUTDIR} ];then
        mkdir ${OUTDIR}
    fi

    echo -e "params {" > nextflow.config
    echo -e '\ttmpdir = '"'${TMPDIR}'" >> nextflow.config
    echo -e "}" >> nextflow.config
    echo -e "env {" >> nextflow.config
    echo -e '\tTMPDIR = params.tmpdir' >> nextflow.config
    echo -e "}" >> nextflow.config
    echo -e "process {" >> nextflow.config
    echo -e "\texecutor = 'slurm'" >> nextflow.config
    echo -e "\tscratch = true" >> nextflow.config
    echo -e '\tclusterOptions = "--partition=qib-long,ei-long --export=ALL,TMPDIR='\${params.tmpdir}'"' >> nextflow.config
    echo -e "\tbeforeScript = 'mkdir -p \"\$TMPDIR\"'" >> nextflow.config
    echo -e "}" >> nextflow.config

    chmod 770 ${TMPDIR} #To referer from the calculating nodes
    
    echo -e ">>>LOG\tTEMPDIR\t$TMPDIR"
    export NXF_TEMP=${TMPDIR}
    export NXF_DISABLE_CHECK_LATEST=true #Deactivate
    export NXF_OPTS='-Xms256m -Xmx1g -XX:MaxRAMPercentage=40 -Djava.io.tmpdir='"$TMPDIR"' -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -XX:ErrorFile='"$TMPDIR"'/hs_err_pid%p.log'
    nextflow run main.nf --params-file config.json

else
    echo -e ">>>LOG\tPrepare yout Config File Exit"
fi

