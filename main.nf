#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import groovy.json.JsonSlurper

def cfg_path = params.config ?: "config.json"
def cfg = new JsonSlurper().parseText(file(cfg_path).text)

println "GeneCat Directory:\t${cfg.genecat}"
println "Out Directory:\t${cfg.outdir}"
println "Scratch Directory :\t${cfg.tempdir}"
println "MGS List (if not specified, ALL):\t${cfg.mgslist}"
println "Parameter to be used:\t${cfg.param}"
println "Working Environment:\t${cfg.env}"

process SampleListGet {
    tag "SampleListGet"
    publishDir "${cfg.outdir}/MAGQuality", mode: 'copy' 

    output:
        path "*.txt", emit: tsv
    script:
        """
        Mapfile=${cfg.genecat}/LOGandSUB/map.0.txt

        IFSback=\${IFS}
        IFS=\$'\n'
        echo "Sample" > Cross_sampeId_temp.txt
        echo "Sample" > Long_sampeId_temp.txt
        for line in `cat \${Mapfile}| grep -v "#"`;do
            Judge=`echo \${line}| cut -f 3`
            if [ \${#Judge} == "0" ];then
                echo -e "Listed\t\${line}| cut -f 1"
                echo \${line}| cut -f 1 >> Cross_sampeId_temp.txt
            else
                echo -e "Listed\t\${line}| cut -f 3"
                echo \${line}| cut -f 3 >> Long_sampeId_temp.txt

            fi
        done

        cat Cross_sampeId_temp.txt | tail -n +2 > Cross_sampeId.txt
        rm Cross_sampeId_temp.txt

        cat Long_sampeId_temp.txt | tail -n +2| sort|uniq > Long_sampeId.txt
        rm Long_sampeId_temp.txt

        IFS=\${IFSback}



        if [ ${cfg.mgslist} == "ALL" ];then
            cat ${cfg.genecat}/Bin_SB/SB.clusters.obs | tail -n +2| cut -f 1 > MGS.txt
        else
            cp ${cfg.mgslist} MGS.txt
        fi
        
        """
}
process Cross_CheckMProfile {
    tag "Cross_CheckMProfile"
    publishDir "${cfg.tempdir}/Cross_CheckMProfile", mode: 'copy' 

    input: 
        val(filename)
    output:
        path "*.tsv"
        val(true), emit: cross_done
    script:
        """  
        Dirname=`dirname ${cfg.genecat}` 
        tail -n +2 \${Dirname}/${filename}/assemblies/metag/Binning/SB/${filename}.cm2 > "${filename}.tsv"
        """
}
process Long_CheckMProfile {
    tag "Long_CheckMProfile"
    publishDir "${cfg.tempdir}/Long_CheckMProfile", mode: 'copy' 

    input: 
        val(filename)
    output:
        path "*.tsv", emit: tsv
        val(true), emit: long_done
    script:
        """  
        echo "hello" > temp.tsv
        Mapfile=${cfg.genecat}/LOGandSUB/map.0.txt
        SampleNum=`cat \${Mapfile} | grep -w ${filename}| wc -l`

        if [ \${SampleNum} = 1 ];then
            sample=`cat \${Mapfile} | grep -w ${filename}|cut -f 1`
            Dirname=`dirname ${cfg.genecat}` 
            tail -n +2 \${Dirname}/\${sample}/assemblies/metag/Binning/SB/\${sample}.cm2 > "\${sample}.tsv"
        else
            sampleRepName=`cat \${Mapfile} | grep -w ${filename}|tail -n 1|cut -f 1`
            Dirname=`dirname ${cfg.genecat}` 
            tail -n +2 \${Dirname}/AssmblGrp_${filename}/metag/Binning/SB/\${sampleRepName}.cm2 > "\${sampleRepName}M\${SampleNum}.tsv"
        fi
        
        cat \${Mapfile} | grep -w ${filename}
        """
}

process CheckMFilter {
    tag "CheckMFilter"
    publishDir "${cfg.outdir}/${filename}", mode: "copy"

    input: 
        tuple val(filename), val(flag)
    output: 
        path "*.tsv", emit: tsv
        val(true), emit: checkMdone
    script:

        """
        #!/usr/bin/env python

        import pandas as pd
        from sys import argv
        import os
        import subprocess
        import json

        MGS="${filename}"
        CrossTEMPDIR=f"${cfg.tempdir}/Cross_CheckMProfile"
        LongTEMPDIR=f"${cfg.tempdir}/Long_CheckMProfile"
        checkMCol=["Name","Completeness","Contamination","Completeness_Model_Used","Translation_Table_Used","Coding_Density","Contig_N50","Average_Gene_Length","Genome_Size","GC_Content","Total_Coding_Sequences","Total_Contigs","Max_Contig_Length","Additional_Notes"]

        CrossSamples=pd.read_table(f"${cfg.outdir}/MAGQuality/Cross_sampeId.txt", names=['sample'])['sample'].tolist()

        GTDBTK= pd.read_table("${cfg.genecat}" + "/Bin_SB/Annotation/GTDBTK.tax",index_col=0, low_memory=False).loc[MGS]['classification']
        Cluster= pd.read_table("${cfg.genecat}" + "/Bin_SB/SB.clusters.obs",index_col=0, low_memory=False).loc[[MGS]]
        Cluster.at[MGS, 'taxon'] = GTDBTK
        
        Cluster.to_csv(f"{MGS}_metadata.tsv", sep='\t')

        binDF = pd.DataFrame(index=Cluster.at[MGS, 'Members'].split(","))
        binDF.index.name = 'bin'
        binDF['sample'] = [i.split("__")[0] for i in binDF.index]

        binDF = binDF.loc[[i for i in binDF.index if "Cano__" not in i]]
        
        with open("${cfg.param}", "r") as f:
            param = json.load(f)
        
        Completeness=param['comp']
        Contamination=param['cont']

        CheckMSummary = pd.DataFrame()
        for sample, binDFTemp in binDF.groupby('sample'):
            if sample in CrossSamples:
                CheckM = pd.read_table(f"{CrossTEMPDIR}/{sample}.tsv", names=checkMCol, index_col=0)
            else:
                CheckM = pd.read_table(f"{LongTEMPDIR}/{sample}.tsv", names=checkMCol, index_col=0)

            CheckM = CheckM.loc[[i.split("__")[-1] for i in binDFTemp.index]].reset_index()
            CheckM['Sample'] = sample
            CheckM = CheckM[['Sample']+checkMCol]
            CheckMSummary = pd.concat([CheckMSummary, CheckM] , axis=0)
        print(CheckMSummary)
        CheckMSummary.to_csv(f"{MGS}_cm_raw.tsv", index=False, sep='\t')
        CheckMSummary = CheckMSummary[CheckMSummary['Completeness'] > Completeness]
        CheckMSummary = CheckMSummary[CheckMSummary['Contamination'] < Contamination]
        CheckMSummary.to_csv(f"{MGS}_cm_filtered.tsv", index=False, sep='\t')
        """
}

process MAGConstruction {
    tag "MAGConstruction"
    publishDir "${cfg.outdir}/${filename}/Link"     , pattern:"*.tsv", mode: "copy"
    publishDir "${cfg.outdir}/${filename}/Fasta"     , pattern:"*.fa.gz", mode: "copy"


    input: 
        tuple val(filename), val(flag)
    output: 
        path "*.tsv", emit: tsv
        path "*.fa.gz", emit: fa_gz

        stdout emit: log
    script:

        """
        Dirname=`dirname ${cfg.genecat}` 

        Mapfile=${cfg.genecat}/LOGandSUB/map.0.txt
        checkMFiltered=${cfg.outdir}/${filename}/${filename}_cm_filtered.tsv
        IFSback=\${IFS}
        IFS=\$'\n'

        for line in `cat \${checkMFiltered}| tail -n +2| cut -f 1,2`;do
            sample=`echo \${line}| cut -f 1`
            sample2=\${sample%M[0-9]*}
            bin=`echo \${line}| cut -f 2|cut -f 1 -d "."`
            binComplete="\${sample2}__\${bin}"
            mapLine=`cat \${Mapfile}| grep -w \${sample2}`

            if [ `echo \${mapLine}| cut -f 3| grep "Grp"` ];then
                group=`echo \${mapLine}| cut -f 3` #group
                flag1=Longitudinal

                #Longitudinal or Single
                numSampleinGrp=`cat \${Mapfile}|grep \${group}|wc -l`

                if [ \${numSampleinGrp} = 1 ];then
                    flag2=Single
                    AssmblGrp=""
                    Dirname2=\${Dirname}/\${sample2}/assemblies/metag
                else
                    flag2=Multiple
                    AssmblGrp=\${group}
                    Dirname2=\${Dirname}/AssmblGrp_\${group}/metag
                fi
            else
                flag1=CrossSectional
                flag2=""
                AssmblGrp=""
                Dirname2=\${Dirname}/\${sample2}/assemblies/metag
            fi
            assemble=\${Dirname2}/scaffolds.fasta.filt
            connect=\${Dirname2}/Binning/SB/\${sample2}
            gfffile=\${Dirname2}/genePred/genes.gff

            cat \${connect}| grep -w "\${bin}.fa.gz"| cut -f 1 > \${binComplete}.tsv


            seqkit grep -n -f \${binComplete}.tsv \${assemble} > \${binComplete}.fa
            gzip -f \${binComplete}.fa


        done

        IFS=\${IFSback}
        """
}

workflow {

    def sampleList = SampleListGet()
    def (crossSample, longSample, mgsList) = sampleList.multiMap {
        first: it[0]
        second: it[1]
        third: it[2]
    }
    def mgsListpath = mgsList.splitText().map { it.trim() }

    def crossLines = crossSample.flatten().splitText().map { it.trim() }.filter { it }
    def (cross_tsv_ch, cross_log) = Cross_CheckMProfile(crossLines)
    
    def longLines = longSample.flatten().splitText().map { it.trim() }.filter { it }
    def (long_tsv_ch, long_log) = Long_CheckMProfile(longLines)

    def all_done = long_log.mix(cross_log).count().map { true }
    def combined = mgsListpath.combine(all_done)
    def (checkM_tsv_ch, checkMdone) = CheckMFilter(combined)
    
    def combined_checkM = mgsListpath.combine(checkMdone)
    def (fasta_tsv_ch, fa_ch, log_ch) = MAGConstruction(combined_checkM)


}
