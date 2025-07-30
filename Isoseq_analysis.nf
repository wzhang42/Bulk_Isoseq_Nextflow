#!/usr/bin/env nextflow
/*******************************************************************************
 This Nextflow pipeline, is the downstream Pacbio isoseq analysis inlcuding Isoseq mapping, 
 Isoform collapsing & Classifiation and/ or fusion trasncript finding.
 Input: Pacbio isoseq fastq file list 
 Output: Isoseq classification or fusion transcript finding result 
 Written by: Wenchao Zhang
 The Center for Applied Bioinformatics,St. Jude Children Research Hosptial
 Date: Initially witten in 05/30/2022,  added the LongQC and LongGF modules in 02/17/2023 
*********************************************************************************/

REFERENCE=""
REFERENCE_PATH=""
REFERENCE_BASE=""

GENECODE=""
POLYA_SETTING=""

def parse_config_parameters() {
// Parse Nextflow Envrionment Paramters 

if( params.genome_build == 'hg38' )
{
  
  REFERENCE = params.Build_hg38.REFERENCE
  GENECODE =  params.Build_hg38.GENECODE
    
  params.ANNOVAR.ANNOVAR_DB = params.Build_hg38.ANNOVAR_DB
  params.ANNOVAR.CONFIG_ANNOVAR = params.Build_hg38.ANNOVAR_CONFIG 

  params.ANNOVAR.PROTOCOL = params.Build_hg38.ANNOVAR_PROTOCOL
  params.ANNOVAR.OPERATION = params.Build_hg38.ANNOVAR_OPERATION 

}
else if( params.genome_build == 'hg19' )
{

  REFERENCE = params.Build_hg19.REFERENCE
  GENECODE =  params.Build_hg19.GENECODE

  params.ANNOVAR.ANNOVAR_DB = params.Build_hg19.ANNOVAR_DB
  params.ANNOVAR.CONFIG_ANNOVAR = params.Build_hg19.ANNOVAR_CONFIG

  params.ANNOVAR.PROTOCOL = params.Build_hg19.ANNOVAR_PROTOCOL
  params.ANNOVAR.OPERATION = params.Build_hg19.ANNOVAR_OPERATION

   
}
else if( params.genome_build == 'b37' ) 
{
  
  REFERENCE = params.Build_b37.REFERENCE
  GENECODE =  params.Build_b37.GENECODE 

  params.ANNOVAR.ANNOVAR_DB = params.Build_b37.ANNOVAR_DB
  params.ANNOVAR.CONFIG_ANNOVAR = params.Build_b37.ANNOVAR_CONFIG

  params.ANNOVAR.PROTOCOL = params.Build_b37.ANNOVAR_PROTOCOL
  params.ANNOVAR.OPERATION = params.Build_b37.ANNOVAR_OPERATION

}
else if( params.genome_build == 'vMg27' )
{
  REFERENCE = params.Build_mg27.REFERENCE
  GENECODE =  params.Build_mg27.GENECODE 

  params.ANNOVAR.ANNOVAR_DB = params.Build_mg27.ANNOVAR_DB
  params.ANNOVAR.CONFIG_ANNOVAR = params.Build_mg27.ANNOVAR_CONFIG

  params.ANNOVAR.PROTOCOL = params.Build_mg27.ANNOVAR_PROTOCOL
  params.ANNOVAR.OPERATION = params.Build_mg27.ANNOVAR_OPERATION  
}
else
{
  error "Invalid geome version: ${params.genome_build}"
  exit 1
}

REFERENCE_PATH = file(REFERENCE).getParent()
REFERENCE_BASE = file(REFERENCE).getBaseName()

if (params.Protocol_Polya == 'Y') 
{
   POLYA_SETTING = "--require-polya"
}

}
//mapper                  : ${params.mapper}
def DispConfig() {
 log.info """
Welocme to run Nextflow Pipeline Isoseq_analysis.nf  
Your configuration are the following:
  project                   : ${params.project}
  isoseq_filelist           : ${params.isoseq_filelist}
  outdir                    : ${params.outdir}
    
  genome_build              : ${params.genome_build}
  Protocol_Polya            : ${params.Protocol_Polya}
  PRIMERS_FA                : ${params.PRIMERS_FA}
 
  Select_Isoseq_Lima        : ${params.Select_Isoseq_Lima}
  Select_Isoseq_Refine      : ${params.Select_Isoseq_Refine}
  Select_Isoseq_Cluster     : ${params.Select_Isoseq_Cluster}
  Cluster_Singleton         : ${params.Isoseq_Cluster.Cluster_Singleton}
  
  Select_minimap2           : ${params.Select_minimap2}
  Select_PBMM2              : ${params.Select_PBMM2}

  Select_SplitNCigarReads   : ${params.Select_SplitNCigarReads}
  Select_Deepvariant        : ${params.Select_Deepvariant}
  Select_SNV_ANNOVAR        : ${params.Select_SNV_ANNOVAR}
  vsc_min_fraction_snp:     : ${params.Deepvariant.vsc_min_fraction_snp}
  vsc_min_fraction_indel    : ${params.Deepvariant.vsc_min_fraction_indel}

  Select_LongQC_Isoseq      : ${params.Select_LongQC_Isoseq}
  Sequencing_platform       : ${params.LongQC.Sequencing_platform}
  LongQC_NSAMPLE            : ${params.LongQC.NSAMPLE}

  Select_LongGF             : ${params.Select_LongGF}
  LongGF_min_overlap_len    : ${params.LongGF.min_overlap_len}
  LongGF_bin_size           : ${params.LongGF.bin_size}
  LongGF_min_map_len        : ${params.LongGF.min_map_len}

  Select_PBFusion           : ${params.Select_PBFusion}
  PBFusion_sorted_gtf       : ${params.PBFusion.gencode_annotation_sorted_gtf}
   
  Select_Isoseq_Collapse    : ${params.Select_Isoseq_Collapse}
  Collaspe_Min_Coverage     : ${params.Isoseq_Collapse.Collaspe_Min_Coverage}
  Collaspe_Min_Identity     : ${params.Isoseq_Collapse.Collaspe_Min_Identity}
 
  Select_Transcript_Filtering:${params.Select_Transcript_Filtering}

  Select_Transcript_Classify: ${params.Select_Transcript_Classify}
  Classify_cage_refTSS_bed  : ${params.Transcript_Classify.cage_refTSS_bed}
  Classify_polyA_list       : ${params.Transcript_Classify.polyA_list}

  Select_IsoQuant           : ${params.Select_IsoQuant}
  GENECOD_gtf_db            : ${params.IsoQuant.GENECOD_gtf_db}

  Python_Script_Path        : ${params.Python_Script_Path}
  Select_Pigeon_Classify_Report : ${params.Select_Pigeon_Classify_Report}
  gencode_annotation_sorted_gtf : ${params.Pigeon_Classify_Report.gencode_annotation_sorted_gtf}
  cage_peak_sorted_bed          : ${params.Pigeon_Classify_Report.cage_peak_sorted_bed}
  polyA_list_txt                : ${params.Pigeon_Classify_Report.polyA_list_txt}
  covearge_min_count_10_modified2_sorted_tsv: ${params.Pigeon_Classify_Report.covearge_min_count_10_modified2_sorted_tsv} 
  
The parsing information based your configuration are the following:
  Reference Sequence        : $REFERENCE
  Reference GENECODE        : $GENECODE
 """  
//  exit 0    
}

// --mapper                          Isoseq Mapping tool. Minimap2 (default)
def helpMessage() {
  log.info """
        Welocme to run Nextflow Pipeline Isoseq_analysis.nf 
        Usage:
        A typical command for running the pipeline is as follows:
        nextflow run Isoseq_analysis.nf -profile cluster -w ~/Nextflow_work --outdir ./
        A command with more configurable arguments can be           
        nextflow run Isoseq_analysis.nf -profile cluster -w ~/Nextflow_work --isoseq_fastq_filelist ./Isoseq_fastq_list.txt --project Ptest --genome_build hg38 --Isoform_classification Y --Fusion_transcript_find N --outdir ./

        Configurable arguments:
        --project                         Project name/atlas, a folder that used to distinguish other analysis 
        --isoseq_filelist                 Fastq file list that you want to do Isoseq analysis     
        --outdir                          The directory for the pipeline output           
        --genome_build                    Reference genome version. hg38(default) | hg19 | b37 
        --PRIMERS_FA                      PRIMA_FASTA file for Isoseq_Lima(to remove cDNA primer) and Isoseq_Refine to PolyA tail and artifial concatemers       
        --Select_Isoseq_Lima              Whether select to run Isoseq_lima to remove cDNA Primer. Y(default) | N
        --Select_Isoseq_Refine            Whether select to run Isoseq_Refine to remove PolyA tail and artifical concatemers. Y(default) | N
        --Protocol_Polya                  Whether the FL reads have a poly(A) tail and need to be remove. Y(default) | N
        --Select_Isoseq_Cluster           Whether select to run Isoseq_Cluster. Y(default) | N
        --Select_minimap2                 Whether select to run minimap2 for isoseq mapping and the downstream function modules. Y(default) | N 
        --Select_PBMM2                    Whether select to run PBMM2 for isoseq mapping and the downstream function modules. Y(default) | N 
        --Select_SplitNCigarReads         Whether select to run GATK_SplitNCigarReads and the downstream function modules. Y(default) | N 
        --Select_Deepvariant              Whether select to run Deepvariant SNV calling in the GPU nodes. Y(default) | N 
        --Select_SNV_ANNOVAR              Whether select to run ANNOVAR for annotating the SNV called by Deepvariant. Y(default) | N  
        
        --Select_LongQC_Isoseq            Whether select LongQC. Y(default) | N  

        --Select_LongGF                   Whether select LongGF for Fusion detection. Y(default) | N 
        
        --Select_Isoseq_Collapse          Whether select to run Isoseq_Collapse to use exonic structure to collapse redundant transcripts, and following function modules. Y(default) | N
        --Select_Transcript_Classify      Whether select to run Transcript Classify for transcript classification. Y(default) | N
        --Select_Transcript_Filtering     Whether select to run Transcript filtering after transcript classification. Y(default) | N   
   
        --Select_IsoQuant                 Whether select to run IsoQuant for longread RNAseq Quantification analysis. Y | N
        --Select_Pigeon_Classify_Report   Whether select to run Pigeon for Classification and Report generation after run Isoquant. Y | N

        --help | h                        This usage statement.
       Note:
         All the above arguments can be configured by the command line interface or in the nextflow.config (default)  
        """
}


// Show help message
if ((params.help) || (params.h)) {
    helpMessage()
    exit 0
}

//Parse the input Parameters and configs. 
parse_config_parameters() 
// Display the configuration
DispConfig() 


// Set up the two external input channels for samples (fastq files) and the reference genome 
Isoseq_Sample_Ch    = Channel
            .fromPath(params.isoseq_filelist)
            .splitText()
            .splitCsv(sep: '\t')
//          .view()

//Set up the two related geneome value channels 

GENECODE_GTF       = file(GENECODE)
REFERENCE          = file(REFERENCE)

Isoseq_Sample_Ch.into { Isoseq_Sample_Ch1; Isoseq_Sample_Ch2; Isoseq_Sample_Ch3; Isoseq_Sample_Ch4; Isoseq_Sample_Ch5 } 
/***************************************************************************************************************************************
* This Process is be responsible for iso-seq lima to remove cDNA primers
***************************************************************************************************************************************/
Isoseq_Lima_In_Ch0 = (params.Select_Isoseq_Lima == "Y" )? Isoseq_Sample_Ch1 : Channel.empty()
Isoseq_Lima_In_Ch = (params.Isoseq_DataType== "BAM" )? Isoseq_Lima_In_Ch0 : Channel.empty() 
process Isoseq_Lima {
   publishDir "${params.outdir}/${params.project}/Isoseq_Lima/$SampleName", mode: 'copy', overwrite: true
   input:
   set SampleName, CCS_bam from Isoseq_Lima_In_Ch

   output:
   tuple SampleName, file("${SampleName}.*.bam"), file("${SampleName}.*.bam.pbi"), file("${SampleName}.*.xml"), file("${SampleName}.json"), file("${SampleName}.lima.*")
   set SampleName, file("${SampleName}.*.bam") into Isoseq_Lima_Ch1
   set SampleName, file("${SampleName}.*.bam") into Isoseq_Lima_Ch2
   set SampleName, file("${SampleName}.*.bam") into Isoseq_Lima_Ch3
   set SampleName, file("${SampleName}.*.bam") into Isoseq_Lima_Ch4

   script:
   """
     lima \
     --isoseq \
     --peek-guess \
     --num-threads ${task.cpus} \
     ${CCS_bam} \
     ${params.PRIMERS_FA} \
     ${SampleName}.bam
   """
}

/***************************************************************************************************************************************
* This Process is be responsible for removing polyA tail (if required) and artifical concatemers
***************************************************************************************************************************************/
Isoseq_Refine_In_Ch0 = (params.Select_Isoseq_Lima == "Y" ) ? Isoseq_Lima_Ch1 : Isoseq_Sample_Ch2
Isoseq_Refine_In_Ch1 = (params.Select_Isoseq_Refine == "Y" ) ? Isoseq_Refine_In_Ch0 : Channel.empty()
Isoseq_Refine_In_Ch = (params.Isoseq_DataType == "BAM" ) ? Isoseq_Refine_In_Ch1 : Channel.empty() 
process Isoseq_Refine {
    publishDir "${params.outdir}/${params.project}/Isoseq_Refine/$SampleName", mode: 'copy', overwrite: true
    input:
    set SampleName, flt_bam from Isoseq_Refine_In_Ch  

    output:
    tuple SampleName, file("${SampleName}.refined.bam"), file("${SampleName}.refined.bam.pbi"), file("${SampleName}.refined.consensusreadset.xml"), file("${SampleName}.refined.filter_summary.report.json"), file("${SampleName}.refined.report.csv")
    set SampleName, file("${SampleName}.refined.bam") into Isoseq_Refine_Ch1
    set SampleName, file("${SampleName}.refined.bam") into Isoseq_Refine_Ch2 
    set SampleName, file("${SampleName}.refined.bam") into Isoseq_Refine_Ch3 
    set SampleName, file("${SampleName}.refined.bam"), file("${SampleName}.refined.bam.pbi") into Isoseq_Refine_Ch4  

    script:
    """
      isoseq refine \
      ${POLYA_SETTING} \
      --num-threads ${task.cpus} \
      ${flt_bam} \
      ${params.PRIMERS_FA} \
      ${SampleName}.refined.bam
    """
}

Isoseq_Cluster_In_Ch0 = (params.Select_Isoseq_Lima == "Y" ) ? Isoseq_Lima_Ch2 : Isoseq_Sample_Ch3
Isoseq_Cluster_In_Ch1 = (params.Select_Isoseq_Refine == "Y" ) ? Isoseq_Refine_Ch1 : Isoseq_Cluster_In_Ch0
Isoseq_Cluster_In_Ch2 = (params.Select_Isoseq_Cluster == "Y" ) ? Isoseq_Cluster_In_Ch1 : Channel.empty()
Isoseq_Cluster_In_Ch = (params.Isoseq_DataType == "BAM" ) ? Isoseq_Cluster_In_Ch2 : Channel.empty()
process Isoseq_Cluster {
    publishDir "${params.outdir}/${params.project}/Isoseq_Cluster/$SampleName", mode: 'copy', overwrite: true
    input:
    set SampleName, flnc_bam, flnc_pbi from Isoseq_Cluster_In_Ch

    output:
    tuple SampleName, file("${SampleName}_transcripts.bam"), file("${SampleName}_transcripts.bam.pbi"),file("${SampleName}_transcripts.*.csv")
    set SampleName, file("${SampleName}_transcripts.bam") into Isoseq_Cluster_Ch1
    set SampleName, file("${SampleName}_transcripts.bam") into Isoseq_Cluster_Ch2    

    script:
    """
       isoseq cluster2 \
       -j ${task.cpus} \
        ${params.Isoseq_Cluster.Cluster_Singleton} \
        ${flnc_bam} \
        ${SampleName}_transcripts.bam   
    """       
}

// Isoseq_Refine_Ch.into {Isoseq_Refine_Ch1; Isoseq_Refine_Ch2}

/***************************************************************************************************************************************
* This Process is be responsible for converting the refined.bam into Fastq file 
***************************************************************************************************************************************/
BAM_To_Fasta_In_Ch0= (params.Select_Isoseq_Lima =="Y")? Isoseq_Lima_Ch3 : Isoseq_Sample_Ch4
BAM_To_Fasta_In_Ch1= (params.Select_Isoseq_Refine =="Y")? Isoseq_Refine_Ch2 : BAM_To_Fasta_In_Ch0
BAM_To_Fasta_In_Ch= (params.Select_Isoseq_Cluster =="Y")? Isoseq_Cluster_Ch1 :BAM_To_Fasta_In_Ch1 
process Pacbio_BAM_To_Fasta {
  publishDir "${params.outdir}/${params.project}/Pacbio_BAM_To_Fasta/$SampleName", mode: 'copy', overwrite: true
    input:
    set SampleName, bam_fastq from BAM_To_Fasta_In_Ch
    
    output:
    set SampleName, file("${SampleName}.fast*") into isoseq_fastq_ch

    // make sure that "bam2fastq $bam_fastq -c 9 -u -o ${SampleName}"  will output ${SampleName}.fastq
    // make sure that "bam2fastq $bam_fastq -c 9 -o ${SampleName}"  will output ${SampleName}.fastq.gz
    // bam2fastq $bam_fastq -c 9 -u -o ${SampleName} or bam2fasta $bam_fastq -c 9 -u -o ${SampleName}.  Isoseq-cluster output' ubam only can be converted to .fasta not .fastq 

    script:
    if (file(bam_fastq).getExtension() == "bam") //treat it as Pacbio uBam // Isoseq-cluster output' ubam only can be converted to .fasta not .fastq  
    """
      bam2fasta $bam_fastq -c 9 -u -o ${SampleName}    
    """
    else if (file(bam_fastq).getExtension() == "fastq")   //treat it as .fastq
    """
      ln -s $bam_fastq ${SampleName}.fastq
    """
    else if (file(bam_fastq).getExtension() == ".gz")   //treat it as .fastq.gz
    """
      gzip -d $bam_fastq ${SampleName}.fastq
    else
      error " $bam_fastq is not .unaligned.bam or .fastq or .fastq.gz   "
    """   
}

isoseq_fastq_ch.into { isoseq_fastq_ch1; isoseq_fastq_ch2 }
/*This Process is responsible for QC analysis of the Pacbio Isoseq data*/
LongQC_Isoseq_In_Ch= (params.Select_LongQC_Isoseq=="Y")? isoseq_fastq_ch1 : Channel.empty()
process LongQC_Isoseq {
    publishDir "${params.outdir}/${params.project}/LongQC_Isoseq/$SampleName", mode: 'copy', overwrite: true

    input:
    set SampleName, isoform_fastq from LongQC_Isoseq_In_Ch
    output:
    set SampleName, file("${SampleName}.tar.gz") into LongQC_Isoseq_ch

    script:
    """
     longqc sampleqc \
     -x ${params.LongQC.Sequencing_platform} \
     -t \
     -s ${SampleName} \
     -n ${params.LongQC.NSAMPLE} \
     -p 4 \
     -o ${SampleName} \
     $isoform_fastq
   
     tar -zcvf ${SampleName}.tar.gz ${SampleName}
    """ 
}


/***************************************************************************************************************************************
* This Process is be responsible for iso-seq mapping by minimap2     
***************************************************************************************************************************************/
minimap2_isoseq_fastq_ch= (params.Select_minimap2=="Y")? isoseq_fastq_ch2 : Channel.empty()
process minimap2_isoseq_mapping {
	publishDir "${params.outdir}/${params.project}/minimap2_isoseq_mapping/$SampleName", mode: 'copy', overwrite: true
    input:
    val ref_fa from REFERENCE
    set SampleName, isoform_fastq from minimap2_isoseq_fastq_ch

    output:
//  set SampleName, isoform_fastq, file("${SampleName}.aln.*") into minimap2_mapping_ch
    set SampleName, isoform_fastq, file("${SampleName}.aln.sort.bam"), file("${SampleName}.aln.sort.bam.bai"), file("${SampleName}.aln.sort.sam") into minimap2_mapping_ch

    script:
    """
    minimap2 -t ${task.cpus} -ax splice -uf --secondary=no -C5 $ref_fa $isoform_fastq > ${SampleName}.aln.sam 2> ${SampleName}.aln.log 
    samtools view -Sb ${SampleName}.aln.sam > ${SampleName}.aln.bam 
    samtools sort ${SampleName}.aln.bam > ${SampleName}.aln.sort.bam 
    samtools index ${SampleName}.aln.sort.bam
    samtools view -h ${SampleName}.aln.sort.bam > ${SampleName}.aln.sort.sam
    
    """
}


mapping_ch_LongGF = (params.Select_LongGF =="Y")? minimap2_mapping_ch : Channel.empty()
/*****************************************************************************************************
* This Process is be responsible for detect inter-chrosomoe Fusion by LongGF
***************************************************************************************************/
process LongGF {
  publishDir "${params.outdir}/${params.project}/LongGF/$SampleName", mode: 'copy', overwrite: true
      input:
      set SampleName, isoform_fastq, align_bam, align_bam_bai, align_sam from mapping_ch_LongGF
      val ref_gtf from GENECODE

      output:
      set SampleName, file("${SampleName}.LongGF_Result.txt"), file("${SampleName}.LongGF_Summary.txt") into LongGF_Ch

      script:
      """
       samtools sort \
       -n ${align_bam} \
       -o ${SampleName}.aln.sort_by_name.bam

       LongGF \
       ${SampleName}.aln.sort_by_name.bam \
       ${ref_gtf} \
       ${params.LongGF.min_overlap_len} \
       ${params.LongGF.bin_size} \
       ${params.LongGF.min_map_len} > ${SampleName}.LongGF_Result.txt

       grep "SumGF" ${SampleName}.LongGF_Result.txt > ${SampleName}.LongGF_Summary.txt

      """
}

/*************************************************************************************************************************************
* This Process is be responsible for Iso-Seq mapping by PBMM2  
*  --preset ISOSEQ \
*   --sample $SampleName \
*   --rg '@RG\tID:${SampleName}' \
*   -j ${params.PBMM2_ThreadN} \
**************************************************************************************************************************************/
PBMM2_Isoseq_In_Ch0= (params.Select_Isoseq_Lima == "Y")? Isoseq_Lima_Ch4  : Isoseq_Sample_Ch5 
PBMM2_Isoseq_In_Ch1= (params.Select_Isoseq_Refine == "Y")? Isoseq_Refine_Ch3 : PBMM2_Isoseq_In_Ch0
PBMM2_Isoseq_In_Ch= (params.Select_Isoseq_Cluster == "Y")? Isoseq_Cluster_Ch2 : PBMM2_Isoseq_In_Ch1
process PBMM2_Isoseq_mapping {
   publishDir "${params.outdir}/${params.project}/PBMM2_IsoSeq_mapping", mode: 'copy', overwrite: true

   input:
   val ref_fa from REFERENCE
   set SampleName, isoseq_bam from PBMM2_Isoseq_In_Ch

   output:
   set val(SampleName), file("${SampleName}.aln.bam"), file("${SampleName}.aln.bam.bai") into PBMM2_Isoseq_Ch

   script:
   """
   export TMPDIR=${params.TMP_DIR}
   pbmm2 align \
   $ref_fa \
   $isoseq_bam \
   ${SampleName}.aln.bam \
   --sort \
   --preset ISOSEQ \
   --sample $SampleName \
   -j ${task.cpus} \
   -J 4
   """
}

PBMM2_Isoseq_Ch.into {PBMM2_Isoseq_Ch1; PBMM2_Isoseq_Ch2; PBMM2_Isoseq_Ch3; PBMM2_Isoseq_Ch4}

IsoQuant_Adaptor_In_Ch = (params.Select_IsoQuant =="Y")? PBMM2_Isoseq_Ch1 : Channel.empty()
/*This process function as an adaptor between the upstream PBMM2 and downstream Isoquant, which
 just provide a buffer to collect all of the aligned bam across samples */
process PBMM2_Adaptor_IsoQuant {
  input:  
  set SampleName, file(align_bam), file(align_bam_bai) from  IsoQuant_Adaptor_In_Ch
  
  output:
  val(SampleName) into IsoQuant_SampleName_Ch
  file(align_bam)  into IsoQuant_align_bam_Ch
  file(align_bam_bai) into IsoQuant_align_bam_bai_Ch

  script:
  """
  """
}

/* This proces is for Calling Isoquant ( https://github.com/ablab/IsoQuant ) for longread quantification analysis across samples */
process IsoQuant {
  publishDir "${params.outdir}/${params.project}/IsoQuant", mode: 'copy', overwrite: true 
  input:
  val ref_fa from REFERENCE
  val ref_gtf from GENECODE
  
  val(SampleNames) from IsoQuant_SampleName_Ch.collect()  
  file(align_bams) from IsoQuant_align_bam_Ch.collect()
  file(align_bam_bais) from IsoQuant_align_bam_bai_Ch.collect()

  output:
  path("${params.project}_IsoQuant")
  set file("${params.project}_IsoQuant/*.transcript_models.gtf"), file("${params.project}_IsoQuant/*.discovered_transcript_grouped_counts.tsv") into IsoQuant_Ch  

  // --genedb $ref_gtf
  script:
  """
      isoquant.py \
      --threads ${task.cpus} \
      --reference $ref_fa \
      --genedb ${params.IsoQuant.GENECOD_gtf_db} \
      --bam ${align_bams.collect {" $it"}.join()} \
      --data_type pacbio_ccs \
      -o .  \
      --prefix ${params.project}_IsoQuant \
      --labels ${SampleNames.collect {" $it"}.join()}  
      
      sort -k1,1 -k4,4n ${params.project}_IsoQuant.extended_annotation.gtf > ${params.project}_IsoQuant.extended_annotation.sorted.gtf
      bgzip ${params.project}_IsoQuant.extended_annotation.sorted.gtf
      tabix ${params.project}_IsoQuant.extended_annotation.sorted.gtf.gz 
   """

}

/* This process is for calling a python script to convert the IsoQuant ouput results for the downstrea  */
Isoquant2Pigeon_In_Ch= (params.Select_Pigeon_Classify_Report =="Y")? IsoQuant_Ch : Channel.empty()
process Isoquant2Pigeon {
    publishDir "${params.outdir}/${params.project}/Isoquant2Pigeon", mode: 'copy', overwrite: true

    input:
    set file(transcript_models_gtf), file(transcript_grouped_counts_tsv) from Isoquant2Pigeon_In_Ch         
    output:
    set file(transcript_models_gtf), file("${params.project}.transcript_model_grouped_counts.csv") into Isoquant2Pigeon_Ch
   
    // python /research/groups/cab/projects/Control/common/Bulk_Isoseq_Nextflow/Python_script/isoquant2pigeon.py 
    script:
    """
      python ${params.Python_Script_Path}/isoquant2pigeon.py \
      --gtf ${transcript_models_gtf} \
      --tsv ${transcript_grouped_counts_tsv} \
      --output ${params.project}.transcript_model_grouped_counts.csv
    """   
}

Pigeon_Classify_Report_In_Ch= (params.Select_Pigeon_Classify_Report =="Y")? Isoquant2Pigeon_Ch : Channel.empty()
process Pigeon_Classify_Report {
   publishDir "${params.outdir}/${params.project}/Pigeon_Classify_Report", mode: 'copy', overwrite: true

   input:
   val ref_fa from REFERENCE
   set file(transcript_models_gtf), file(transcript_model_grouped_counts_csv) from Pigeon_Classify_Report_In_Ch   

   output:
   set file("${params.project}.pigeon_junctions.txt"), file("${params.project}.pigeon.report.json"), file("${params.project}.pigeon.summary.txt"), file("${params.project}.pigeon_classification.txt"), file("${params.project}.pigeon_classification.report.txt") into Pigeon_Classify_Report_Ch      

   /*
    --cage-peak $HG38/refTSS_v3.3_human_coordinate.hg38.sorted.bed \
      --poly-a $HG38/polyA.list.txt \
      --coverage $HG38/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified2.sorted.tsv
   */
   script:
   """
     pigeon classify \
      -j ${task.cpus} \
      -o ${params.project}.pigeon \
      ${transcript_models_gtf} \
      ${params.Pigeon_Classify_Report.gencode_annotation_sorted_gtf} \
      ${ref_fa} \
      --flnc ${transcript_model_grouped_counts_csv} \
      --cage-peak ${params.Pigeon_Classify_Report.cage_peak_sorted_bed} \
      --poly-a ${params.Pigeon_Classify_Report.polyA_list_txt} \
      --coverage ${params.Pigeon_Classify_Report.covearge_min_count_10_modified2_sorted_tsv}  

      pigeon report \
      -j 2 \
      ${params.project}.pigeon_classification.txt \
      ${params.project}.pigeon_classification.report.txt
   """
}

/*
process DESeq {
  
}
*/

/***********************************************************************************************
* This Processe is for RNA Variant Calling, is responsible for SplitNCigarReads 
************************************************************************************************/
GATK_SplitNCigarReads_In_Ch= (params.Select_SplitNCigarReads =="Y")? PBMM2_Isoseq_Ch2 : Channel.empty()
process GATK_SplitNCigarReads {
   publishDir "${params.outdir}/${params.project}/GATK_SplitNCigarReads", mode: 'copy', overwrite: true

   input:
   val ref_fa from REFERENCE
   set SampleName, file(PBMM2_Aligned_out_bam), file(PBMM2_Aligned_out_bam_bai) from GATK_SplitNCigarReads_In_Ch

   output:
   set SampleName, file("${SampleName}.splitN.bam"), file("${SampleName}.splitN.bai") into GATK_SplitNCigarRead_Ch

   script:
   """
   export TMPDIR=${params.TMP_DIR}
   gatk SplitNCigarReads \
    --java-options "-Xms8g" \
    --reference ${ref_fa} \
    --input ${PBMM2_Aligned_out_bam} \
    --output ${SampleName}.splitN.bam
   """
}

/*************************************************************************************************************************************
* This Process is be responsible for Isoseq SNV calling
**************************************************************************************************************************************/
Deepvariant_Caller_In_Ch= (params.Select_Deepvariant =="Y")? GATK_SplitNCigarRead_Ch : Channel.empty()
process SNV_Deepvariant {
      publishDir "${params.outdir}/${params.project}/SNV_Deepvariant", mode: 'copy', overwrite: true
      input: 
      set SampleName, file(aln_bam), file(aln_bam_bai) from Deepvariant_Caller_In_Ch     
      val ref_genome from REFERENCE

      output:
      set SampleName, file("${SampleName}.deepvariant.vcf") into SNV_Deepvariant_Caller_Ch
      
      script:
      """
      pbrun deepvariant \
      --num-gpus 2 \
      --tmp-dir ${params.TMP_DIR} \
      --ref ${ref_genome} \
      --in-bam ${aln_bam} \
      --mode pacbio \
      --vsc-min-fraction-snps ${params.Deepvariant.vsc_min_fraction_snp} \
      --vsc-min-fraction-indels ${params.Deepvariant.vsc_min_fraction_indel} \
      --out-variants ${SampleName}.deepvariant.vcf                  
      """  
}

/********************************************************************************************
*  This process is responsible for filtering out Deepvariant called varaints
*********************************************************************************************/
SNV_Deepvariant_Filter_In_Ch= (params.Select_SNV_ANNOVAR=="Y")? SNV_Deepvariant_Caller_Ch: Channel.empty()
process SNV_Variant_Filtering{
      publishDir "${params.outdir}/${params.project}/SNV_Variant_Filtering", mode: 'copy', overwrite: true
      input: 
      set SampleName, file(deepvariant_vcf) from SNV_Deepvariant_Filter_In_Ch     

      output:
      set SampleName, file("${SampleName}.deepvariant.filter.vcf") into SNV_Deepvariant_Filter_Ch
      
      script:
      """
      cat $deepvariant_vcf| grep -v "RefCall" > ${SampleName}.deepvariant.filter.vcf                 
      """  
}

/********************************************************************************************
*  This process is responsible for annotating the RNA-seq SNV varaints
*********************************************************************************************/
process Variant_SNV_ANNOVAR {
   publishDir "${params.outdir}/${params.project}/Variant_SNV_ANNOVAR/${SampleName}", mode: 'copy', overwrite: true
   input:
   set SampleName, file(deepvariant_filter_vcf) from SNV_Deepvariant_Filter_Ch 
   
   output:
   set file("${SampleName}.vcf.ann.avinput"), file("${SampleName}.vcf.ann.*.txt"), file("${SampleName}.vcf.ann.*.vcf.gz"), file("${SampleName}.vcf.ann.*.vcf.gz.tbi") into Variant_ANNOVAR_Ch
   
   script:
   """
  ${params.ANNOVAR.ANNOVAR_PATH}/table_annovar.pl \
   ${deepvariant_filter_vcf} \
   ${params.ANNOVAR.ANNOVAR_DB} \
   -buildver ${params.ANNOVAR.REF_BUILD} \
   -out ${SampleName}.vcf.ann \
   -remove \
   -protocol ${params.ANNOVAR.PROTOCOL} \
   -operation ${params.ANNOVAR.OPERATION} \
   -nastring . -vcfinput

  bgzip ${SampleName}.vcf.ann.${params.ANNOVAR.REF_BUILD}_multianno.vcf
  tabix ${SampleName}.vcf.ann.${params.ANNOVAR.REF_BUILD}_multianno.vcf.gz

  gatk VariantsToTable \
   -V ${SampleName}.vcf.ann.${params.ANNOVAR.REF_BUILD}_multianno.vcf.gz \
   -O ${SampleName}.GATK.annovar.${params.ANNOVAR.REF_BUILD}.tab \
   ${params.ANNOVAR.CONFIG_STANDARD} ${params.ANNOVAR.CONFIG_INFO} ${params.ANNOVAR.CONFIG_FORMAT} ${params.ANNOVAR.CONFIG_ANNOVAR}

   """
}

/***********************************************************************************************
* This Processe is for collapsing the redundant transcripts based on exonic structures 
************************************************************************************************/
PBMM2_Isoseq_Combine_Ch= PBMM2_Isoseq_Ch3.combine(Isoseq_Refine_Ch4, by:0)
Isoseq_Collapse_In_Ch= (params.Select_Isoseq_Collapse =="Y")? PBMM2_Isoseq_Combine_Ch : Channel.empty()
process Isoseq_Collapse {
   publishDir "${params.outdir}/${params.project}/Isoseq_Collapse", mode: 'copy', overwrite: true

   input:
   set SampleName, file(PBMM2_Aligned_out_bam), file(PBMM2_Aligned_out_bam_bai), file(flnt_bam), file(flnt_bam_pbi) from Isoseq_Collapse_In_Ch

   output: 
   set SampleName, file("${SampleName}.gff"),file("${SampleName}.fasta"), file("${SampleName}.flnc_count.txt"), file("${SampleName}.abundance.txt"), file("${SampleName}.read_stat.txt"),file("${SampleName}.group.txt"),file("${SampleName}.report.json" ) into Isoseq_Collapse_Ch
     
   /* export TMPDIR=${params.TMP_DIR}
   isoseq collapse \
     --num-threads ${params.Isoseq_Collapse_threadsN} \
     --min-aln-coverage ${params.Collaspe_Min_Coverage} \
     --min-aln-identity ${params.Collaspe_Min_Identity} \
     --do-not-collapse-extra-5exons \
     ${PBMM2_Aligned_out_bam} \
     ${SampleName}.gff
   */
 
   script: 
   """
   export TMPDIR=${params.TMP_DIR}
   isoseq collapse \
   --num-threads ${task.cpus} \
   --min-aln-coverage ${params.Isoseq_Collapse.Collaspe_Min_Coverage} \
   --min-aln-identity ${params.Isoseq_Collapse.Collaspe_Min_Identity} \
   --do-not-collapse-extra-5exons \
   ${PBMM2_Aligned_out_bam} \
   ${flnt_bam} \
   ${SampleName}.gff
   """
}

/*************************************************************************************************************************************
* This Process is be responsible for using pigeon to classify Isoseq trasncript based on annotation
**************************************************************************************************************************************/
Transcript_Classify_In_Ch= (params.Select_Transcript_Classify =="Y")? Isoseq_Collapse_Ch : Channel.empty()
process Transcript_Classify {
      publishDir "${params.outdir}/${params.project}/Transcript_Classify", mode: 'copy', overwrite: true
      input: 
      set SampleName, file(collapse_gff), file(collapse_fasta), file(flnc_count_txt), file(abundance_txt), file(read_stat_txt), file(group_txt), file(report_json) from Transcript_Classify_In_Ch     
      val ref_fa from REFERENCE
      val genecode_gtf from GENECODE_GTF

      output:
      set SampleName, file("${SampleName}_classification.txt"), file("${SampleName}_junctions.txt"),  file("${SampleName}.report.json"), file("${SampleName}.summary.txt") into Transcript_Classify_Ch
      
      script:
      """
      pigeon sort ${collapse_gff} 
      pigeon classify \
      --num-threads ${task.cpus} \
      ${collapse_gff}.sorted \
      ${genecode_gtf} \
      ${ref_fa} \
      --cage-peak ${params.Transcript_Classify.cage_refTSS_bed} \
      --poly-a ${params.Transcript_Classify.polyA_list} \
      --fl ${flnc_count_txt} \
      --out-prefix ${SampleName}
      """  
}

/********************************************************************************************
*  This process is responsible for filtering out the classified transcripts
*********************************************************************************************/
Transcript_Filtering_In_Ch= (params.Select_Transcript_Filtering=="Y")? Transcript_Classify_Ch: Channel.empty()
process Transcript_Filtering{
      publishDir "${params.outdir}/${params.project}/Transcript_Filtering", mode: 'copy', overwrite: true
      input: 
      set SampleName, file("${SampleName}_classification.txt"), file(junction_txt), file(report_json), file(summary_txt) from Transcript_Filtering_In_Ch     

      output:
      set SampleName, file("${SampleName}_classification.filtered_lite_classification.txt"), file("${SampleName}_classification.filtered_lite_junctions.txt"), file("${SampleName}_classification.filtered_lite_reasons.txt"), file("${SampleName}_classification.filtered.report.json"),file("${SampleName}_classification.filtered.summary.txt") into Transcript_Filtering_Ch
      // set SampleName, file("${SampleName}.filtered_lite_classification.txt"), file("${SampleName}.filtered_lite_junctions.txt"), file("${SampleName}.filtered_lite_reasons.txt"), file("${SampleName}.filtered.report.json"), file("${SampleName}.filtered.summary.txt") into Transcript_Filtering_Ch     
      // ln -s ${classification_txt} ${SampleName}.txt
      // --mono-exon  --skip-junctions \

      script:
      """
      pigeon filter \
      --mono-exon \
      ${SampleName}_classification.txt
      """  
}

/********************************************************************************************
*  This process is responsible for use pbfusion(Pacbio recommended fusion calling method)  
to detect fusion genes.  https://github.com/pacificbiosciences/pbfusion/
*********************************************************************************************/
PBFusion_In_Ch= (params.Select_PBFusion=="Y")? PBMM2_Isoseq_Ch4 : Channel.empty()
process PBFusion {
   publishDir "${params.outdir}/${params.project}/PBFusion/", mode: 'copy', overwrite: true
   input:
   set SampleName, file(Aligned_bam), file(Aligned_bam_bai) from PBFusion_In_Ch

   output:
   set SampleName, file("${SampleName}.pbfusion.breakpoints.groups.bed") into PBFusion_Ch

   script:
   """
    pbfusion discover \
    -t ${task.cpus} \
    --gtf ${params.PBFusion.gencode_annotation_sorted_gtf} \
    --output-prefix ${SampleName}.pbfusion \
    ${Aligned_bam}
   """
}
