params.reference = null
params.outBAM = "$projectDir/BAMs"
params.outVCF = "$projectDir/VCFs"
//A CSV of reads to process is required. This is provided in the GitHub repo.
params.index = "$projectDir/Reads.csv"
//Threads is hardcoded. Changing the threads after a run will prevent the resume from working.
params.threads = 8

//Reference is required on the fly
if (params.reference == null) error "Please specify a reference genome with --reference"

//Scale the number of alignment jobs to the computer being used.
maxcpus = Runtime.runtime.availableProcessors()
maxJobs = maxcpus / params.threads

//Print some log information
log.info """\
	ReadAlign-NF Pipeline
	=====================
	assembly: ${params.reference}
	BAMs output: ${params.outBAM}
	VCFs output: ${params.outVCF}
	"""
	.stripIndent()

//Index the reference assembly wtih BWA
process Index {
	conda 'bwa=0.7.15'
	input:
	path reference

	output:
	path "${reference}.{'',amb,ann,bwt,pac,sa}"

	"""
	bwa index $reference
	"""
}

//Align paired reads to the indexed reference assembly
process Align {
	//Modify maxForks if you want to manually set the number of BWA jobs to run.
	maxForks maxJobs
	conda 'bwa=0.7.15 samtools=1.17'
	input:
	path reference
	tuple val(SampleID), val(FCID), file(Read1), file(Read2)

	output:
	tuple val(SampleID), path("${SampleID}.${FCID}.bam")

	script:
	def idxbase = reference[0].baseName

	"""
	bwa mem -t ${params.threads} $idxbase $Read1 $Read2 | samtools sort -o ${SampleID}.${FCID}.bam -
	"""
}

//Check the sample ID for duplicates and merge them. If only a single entry found, the merge command in replaced by a sym link.
//The process will output the BAM files into the BAM directory.
process Merge {
	debug true
	tag { SampleID }
	conda 'samtools=1.17 bedtools=2.31.0 freebayes=1.3.6'
	input:
	tuple val(SampleID), val(FCID), path(SampleBAM)

//If copy is preferred add, "mode: 'copy', overwrite: true" below
	publishDir params.outBAM, saveAs: { filename -> filename.equals('cov.txt') ? null : filename }
	output:
	path "${SampleID}*merge.bam"
	tuple path("${SampleID}*merge.bam"), val(SampleID), path("cov.txt"), emit: tuple

//The code is shell so that environmentatl variable can be set
	shell:
	'''
	if [[ ! -f !{PWD}/coverage.txt ]] ; then echo -e "SampleID\tBAMs\tCoverage\tVCF\tBAM\tFlowCellIDs" > !{PWD}/coverage.txt ; fi
	length=$(echo !{SampleBAM} | tr ' ' '\n' | wc -l)
	if [[ $length == 1 ]]; then
		ln -s !{SampleBAM} !{SampleID}.1FCs_merge.bam
	else
		samtools merge !{SampleID}.${length}FCs_merge.bam !{SampleBAM}
	fi
	bedtools genomecov -ibam !{SampleID}.${length}FCs_merge.bam > !{SampleID}.${length}FCs_merge.GenCov
	cov=$(sort -grk5,5 !{SampleID}.${length}FCs_merge.GenCov | head -n 1 | cut -f 2)
	echo $cov > cov.txt
	echo -e "!{SampleID}\t$length\t$cov\t!{params.outVCF}/!{SampleID}.vcf\t!{params.outBAM}/!{SampleID}.${length}FCs_merge.bam\t!{FCID}" >> !{PWD}/coverage.txt
	'''
}

//Run FreeBayes. This process will output the VCF into the VCF directory
process FreeBayes {
	debug true
	conda 'freebayes=1.3.6'
	input:
	tuple path(SampleBAM), val(SampleID), path(cov)
	path reference

	publishDir params.outVCF
	output:
	path "${SampleID}*.vcf"

	"""
	cov="\$(cat ${cov})"
	freebayes -f $reference -C2 $SampleBAM > ${SampleID}_\${cov}x.vcf
	"""
}

//Run the workflow
workflow {
	ref = file(params.reference)
	index_ch = Index(ref)
//Code to split the CSV
	Channel.fromPath(params.index) \
	| splitCsv(header:true) \
	| map { row -> tuple(row.SampleID, row.FCID, file(row.Read1), file(row.Read2)) } \
	| set { read_ch }

	align_ch = Align(index_ch, read_ch)
//Code to merge channels with the same SampleID
	merge_ch = read_ch \
	| join(align_ch, by: [0] ) \
	| map { SampleID, FCID, Read1, Read2, SampleBAM -> tuple( SampleID, FCID, SampleBAM ) } \
	| groupTuple() \
	| map { SampleID, FCID, SampleBAM -> tuple( SampleID.toString(), FCID.flatten(), SampleBAM.flatten() ) } \
	| set { merge_in }
	Merge(merge_in)
//Code to rename tuple headers
	FB_ch = Merge.out.tuple \
	| map { SampleBAM, SampleID, cov -> tuple( SampleBAM, SampleID, cov ) }
	FreeBayes(FB_ch, ref)
}
