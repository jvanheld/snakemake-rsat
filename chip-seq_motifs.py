"""

Workflow applying various motif discovery algorithms in a peak set from ChIP-seq.

"""

#QSUB_PARAM = " -V -m a "
GENOME="mm9"

FACTORS="CEBPA HNF4A".split()
SWEMBL_R="0.05 0.02 0.01".split()

rule all:
  """Run the workflow on each peak set.
  """
  input:expand("test/fg-chip-seq/SWEMBL_mmus_{factor}_vs_mmus_Input_peaks_R{swembl_r}_nof.fasta", factor=FACTORS, swembl_r=SWEMBL_R)

rule fasta_from_bed:
    """Fetch sequences from UCSC to obtain a fasta file from a set of
    genomic coordinates described in a bed file.

    Example: 
    mkdir -p test/fetch_seq; cd test/fetch_seq; wget http://pedagogix-tagc.univ-mrs.fr/rsat/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed; cd -
    snakemake --snakefile ${RSAT}/snakemake_files/chip-seq_motifs.py test/fetch_seq/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.fasta

    """
    input: "{peaks}.bed"
    output: "{peaks}.fasta"
    log: "{peaks}_fetch-seq.log"
    benchmark: "{peaks}_fetch-seq_benchmark.json"
    params:
#    params: qsub = QSUB_PARAM + " -q short -e {reads}_bam_to_bed_qsub.err -o {reads}_bam_to_bed_qsub.out"
    shell:"fetch-sequences -i {input} -header_format galaxy -genome {GENOME} -o {output} 2> {log} "
