"""Extract all the sequences surrounding a given genomic reference
point (TSS, TTS, start codon, ...) for a given genome, and run
position-analysis to discover k-mers characterized by specific
positional profiles relative to this reference point.

"""



ORG=chlamydomonas_reinhardtii

rule all_seq:
    	"""Retrieve all 
	Required parameters:
		config['qsub']
		config["genome"]["size"]

	Usage: 
		PEAKS_MACS2 = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/macs2/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_macs2_peaks.narrowPeak",
	               zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])
	"""
        input: 	treatment="{result_dir}/{treatment}/{treatment}_{aligner}.bed", \
            	control="{result_dir}/{control}/{control}_{aligner}.bed"
        params: outdir="{result_dir}/{treatment}_vs_{control}/macs2/", \
            name="{treatment}_vs_{control}_{aligner}_qval{qval}_macs2", \
            qval = config["macs2"]["qval"], \
            call_summits = config["macs2"]["call_summits"], \
            keep_dup = config["macs2"]["keep_dup"], \
            band_width = config["macs2"]["band_width"], \
            mfold_min = config["macs2"]["mfold_min"], \
            mfold_max = config["macs2"]["mfold_max"], \
            other_options = config["macs2"]["other_options"], \
            genome_size = config["genome"]["size"], \
            qsub=config["qsub"] \
                  + " -q long -e {result_dir}/{treatment}_vs_{control}/macs2/{treatment}_vs_{control}_{aligner}_qval{qval}_macs2_qsub.err" \
                  + " -o {result_dir}/{treatment}_vs_{control}/macs2/{treatment}_vs_{control}_{aligner}_qval{qval}_macs2_qsub.out" 
	output: narrowpeaks = "{result_dir}/{treatment}_vs_{control}/macs2/{treatment}_vs_{control}_{aligner}_qval{qval}_macs2_peaks.narrowPeak", \
                peaks_bed = "{result_dir}/{treatment}_vs_{control}/macs2/{treatment}_vs_{control}_{aligner}_qval{qval}_macs2_peaks.bed", \
                peak_len = "{result_dir}/{treatment}_vs_{control}/macs2/{treatment}_vs_{control}_{aligner}_qval{qval}_macs2_peaklen.tab", \
                summits = "{result_dir}/{treatment}_vs_{control}/macs2/{treatment}_vs_{control}_{aligner}_qval{qval}_macs2_summits.bed"
        log: "{result_dir}/{treatment}_vs_{control}/macs2/{treatment}_vs_{control}_{aligner}_qval{qval}_macs2.log"
        benchmark: "{result_dir}/{treatment}_vs_{control}/macs2/{treatment}_vs_{control}_{aligner}_qval{qval}_macs2_benchmark.json"
        
	shell: "(macs2 \
callpeak -t {input.treatment} -c {input.control} --gsize {params.genome_size} \
--qvalue {params.qval} \
--keep-dup {params.keep_dup} \
--bdg \
--bw {params.band_width} \
--mfold {params.mfold_min} {params.mfold_max} \
--outdir {params.outdir} --name {params.name} \
{params.call_summits} {params.other_options} ; \
convert-features -from bed3col -to bed3col -i {output.narrowpeaks} -o {output.peaks_bed}; \
sequence-lengths -i {output.peaks_bed} -in_format bed | classfreq -ci 10 -v 1 -o {output.peak_len}) &> {log}"

