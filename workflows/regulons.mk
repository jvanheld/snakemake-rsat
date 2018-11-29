################################################################
## This makefile collects upstream sequences of regulons and also
## random clusters of the same sime.
##
##  Add -j #threads to run in parallel.
##
## Tested with some maize clusters from http://www.ncbi.nlm.nih.gov/pubmed/25918418
## Bruno Contreras Moreira, Jaime Castro Mondragon, Jacques Van Helden

include ${RSAT}/makefiles/util.mk
MAKEFILE=regulons.mk

## Define parameters
SPECIES=zea_mays
ORTH_SPECIES=sorghum_bicolor
ORTH_SPECIES_FULL=Sorghum_bicolor.Sorbi1.29
ORTHIDENT=70
FROM=-1000
TO=+200
RNDSAMPLES=50#how many random replicates to be generated / regulon
MEMEMKORD=2#for MEME bg

# input folders, one per regulons
REGULONS=ABI4 E2F ERF11 MYB59 WRI1

#GBF6 GLK2 RS2 CIB1 ERF11 GLK1 HOX16 MYB59 WRI1

list_param:
	@echo
	@echo "Parameters"
	@echo " SPECIES         ${SPECIES} "
	@echo " ORTH_SPECIES    ${ORTH_SPECIES}"
	@echo " ORTHIDENT       ${ORTHIDENT}"
	@echo " FROM            ${FROM}"
	@echo " TO              ${TO}"
	@echo " RNDSAMPLES      ${RNDSAMPLES}"
	@echo " MEMEMKORD       ${MEMEMKORD}"


## Obtain upstream sequences of all regulons, convert consensus and get orths
## uses iupac2meme from MEME suite
.PHONY: get-seqs $(REGULONS)
get-seqs: $(REGULONS) 
$(REGULONS): list_param 
	@echo Retrieving upstream sequences of $@
	@retrieve-seq -org ${SPECIES} -feattype gene -from ${FROM} -to ${TO} \
		-noorf -i $@/regulon$@.txt -label id -o $@/regulon$@.raw.fna
	@retrieve-seq -org ${SPECIES} -feattype gene -from ${FROM} -to ${TO} \
		-noorf -rm -i $@/regulon$@.txt -label id -o $@/regulon$@.rm.fna
	@echo
	@echo Converting cognate consensus of $@ 
	@perl -ne 'print `iupac2meme $$_`' $@/consensus$@.txt > $@/consensus$@.meme
	@convert-matrix -i $@/consensus$@.meme -from meme -to tf -o $@/consensus$@.tf
	@echo
	@echo Retrieving orthologues of $@
	@get-orthologs-compara -i $@/regulon$@.txt -ref_org ${ORTH_SPECIES} \
		-ident_query ${ORTHIDENT} -o $@/orthologues.txt
	@retrieve-seq -org ${ORTH_SPECIES_FULL} -feattype gene -from ${FROM} -to ${TO} \
		-noorf -i $@/orthologues.txt -label id -o $@/orthologues.raw.fna
	@retrieve-seq -org ${ORTH_SPECIES_FULL} -feattype gene -from ${FROM} -to ${TO} \
		-noorf -rm -i $@/orthologues.txt -label id -o $@/orthologues.rm.fna
	@cat $@/regulon$@.raw.fna $@/orthologues.raw.fna > \
		$@/regulon$@.orthologues.raw.fna
	@echo
	@echo Creating ${RNDSAMPLES} random clusters based on $@
	@for r in `seq 1 ${RNDSAMPLES}`; do \
		echo Replicate $$r; \
		NSEQS=`wc -l $@/regulon$@.txt`; \
		RNDREGULON=$@/random$@-$$r.txt; \
        RNDRAW=$@/random$@-$$r.raw.fna; \
        RNDMSK=$@/random$@-$$r.rm.fna; \
        random-genes -n $$NSEQS} -org ${SPECIES} -feattype gene -g 1 -o $$RNDREGULON; \
        retrieve-seq -org ${SPECIES} -feattype gene -from ${FROM} -to ${TO} \
        	-noorf -i $$RNDREGULON -label id -o $$RNDRAW; \
        retrieve-seq -org ${SPECIES} -feattype gene -from ${FROM} -to ${TO} \
        	-noorf -rm -i $$RNDREGULON -label id -o $$RNDMSK; \
	done

## Obtain upstream sequences of all genes in a genome 
## To be used by peak-motifs ans MEME
BGMSK=${SPECIES}-upstream.rm.fna
BGMSKMDL=${SPECIES}-upstream.rm.markov${MEMEMKORD}
background: get-seqs 
	@echo
	@retrieve-seq -org ${SPECIES} -feattype gene -from ${FROM} -to ${TO} \
        -noorf -rm -all -label id -o ${BGMSK}
	@fasta-get-markov -m ${MEMEMKORD} < ${BGMSK} > ${BGMSKMDL}

## Check clusters for missing sequences
check-seqs:
	grep WARNING */*.fna || true
	@echo
	@echo Run \'make -f ${MAKEFILE} fix-seq\' to remove those WARNINGs;\

all: get-seqs background

fix-seqs:
	@perl -i -ne "print unless /^;WARNING/" */*.fna

clean:
	rm -f */random* */*orthologues* */*fna */*meme *upstream.rm* 

