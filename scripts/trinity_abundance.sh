#!/usr/bin/env bash

set -eu
bowtie_path="$(readlink -f bin/bowtie1)"
PATH="${bowtie_path}:${PATH}"


 bin/trinity/util/align_and_estimate_abundance.pl\
	--transcripts output/trinity/Trinity.fasta \
	--seqType fq \
	--est_method RSEM \
	--aln_method bowtie \
	--output_dir output/trinity_abundance \
	--prep_reference \
	--SS_lib_type RF \
	--thread_count 1 \
	--trinity_mode \
	--left output/bbduk_trim/abdo_r1.fq.gz,output/bbduk_trim/thorax_r1.fq.gz \
    --right output/bbduk_trim/abdo_r2.fq.gz,output/bbduk_trim/thorax_r2.fq.gz \