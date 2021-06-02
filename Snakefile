configfile: "config.json"

rule all:
	input:
		expand("vcf_files/{genome}_{sample}_done.txt", genome=config["genomes"], sample=config["samples"])

rule bcf_filter:
	input:
		"vcf_files/{genome}_{sample}_filtered_multiallelic.bcf"
	output:
		"vcf_files/{genome}_{sample}_filtered_multiallelic_q20.vcf"
	threads: 1
	conda:
		"envs/whatshap.yaml"

	shell:
		"bcftools view --threads 1 --apply-filters PASS -Ov {input} > {output}"


rule bedtools_cds_only:
	input:
		vcf="vcf_files/{genome}_{sample}_filtered_multiallelic_q20.vcf",
		cds="/data/coding/annotations/ITAG4.0_cds_models.gff"
	output:
		"vcf_files/{genome}_{sample}_filtered_multiallelic_q20_ITAG40_cds.vcf"
	threads: 1
	conda:
		"envs/whatshap.yaml"

	shell:
		"grep $'#' {input.vcf} > {output} && bedtools intersect -a {input.vcf} -b {input.cds} >> {output}"

rule bcftools_rm_dupl:
	input:
		"vcf_files/{genome}_{sample}_filtered_multiallelic_q20_ITAG40_cds.vcf"
	output:
		"vcf_files/{genome}_{sample}_filtered_multiallelic_q20_ITAG40_cds.bcf"
	threads: 1
	conda:
		"envs/whatshap.yaml"
	shell:
		"bcftools norm --rm-dup both --threads 1 -Ob {input} > {output}"
	

rule whatshap_phasing:
	input:
		ref="ref_genome/{genome}.fa",
		bcf="vcf_files/{genome}_{sample}_filtered_multiallelic_q20_ITAG40_cds.bcf",
		bam="bam_files/{genome}_{sample}_fixmate_sorted.bam"
	output:
		"phased/{genome}_{sample}_phased.vcf"
	threads: 1
	conda:
		"envs/whatshap.yaml"
	shell:
		"whatshap phase --indels -o {output} --reference {input.ref} {input.bcf} {input.bam}"

rule bcftools_compress:
	input:
		"phased/{genome}_{sample}_phased.vcf"
	output:
		"phased/{genome}_{sample}_phased.bcf"
	threads: 1
	conda:
		"envs/whatshap.yaml"
	shell:
		"bcftools view -Ob {input} > {output}"

rule clear_up:
	input:
		q20="vcf_files/{genome}_{sample}_filtered_multiallelic_q20.vcf",
		shake="phased/{genome}_{sample}_phased.bcf",
		vcf="phased/{genome}_{sample}_phased.vcf",
		itag_vcf="vcf_files/{genome}_{sample}_filtered_multiallelic_q20_ITAG40_cds.vcf"

	output:
		"vcf_files/{genome}_{sample}_done.txt"
	threads: 1
	shell:
		"rm {input.q20} && rm {input.vcf} && rm {input.itag_vcf} && touch {output}"

