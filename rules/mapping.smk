rule align_bwa:
    output:
        bam=temp(paths.mapping.coord_sorted_bam),
        bai=temp(paths.mapping.coord_sorted_bai),
    input:
        fastq=paths.basecall.fastq,
        refgenome=paths.refgenome.fasta,
        bwa_index=paths.refgenome.bwt,
    params:
        cli_opts=config["software"]["bwa"]["cli_opts"],
        memory=config["software"]["sort"]["memory_per_thread"],
        sort_threads=config["software"]["sort"]["threads"],
    threads: config["software"]["bwa"]["threads"]
    # conda:
    #     "../envs/bwa.yml"
    container: "docker://gerlichlab/sister-pore-c-docker:bwa"
    log:
        to_log(paths.mapping.coord_sorted_bam),
    benchmark:
        to_benchmark(paths.mapping.coord_sorted_bam)
    shell:
        "( bwa {params.cli_opts} -t {threads} "
        "{input.refgenome}  {input.fastq} "
        " | pore_c alignments reformat-bam - - "
        " | samtools sort -O bam -m {params.memory} -@ {params.sort_threads} -o {output.bam} -) 2>{log} ;"
        " samtools index {output.bam} 2>{log} "


def is_phased(wildcards):
    vcf = lookup_value("vcf_path", mapping_df)(wildcards).strip()
    if vcf:
        return True
    else:
        return False


rule haplotag:
    output:
        paths.mapping.haplotagged_aligns,
    input:
        bam=paths.mapping.coord_sorted_bam,
        bai=paths.mapping.coord_sorted_bai,
        refgenome=paths.refgenome.fasta_unzipped,
    params:
        vcf=lookup_value("vcf_path", mapping_df),
        is_phased=is_phased,  #conda: "../envs/whatshap.yml"
    log:
        to_log(paths.mapping.haplotagged_aligns),
    benchmark:
        to_benchmark(paths.mapping.haplotagged_aligns)
    container: "docker://gerlichlab/sister-pore-c-docker:latest"
    wrapper:
        WRAPPER_PREFIX + "/whatshap/haplotag"


rule create_alignment_table:
    input:
        bam=paths.mapping.coord_sorted_bam,
        alignment_haplotypes=paths.mapping.haplotagged_aligns,
    output:
        pq=paths.align_table.alignment,
    log:
        to_log(paths.align_table.alignment),
    benchmark:
        to_benchmark(paths.align_table.alignment)
    threads: config["software"]["pore_c"]["create_alignment_table"]["threads"]
    container: "docker://gerlichlab/sister-pore-c-docker:pore-c"
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        "alignments create-table {input.bam} {output.pq} --alignment-haplotypes {input.alignment_haplotypes} 2>{log}"


rule assign_fragments:
    input:
        align_table=paths.align_table.alignment,
        fragments_table=paths.virtual_digest.fragments,
    output:
        pq=paths.align_table.pore_c,
    log:
        to_log(paths.align_table.pore_c),
    benchmark:
        to_benchmark(paths.align_table.pore_c)
    threads: config["software"]["pore_c"]["create_alignment_table"]["threads"]
    container: "docker://gerlichlab/sister-pore-c-docker:pore-c"
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        "alignments assign-fragments {input.align_table} {input.fragments_table} {output} 2>{log}"


rule remove_index_from_readname:
    output:
        bam=paths.mapping.coord_sorted_bam_wo_index,
        index=paths.mapping.coord_sorted_bai_wo_index
    input:
        mapping=paths.mapping.coord_sorted_bam
    container: "docker://gerlichlab/sister-pore-c-docker:latest"
    shell:
        "python scratch/remove_idx_from_read_name.py -i {input.mapping} > {output.bam}; samtools index {output.bam}"

rule to_contacts:
    input:
        align_table=paths.align_table.pore_c,
    output:
        contacts=paths.contacts.contacts,
    params:
        prefix=to_prefix(paths.contacts.contacts),
    log:
        to_log(paths.contacts.contacts),
    benchmark:
        to_benchmark(paths.contacts.contacts)
    container: "docker://gerlichlab/sister-pore-c-docker:pore-c"
    threads: 1
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        "alignments to-contacts {input} {output.contacts} 2>{log}"


def expand_basecall_batches(path):
    def inner(wildcards):
        # here be magic
        _ = checkpoints.import_basecalls.get(**wildcards).output[0]
        glob_path = expand(paths.basecall.fastq, **wildcards, allow_missing=True)
        _ = glob_wildcards(glob_path[0])
        batch_id = sorted([b for b in _.batch_id if b != "fail"], key=lambda x: int(x.replace("batch", "")))
        res = expand(
            str(path),
            enzyme=[wildcards.enzyme],
            refgenome_id=[wildcards.refgenome_id],
            phase_set_id=[wildcards.phase_set_id],
            run_id=[wildcards.run_id],
            batch_id=batch_id,
        )
        return res

    return inner


rule create_contact_fofn:
    output:
        paths.merged_contacts.fofn,
    input:
        expand_basecall_batches(paths.contacts.contacts),
    run:
        with open(output[0], "w") as fh:
            fh.write("{}\n".format("\n".join(input)))


rule merge_contact_files:
    output:
        directory(paths.merged_contacts.contacts),
    input:
        paths.merged_contacts.fofn,
    log:
        to_log(paths.merged_contacts.contacts),
    benchmark:
        to_benchmark(paths.merged_contacts.contacts)
    threads: 4
    container: "docker://gerlichlab/sister-pore-c-docker:pore-c"
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        "contacts merge {input} {output} --fofn"


rule summarise_contacts:
    output:
        pq=paths.merged_contacts.concatemers,
        csv=paths.merged_contacts.concatemer_summary,
    input:
        contacts=ancient(paths.merged_contacts.contacts),
        read_summary=ancient(paths.basecall.summary),
    params:
        metadata=lookup_json(["run_id", "enzyme", "biospecimen", "refgenome_id", "phase_set_id"], mapping_df),
    log:
        to_log(paths.merged_contacts.concatemers),
    benchmark:
        to_benchmark(paths.merged_contacts.concatemers)
    threads: 10
    container: "docker://gerlichlab/sister-pore-c-docker:pore-c"
    shell:
        "pore_c {DASK_SETTINGS} --dask-num-workers {threads} "
        "contacts summarize {input.contacts} {input.read_summary} {output.pq} {output.csv} "
        "--user-metadata '{params.metadata}' 2>{log}"
