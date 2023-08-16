rule generate_brdu_index:
    output:
        paths.brdu_calling.index
    input:
        fast5=paths.fast5.fast5,
        seq_summary=paths.fast5.seq_summary
    container: "docker://gerlichlab/dnascent-docker:latest"
    log:
        to_log(paths.brdu_calling.index)
    shell:
        "/dnascent/DNAscent/bin/DNAscent index -f {input.fast5} -s {input.seq_summary} -o {output}"

rule call_brdu:
    output:
        paths.brdu_calling.detect
    input:
        index=paths.brdu_calling.index,
        mapping=paths.mapping.coord_sorted_bam_wo_index,
        refgenome=paths.refgenome.fasta
    container: "docker://gerlichlab/dnascent-docker:latest"
    threads: 20
    shell:
        "/dnascent/DNAscent/bin/DNAscent detect -b {input.mapping} -r {input.refgenome} -i {input.index} -o {output} -q 20 -l 25 -t {threads} || touch {output}"