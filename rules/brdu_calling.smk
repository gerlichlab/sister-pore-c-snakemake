container: "docker://gerlichlab/dnascent-docker:latest"

rule generate_brdu_index:
    output:
        paths.brdu_calling.index
    input:
        fast5=paths.fast5.fast5,
        seq_summary=paths.fast5.seq_summary
    shell:
        "/dnascent/DNAscent/bin/DNAscent index -f {input.fast5} -s {input.seq_summary} -o {output}"

rule call_brdu:
    output:
        paths.brdu_calling.detect
    input:
        index=paths.brdu_calling.index,
        mapping=paths.mapping.coord_sorted_bam,
        refgenome=paths.refgenome.fasta
    shell:
        "/dnascent/DNAscent/bin/DNAscent detect -b {input.mapping} -r {input.refgenome} -i {input.index} -o {output} -q 30 -l 500"