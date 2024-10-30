rule generate_brdu_index:
    output:
        paths.brdu_calling.index
    input:
        fast5=paths.fast5.fast5,
        seq_summary=paths.fast5.seq_summary
    # container: "library://mboemo/dnascent/dnascent:4.0.3"
    log:
        to_log(paths.brdu_calling.index)
    shell:
        "singularity run library://mboemo/dnascent/dnascent:4.0.3 index -f {input.fast5} -s {input.seq_summary} -o {output}"
        #TODO: rename fast5 to pod5

rule call_brdu:
    output:
        paths.brdu_calling.detect
    input:
        index=paths.brdu_calling.index,
        mapping=paths.mapping.coord_sorted_bam_wo_index,
        refgenome=paths.refgenome.fasta
    # container: "library://mboemo/dnascent/dnascent:4.0.3"
    threads: 20
    shell:
        "singularity run library://mboemo/dnascent/dnascent:4.0.2 detect -b {input.mapping} -r {input.refgenome} -i {input.index} -o {output} -q 30 -l 500 -t {threads} || touch {output} && touch {output}.err"