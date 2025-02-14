rule generate_brdu_index:
    output:
        paths.brdu_calling.index
    input:
        pod5=paths.pod5.pod5,
    # container: "library://mboemo/dnascent/dnascent:4.0.3"
    log:
        to_log(paths.brdu_calling.index)
    shell:
        "/groups/gerlich/sequencing_data/159/DNAscent/bin/DNAscent index -f {input.pod5} -o {output}"
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
        "/groups/gerlich/sequencing_data/159/DNAscent/bin/DNAscent detect -b {input.mapping} -r {input.refgenome} -i {input.index} -o {output} -q 10 -l 200 -t {threads} || touch {output} && touch {output}.err"