rule to_label_coolers:
    output:
        paths.matrix.cool_label_split,
    input:
        pairs=paths.pairs.label_split,
    params:
        chromsizes=config['chrom_sizes']
    conda:
        PORE_C_CONDA_FILE
    shell:
        "cooler cload pairs {params.chromsizes}:1000 {input.pairs} {output} --chrom1 2 --chrom2 4 --pos1 3 --pos2 5 "

rule zoomify_and_balance_label_coolers:
    input:
        paths.matrix.cool_label_split,
    output:
        paths.matrix.mcool_label_split,
    params:
        resolutions=",".join(map(str, config["matrix_resolutions"]["zoomify"])),
    conda:
        "../envs/cooler.yml"
    threads: 10
    shell:
        "cooler zoomify -n {threads} -r {params.resolutions} -o {output} {input} --balance"

# rule transfer_weights:
#     input:
#         cis_and_trans=expand(paths.matrix.mcool_label_split, contact_type=["cis_and_trans"], allow_incomplete=True),
#         cis=expand(paths.matrix.mcool_label_split, contact_type=["cis"], allow_incomplete=True),
#         trans=expand(paths.matrix.mcool_label_split, contact_type=["trans"], allow_incomplete=True)
#     output:
#         cis_and_trans=expand(paths.matrix.mcool_split_weights_transfered, contact_type=["cis_and_trans"], allow_incomplete=True),
#         cis=expand(paths.matrix.mcool_split_weights_transfered, contact_type=["cis"], allow_incomplete=True),
#         trans=expand(paths.matrix.mcool_split_weights_transfered, contact_type=["trans"], allow_incomplete=True)
#     conda:
#     shell:
#         "python bin/transfer_weights"
    