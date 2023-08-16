rule to_label_coolers:
    output:
        paths.matrix.cool_label_split,
    input:
        pairs=paths.pairs.label_split,
    params:
        chromsizes=config['chrom_sizes']
    container: "docker://gerlichlab/sister-pore-c-docker:latest"
    shell:
        "cooler cload pairs {params.chromsizes}:1000 {input.pairs} {output} --chrom1 2 --chrom2 4 --pos1 3 --pos2 5 "

rule zoomify_and_balance_label_coolers:
    input:
        paths.matrix.cool_label_split,
    output:
        paths.matrix.mcool_label_split,
    params:
        resolutions=",".join(map(str, config["matrix_resolutions"]["zoomify"])),
    container: "docker://gerlichlab/sister-pore-c-docker:latest"
    threads: 10
    shell:
        "cooler zoomify -n {threads} -r {params.resolutions} -o {output} {input} --balance"

rule transfer_weights:
    input:
        cis_and_trans=expand(paths.matrix.mcool_label_split, contact_type=["cis_and_trans"], allow_missing=True),
        cis=expand(paths.matrix.mcool_label_split, contact_type=["cis"], allow_missing=True),
        trans=expand(paths.matrix.mcool_label_split, contact_type=["trans"], allow_missing=True)
    output:
        cis_and_trans=expand(paths.matrix.mcool_split_weights_transfered, contact_type=["cis_and_trans"], allow_missing=True),
        cis=expand(paths.matrix.mcool_split_weights_transfered, contact_type=["cis"], allow_missing=True),
        trans=expand(paths.matrix.mcool_split_weights_transfered, contact_type=["trans"], allow_missing=True)
    params:
        resolutions=",".join(map(str, config["matrix_resolutions"]["zoomify"])),
    container: "docker://gerlichlab/sister-pore-c-docker:latest"
    shell:
        """
        python bin/transfer_weights.py --input_cis_and_trans {input.cis_and_trans}\
                                    --input_cis {input.cis}\
                                    --input_trans {input.trans}\
                                    --resolutions {params.resolutions}\
                                    --output_cis_and_trans {output.cis_and_trans}\
                                    --output_cis {output.cis}\
                                    --output_trans {output.trans}
        """
    