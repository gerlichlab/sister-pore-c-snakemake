rule make_label_library:
    input:
        expand_basecall_batches(paths.brdu_calling.detect)
    output:
        paths.brdu_calling.label_library
    conda:
        PORE_C_CONDA_FILE
    shell:
        "python bin/make_label_library.py --input {input} --output {output}"


rule assign_porec_fragments:
    input:
        fragments=paths.align_table.pore_c,
        label_library=paths.brdu_calling.label_library
    output:
        paths.align_table.annotated_pore_c
    log:
        to_log(paths.align_table.annotated_pore_c)
    conda:
        "../envs/spoc.yml"
    shell:
        "spoc annotate {input.fragments} {input.label_library} {output}"


rule assign_pairs:
    input:
        contacts=paths.merged_contacts.contacts,
        pairs=paths.pairs.pairs,
        label_library=paths.brdu_calling.label_library
    output:
        paths.pairs.assigned_pairs
    conda:
        "../envs/assign_pairs.yml"
    shell:
        "python bin/assign_brdu_pairs.py --pairs {input.pairs} --contacts {input.contacts} --label_lib {input.label_library} --output {output}"

rule split_pairs:
    input:
        paths.pairs.assigned_pairs
    output:
        expand(paths.pairs.label_split, label_type=["labelled_only", "all_reads"], contact_type=["cis", "trans", "cis_and_trans"], allow_missing=True)
    conda:
        "../envs/assign_pairs.yml"
    shell:
        "python bin/split_assigned_pairs.py --pairs {input}"
    