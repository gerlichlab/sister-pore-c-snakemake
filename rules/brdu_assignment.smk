rule make_label_library:
    input:
        expand_basecall_batches(paths.brdu_calling.detect)
    output:
        paths.brdu_calling.label_library
    conda:
        PORE_C_CONDA_FILE
    shell:
        "python bin/make_label_library.py --input {input} --output {output}"


rule assign_pairs:
    input:
        contacts=paths.merged_contacts.contacts,
        pairs=paths.pairs.pairs,
        label_library=paths.brdu_calling.label_library
    output:
        paths.brdu_calling.assigned_pairs
    conda:
        "../envs/assign_pairs.yml"
    shell:
        "python bin/assign_brdu_pairs.py --pairs {input.pairs} --label_lib {input.label_library} --output {output}"
    