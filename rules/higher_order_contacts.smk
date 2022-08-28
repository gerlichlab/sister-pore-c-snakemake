rule expand_triplets:
    input:
        fragments=paths.align_table.annotated_pore_c
    output:
        paths.contacts.triplet_contacts
    log:
        to_log(paths.contacts.triplet_contacts)
    conda:
        "../envs/spoc.yml"
    shell:
        "spoc expand {input.fragments} {output}"

rule merge_triplets:
    input:
        expand_basecall_batches(paths.contacts.triplet_contacts)
    output:
        directory(paths.merged_contacts.triplet_contacts)
    log:
        to_log(paths.merged_contacts.triplet_contacts)
    conda:
        "../envs/spoc.yml"
    shell:
        "spoc merge contacts {input} -o {output}"