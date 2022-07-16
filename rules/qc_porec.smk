rule create_qc_pore_c:
    output:
        paths.qc.pore_c
    input:
        csv=expand_rows(paths.merged_contacts.concatemer_summary, mapping_df),
        pq=expand_rows(paths.merged_contacts.contacts, mapping_df)
    log:
        to_log(paths.qc.pore_c)
    conda:
        "../envs/qc.yml"
    shell:
        "python bin/create_pore_c_qc_plots.py --concatamer_tables {input.pq} --concatamer_summaries {input.csv} --output {output}"