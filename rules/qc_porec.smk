rule create_qc_pore_c:
    output:
        paths.qc.pore_c
    input:
        csv=paths.merged_contacts.concatemer_summary,
        pq=paths.merged_contacts.contacts
    conda:
        "../envs/qc.yml"
    shell:
        "python bin/create_pore_c_qc_plots.py --concatamer_table {input.pq} --concatamer_summary {input.csv} --output {output}"