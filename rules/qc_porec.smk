rule create_qc_pore_c:
    output:
        paths.qc.pore_c
    input:
        csv=paths.merged_contacts.concatemer_summary,
        pq=paths.align_table.alignment