rule create_qc_sister_pore_c:
    output:
        paths.qc.sister_pore_c
    input:
        detect=get_all_detect_files(paths),
        label_library=expand_rows(paths.brdu_calling.label_library, mapping_df),
        all_pairs_paths=expand_rows(paths.pairs.unsorted_pairs, mapping_df),
        all_reads_cis_trans=expand_rows_w_label_types(paths.pairs.label_split, mapping_df, label_type=["all_reads"], contact_type=["cis_and_trans"]),
        all_reads_cis=expand_rows_w_label_types(paths.pairs.label_split, mapping_df, label_type=["all_reads"], contact_type=["cis"]),
        all_reads_trans=expand_rows_w_label_types(paths.pairs.label_split, mapping_df, label_type=["all_reads"], contact_type=["trans"]),
        labelled_cis_trans=expand_rows_w_label_types(paths.pairs.label_split, mapping_df, label_type=["labelled_only"], contact_type=["cis_and_trans"]),
        labelled_cis=expand_rows_w_label_types(paths.pairs.label_split, mapping_df, label_type=["labelled_only"], contact_type=["cis"]),
        labelled_trans=expand_rows_w_label_types(paths.pairs.label_split, mapping_df, label_type=["labelled_only"], contact_type=["trans"])
    conda:
        "envs/qc.yml"
    shell:
        """python bin/create_sister_pore_c_report.py --detect_paths {input.detect} \
                                                   --label_library_paths {input.label_library} \
                                                   --all_pairs_paths {input.all_pairs_paths} \
                                                   --all_reads_cis_trans {input.all_reads_cis_trans} \
                                                   --all_reads_cis {input.all_reads_cis} \
                                                   --all_reads_trans {input.all_reads_trans} \
                                                   --labelled_cis_trans {input.labelled_cis_trans} \
                                                   --labelled_cis {input.labelled_cis} \
                                                   --labelled_trans {input.labelled_trans} \
                                                   --output {output}
        """
