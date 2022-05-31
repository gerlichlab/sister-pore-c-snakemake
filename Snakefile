import pandas as pd
import yaml
from pathlib import Path
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.5.4")


wildcard_constraints:
    enzyme="[^_]+",
    run_id="[^_]+",
    batch_id="batch\d+",


BASE_DIR = Path(workflow.basedir)


### Validation of schemas ###
##### load config and sample sheets ##
configfile: BASE_DIR / "config/config.yaml"


##### load rules #####
include: "rules/common.smk"  # python helper functions

basecall_df, reference_df, mapping_df = create_config_dataframes()
paths = create_path_accessor()

include: "rules/refgenome.smk"  # prepare the reference genome
include: "rules/reads.smk"  # import fastqs, fast5s, sequencing summary
include: "rules/mapping.smk"  # map and process resulting alignments
include: "rules/exports.smk"  # export to alternative formats
include: "rules/methylation.smk"  # use f5c to call cpg methylation
include: "rules/qc_porec.smk" # qc stuff
include: "rules/brdu_calling.smk"
include: "rules/brdu_assignment.smk"
include: "rules/generate_coolers_sister_specific.smk"

##### output paths #####

rule all:
    input:
        basecalls=expand_rows(paths.basecall.catalog, basecall_df),
        refgenome=expand_rows(paths.refgenome.bwt, reference_df),
        contacts=expand_rows(paths.merged_contacts.concatemers, mapping_df),

rule qc:
    input:
        qc=expand_rows(paths.qc.pore_c, mapping_df)

rule create_mcoolers_weight_transferred:
    input:
        mcoolers_weight_transferred=expand_rows_w_label_types(paths.matrix.mcool_split_weights_transfered, mapping_df)

rule create_mcoolers_label_specific:
    input:
        mcoolers=expand_rows_w_label_types(paths.matrix.mcool_label_split, mapping_df)

rule assign_brdu:
    input:
        assigned_pairs=expand_rows(paths.pairs.assigned_pairs, mapping_df)

rule split_assigned_pairs:
    input:
        split_pairs=expand_rows_w_label_types(paths.pairs.label_split, mapping_df)

rule make_library:
    input:
        library=expand_rows(paths.brdu_calling.label_library, mapping_df)

rule brdu_index:
    input:
        indices=expand_rows(paths.brdu_calling.index, mapping_df)

rule cooler:
    input:
        expand_rows(paths.matrix.mcool, mapping_df),


rule pairs:
    input:
        expand_rows(paths.pairs.index, mapping_df),


rule salsa:
    input:
        expand_rows(paths.assembly.salsa_bed, mapping_df),


rule juicer:
    input:
        expand_rows(paths.juicebox.hic, mapping_df),


rule mnd:
    input:
        expand_rows(paths.juicebox.mnd, mapping_df),


rule f5c_indexes:
    # Create f5c indexes, but do not call methylation
    input:
        expand_rows(paths.basecall.f5c_index_all, mapping_df),


rule methylation:
    input:
        expand_rows(paths.methylation.per_locus_methylation, mapping_df),


rule test:
    input:
        expand_rows(paths.merged_contacts.concatemers, mapping_df),
        expand_rows(paths.matrix.mcool, mapping_df),  #expand_rows(paths.matrix.haplotyped_cools, mapping_df),
        expand_rows(paths.pairs.index, mapping_df),
        expand_rows(paths.assembly.salsa_bed, mapping_df),
        expand_rows(paths.juicebox.hic, mapping_df),
        expand_rows(paths.juicebox.mnd, mapping_df),
