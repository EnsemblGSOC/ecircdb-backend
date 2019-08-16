#!/bin/bash
addcore[0]="\`analysis\`"
addcore[1]="\`assembly\`"
addcore[2]="\`backsplice_junction\`"
addcore[3]="\`canonical_junction\`"
addcore[4]="\`exon\`"
addcore[5]="\`locus\`"
addcore[6]="\`locus_expression\`"
addcore[7]="\`sample\`"
addcore[8]="\`species\`"
addcore[9]="\`ensembl_gene\`"

addedcore[0]="\`core_analysis\`"
addedcore[1]="\`core_assembly\`"
addedcore[2]="\`core_backsplicejunction\`"
addedcore[3]="\`core_canonicaljunction\`"
addedcore[4]="\`core_exon\`"
addedcore[5]="\`core_locus\`"
addedcore[6]="\`core_locusexpression\`"
addedcore[7]="\`core_sample\`"
addedcore[8]="\`core_species\`"
addedcore[9]="\`core_ensemblgene\`"

tablemissingfields[0]="INSERT INTO \`core_analysis\` VALUES"
tablemissingfields[1]="INSERT INTO \`core_assembly\` VALUES"
tablemissingfields[2]="INSERT INTO \`core_backsplicejunction\` VALUES"
tablemissingfields[3]="INSERT INTO \`core_canonicaljunction\` VALUES"
tablemissingfields[4]="INSERT INTO \`core_exon\` VALUES"
tablemissingfields[5]="INSERT INTO \`core_locus\` VALUES"
tablemissingfields[6]="INSERT INTO \`core_locusexpression\` VALUES"
tablemissingfields[7]="INSERT INTO \`core_sample\` VALUES"
tablemissingfields[8]="INSERT INTO \`core_species\` VALUES"
tablemissingfields[9]="INSERT INTO \`core_ensemblgene\` VALUES"

tableaddedfields[0]="INSERT INTO \`core_analysis\` (analysis_id, run_date, logic_name, parameters, assembly_id_id, sample_id_id) VALUES"
tableaddedfields[1]="INSERT INTO \`core_assembly\` (assembly_id, assembly_accession, assembly_name, species_id_id) VALUES"
tableaddedfields[2]="INSERT INTO \`core_backsplicejunction\` (junction_id, browser_string, coord_id, genomic_size, seq_region_end, seq_region_name, seq_region_start, seq_region_strand, jpm, abundance_ratio, predicted_exons, raw_count, splice_signal, analysis_id_id, locus_id_id, classification, in_platelets, is_published, junction_repeat_coverage, n_methods, splice_type, spliced_size, splice_3_support, splice_5_support, splice_site_count, tpm, transcript_count, gc_perc) VALUES"
tableaddedfields[3]="INSERT INTO \`core_canonicaljunction\` (junction_id, browser_string, coord_id, genomic_size, seq_region_end, seq_region_name, seq_region_start, seq_region_strand, jpm, predicted_exons, raw_count, splice_signal, analysis_id_id, locus_id_id, splice_type) VALUES"
tableaddedfields[4]="INSERT INTO \`core_exon\` (locus_id, coord_id, genomic_size, seq_region_end, seq_region_name, seq_region_start, seq_region_strand, rank, stable_id) VALUES"
tableaddedfields[5]="INSERT INTO \`core_locus\` (locus_id, browser_string, coord_id, genomic_size, seq_region_end, seq_region_name, seq_region_start, seq_region_strand, gene_name, is_circrna_host, nexons, source, spliced_size, stable_id, total_backsplice_reads, total_splice_reads, assembly_id_id, circrna_abundance_ratio, biotype) VALUES"
tableaddedfields[6]="INSERT INTO \`core_locusexpression\` (expression_id, raw_count, rpkm, rpkm_external, rpkm_internal, rpkm_ratio, tpm, analysis_id_id, locus_id_id) VALUES"
tableaddedfields[7]="INSERT INTO \`core_sample\` (sample_id, accession, backspliced_reads, chimeric_reads, circrna_count, description, fastqc_path, bigwig_path, ftp_path, instrument, library_size, library_strategy, project, gtag_cjuncs_reads, mapped_reads, read_length, fraction, submitter, source, total_spliced_reads, TIN, species_id_id, fastq_path, source_id) VALUES"
tableaddedfields[8]="INSERT INTO \`core_species\` (taxon_id, scientific_name) VALUES"
tableaddedfields[9]="INSERT INTO \`core_ensemblgene\` (gene_id, coord_id, seq_region_end, seq_region_name, seq_region_start, seq_region_strand, gene_name, description, is_circrna_host, biotype, stable_id, assembly_id_id) VALUES"

sed -ri "s|\) ENGINE=InnoDB AUTO_INCREMENT=[0-9]+ DEFAULT CHARSET=latin1;||g" db.sql

# Change table names
for index in ${!addcore[*]}
do
    sed -i "s|${addcore[$index]}|${addedcore[$index]}|" db.sql
    echo $((index))
done

# Change table field names
for index in ${!tablemissingfields[*]}
do
    sed -i "s|${tablemissingfields[$index]}|${tableaddedfields[$index]}|" db.sql
    echo $((index + 10))
done
