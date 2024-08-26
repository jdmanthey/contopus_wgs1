# contopus_wgs1
divergence history of two contopus species and characterization of a few hybrids

Descriptions of scripts to reproduce whole genome sequencing analyses of wood-pewees in Nebraska.

00_setup.sh – prep directories

00_mutations_rate – directory with scripts to estimate mutation rate from four bird genomes and a time-calibrated phylogeny

01_satsuma.sh – rearrange flycatcher genome to a chromosome scale genome

02_rename_filter_genome.r – rename scaffolds and reorder the scaffolds of the reference

03_index_reference.sh – index the new reference genome

04_filter_align.sh – filter the fasq data and align to the reference

05_extract_filtering_info.sh – get info from filtering log files

06_create_genotyping_jobs.r – R script to make SLURM array scripts for use with GATK

07_cat_vcf_files.sh – concatenate the output VCF files after running GATK

08_filter1.sh – intitial filtering of VCF files

09_structure.sh – run ADMIXTURE and PCA

10_filter2.sh – filtering the VCF files for additional analyses

11_cat_vcf_files.sh – concatenate the newly filtered VCF files

12_make_window_stats_script.r – R script to make SLURM array job to calculate stats in sliding windows

13_cat_window_output.sh – concatnetate the stats outputs from all windows

14{various_names} – scripts for runnind ldhat

15_filter_for_window_admixture.sh – filtering VCF files to run ADMIXTURE in sliding windows

16_create_window_admixture_scripts.r – R script to make SLURM array job for running ADMIXTURE for all sliding windows

17_combine_admixture_output.r – R script to combine all the ADMIXTURE outputs

18_filter3_demography.sh – filter VCFs for Stairway analyses

19_demography_sfs_prep.sh – get the SFS for use in Stairway

20_demography_stairway.sh – run Stairway
21_{various} – estimate and plot coverage for all samples
22_another_PCA.sh – run another PCA without the outgroup (during revisions)


