# vtamR 1.0.4 (prerelease) 2026/07/10

- Add pattern argument to history_by and summarize_by
- Modify history_by to accept a vector of values, not just a singe value
- Correct filter_asv_global =< cutoff instead of < cutoff
- Add concatenate_files function: Reads files from multiple directories 
  and concatenates the contents of files sharing the same filename
- Modify filter_pcr_error: Add min_read_count argument. If less than min_read_count
  in the sample/run, the variant is not checked.
- Correct download_osf to work both of windows and linux
- Add demultiplex_fastq_pairs
- Correct count_reads when using bash commands
- Modify check_dir: delete slash at the end of dir and return invisible path


# vtamR 1.0.3 (prerelease) 2026/04/16

- Short tutoriel added
- Initialize stat_df from get_stat
- Make dir, outdir obligatory

# vtamR 1.0.2 (prerelease) 2026/04/16

Documentation updated

# vtamR 1.0.1 (prerelease) 2026/04/13

Tutorial updated

# vtamR 1.0.0 (prerelease)

- Function Names Harmonized: Correspondences between old and new function 
names can be found in `vtamR/R_divers/vtam_functions_table.csv`.


# vtamR 0.3.2 (prerelease)

- PoolDatasets is split to pool_datasets and pool_markers
- Add mean, sum, min, max as aggregation methods to pool_datasets, pool_markers
PoolReplicates, WriteASVtable, 
- Correct the read_counts returned by RandomSeq
- Add denoise_by_swarm
- Complete tutorial with pooling two plates or two markers
- Troubleshooting section to pool runs before filtering
- Correct output file stucture

# vtamR 0.3.1 (prerelease)

- Fixed version number in `DESCRIPTION` to match the prerelease tag.
- Updated package metadata for consistency with GitHub prerelease.
- Minor internal cleanups.

# vtamR 0.3.0

- Added functions for clustering and to optimize clustering thresholds
	- ClusterASV
	- ClassifyClusters
	- PairwiseIdentity
	- PairwiseIdentityPlotPerSwarmD
	- PairwiseIdentityPlotPerClusterIdentityThreshold
	- ClassifyClusters
	- PlotClusterClasstification	
- Added functions for taxassign	
	- TaxAssigRDP
	- TaxAssig named to TaxAssigLTG	
- Other new functions
	- MakeMockCompositionLTG
	- ASVspecificCutoff
- Improved handling of path variable, num_threads, file compression, third partie programs
- Improved documentation

# vtamR 0.2.0 (prerelease)

- Improoved tutorial
- Improoved install instructions from github
- Improoved documentation of the functions

# vtamR 0.1.0 (prerelease)

- Initial stable prerelease of vtamR.
- Core functions implemented for:
	- fastq prepocesssing (Merge, SortReads, RandomSeq, TrimPrimer)
	- filtering steps
	- optimize functions
	- Taxonomic assignment
	- Reporting (Variant table, read counts)
- Added vignettes and introductory examples.
