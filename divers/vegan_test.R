
library(vegan)

data(package = "vegan") ## names of data sets in the package

# species in col, site in line, no name for sites
data(dune) # Vegetation and Environment in Dutch Dune Meadows
str(dune) #a data frame of observations of 30 species at 20 sites

# species in line taxlevel in col. Informative rownames et colnames
data(dune.taxon)
str(dune.taxon)

data(dune.env)
str(dune.env)


x = otu

otu <- "/home/meglecz/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/5_cluster/16_mOTU_vsearch.csv"
tax <- "/home/meglecz/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/4_filter/assign_taxonomy_ltg.csv"
sample_type <- "/home/meglecz/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/3_demultiplexed/sampleinfo.csv"
rm_coltrol = TRUE
sep = ","
outfile_motu = "/home/meglecz/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/vegan/motu_table_vegan.csv"
outfile_taxa = "/home/meglecz/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/vegan/tax_table_vegan.csv"

vegan_dfs <- format_for_vegan(otu, tax = tax, rm_coltrol=TRUE, sample_type = sample_type, sep = ",")
vegan_dfs <- format_for_vegan(otu, tax = tax, rm_coltrol=TRUE, sample_type = sample_type, outfile_motu = outfile_motu, outfile_taxa = outfile_taxa)
vegan_dfs <- format_for_vegan(otu, tax = tax, rm_coltrol=FALSE, sep = ",")
vegan_motu_df <- vegan_dfs[[1]]
vegan_tax_df <- vegan_dfs[[2]]



simpson <- diversity(vegan_motu_df, "simpson") # or assign to var.
simpson 


raremin <- min(rowSums(vegan_motu_df)) 
raremin
sRare <- rarefy(vegan_motu_df, raremin) # now use function rarefy
sRare #gives an "expected"rarefied" number of species (not obs) if only 15 individuals were present
rarecurve(vegan_motu_df, col = "blue")


set.seed(2) # random no. generator / way to specify seeds, 2=no. of integers?
community_matrix=matrix(
  sample(1:100,300,replace=T),nrow=10, # counts up to 100, 300 cells
  dimnames=list(paste("community",1:10,sep=""),paste("sp",1:30,sep="")))
head(community_matrix)
head(vegan_motu_df)

example_NMDS=metaMDS(vegan_motu_df, # Our community-by-species matrix
                     k=2) # The number of reduced dimensions. Increase if high stress is problem. 
plot(example_NMDS)
ordiplot(example_NMDS,type="n") #Ordination plot function especially for congested plots
orditorp(example_NMDS,display="species",col="red",air=0.01) #The function adds text or points to ordination plots
orditorp(example_NMDS,display="sites",cex=1.25,air=0.01)


###################
library(dplyr)
library(tidyr)
library(vegan)
library(tibble)

# Example taxonomy table: one row per taxon (e.g. OTU/ASV)
taxonomy <- data.frame(
  taxon_id = paste0("OTU", 1:6),
  Kingdom  = "Bacteria",
  Phylum   = c("Firmicutes", "Firmicutes", "Proteobacteria",
               "Proteobacteria", "Bacteroidetes", "Bacteroidetes"),
  Genus    = c("Lactobacillus", "Clostridium", "Escherichia",
               "Pseudomonas", "Bacteroides", "Prevotella")
)

# Example abundance table: long format (sample, taxon_id, count)
abundance_long <- data.frame(
  sample_id = rep(c("S1", "S2", "S3"), each = 6),
  taxon_id  = rep(paste0("OTU", 1:6), times = 3),
  count     = c(10, 0, 5, 2, 0, 3,
                4, 8, 0, 0, 6, 1,
                0, 2, 7, 5, 3, 0)
)


genus_long <- abundance_long %>%
  left_join(taxonomy, by = "taxon_id") %>%
  group_by(sample_id, Genus) %>%
  summarise(count = sum(count), .groups = "drop")


genus_wide <- genus_long %>%
  pivot_wider(names_from = Genus, values_from = count, values_fill = 0) %>%
  column_to_rownames("sample_id")

genus_wide