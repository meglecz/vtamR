library(phyloseq)
library(dplyr)
library(vtamR)


otu <- "/home/meglecz/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/5_cluster/16_mOTU_vsearch.csv"
tax <- "/home/meglecz/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/4_filter/assign_taxonomy_ltg.csv"
samples <- "/home/meglecz/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/3_demultiplexed/sampleinfo.csv"
# Create a data frame containing sample metadata
samples <- read_input(samples, sep=",") %>%
  select(sample, sample_type, habitat) %>%
  distinct()

samples = NULL
rm_control = TRUE
sep = ","

phyloseq_obj <- format_for_phyloseq(otu = otu,  samples = NULL, rm_control = FALSE)
phyloseq_obj <- format_for_phyloseq(otu = otu, tax = tax, samples = NULL, rm_control = TRUE)






samples <- sampleinfo_df %>%
  select(sample, sample_type, habitat) %>%
  distinct()

phy_object <- make_phyloseq(
  otu = otu,
  tax = tax,
  sep = ","
)

sd = as.matrix(sample_data(phy_object))
  
phy_object
OTU1 <- as(otu_table(phy_object), "matrix")
class(phy_object)
sample_names(phy_object)
rank_names(phy_object)
sample_variables(phy_object)
phy_object <- subset_samples(phy_object, sample_type =="real")
SAMPLE1 <- as(sample_data(phy_object), "matrix")
OTU2 <- as(otu_table(phy_object), "matrix")
phy_object

total = median(sample_sums(phy_object))
standf = function(x, t=total) round(t * (x / sum(x)))
phy_object = transform_sample_counts(phy_object, standf)
OTU3 <- as(otu_table(phy_object), "matrix")

plot_bar(phy_object, fill = "class") +
  geom_bar(aes(color=class, fill=class), stat="identity", position="stack")

plot_richness(phy_object, measures=c("Chao1", "Shannon"))

phy_object.ord <- ordinate(phy_object, "NMDS", "bray")
plot_ordination(phy_object, phy_object.ord, type="taxa", color="class", shape= "phylum", 
                title="OTUs")

plot_ordination(phy_object, phy_object.ord, type="samples", color="sample_type", 
                shape="habitat", title="Samples") + geom_point(size=3)
