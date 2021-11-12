#############################################################################################################################
######### SARS-CoV-2 XB lineage investigation GISAID wrangler ###############################################################
#############################################################################################################################

######### Bernardo Gutierrez ################################################################################################

library(seqinr)
library(ape)
library(tidyverse)
library(data.table)
library(lubridate)
library(plyr)

######### Read metadata and sequence data ###################################################################################

# Read metadata files
xb_meta_full_store <- read.csv("XB_GISAID/XB_analysis_metadata_full.tsv", sep = "\t", strip.white = TRUE) %>% 
  select(c(strain, gisaid_epi_isl, date, region, country, division, location, region_exposure, country_exposure,
           division_exposure, pangolin_lineage, date_submitted))

xb_meta_filt_store <- read.csv("XB_GISAID/XB_analysis_metadata_filtered.tsv", sep = "\t", strip.white = TRUE) %>% 
  select(c(strain, gisaid_epi_isl, date, region, country, division, location, region_exposure, country_exposure,
           division_exposure, pangolin_lineage, date_submitted))

# Read aligned sequence files
#xb_aln_full <- read.fasta(file = "XB_GISAID/xxxxxx.fasta", seqtype = "DNA", as.string = TRUE,
#                   set.attributes = FALSE)

xb_aln_filt <- read.fasta(file = "XB_GISAID/XB_analysis_sequences_filtered_aln_masked.fasta", seqtype = "DNA",
                          as.string = TRUE, set.attributes = FALSE)

# Read updated Pangolin output (Pangolin v3.1.11)
xb_pango <- read.csv("Pangolin_XB_analysis_20210903.csv", sep = ",", strip.white = TRUE) %>% 
  select(c(X...Sequence.name, Lineage)) # <-- ID column name may change depending on where the code is run from
colnames(xb_pango) <- c("id", "lineage")

# Read updated phylogenetic lineage assignment
xb_phylo_new_lin <- read.csv("XB_20210903_phylogenetic_reassignments.csv", sep = "\t", strip.white = TRUE)


######### Clean and filter metadata #########################################################################################

# Remove white spaces
xb_meta_full <- as.data.frame(apply(xb_meta_full_store, 2, function(x) gsub('\\s+', '_',x)))
xb_meta_filt <- as.data.frame(apply(xb_meta_filt_store, 2, function(x) gsub('\\s+', '_',x)))

# Clean date formats
xb_meta_full$date <- ymd(xb_meta_full$date)
xb_meta_full$date_submitted <- ymd(xb_meta_full$date_submitted)
xb_meta_filt$date <- ymd(xb_meta_filt$date)
xb_meta_filt$date_submitted <- ymd(xb_meta_filt$date_submitted)

# Clean up region column for Central American countries
xb_meta_full$region[xb_meta_full$country=="Costa_Rica"] <- "Central_America"
xb_meta_full$region[xb_meta_full$country=="El_Salvador"] <- "Central_America"
xb_meta_full$region[xb_meta_full$country=="Guatemala"] <- "Central_America"
xb_meta_full$region[xb_meta_full$country=="Honduras"] <- "Central_America"

xb_meta_filt$region[xb_meta_filt$country=="Costa_Rica"] <- "Central_America"
xb_meta_filt$region[xb_meta_filt$country=="El_Salvador"] <- "Central_America"
xb_meta_filt$region[xb_meta_filt$country=="Guatemala"] <- "Central_America"
xb_meta_filt$region[xb_meta_filt$country=="Honduras"] <- "Central_America"

# Clean up region column for Caribbean countries
xb_meta_full$region[xb_meta_full$country=="Cayman_Islands"] <- "Caribbean"
xb_meta_filt$region[xb_meta_filt$country=="Cayman_Islands"] <- "Caribbean"

# Set variables as factors
x <- xb_meta_full %>% select(c(strain, gisaid_epi_isl, date, date_submitted))
y <- xb_meta_full %>% select(c(region, country, division, location, region_exposure, country_exposure, division_exposure,
                               pangolin_lineage)) %>% lapply(as.factor) %>% as.data.frame()
xb_meta_full <- bind_cols(x, y)

x <- xb_meta_filt %>% select(c(strain, gisaid_epi_isl, date, date_submitted))
y <- xb_meta_filt %>% select(c(region, country, division, location, region_exposure, country_exposure, division_exposure,
                               pangolin_lineage)) %>% lapply(as.factor) %>% as.data.frame()
xb_meta_filt <- bind_cols(x, y)

# Add second Pango lineage column from updated Pangolin run
xb_meta_full$pangolin_lineage_update <- xb_pango$lineage[xb_pango$id %in% xb_meta_full$strain]
xb_meta_filt$pangolin_lineage_update <- xb_pango$lineage[xb_pango$id %in% xb_meta_filt$strain]

# Add third phylogenetically assigned lineage column; i.e. identify Minor lineages
xb_meta_filt$pangolin_lineage_phylo <- as.character(xb_meta_filt$pangolin_lineage)
xb_meta_filt$pangolin_lineage_phylo[xb_meta_filt$strain %in% xb_phylo_new_lin$strain] <- "B.1.628minor"
xb_meta_filt$pangolin_lineage_phylo <- as.factor(xb_meta_filt$pangolin_lineage_phylo)

# Add marker for sequences before and after April 2021 for B.1.628
# Rationale: long tail of few B.1.628 sequences during this time
xb_meta_filt$april_threshold <- rep("NA", nrow(xb_meta_filt))
xb_meta_filt$april_threshold[xb_meta_filt$pangolin_lineage_phylo=="B.1.628" &
                               xb_meta_filt$date>"2021-04-01"] <- "post_April"
xb_meta_filt$april_threshold[xb_meta_filt$pangolin_lineage_phylo=="B.1.628" &
                               xb_meta_filt$date<="2021-04-01"] <- "pre_April"

# Add marker for sequences before and after May 15 2021 for B.1.628
# Rationale: dip in new sequences (for the Americas) at this time
xb_meta_filt$may_threshold <- rep("NA", nrow(xb_meta_filt))
xb_meta_filt$may_threshold[xb_meta_filt$pangolin_lineage_phylo=="B.1.628" &
                               xb_meta_filt$date>"2021-05-13"] <- "post_May_13"
xb_meta_filt$may_threshold[xb_meta_filt$pangolin_lineage_phylo=="B.1.628" &
                               xb_meta_filt$date<="2021-05-13"] <- "pre_May_13"

######### Attach geocoding info ###########################################################################################

z <- xb_geocoded_dt %>% select(lon, lat)
xb_meta_full <- bind_cols(xb_meta_full, z)


######### Downsampling metadata files #####################################################################################

# American downsample, all sequences excluding Africa, Asia and Europe
xb_meta_filt_Americas <- xb_meta_filt[xb_meta_filt$region!="Africa" & xb_meta_filt$region!="Asia" & xb_meta_filt$region!="Europe",]

# Systematic downsample, one sample per administrative division per day
xb_meta_filt_down <-  plyr::ddply(xb_meta_filt,.(division, date), function(x) x[sample(nrow(x),1),]) # <-- systematic downsampling

# Double downsample for GARD, 200 random sequences from the systematic downsample
xb_meta_filt_gard <-  xb_meta_filt_down[sample(nrow(xb_meta_filt_down), 200), ] # <-- GARD subsampling

# Random downsampled data sets - 10 independent random downsamples (same size as the systematic downsample)
# [!] Note that new random datasets are created every time these lines are run [!]
xb_meta_filt_rand1 <-  xb_meta_filt[sample(nrow(xb_meta_filt), nrow(xb_meta_filt_down)), ] # <-- random subsampling no.1
xb_meta_filt_rand2 <-  xb_meta_filt[sample(nrow(xb_meta_filt), nrow(xb_meta_filt_down)), ] # <-- random subsampling no.2
xb_meta_filt_rand3 <-  xb_meta_filt[sample(nrow(xb_meta_filt), nrow(xb_meta_filt_down)), ] # <-- random subsampling no.3
xb_meta_filt_rand4 <-  xb_meta_filt[sample(nrow(xb_meta_filt), nrow(xb_meta_filt_down)), ] # <-- random subsampling no.4
xb_meta_filt_rand5 <-  xb_meta_filt[sample(nrow(xb_meta_filt), nrow(xb_meta_filt_down)), ] # <-- random subsampling no.5
xb_meta_filt_rand6 <-  xb_meta_filt[sample(nrow(xb_meta_filt), nrow(xb_meta_filt_down)), ] # <-- random subsampling no.6
xb_meta_filt_rand7 <-  xb_meta_filt[sample(nrow(xb_meta_filt), nrow(xb_meta_filt_down)), ] # <-- random subsampling no.7
xb_meta_filt_rand8 <-  xb_meta_filt[sample(nrow(xb_meta_filt), nrow(xb_meta_filt_down)), ] # <-- random subsampling no.8
xb_meta_filt_rand9 <-  xb_meta_filt[sample(nrow(xb_meta_filt), nrow(xb_meta_filt_down)), ] # <-- random subsampling no.9
xb_meta_filt_rand10 <-  xb_meta_filt[sample(nrow(xb_meta_filt), nrow(xb_meta_filt_down)), ] # <-- random subsampling no.10

# Recombination test downsample, five random sequences per lineage 
xb_meta_filt_rec <- plyr::ddply(xb_meta_filt_Americas,.(pangolin_lineage_phylo), function(x) x[sample(nrow(x),5),])

######### Cleanup and data export #########################################################################################

# Write FASTA files from downsampled datasets
write.fasta(sequences = xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_Americas$strain],
            names = names(xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_Americas$strain]),
            nbchar = 300, file.out = "XB_analysis_20210831_Americas.fasta")

write.fasta(sequences = xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_gard$strain],
            names = names(xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_gard$strain]),
            nbchar = 300, file.out = "XB_analysis_20210831_subset_GARD.fasta")

write.fasta(sequences = xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_down$strain],
            names = names(xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_down$strain]),
            nbchar = 300, file.out = "XB_analysis_20210831_syst_subset.fasta")

write.fasta(sequences = xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand1$strain],
            names = names(xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand1$strain]),
            nbchar = 300, file.out = "XB_analysis_20210831_subset_rand1.fasta")

write.fasta(sequences = xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand2$strain],
            names = names(xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand2$strain]),
            nbchar = 300, file.out = "XB_analysis_20210831_subset_rand2.fasta")

write.fasta(sequences = xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand3$strain],
            names = names(xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand3$strain]),
            nbchar = 300, file.out = "XB_analysis_20210831_subset_rand3.fasta")

write.fasta(sequences = xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand4$strain],
            names = names(xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand4$strain]),
            nbchar = 300, file.out = "XB_analysis_20210831_subset_rand4.fasta")

write.fasta(sequences = xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand5$strain],
            names = names(xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand5$strain]),
            nbchar = 300, file.out = "XB_analysis_20210831_subset_rand5.fasta")

write.fasta(sequences = xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand6$strain],
            names = names(xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand6$strain]),
            nbchar = 300, file.out = "XB_analysis_20210831_subset_rand6.fasta")

write.fasta(sequences = xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand7$strain],
            names = names(xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand7$strain]),
            nbchar = 300, file.out = "XB_analysis_20210831_subset_rand7.fasta")

write.fasta(sequences = xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand8$strain],
            names = names(xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand8$strain]),
            nbchar = 300, file.out = "XB_analysis_20210831_subset_rand8.fasta")

write.fasta(sequences = xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand9$strain],
            names = names(xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand9$strain]),
            nbchar = 300, file.out = "XB_analysis_20210831_subset_rand9.fasta")

write.fasta(sequences = xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand10$strain],
            names = names(xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rand10$strain]),
            nbchar = 300, file.out = "XB_analysis_20210831_subset_rand10.fasta")

# Sequences for recombination data set (manually edit alignment after export)
write.fasta(sequences = xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rec$strain],
            names = names(xb_aln_filt[names(xb_aln_filt) %in% xb_meta_filt_rec$strain]),
            nbchar = 300, file.out = "XB_analysis_20210831_recombination_test.fasta")

# Write date metadata file
write.table(data.frame(id = xb_meta_filt$strain, date = xb_meta_filt$date), quote = FALSE,
          file = "xb_meta_filt_dates.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

# Write date threshold metadata file
write.table(data.frame(id = xb_meta_filt$strain,
                       April_threshold = xb_meta_filt$april_threshold,
                       May_13_threshold = xb_meta_filt$may_threshold), quote = FALSE,
            file = "xb_meta_filt_thresholds.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

rm(x, y, z, xb_meta_filt_store, xb_meta_full_store)
