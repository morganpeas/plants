library(dplyr) # data manipulation#
library(tidyr) # ditto
library(tibble) # data structure
library(readxl) # you can guess... reading excel files
library(traits) # importing trait data
library(missForest) # trait imputation
library(ape) # working with phylogenies
library(purrr) # looping over a number of objects
library(phytools) # manipulating phylogenies and running comparative models
library(nlme) # generalised linear and mixed effects models with correlated data
library(caper) # working with phylogenetic data
library(caret) # creating data partitions and running cross validation/training-testing

# S.Phylomaker function that was used to produce the Qian tree
source("https://raw.githubusercontent.com/jinyizju/S.PhyloMaker/master/R_codes%20for%20S.PhyloMaker")
# devtools::install_github("ropenscilabs/datastorr")
# devtools::install_github("wcornwell/taxonlookup")
library(taxonlookup) # taxonomic lookup for land plants (useful for phylogeny generation)
#Going to skip this bit for now. Not sure that I want to generate the missing bits of phylogeny.
# devtools::install_github("traitecoevo/phyndr")
library(phyndr) # maximising overlap between traits and phylogeny 


# read in the monocot test dataset
library(readxl)
monocot_rl <- read_excel("~/Desktop/Project_2/monocots/monocots_test.xlsx")

#Qian phylogeny
# These appear to be the same species, so just keep one
monocot_rl <- distinct(monocot_rl, binomial, .keep_all = TRUE)
# read in the Qian phylogeny
qian_tree <- read.tree("~/Desktop/Project_2/QianTree.tre")
# read in the nodes information of the megaphylogeny
qian_clades <- read.delim("~/Desktop/Project_2/QianNodes.txt",
                          header = TRUE,
                          sep = "\t")
#make sure that the tree is ultrametric
plot(qian_tree_ultrametric$edge.length, qian_tree$edge.length)
# get the root-to-tip distances
diffs <- diag(vcv.phylo(qian_tree))
hist(diffs)
min(diffs)-max(diffs) / mean(diffs)
#safe to say that the tree really is ultrametric and it is a rounding error that 
#is causing it to appear non-ultrametric#

#Explore the overlap between datasets#
head(qian_tree$tip.label)
head(monocot_rl$binomial)
#check how many of the species names match
sum(qian_tree$tip.label %in% monocot_rl$binomial)

#how many of the genus names match# #CONFUSED HERE#
sum(unique(monocot_rl$genus) %in% unique(qian_tree$Genus))

#Add in missing species to the phylogeny
qian_tree$node.label <- 1:length(qian_tree$node.label)
qian_full_tree_ultrametric <- force.ultrametric(qian_tree, method = "extend")

frm <- ~family/genus/binomial
monocot_rl_tr <- monocot_rl %>%
  dplyr::select(family, genus, binomial) %>%
  
  # add in an outgroup
  bind_rows(data.frame(family = "OutgroupFam", 
                       genus = "OutgroupGen", 
                       binomial = "OutgroupSp")) %>%
  mutate(family = factor(family),
         genus = factor(genus),
         binomial = factor(binomial))

tr <- as.phylo(frm, data = monocot_rl_tr)
plot(tr, show.tip.label = FALSE)

tree <- root(tr, outgroup = "OutgroupSp", resolve.root = TRUE)

#Analysis
# remove NA red list assessment status as these are duplicates
  filter(!is.na(iucn_rating))

#Prepare the data for analysis#
#Red list status as a continuous variable#
monocot_rl <- monocot_rl %>%,
# set up a column to contain the continous red list status#
  RL_cont = monocot_rl$iucn_rating

# if the rating is DD, turn it ino NA#
  RL_cont = na_if(RL_cont, "DD")

# recode red list status into a continuous variable#
  RL_cont = 
    recode(RL_cont, 
         LC = 1,
         NT = 2,
         VU = 3,
         EN = 4, 
         CR = 5)
  RL_cont = as.numeric(RL_cont)

# set up a column to contain the binary red list status#
  RL_binary = monocot_rl$iucn_rating

# if the rating is DD, turn it ino NA#
  RL_binary = na_if(RL_binary, "DD")

# recode red list status into a binary (0|1) variable#
  RL_binary = 
  recode(RL_binary,
         LC = 0,
         NT = 0,
         VU = 1,
         EN = 1, 
         CR = 1)
  RL_binary = as.numeric(RL_binary) 
  
#add columns to monocot dataset#
monocot_rl$RL_binary <- RL_binary 
monocot_rl$RL_cont <- RL_cont

# the DD data are the data we're going to predict once we have a model
# so let's pull these data out
monocot_noDD <- monocot_rl %>%
  filter(iucn_rating != "DD")
monocot_DD <- monocot_rl %>%
  filter(iucn_rating == "DD")

# match the dataset to the phylogeny
monocot_comp <- comparative.data(phy = qian_full_tree_ultrametric, 
                                 data = as.data.frame(monocot_noDD), 
                                 names.col = 'binomial',
                                 na.omit = FALSE)

#Multicollinearity#
monocot_multi <- read_excel("~/Desktop/Project_2/monocots/monocots_collinearitytest.xlsx")
library(car)
all(is.na(monocot_multi))

source("http://highstat.com/Books/BGS/GAMM/RCodeP2/HighstatLibV6.R")
test <- corvif(monocot_multi) #dataframe with columns with only explanatory variables#
test
plot(monocot_multi)
scale(monocot_multi$maximum_height_m)


# A gls without accounting for phylogeny#
#Categorical variables#
#growth form#
fit <- gls(RL_cont ~ growth_form,
           monocot_comp$data,
           method = "ML")
phylosig(tree = monocot_comp$phy,
         x = resid(fit),
         method = "lambda",
         test = TRUE)
#habitat#
fit <- gls(RL_cont ~ habitat,
           monocot_comp$data,
           method = "ML")
phylosig(tree = monocot_comp$phy,
         x = resid(fit),
         method = "lambda",
         test = TRUE)
#perennial_annual#
fit <- gls(RL_cont ~ perennial_annual,
           monocot_comp$data,
           method = "ML")
phylosig(tree = monocot_comp$phy,
         x = resid(fit),
         method = "lambda",
         test = TRUE)
#leaf_shape#
fit <- gls(RL_cont ~ leaf_shape,
           monocot_comp$data,
           method = "ML")
phylosig(tree = monocot_comp$phy,
         x = resid(fit),
         method = "lambda",
         test = TRUE)
#fruit_type#
fit <- gls(RL_cont ~ fruit_type,
           monocot_comp$data,
           method = "ML")
phylosig(tree = monocot_comp$phy,
         x = resid(fit),
         method = "lambda",
         test = TRUE)

#Numeric Variables#
#TO DO: NORMALIZE DATA AND FIGURE OUT HOW TO GET RID OF THE 'NA' IN EACH#
#max_height#
hist(monocot_comp$maximum_height_m)
fit <- gls(RL_cont ~ max_height,
           monocot_comp$data,
           method = "ML")
phylosig(tree = monocot_comp$phy,
         x = resid(fit),
         method = "lambda",
         test = TRUE)
#leaf_maximum_width_cm#
fit <- gls(RL_cont ~ leaf_maximum_width_cm,
           monocot_comp$data,
           method = "ML")
phylosig(tree = monocot_comp$phy,
         x = resid(fit),
         method = "lambda",
         test = TRUE)
#leaf_maximum_length_cm#
fit <- gls(RL_cont ~ leaf_maximum_length_cm,
           monocot_comp$data,
           method = "ML")
phylosig(tree = monocot_comp$phy,
         x = resid(fit),
         method = "lambda",
         test = TRUE)

#A gls accounting for phylogenetic signal - pgls#
#Categorical variables#
#growth_form
fit <- gls(RL_cont ~ growth_form,
           monocot_comp$data, 
           correlation = corPagel(1, 
                                  monocot_comp$phy,
                                  fixed = FALSE),
           method = "ML")
summary(fit)
#habitat
fit <- gls(RL_cont ~ habitat,
           monocot_comp$data, 
           correlation = corPagel(1, 
                                  monocot_comp$phy,
                                  fixed = FALSE),
           method = "ML")
summary(fit)
#perennial_annual
fit <- gls(RL_cont ~ perennial_annual,
           monocot_comp$data, 
           correlation = corPagel(1, 
                                  monocot_comp$phy,
                                  fixed = FALSE),
           method = "ML")
summary(fit)
#leaf_shape
fit <- gls(RL_cont ~ leaf_shape,
           monocot_comp$data, 
           correlation = corPagel(1, 
                                  monocot_comp$phy,
                                  fixed = FALSE),
           method = "ML")
summary(fit)
#fruit_type
fit <- gls(RL_cont ~ fruit_type,
           monocot_comp$data, 
           correlation = corPagel(1, 
                                  monocot_comp$phy,
                                  fixed = FALSE),
           method = "ML")
summary(fit)
#Numeric Variables#
#max_height#
fit <- gls(RL_cont ~ max_height,
           monocot_comp$data, 
           correlation = corPagel(1, 
                                  monocot_comp$phy,
                                  fixed = FALSE),
           method = "ML")
summary(fit)
#leaf_maximum_width_cm
fit <- gls(RL_cont ~ leaf_maximum_width_cm,
           monocot_comp$data, 
           correlation = corPagel(1, 
                                  monocot_comp$phy,
                                  fixed = FALSE),
           method = "ML")
summary(fit)
#leaf_maximum_length_cm
fit <- gls(RL_cont ~ leaf_maximum_length_cm,
           monocot_comp$data, 
           correlation = corPagel(1, 
                                  monocot_comp$phy,
                                  fixed = FALSE),
           method = "ML")
summary(fit)

