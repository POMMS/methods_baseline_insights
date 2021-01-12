library(philr); packageVersion("philr")
library(phyloseq); packageVersion("phyloseq")
library(ape); packageVersion("ape")
library(ggplot2); packageVersion("ggplot2")
library(glmnet); packageVersion('glmnet')
library(ggtree); packageVersion("ggtree")
library(ggimage); packageVersion("ggimage")
library(dplyr); packageVersion('dplyr')
library(stringr); packageVersion('stringr')

setwd("C:/Users/kim/Documents/R24_interim_paper")

# this code largely duplicates the PhILR vignette here:
# https://github.com/jsilve24/philr

# note: this workflow involves filtering out low counts entirely; we might prefer to
# agglomerate these into a low-count group, i.e. via
#   counts.genus <- otu_table(data.genus)@.Data
#   N <- ncol(counts)
#   include_vector <- apply(counts.genus, 1, function(x) {
#     sum(x >= 5) / N >= 0.1
#   })
#   other_label <- taxa_names(data.genus)[which(include_vector == FALSE)[1]]
#   data.merged <- merge_taxa(data.genus, taxa_names(data.genus)[!include_vector], 1)

# read in data from a phyloseq object
data <- readRDS("data/jj_r24_exp_v4tree_noOutliersV2.rds")

# optional: agglomerate to species level
# data.species <- tax_glom(data, taxrank="Species", NArm=FALSE)

# filter low-abundance taxa taxa into an "Other" group
# inclusion criterion: > 3-count in more than 10% of samples (across patients)
data.filtered <-  filter_taxa(data, function(x) sum(x > 3) > (0.1*length(x)), TRUE)

counts <- otu_table(data.filtered)@.Data

cat("Proportion zeros after filtering:",round(sum(counts==0)/length(c(counts)), 3),"\n")

# optional: filter taxa with a very low coefficient of variation
# data.filtered <-  filter_taxa(data.filtered, function(x) sd(x)/mean(x) > 2.0, TRUE)

# add a pseudocount
data.final <- transform_sample_counts(data.filtered, function(x) x+1)

cat("Remaining taxa:",ntaxa(data.final),"\n")

# check tree is rooted and binary
is.rooted(phy_tree(data.final))
is.binary.tree(phy_tree(data.final))

# label nodes n###
phy_tree(data.final) <- makeNodeLabel(phy_tree(data.final), method="number", prefix='n')
name.balance(phy_tree(data.final), tax_table(data.final), 'n1')

# un-agglomerated root is: [1] "Order_Bacteroidales/Kingdom_Bacteria"
# species-agglomerated root is: [1] "Order_Bacteroidales/Unclear_Lineage_Identity"

# fix level naming in accordance with other figures
levels(sample_data(data.final)$Health) <- c("HWC", "OB")
metadata <- sample_data(data.final)

# pull stuff out of the phyloseq object
otu.table <- t(otu_table(data.final)) # dimensions: samples (patients) x taxa
tree <- phy_tree(data.final)
tax <- tax_table(data.final)

# to translate an SV idenfitier to a parsable taxonomic label:
# as.vector(tax[rownames(tax) == colnames(otu.table[1,1])[[1]]])

data.philr <- philr(otu.table, tree, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt') # dimensions: samples (patients) x taxa-1

# ordinate
data.dist <- dist(data.philr, method="euclidean")
data.pcoa <- ordinate(data.final, 'PCoA', distance=data.dist)
plot_ordination(data.final, data.pcoa, color='Health') + geom_point(size=4)

# fit a sparse logistic regression model to identify distinguishing balances between HWC and OB cohorts
# note: alpha=1 denotes LASSO
# this should internally use CV to choose a value for lambda (penalization weight)
glmmod <- glmnet(data.philr, sample_data(data.final)$Health, alpha=1, family="binomial")

# show betas
# glmmod$lambda

# threshold balances by coefficient magnitude (this is pretty arbitrary, FYI)
top.coords <- as.matrix(coefficients(glmmod, s=0.15))

# save the full data set: names + beta values of all
top.values <- top.coords

# get the names of the thresholded balances
top.coords <- rownames(top.coords)[which(top.coords != 0)]

# remove the intercept
(top.coords <- top.coords[2:length(top.coords)])
top.values <- sort(top.values[2:length(top.values)])

# sanity check: plot the above/below threshold coefficients 
# df <- data.frame(x=1:length(top.values), y=top.values, threshold=as.factor(top.values != 0))
# levels(df$threshold) <- c("below", "above")
# p <- ggplot(df) +
#   geom_point(aes(x=x, y=y, color=threshold))
# show(p)

# get names associated with the balances; this is done by finding the finest level
# taxonomic identifier at which 95% (default) of components of the subtree agree
# setting thresh=1 requires all components of the subtree agree, so this will assign
# the finest taxonomic level to each 1/2 of the balance at which all components agree
tc.names <- sapply(top.coords, function(x) name.balance(tree, tax, x, thresh=1))

# let's build full list of components for disambiguation
# function `get.ud.tips` stolen from PhILR: use phangorn to get the tip names of children of this node

# Returns a list of the 'up' and 'down' subtree's values as a vector of tip ids (corresponds
# to up and down used for sbp creation)
# Each value is the ID of a tip
get.ud.tips <- function(tr,coord){
  l.tips <- list()
  child <- phangorn::Children(tr, name.to.nn(tr,coord))
  if (length(child) < 2) stop(paste0(coord,' is a tip'))
  if (length(child) > 2) stop("Tree is not soley binary.") #TODO: Bit of validation - consider better location
  l.tips[['up']] <- sapply(unlist(phangorn::Descendants(tr,child[1],type='tips')), function(x) nn.to.name(tr, x))
  l.tips[['down']] <- sapply(unlist(phangorn::Descendants(tr,child[2],type='tips')), function(x) nn.to.name(tr, x))
  return(l.tips)
}

full.labels <- list()
for(balance_id in top.coords) {
  all.tips <- get.ud.tips(tree, balance_id)
  full.labels[[balance_id]]$up <- c()
  for(tip_id in all.tips$up) {
    idx <- which(rownames(tax) == tip_id)
    label <- paste0("Family ",tax[idx,"Family"],
               ", Genus ",tax[idx,"Genus"],
               ", Species ",tax[idx,"Species"])
    full.labels[[balance_id]]$up <- c(full.labels[[balance_id]]$up, label)
  }
  full.labels[[balance_id]]$down <- c()
  for(tip_id in all.tips$down) {
    idx <- which(rownames(tax) == tip_id)
    label <- paste0("Family ",tax[idx,"Family"],
                    ", Genus ",tax[idx,"Genus"],
                    ", Species ",tax[idx,"Species"])
    full.labels[[balance_id]]$down <- c(full.labels[[balance_id]]$down, label)
  }
}

# optional: we can examine the 'votes' cast by child nodes in the determination of `tc.names``
# balance_id <- "n151"
# votes <- name.balance(tree, tax, balance_id, return.votes=c('up', 'down'), thresh=1)
# str(votes)
# tc.names[which(names(tc.names) == balance_id)] <- "Genus Eubacterium brachy group, species uncultured bacterium/Family XIII, genus AD3011 group, species uncultured bacterium"

# visualize differential balances
data.philr.long <- convert_to_long(data.philr, get_variable(data.final, 'Health')) %>%
  filter(coord %in% top.coords) %>%
  #mutate(readable_coord = str_replace(tc.names[coord], "/", "/\n"))
  mutate(readable_coord = tc.names[coord])
names(data.philr.long) <- c("sample", "cohort", "coord", "value", "readable_coord")
head(data.philr.long)

data.philr.stats <- data.philr.long %>%
  group_by(cohort, coord) %>%
  mutate(mean=mean(value), sd=sd(value)) %>%
  select(cohort, coord, sd, mean, readable_coord)

# bar plot differences and dump constituent taxa to a file
sink("out.txt")
for(coord in unique(data.philr.long$coord)) {
  ggplot(data.philr.long[data.philr.long$coord == coord,], aes(x=cohort, y=value)) +
    geom_boxplot(fill='lightgrey') +
    #facet_wrap(.~ coord, scales='free_y', nrow=1, ncol=5) +
    xlab('Health') + ylab('PhILR balance Value') +
    theme_bw()
  
  # dump the UP/DOWN constituent microbes for this balance/node
  cat("Node:",coord,"\n\n")
  cat("Short name:",tc.names[[coord]],"\n")
  cat("Constituents (top):\n")
  for(x in full.labels[[coord]]$up) {
    cat(paste0("\t",x,"\n"))
  }
  cat("Constituents (bottom):\n")
  for(x in full.labels[[coord]]$down) {
    cat(paste0("\t",x,"\n"))
  }
  cat("\n\n")
}
sink()

# save data to a file for Jess; Jess has been making the final version of this figure in PRISM
data.philr.out <- data.philr.long
data.philr.out$readable_coord <- sapply(data.philr.out$readable_coord, function(x) gsub('\n', '', x))
write.table(data.philr.out, file="F6_data.tsv", quote=FALSE, sep='\t', col.names = NA)

# alternative: crude version of dot plots with error bars
#  ggplot() +
#   geom_jitter(data=data.philr.long, aes(x=cohort, y=value),
#               position=position_jitter(0.2), size=2) +
#   geom_errorbar(data=data.philr.stats, aes(x=cohort, ymin=mean, ymax=mean), width = 0.5) +
#   geom_errorbar(data=data.philr.stats, aes(x=cohort, ymin=mean-sd, ymax=mean+sd), width = 0.3) +
#   facet_wrap(.~ readable_coord) +
#   theme_minimal() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme( axis.line = element_line(color="black", size=1, linetype="solid")) +
#   xlab("Health Status") +
#   ylab("PhILR value")

# visualize the position of these 'significant' balances on the phylogenetic tree
# label nodes of tree by (1) effect size in logistic regression (2) balance ratio

# get nodes at such-and-such depth from root (n1); we'll ignore tips at this level
if(FALSE) {
  depth <- 4
  nodelist <- c("n1")
  for(dd in 1:depth) {
    for(node in nodelist) {
      new_nodelist <- c()
      edges <- phangorn::Children(tree, name.to.nn(tree,node))
      if(length(edges) < 2) {
        # tip; ignore
      } else {
        # what to do about childless nodes??? (e.g. n44)
        left.children <- sapply(unlist(phangorn::Descendants(tree,edges[1],type='children')), function(x) nn.to.name(tree, x))
        right.children <- sapply(unlist(phangorn::Descendants(tree,edges[2],type='children')), function(x) nn.to.name(tree, x))
        new_nodelist <- c(new_nodelist, left.children)
        new_nodelist <- c(new_nodelist, right.children)
      }
    }
    nodelist <- new_nodelist
  }
}

# having some trouble (above) traversing the tree to automatically grab (there are some leaves and childless nodes [n44?] at
# shallow depth that are very time-consuming to work through parsing); it's easy enough to identify nodes by
# eye; we'll get descendants of these and find a consensus order from which to label clades
nodelist <- c("n5", "n34")

# get tips associated with these nodes
clade_labels <- list()
for(node in nodelist) {
  edges <- phangorn::Children(tree, name.to.nn(tree,node))
  left.tips <- sapply(unlist(phangorn::Descendants(tree,edges[1],type='tips')), function(x) nn.to.name(tree, x))
  right.tips <- sapply(unlist(phangorn::Descendants(tree,edges[2],type='tips')), function(x) nn.to.name(tree, x))
  tips <- c(left.tips, right.tips)
  order_labels <- c()
  for(tip in tips) {
    idx <- which(rownames(tax) == tip)
    order_labels <- c(order_labels, tax[idx,"Family"])
  }
  clade_labels[[node]] <- unique(order_labels)
}

# get +/- status of balance indexed by node label (for arrow image choice)
avg_balance <- data.philr.long %>%
  group_by(coord) %>%
  summarize(mean_coord = mean(value)) %>%
  arrange(desc(coord))

# temp copy of tree
plot_tree <- tree

# eliminate hexidecimal tax IDs from tree tips to prevent these being rendered
plot_tree$tip.label <- rep("", length(plot_tree$tip.label))

pg <- ggtree(plot_tree) # optional: branch.length='none'

# add up/down images to tree
node_idx <- which(tree$node.label %in% names(tc.names)) + (length(tree$edge.length) - tree$Nnode)
d <- data.frame(node=node_idx, images=rep(0, 5))
for(i in 1:length(node_idx)) {
  if(avg_balance[avg_balance$coord==names(tc.names)[i],]$mean_coord > 0) {
    d[d$node==node_idx[i],"images"] <- "up_arrow.png"
  } else {
    d[d$node==node_idx[i],"images"] <- "down_arrow.png"
  }
}
pg <- pg +
  theme_tree2() + 
  geom_tiplab(align=TRUE, linesize=0.01) +
  xlim(-0.1, 1.5)
pg <- pg %<+% d + geom_nodelab(aes(image=images), geom="image", size=0.02, alpha=0.5)

# label the tree for identified balances
if(FALSE) {
  for(i in 1:length(node_idx)) {
    pg <- pg +
          geom_cladelabel(node=node_idx[i], label=tc.names[i], color="black", offset=0, align=TRUE, fontsize=3)
  }
}

# label the tree by clade (family level)
for(node in nodelist) {
  node_label <- gsub("n", "", node)
  pg <- pg +
    #geom_hilight(node=as.numeric(node_label), fill="steelblue", alpha=0.5) +
    geom_cladelabel(node=node_label, label=clade_labels[[node]], color="black", offset=0, barsize=2, align=TRUE, fontsize=3)
}

plot(pg)










