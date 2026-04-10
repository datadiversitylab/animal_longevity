# Load packages
library(ape)
library(caper)
library(plotrix)

# Datasets
tree <- read.nexus((here("NIH_R01", "agemammals.trees")))

data <- read_csv((here("NIH_R01", "mammal_agediet.csv"))) 

# Match data to tree tips
rownames(data) <- data$species
data <- data[tree$tip.label, ]

# Deduplicate and clean for caper (keeps first occurrence)
data_caper <- data[!duplicated(data$species), ]
data_caper <- data_caper[complete.cases(data_caper[, c("mean_long", "logmass")]), ]

# Prune tree to match cleaned data
tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% data_caper$species])

# Build comparative data object for caper
comp_data <- comparative.data(
  phy = tree_pruned,
  data = data_caper,
  names.col = species,
  vcv = TRUE,
  na.omit = FALSE,
  warn.dropped = TRUE
)

# Fit PGLS model (log_longevity ~ body mass)
pgls_mod <- pgls(log(mean_long) ~ logmass, 
                 data = comp_data,
                 lambda = "ML")

# Extract residuals matched to tips
pgls_resid <- residuals(pgls_mod)

# Ancestral state reconstruction on log-scaled residuals
ace_resid <- ace(pgls_resid, tree_pruned, type = "continuous", method = "ML")

# Combine tip + node values in ape's node order
all_resid <- c(pgls_resid, ace_resid$ace)

# Mean parent-child residual per edge (for colors)
edge_vals <- sapply(1:nrow(tree_pruned$edge), function(i) {
  mean(all_resid[tree_pruned$edge[i, ]])
})

# Map edge values to palette
resid_pal <- colorRampPalette(c("#2166AC", "#B2182B"))
n_cols <- 20
edge_cols_idx <- cut(edge_vals, breaks = n_cols, labels = FALSE)
edge_cols <- resid_pal(n_cols)[edge_cols_idx]

# Matched data to tips
data_plot <- comp_data$data

# Define colors for diet
diet_cols <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A"),
  c("Carnivore", "Herbivore", "Omnivore")
)
tip_colors <- diet_cols[data_plot$diet.new]

# Create size vector for points
point_size <- scales::rescale(data_plot$mean_long, to = c(0.5, 2.5))

# Save as PDF
pdf("Figure2_Mammal_Phylogeny.pdf", width = 12, height = 12)

# Plot circular tree with residual-colored branches
plot(tree_pruned, 
     type = "fan", 
     show.tip.label = FALSE,
     edge.width = 1.5,
     edge.color = edge_cols,
     no.margin = FALSE)

# Add tip points
tiplabels(pch = 21, 
          col = "white",
          bg = tip_colors,
          cex = point_size,
          offset = 0.5)

# Add legend for diet
legend("topleft",
       legend = c("Carnivore", "Herbivore", "Omnivore"),
       pch = 21,
       pt.bg = diet_cols,
       col = "white",
       pt.cex = 2,
       bty = "n",
       title = "Diet",
       cex = 1.2)

# Add size reference
legend("topright",
       legend = c("Short", "Medium", "Long"),
       pch = 21,
       pt.bg = "grey50",
       col = "white",
       pt.cex = c(0.8, 1.5, 2.5),
       bty = "n",
       title = "Longevity",
       cex = 1.2)

# Add continuous color bar for PGLS residuals (log scale)
color.legend(
  xl = -12, yb = 10, xr = 12, yt = 15,
  legend = round(range(edge_vals), 2),
  rect.col = resid_pal(n_cols),
  align = "lt",
  gradient = "x",
  cex = 0.5
)
text(0, 35, "PGLS Residuals", cex = 1)

dev.off()
