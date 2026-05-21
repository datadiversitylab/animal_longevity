library(ganttrify)
library(ggplot2)

project <- data.frame(
  wp = c(
    "Aim 1", "Aim 1", "Aim 1",
    "Aim 2", "Aim 2", "Aim 2",
    "Aim 3", "Aim 3", "Aim 3"
  ),
  activity = c(
    "Phylogenetic signal, body size regression & outlier identification",
    "Ancestral state reconstruction",
    "Focal species + matched controls confirmed",
    "phyca gene content extraction + dN/dS estimation",
    "Ensemble classifiers + SHAP interpretation",
    "Hierarchical clustering and convergence tests",
    "Bayesian predictive models (brms)",
    "Longevity predictions for new species",
    "Public database and web interface"
  ),
  start_date = c(
    1,  7, 10,
    13, 19, 25,
    25, 37, 37
  ),
  end_date = c(
    9,  18, 21,
    24, 27, 33,
    42, 48, 60
  ),
  stringsAsFactors = FALSE
)


spots <- data.frame(
  activity   = c(
    "Focal species + matched controls confirmed",
    "Hierarchical clustering and convergence tests",
    "Longevity predictions for new species"
  ),
  spot_type  = c("M", "M", "M"),
  spot_date  = c(22, 34, 48),
  stringsAsFactors = FALSE
)

p <- ganttrify(
  project            = project,
  spots              = spots,
  project_start_date = "2026-01",
  mark_quarters      = TRUE,
  month_breaks       = 3,
  size_text_relative = 0.9,
  alpha_wp           = 0.9,
  alpha_activity     = 0.7,
  line_end_wp        = "round",
  line_end_activity  = "round",
  axis_text_align    = "left",
  label_wrap         = 45,
  colour_palette     = c(
    "#4A90D9",
    "#5BAD72",  
    "#E09B3D"  
  ),
  font_family        = "sans"
) +
  theme(
    plot.title    = element_text(size = 13, face = "bold",  color = "grey20"),
    plot.subtitle = element_text(size =  9, color = "grey40"),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  filename = "Timeline_R01.pdf",
  plot     = p,
  width    = 10,
  height   = 4
)

ggsave(
  filename = "Timeline_R01.jpg",
  plot     = p,
  width    = 10,
  height   = 4
)
