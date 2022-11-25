# MULTISPECIES RESULTS -----
# Read Jaccard metrics ----
res <- read.csv("results/jaccard_pa_lgcp.csv")
res[res[,1] == "current", ]

# Read thresholded differences ----
diffs <- read.csv("results/thresh_diffs_pa_lgcp.csv")
diffs[diffs$species == "acch", ]

# Croped version (to species just from the south)
diffs.cr <- read.csv("results/thresh_diffs_south_pa_lgcp.csv")
diffs.cr[diffs.cr$species == "spam", ]

# INDIVIDUAL SPECIES RESULTS ----
# Define species
sp <- "mybo"

# Read summaries ----
summs <- read.csv(paste0("results/", sp, "/", sp, "_model_summary.csv"),
                  row.names = 1)
round(summs, 2)

# Read CV metrics ----
cv <- read.csv(paste0("results/", sp, "/", sp, "_cv_metrics.csv"),
               row.names = 1)
cv

cv.pa <- read.csv(paste0("results/", sp, "/", sp, "_pa_cv_metrics.csv"),
               row.names = 1)
cv.pa