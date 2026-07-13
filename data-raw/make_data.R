# Rebuilds data/*.rda from the canonical CSVs in data-raw/ (themselves copied
# from the repo-root testdata/ that all three language implementations share).
# Run from the package root:  Rscript data-raw/make_data.R
# The .rda object name must match the file base name for lazy-data loading.

vdem    <- utils::read.csv("data-raw/eiv_vdem_panel.csv",  stringsAsFactors = FALSE)
psid    <- utils::read.csv("data-raw/psid_wages_panel.csv", stringsAsFactors = FALSE)
castle  <- utils::read.csv("data-raw/castle_panel.csv",    stringsAsFactors = FALSE)
divorce <- utils::read.csv("data-raw/divorce_panel.csv",   stringsAsFactors = FALSE)
fscore  <- utils::read.csv("data-raw/f_score_panel.csv",   stringsAsFactors = FALSE)

if (!dir.exists("data")) dir.create("data")
save(vdem,    file = "data/vdem.rda",    compress = "xz")
save(psid,    file = "data/psid.rda",    compress = "xz")
save(castle,  file = "data/castle.rda",  compress = "xz")
save(divorce, file = "data/divorce.rda", compress = "xz")
save(fscore,  file = "data/fscore.rda",  compress = "xz")

message("built: ", paste(list.files("data"), collapse = ", "))
