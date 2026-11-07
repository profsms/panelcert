# Rebuilds data/*.rda from the canonical CSVs in data-raw/ (themselves copied
# from the repo-root testdata/ that all three language implementations share).
# Run from the package root:  Rscript data-raw/make_data.R
# The .rda object name must match the file base name for lazy-data loading.

vdem    <- utils::read.csv("data-raw/eiv_vdem_panel.csv",  stringsAsFactors = FALSE)
psid    <- utils::read.csv("data-raw/psid_wages_panel.csv", stringsAsFactors = FALSE)
twins   <- utils::read.csv("data-raw/twins.csv", stringsAsFactors = FALSE)
castle  <- utils::read.csv("data-raw/castle_panel.csv",    stringsAsFactors = FALSE)
divorce <- utils::read.csv("data-raw/divorce_panel.csv",   stringsAsFactors = FALSE)
minimum_wage <- utils::read.csv("data-raw/minimum_wage_panel.csv",
                                stringsAsFactors = FALSE)
brazil <- utils::read.csv("data-raw/brazil_property_tax_panel.csv",
                          stringsAsFactors = FALSE)
fscore  <- utils::read.csv("data-raw/f_score_panel.csv",   stringsAsFactors = FALSE)
# Stepped-wedge trial extracts (CC BY 4.0; built by data-raw/make_sw_trials.py)
haines2017 <- utils::read.csv("data-raw/haines2017_trial1.csv", stringsAsFactors = FALSE)
dots99     <- utils::read.csv("data-raw/dots99_uganda.csv", stringsAsFactors = FALSE)

if (!dir.exists("data")) dir.create("data")
save(vdem,    file = "data/vdem.rda",    compress = "xz")
save(psid,    file = "data/psid.rda",    compress = "xz")
save(twins,   file = "data/twins.rda",   compress = "xz")
save(castle,  file = "data/castle.rda",  compress = "xz")
save(divorce, file = "data/divorce.rda", compress = "xz")
save(minimum_wage, file = "data/minimum_wage.rda", compress = "xz")
save(brazil, file = "data/brazil.rda", compress = "xz")
save(fscore,  file = "data/fscore.rda",  compress = "xz")
save(haines2017, file = "data/haines2017.rda", compress = "xz")
save(dots99,     file = "data/dots99.rda",     compress = "xz")

message("built: ", paste(list.files("data"), collapse = ", "))
