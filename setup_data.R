## ============================================================================
## setup_data.R -- STEP 1 of every route through this repository.
##
## The data are shipped as CSV (data/*.csv, data/*.csv.gz). The R scripts read
## binary .rds files. This script builds one data/<name>.rds from the CSV file(s)
## of each table, with every number restored bit-for-bit (lib/data_io.R explains
## how). Nothing is computed; it takes about two minutes.
##
## The .rds files are a local cache and are not under git (.gitignore). A file
## that already exists is left alone, so a table you regenerated with a pipeline
## stage is not overwritten; FORCE=1 rebuilds everything from the CSV.
##
## in : ./data/*.csv, ./data/*.csv.gz, ./data/schema.csv
## out: ./data/<name>.rds for every table in data/schema.csv
## RUN: Rscript setup_data.R            (from the repository root)
##      FORCE=1 Rscript setup_data.R    (rebuild all)
##      ONLY=selected_big,rescue_sterile Rscript setup_data.R
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
source("./lib/data_io.R")
FORCE = nzchar(Sys.getenv("FORCE", ""))
ONLY  = strsplit(Sys.getenv("ONLY", ""), ",")[[1]]

if (!file.exists(DATA_SCHEMA)) stop(DATA_SCHEMA, " not found -- run this from the repository root")
## one object per table; the members of a list (<name>__<member>.csv) form one object
OBJ = unique(sub("__.*$", "", sub("\\.csv$", "", .schema_read()$file)))
if (length(ONLY)) OBJ = intersect(OBJ, ONLY)

t0 = Sys.time(); n_made = 0L
for (o in OBJ) {
  out = sprintf("./data/%s.rds", o)
  if (file.exists(out) && !FORCE) { cat(sprintf("  kept   %s\n", out)); next }
  stem = file.path("./data", o)
  if (!length(data_csv_files(stem))) stop("the CSV of ", o, " is missing from ./data")
  t1 = Sys.time()
  saveRDS(read_data_csv(stem), out)
  n_made = n_made + 1L
  cat(sprintf("  wrote  %-58s %5.1f s\n", out, as.numeric(difftime(Sys.time(), t1, units = "secs"))))
}
cat(sprintf("\n%d of %d tables written to ./data as .rds in %.1f min. Next: bash make_figures.sh\n",
            n_made, length(OBJ), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
