## ============================================================================
## lib/data_io.R -- the CSV <-> R table layer of the repository.
##
## The repository ships its data as CSV (data/*.csv, data/*.csv.gz) because CSV
## is an open format any tool can read. The R scripts read binary .rds files,
## which load fast and hold every number bit-for-bit. This file converts between
## the two WITHOUT LOSING A SINGLE BIT, which ordinary CSV writing does not do:
## write.csv() and data.table::fwrite() keep 15 significant digits, and a
## double needs up to 17. The scripts join tables on double-valued keys
## (beta_ros, r_L_P), so a number that comes back one unit in the last place
## off silently drops rows from a join.
##
## How a number is written (num_to_text):
##   * with 15 significant digits ("%.15g") when that text reads back as the
##     same double -- 0.0015 is written 0.0015, 15 is written 15;
##   * otherwise with exactly 17 significant digits ("%#.17g"), which identifies
##     one double and no other.
## How it is read back (text_to_num): R's own text-to-number conversion can be
## off in the last places for 17-digit text (it is on Apple silicon, where long
## double is double; at 1e-298 it is off in the 11th digit). So every value
## read is printed again and compared with the text in the file; where the two
## differ, the value is corrected until it prints exactly as written.
##
## data/schema.csv records, for every column of every file, its R class
## (integer, numeric, logical, character, or list:<class> for a list column,
## whose cells are written as values joined by "|"), so that an integer-valued
## double column comes back as double and an empty text cell as "" or NA as it
## was. An R list of tables (e.g. list(sm = , thresh = , settings = )) is
## written as one file per member, <name>__<member>.csv. A table whose
## compressed file would exceed PART_MAX bytes is cut by rows into
## <name>.part1of4.csv.gz, ...; GitHub refuses files above 100 MB.
##
##   write_data_csv(x, "./data/name")   table or list  ->  CSV file(s) + schema rows
##   read_data_csv("./data/name")       CSV file(s)    ->  the same R object
##   data_identical(a, b)               same columns, classes, keys and bits
## Used by setup_data.R (CSV -> .rds) and pipeline/90_export_csv.R (.rds -> CSV,
## with a read-back check that the result is identical).
## ============================================================================
suppressMessages(library(data.table))

DATA_SCHEMA = "./data/schema.csv"
GZ_ABOVE    = 8e6        # plain .csv below this many bytes, .csv.gz above
PART_MAX    = 45e6       # a .csv.gz above this many bytes is written in parts
LIST_SEP    = "|"

## ---- numbers <-> text, exactly ---------------------------------------------
num_to_text = function(x) {
  s = sprintf("%.15g", x)
  fin = which(is.finite(x))
  bad = fin[suppressWarnings(as.numeric(s[fin])) != x[fin]]
  if (length(bad)) s[bad] = sprintf("%#.17g", x[bad])
  s[is.na(x) & !is.nan(x)] = NA_character_        # NA -> empty cell; NaN, Inf, -Inf as text
  s
}

## mantissa and decimal exponent of printed numbers ("1.5e-298" -> 1.5, -298)
.mant_exp = function(s) {
  e = regexpr("e", s, fixed = TRUE)
  list(m = as.numeric(ifelse(e > 0, substr(s, 1, e - 1), s)),
       e = ifelse(e > 0, suppressWarnings(as.integer(substr(s, e + 1, nchar(s)))), 0L))
}
text_to_num = function(s) {
  y = suppressWarnings(as.numeric(s))
  fin = which(is.finite(y))
  bad = fin[sprintf("%#.17g", y[fin]) != s[fin]]
  if (length(bad)) bad = bad[sprintf("%.15g", y[bad]) != s[bad]]
  if (length(bad)) {
    ## 17-digit text the parser missed. Far from 1 (1e-298) R can be off by 1e-11
    ## relative, so first correct by the ratio of the written to the printed value ...
    a = .mant_exp(s[bad]); b = .mant_exp(sprintf("%#.17g", y[bad]))
    y[bad] = y[bad] * ((a$m / b$m) * 10^(a$e - b$e))
    bad = bad[sprintf("%#.17g", y[bad]) != s[bad]]
  }
  if (length(bad)) {
    ## ... then walk the neighbouring doubles until one prints exactly as written
    ## (in half units of the last place, so that a value just below a power of
    ## two is not stepped over)
    ulp = 2^(pmax(floor(log2(abs(y[bad]))), -1021) - 53)
    for (k in as.vector(rbind(1:64, -(1:64)))) {                 # 1, -1, 2, -2, ...
      cand = y[bad] + k * ulp
      hit = sprintf("%#.17g", cand) == s[bad]
      y[bad[hit]] = cand[hit]
      bad = bad[!hit]; ulp = ulp[!hit]
      if (!length(bad)) break
    }
    if (length(bad)) stop(length(bad), " numbers could not be read back exactly, e.g. ", s[bad[1]])
  }
  y
}

## ---- list columns <-> "a|b|c" ----------------------------------------------
## a cell holding a zero-length vector is written as the text "()" and a NULL
## cell as an empty cell
.col_class = function(v) {
  if (is.list(v)) {
    ## a one-dimensional array cell is written as the plain vector it holds
    cl = unique(vapply(v, function(e) if (length(dim(e)) == 1L) class(as.vector(e))[1] else class(e)[1], ""))
    cl = setdiff(cl, "NULL")
    if (length(cl) == 0L) cl = "logical"
    if (length(cl) != 1L || !cl %in% c("integer", "numeric", "logical", "character"))
      stop("list column with cells of class ", paste(cl, collapse = "/"), " cannot be written as CSV")
    return(paste0("list:", cl))
  }
  cl = class(v)[1]
  if (!cl %in% c("integer", "numeric", "logical", "character")) stop("column of class ", cl, " cannot be written as CSV")
  cl
}
.atom_to_text = function(v, cl) switch(cl, numeric = num_to_text(v), as.character(v))
.text_to_atom = function(s, cl) switch(cl, numeric = text_to_num(s), integer = as.integer(s),
                                       logical = as.logical(s), character = s)
.list_to_text = function(v, cl) {
  n = lengths(v); isnull = vapply(v, is.null, NA)
  flat = .atom_to_text(unlist(v, use.names = FALSE), cl)
  flat[is.na(flat)] = "NA"
  if (cl == "character" && any(grepl(LIST_SEP, flat, fixed = TRUE))) stop("a list cell contains '", LIST_SEP, "'")
  out = rep("()", length(v)); out[isnull] = NA_character_
  if (length(flat)) {
    g = rep(seq_along(v), n)
    j = data.table(g = g, f = flat)[, .(s = paste(f, collapse = LIST_SEP)), by = g]
    out[j$g] = j$s
  }
  out
}
.text_to_list = function(s, cl) {
  out = vector("list", length(s))
  empty = !is.na(s) & s == "()"
  out[empty] = list(.text_to_atom(character(0), cl))
  i = which(!is.na(s) & !empty)
  if (length(i)) {
    parts = strsplit(s[i], LIST_SEP, fixed = TRUE)
    n = lengths(parts); flat = unlist(parts, use.names = FALSE)
    flat[flat == "NA"] = NA_character_
    val = .text_to_atom(flat, cl)
    out[i] = unname(split(val, factor(rep(seq_along(i), n), levels = seq_along(i))))
  }
  out
}

## ---- schema -----------------------------------------------------------------
.schema_read = function() if (file.exists(DATA_SCHEMA))
  fread(DATA_SCHEMA, colClasses = c("character", "character", "character", "integer", "character"), na.strings = "") else
  data.table(file = character(), column = character(), class = character(), key_pos = integer(), table_class = character())
.schema_put = function(fn, rows) {
  S = .schema_read(); S = S[S$file != fn]
  S = rbind(S, rows)
  fwrite(S[order(file)], DATA_SCHEMA)
}
## the schema name of a CSV file: <name>.csv, whatever its compression or part
.schema_name = function(f) sub("(\\.part[0-9]+of[0-9]+)?\\.csv(\\.gz)?$", ".csv", basename(f))

## ---- gzip without extra packages -------------------------------------------
.gunzip_tmp = function(f) {
  tmp = tempfile(fileext = ".csv"); i = gzfile(f, "rb"); o = file(tmp, "wb")
  on.exit({ close(i); close(o) })
  repeat { b = readBin(i, raw(), 2^24); if (!length(b)) break; writeBin(b, o) }
  tmp
}

## ---- one table ---------------------------------------------------------------
.table_files = function(stem) c(paste0(stem, c(".csv", ".csv.gz")), Sys.glob(paste0(stem, ".part*of*.csv.gz")))
.write_table = function(x, stem) {
  tc = if (is.data.table(x)) "data.table" else "data.frame"
  key = if (is.data.table(x)) key(x) else NULL
  x = as.data.table(x)
  cls = vapply(x, .col_class, "")
  txt = as.data.table(lapply(seq_along(x), function(j) {
    v = x[[j]]; cl = cls[[j]]
    if (startsWith(cl, "list:")) .list_to_text(v, sub("list:", "", cl))
    else if (cl == "numeric") num_to_text(v) else v }))
  setnames(txt, names(x))
  old = .table_files(stem); file.remove(old[file.exists(old)])
  f = paste0(stem, ".csv"); fwrite(txt, f)
  if (file.size(f) > GZ_ABOVE) { file.remove(f); f = paste0(stem, ".csv.gz"); fwrite(txt, f, compress = "gzip") }
  if (file.size(f) > PART_MAX) {
    n = ceiling(file.size(f) / PART_MAX); file.remove(f)
    rows = split(seq_len(nrow(txt)), ceiling(seq_len(nrow(txt)) * n / nrow(txt)))
    f = vapply(seq_len(n), function(i) { g = sprintf("%s.part%dof%d.csv.gz", stem, i, n)
      fwrite(txt[rows[[i]]], g, compress = "gzip"); g }, "")
  }
  fn = paste0(basename(stem), ".csv")
  .schema_put(fn, data.table(file = fn, column = names(x), class = unname(cls),
                             key_pos = match(names(x), key), table_class = tc))
  f
}
## f: one file, or the parts of one table in order
.read_table = function(f) {
  fn = .schema_name(f[1])
  S = .schema_read(); S = S[S$file == fn]
  if (!nrow(S)) stop("no rows for ", fn, " in ", DATA_SCHEMA)
  ## numbers and list cells are read as text and converted here; the rest by fread
  cc = ifelse(S$class %in% c("integer", "logical"), S$class, "character")
  x = rbindlist(lapply(f, function(g) {
    src = if (grepl("\\.gz$", g)) .gunzip_tmp(g) else g
    if (src != g) on.exit(file.remove(src))
    fread(src, colClasses = cc, na.strings = "", sep = ",", header = TRUE) }))
  if (!identical(names(x), S$column)) stop(fn, ": the columns do not match ", DATA_SCHEMA)
  for (j in seq_len(nrow(S))) {
    cl = S$class[j]
    if (cl == "numeric") set(x, j = j, value = text_to_num(x[[j]]))
    else if (startsWith(cl, "list:")) set(x, j = j, value = .text_to_list(x[[j]], sub("list:", "", cl)))
  }
  k = S[!is.na(key_pos)][order(key_pos), column]
  if (length(k)) setattr(x, "sorted", k)
  if (S$table_class[1] == "data.frame") x = as.data.frame(x)
  x
}

## ---- public: a table, or a list of tables and plain lists ---------------------
## A list member that is itself a plain list of vectors (settings, cache keys)
## is written as a one-row table of list columns, <name>__<member>.list.csv.
write_data_csv = function(x, stem) {
  if (is.data.frame(x)) return(.write_table(x, stem))
  if (!is.list(x) || is.null(names(x))) stop("write_data_csv: a table or a named list of tables")
  unlist(lapply(names(x), function(nm) {
    m = x[[nm]]
    if (is.data.frame(m)) .write_table(m, paste0(stem, "__", nm))
    else .write_table(as.data.table(lapply(m, list)), paste0(stem, "__", nm, ".list"))
  }))
}
## the CSV files of one object: a single table (possibly in parts), or the members of a list
.parts_sorted = function(f) f[order(as.integer(sub("^.*\\.part([0-9]+)of[0-9]+\\.csv\\.gz$", "\\1", f)))]
data_csv_files = function(stem) {
  one = paste0(stem, c(".csv", ".csv.gz")); one = one[file.exists(one)]
  if (length(one)) return(one[1])
  p = Sys.glob(paste0(stem, ".part*of*.csv.gz"))
  if (length(p)) return(.parts_sorted(p))
  sort(Sys.glob(paste0(stem, "__*.csv*")))
}
read_data_csv = function(stem) {
  f = data_csv_files(stem)
  if (!length(f)) stop("no CSV for ", stem)
  if (!any(grepl("__", basename(f), fixed = TRUE))) return(.read_table(f))
  ## members: <stem>__<member>[.list][.partIofN].csv[.gz]
  member = sub("(\\.list)?(\\.part[0-9]+of[0-9]+)?\\.csv(\\.gz)?$", "", sub("^.*__", "", basename(f)))
  out = lapply(unique(member), function(m) { g = f[member == m]
    if (length(g) > 1L) g = .parts_sorted(g)
    t = .read_table(g)
    if (grepl("\\.list\\.csv", g[1])) lapply(as.list(t), `[[`, 1L) else t })
  setNames(out, unique(member))
}

## ---- two R objects hold the same data: same columns, classes, keys and bits ----
## list cells that are one-dimensional arrays count as the plain vectors they hold
.plain_cells = function(v) {
  if (!is.list(v)) return(v)
  d = which(vapply(v, function(e) length(dim(e)) == 1L, NA))
  if (length(d)) v[d] = lapply(v[d], as.vector)
  v
}
data_identical = function(a, b) {
  if (is.data.frame(a) != is.data.frame(b)) return(FALSE)
  if (!is.list(a) || !is.list(b)) return(identical(a, b))
  if (!is.data.frame(a))                     # a list of members: by name, in any order
    return(setequal(names(a), names(b)) && all(vapply(names(a), function(n) data_identical(a[[n]], b[[n]]), NA)))
  identical(class(a), class(b)) && identical(names(a), names(b)) && nrow(a) == nrow(b) &&
    identical(attr(a, "sorted"), attr(b, "sorted")) &&
    all(vapply(seq_along(a), function(j) identical(.plain_cells(a[[j]]), .plain_cells(b[[j]])), NA))
}
