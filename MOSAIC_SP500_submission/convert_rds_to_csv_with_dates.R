args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  normalizePath(getwd())
}

rds_path <- file.path(script_dir, "data", "stocks", "sp500-data-2016-2020.rds")
csv_path <- file.path(script_dir, "data", "stocks", "sp500-data-2016-2020_converted.csv")

cat("Reading:\n")
cat("  ", rds_path, "\n\n")

obj <- readRDS(rds_path)

cat("Object class:\n")
print(class(obj))

cat("\nObject dimensions:\n")
print(dim(obj))

idx <- attr(obj, "index")

if (is.null(idx)) {
  stop("No index attribute found. Cannot recover dates.")
}

dates <- as.POSIXct(idx, origin = "1970-01-01", tz = "UTC")
dates_date <- as.Date(dates)

cat("\nRecovered date range:\n")
cat("  first:", as.character(min(dates_date)), "\n")
cat("  last :", as.character(max(dates_date)), "\n")
cat("  number of dates:", length(dates_date), "\n\n")

if (min(dates_date) < as.Date("2015-01-01") || max(dates_date) > as.Date("2021-12-31")) {
  warning("Recovered dates are outside expected 2016-2020 range. Check carefully.")
} else {
  cat("SUCCESS: recovered dates look compatible with 2016-2020 data.\n\n")
}

# ------------------------------------------------------------
# Safe conversion from xts/zoo object to plain numeric matrix.
# This avoids as.data.frame.xts / dimnames issues.
# ------------------------------------------------------------
n <- nrow(obj)
p <- ncol(obj)

mat <- matrix(as.numeric(obj), nrow = n, ncol = p)
colnames(mat) <- colnames(obj)

df <- data.frame(
  Date = as.character(dates_date),
  mat,
  check.names = FALSE
)

cat("CSV dimensions will be:\n")
cat("  rows:", nrow(df), " columns:", ncol(df), "\n")
cat("First columns:\n")
print(colnames(df)[1:min(10, ncol(df))])

write.csv(df, csv_path, row.names = FALSE)

cat("\nSaved CSV with real Date column to:\n")
cat("  ", csv_path, "\n\n")

cat("First few rows:\n")
print(head(df[, 1:min(6, ncol(df))]))
