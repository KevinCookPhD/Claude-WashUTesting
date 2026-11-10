# find_talk_pptx.R
# Search a directory tree for the PowerPoint whose title slide mentions a given date.
#
# Usage (edit the settings below, then source / Rscript this file):
#   Rscript find_talk_pptx.R "C:/path/to/your/directory"
#
# Base R only; no packages needed. Works on .pptx/.pptm/.ppsx (not legacy .ppt).

# ---- Settings ---------------------------------------------------------------
args     <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) >= 1) args[1] else
  "C:/Users/Kevin/OneDrive - Childrens National Hospital/"
target   <- as.Date("2025-12-03")
# Also flag files whose internal created/modified date is within this many days
# of the target (useful if the date isn't printed on the title slide).
meta_window_days <- 14
out_csv  <- "pptx_title_slides.csv"
# -----------------------------------------------------------------------------

# Build regexes for many ways a date can be written
date_patterns <- function(d) {
  y  <- format(d, "%Y");  m <- as.integer(format(d, "%m")); dd <- as.integer(format(d, "%d"))
  mon_full <- month.name[m]; mon_abb <- month.abb[m]
  mon <- sprintf("(%s|%s)", mon_full, mon_abb)
  day <- sprintf("0?%d(st|nd|rd|th)?", dd)
  c(
    sprintf("%s\\.?\\s+%s,?\\s+%s", mon, day, y),               # December 3, 2025 / Dec. 3rd 2025
    sprintf("%s\\s+%s\\.?,?\\s+%s", day, mon, y),               # 3 December 2025
    sprintf("\\b0?%d[/.-]0?%d[/.-](%s|%s)\\b", m, dd, y, substr(y, 3, 4)),  # 12/3/2025, 12-03-25
    sprintf("\\b0?%d[/.-]0?%d[/.-]%s\\b", dd, m, y),            # 3/12/2025 (day-first)
    sprintf("\\b%s-0?%d-0?%d\\b", y, m, dd)                     # 2025-12-03
  )
}
pats <- date_patterns(target)

read_xml_text <- function(zipfile, member) {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  ok <- tryCatch(unzip(zipfile, files = member, exdir = tmp, junkpaths = FALSE),
                 error = function(e) character(0), warning = function(w) character(0))
  if (length(ok) == 0) return(NA_character_)
  paste(readLines(file.path(tmp, member), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

# Pull visible text (<a:t> runs) from slide XML; paragraphs (<a:p>) become line breaks
xml_text <- function(xml) {
  if (is.na(xml)) return("")
  xml <- gsub("</a:p>", "\n", xml, fixed = TRUE)
  runs <- regmatches(xml, gregexpr("<a:t(?: [^>]*)?>[^<]*</a:t>|\n", xml, perl = TRUE))[[1]]
  txt <- ifelse(runs == "\n", "\n", sub("^<a:t[^>]*>", "", sub("</a:t>$", "", runs)))
  txt <- paste(txt, collapse = " ")
  txt <- gsub("&amp;", "&", txt, fixed = TRUE)
  txt <- gsub("&lt;", "<", txt, fixed = TRUE)
  txt <- gsub("&gt;", ">", txt, fixed = TRUE)
  txt <- gsub("&quot;", "\"", txt, fixed = TRUE)
  txt <- gsub("&apos;", "'", txt, fixed = TRUE)
  trimws(gsub("[ \t]*\n[ \t\n]*", " | ", txt))
}

# Find the first slide in presentation order (not just "slide1.xml")
first_slide_member <- function(f) {
  pres <- read_xml_text(f, "ppt/presentation.xml")
  rels <- read_xml_text(f, "ppt/_rels/presentation.xml.rels")
  fallback <- "ppt/slides/slide1.xml"
  if (is.na(pres) || is.na(rels)) return(fallback)
  rid <- regmatches(pres, regexpr("<p:sldId [^>]*r:id=\"[^\"]+\"", pres))
  if (length(rid) == 0) return(fallback)
  rid <- sub(".*r:id=\"([^\"]+)\".*", "\\1", rid)
  rel_tags <- regmatches(rels, gregexpr("<Relationship [^>]*>", rels))[[1]]
  hit <- rel_tags[grepl(sprintf("Id=\"%s\"", rid), rel_tags, fixed = TRUE)]
  if (length(hit) == 0) return(fallback)
  tgt <- sub(".*Target=\"([^\"]+)\".*", "\\1", hit[1])
  tgt <- sub("^/", "", tgt)
  if (!startsWith(tgt, "ppt/")) tgt <- paste0("ppt/", tgt)
  tgt
}

get_meta_dates <- function(f) {
  core <- read_xml_text(f, "docProps/core.xml")
  g <- function(tag) {
    m <- regmatches(core, regexpr(sprintf("<dcterms:%s[^>]*>[^<]+", tag), core))
    if (length(m) == 0) return(as.Date(NA))
    as.Date(substr(sub(".*>", "", m), 1, 10))
  }
  if (is.na(core)) return(list(created = as.Date(NA), modified = as.Date(NA)))
  list(created = g("created"), modified = g("modified"))
}

files <- list.files(root_dir, pattern = "\\.(pptx|pptm|ppsx|potx)$",
                    recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
files <- files[!grepl("(^|/)~\\$", files)]   # skip PowerPoint temp/lock files
cat(sprintf("Scanning %d PowerPoint files under '%s' ...\n", length(files), root_dir))

rows <- vector("list", length(files))
for (i in seq_along(files)) {
  f <- files[i]
  res <- tryCatch({
    member <- first_slide_member(f)
    title_text <- xml_text(read_xml_text(f, member))
    md <- get_meta_dates(f)
    date_in_text <- any(vapply(pats, function(p) grepl(p, title_text, ignore.case = TRUE, perl = TRUE), logical(1)))
    near <- function(d) !is.na(d) && abs(as.numeric(d - target)) <= meta_window_days
    data.frame(file = f, title_slide_text = title_text,
               date_on_title_slide = date_in_text,
               meta_created = md$created, meta_modified = md$modified,
               meta_near_target = near(md$created) || near(md$modified),
               error = "", stringsAsFactors = FALSE)
  }, error = function(e) {
    data.frame(file = f, title_slide_text = "", date_on_title_slide = FALSE,
               meta_created = as.Date(NA), meta_modified = as.Date(NA),
               meta_near_target = FALSE, error = conditionMessage(e), stringsAsFactors = FALSE)
  })
  rows[[i]] <- res
}
results <- do.call(rbind, rows)

if (is.null(results) || nrow(results) == 0) {
  cat("No PowerPoint files found.\n")
} else {
  results <- results[order(-results$date_on_title_slide, -results$meta_near_target, results$file), ]
  write.csv(results, out_csv, row.names = FALSE)

  strong <- results[results$date_on_title_slide, ]
  weak   <- results[!results$date_on_title_slide & results$meta_near_target, ]

  cat(sprintf("\n=== Date %s found on title slide: %d file(s) ===\n", format(target), nrow(strong)))
  for (k in seq_len(nrow(strong)))
    cat(sprintf("\n%s\n   %s\n", strong$file[k], substr(strong$title_slide_text[k], 1, 200)))

  cat(sprintf("\n=== Not on title slide, but created/modified within %d days of target: %d file(s) ===\n",
              meta_window_days, nrow(weak)))
  for (k in seq_len(nrow(weak)))
    cat(sprintf("\n%s\n   created %s, modified %s\n   %s\n", weak$file[k],
                weak$meta_created[k], weak$meta_modified[k], substr(weak$title_slide_text[k], 1, 200)))

  if (sum(results$error != "") > 0)
    cat(sprintf("\n%d file(s) could not be read (see 'error' column in %s).\n", sum(results$error != ""), out_csv))
  cat(sprintf("\nFull table of every title slide saved to: %s\n", normalizePath(out_csv)))
}
