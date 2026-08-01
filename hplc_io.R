# -----------------------------------------------------------------------------
# Shimadzu LabSolutions ASCII 匯出檔解析器（共用模組）
#
# 由 caffeine_calibration.R 與 align_chromatograms.R 共同 source。
# 僅使用 base R。
# -----------------------------------------------------------------------------

read_lines_big5 <- function(path) {
  # 檔案為 Big5 (繁中 Windows)。注意 readLines(encoding=) 只是「標記」而不轉碼，
  # 直接用會讓含中文的列在 strsplit 時變成 NA，必須以 iconv 實際轉成 UTF-8。
  raw <- readLines(path, warn = FALSE, encoding = "bytes")
  out <- iconv(raw, from = "BIG5", to = "UTF-8", sub = "?")
  ifelse(is.na(out), iconv(raw, from = "latin1", to = "UTF-8", sub = "?"), out)
}

split_sections <- function(lines) {
  hits <- grep("^\\[.+\\]\\s*$", lines)
  if (!length(hits)) return(list())
  names_ <- gsub("^\\[|\\]\\s*$", "", lines[hits])
  ends <- c(hits[-1] - 1L, length(lines))
  out <- lapply(seq_along(hits), function(i) {
    if (hits[i] + 1L > ends[i]) character(0) else lines[(hits[i] + 1L):ends[i]]
  })
  names(out) <- names_
  out
}

parse_kv <- function(lines) {
  lines <- lines[grepl("\t", lines)]
  if (!length(lines)) return(list())
  parts <- strsplit(lines, "\t", fixed = TRUE)
  setNames(lapply(parts, function(p) trimws(paste(p[-1], collapse = "\t"))),
           trimws(vapply(parts, `[`, "", 1)))
}

parse_peak_table <- function(lines) {
  if (!length(lines)) return(NULL)
  np <- grep("^# of Peaks", lines)
  n_peaks <- if (length(np)) as.integer(strsplit(lines[np[1]], "\t")[[1]][2]) else 0L
  hdr <- grep("^Peak#", lines)
  if (!length(hdr) || is.na(n_peaks) || n_peaks == 0L) return(NULL)
  cols <- strsplit(lines[hdr[1]], "\t", fixed = TRUE)[[1]]
  body <- lines[(hdr[1] + 1L):(hdr[1] + n_peaks)]
  mat <- do.call(rbind, lapply(body, function(l) {
    p <- strsplit(l, "\t", fixed = TRUE)[[1]]
    length(p) <- length(cols); p
  }))
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(df) <- cols
  for (cn in c("R.Time", "I.Time", "F.Time", "Area", "Height", "Conc."))
    if (cn %in% names(df)) df[[cn]] <- suppressWarnings(as.numeric(df[[cn]]))
  df
}

parse_chrom <- function(lines) {
  if (!length(lines)) return(NULL)
  hdr <- grep("^R.Time", lines)
  if (!length(hdr)) return(NULL)
  body <- lines[(hdr[1] + 1L):length(lines)]
  body <- body[grepl("^[-0-9.]+\t[-0-9.]+\\s*$", body)]
  if (!length(body)) return(NULL)
  m <- do.call(rbind, strsplit(body, "\t", fixed = TRUE))
  data.frame(time_min = as.numeric(m[, 1]), intensity = as.numeric(m[, 2]))
}

#' 讀入一個 LabSolutions 匯出檔
#' @param path 檔案完整路徑
load_run_file <- function(path) {
  if (!file.exists(path)) stop("找不到檔案: ", path)
  sec <- split_sections(read_lines_big5(path))
  list(file  = basename(path),
       info  = parse_kv(sec[["Sample Information"]]),
       peaks = parse_peak_table(sec[["Peak Table(Detector A-Ch1)"]]),
       chrom = parse_chrom(sec[["LC Chromatogram(Detector A-Ch1)"]]))
}

#' 以峰起訖點連線為基線，梯形法積分 (uV*sec)
trapz_baseline <- function(chrom, t0, t1) {
  if (is.null(chrom)) return(NA_real_)
  seg <- chrom[chrom$time_min >= t0 & chrom$time_min <= t1, , drop = FALSE]
  if (nrow(seg) < 3) return(NA_real_)
  t <- seg$time_min * 60
  y <- seg$intensity
  bl <- seq(y[1], y[length(y)], length.out = length(y))
  d <- y - bl
  sum(diff(t) * (head(d, -1) + tail(d, -1)) / 2)
}
