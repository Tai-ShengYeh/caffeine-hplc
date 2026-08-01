# -----------------------------------------------------------------------------
# 用 alignDE 對 HPLC 層析圖做滯留時間對齊 (peak alignment)，再重新定量
#
# 方法 (alignDE; Zhang, Chen & Liang 2011 [1]):
#   1. CWT (Mexican Hat) 偵測峰位 -> Haar 小波估計峰寬          [2]
#   2. 懲罰最小平方 (Whittaker) 擬合基線並扣除                   [3]
#   3. 相鄰過近的峰合併成群 (peakClustering)
#   4. 差分演化 (DEoptim) 搜尋每群峰的位移，最大化與參考圖的相關係數  [4][5]
#      峰內不變形、峰間以線性內插伸縮 -> 峰形/峰高/峰面積保留
#
# 對齊後改用「所有檔案共用的固定積分視窗」重新積分，再建檢量線並定量。
#
# 參考文獻:
#   [1] Zhang Z-M, Chen S, Liang Y-Z (2011) Talanta 83(4):1108-1117.
#       doi:10.1016/j.talanta.2010.08.008
#   [2] Du P, Kibbe WA, Lin SM (2006) Bioinformatics 22(17):2059-2065.
#       doi:10.1093/bioinformatics/btl355
#   [3] Eilers PHC (2003) Anal Chem 75(14):3631-3636. doi:10.1021/ac034173t
#   [4] Storn R, Price K (1997) J Glob Optim 11(4):341-359.
#       doi:10.1023/A:1008202821328
#   [5] Mullen K et al. (2011) J Stat Softw 40(6):1-26. doi:10.18637/jss.v040.i06
#   對照方法 COW: Nielsen N-PV et al. (1998) J Chromatogr A 805(1-2):17-35.
#       doi:10.1016/S0021-9673(98)00021-1
#
# 用法:  Rscript align_chromatograms.R
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(alignDE))
options(stringsAsFactors = FALSE)

base_dir <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
out_dir <- file.path(base_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE)
source(file.path(base_dir, "hplc_io.R"))

# ---- 設定 -------------------------------------------------------------------
standards <- c("std-0.txt" = 0, "std-25.txt" = 25, "std-50.txt" = 50,
               "std-100.txt" = 100, "std-150.txt" = 150, "std-200.txt" = 200)
samples   <- c("sample-1.txt", "sample-2.txt", "sample-3.txt", "sample-4.txt")
all_files <- c(names(standards), samples)

REFERENCE <- "std-100.txt"   # 參考圖：中間濃度、峰形乾淨、單一主峰

DT_SEC  <- 0.5               # 取樣間隔
T_END   <- 6.0               # 跑程長度 (min)
SCALES  <- seq(1, 56, 1)     # CWT 尺度

# alignDE 參數（依 README 的選參指引）
SNR_TH      <- 3      # 峰偵測訊噪比門檻
RIDGE_LEN   <- 5      # ridge 最短長度
BL_LAMBDA   <- 100    # 基線平滑度
BL_THRESH   <- 0.3    # 峰形門檻
CLUSTER_GAP <- 5      # 峰群合併間距
SLACK       <- 45     # 最大可位移點數 (45 點 = 0.375 min，實測最大漂移約 30 點)
GROUP_N     <- 4      # 同時最佳化的峰數
DE_NP       <- 90     # 族群大小，約 2*slack
DE_ITERMAX  <- 200

set.seed(20161017)    # DEoptim 為隨機演算法，固定種子確保可重現

# ---- 讀入 -------------------------------------------------------------------
cat("讀入層析圖 ...\n")
runs <- lapply(all_files, function(fn) load_run_file(file.path(base_dir, fn)))
names(runs) <- all_files
N <- nrow(runs[[1]]$chrom)
tvec <- runs[[1]]$chrom$time_min
stopifnot(all(vapply(runs, function(r) nrow(r$chrom), 0) == N))

# ---- 1~2. 峰偵測 + 基線扣除 -------------------------------------------------
prepare <- function(p) {
  wCoefs    <- cwt(p, scales = SCALES, wavelet = "mexh")
  ridgeList <- getRidge(getLocalMaximumCWT(wCoefs), gapTh = 3, skip = 2)
  info      <- identifyMajorPeaks(p, ridgeList, wCoefs,
                                  SNR.Th = SNR_TH, ridgeLength = RIDGE_LEN)
  peakWidth <- widthEstimationCWT(p, info)
  list(signal    = p - baselineCorrectionCWT(p, peakWidth,
                                             threshold = BL_THRESH, lambda = BL_LAMBDA),
       peakWidth = peakWidth)
}

cat("CWT 峰偵測與基線扣除 ...\n")
prep <- lapply(all_files, function(fn) {
  x <- prepare(runs[[fn]]$chrom$intensity)
  cat(sprintf("  %-13s 偵測到 %d 個峰\n", fn, length(x$peakWidth$peakIndex)))
  x
})
names(prep) <- all_files

REF <- prep[[REFERENCE]]$signal

# ---- 3~4. 對齊 --------------------------------------------------------------
cat(sprintf("\n以 %s 為參考，差分演化對齊 (slack=%d, NP=%d, itermax=%d) ...\n",
            REFERENCE, SLACK, DE_NP, DE_ITERMAX))

fit_len <- function(v, n) {
  # alignDE 的輸出長度可能與輸入不同（峰間線性伸縮），統一裁切/補齊
  if (length(v) >= n) return(v[seq_len(n)])
  c(v, rep(v[length(v)], n - length(v)))
}

aligned <- list()
for (fn in all_files) {
  if (fn == REFERENCE) { aligned[[fn]] <- REF; next }
  pw <- peakClustering(prep[[fn]]$peakWidth, n = CLUSTER_GAP)
  t0 <- Sys.time()
  a <- tryCatch(
    alignDE(prep[[fn]]$signal, pw, REF, slack = SLACK, n = GROUP_N,
            control = list(NP = DE_NP, itermax = DE_ITERMAX, trace = FALSE)),
    error = function(e) { cat("    ! ", conditionMessage(e), "\n"); prep[[fn]]$signal })
  aligned[[fn]] <- fit_len(a, N)
  cat(sprintf("  %-13s %d 個峰群, %.1f 秒\n", fn, length(pw$peakIndex),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# ---- 診斷：相關係數與咖啡因峰頂位置 -----------------------------------------
caf_lo <- 4.0; caf_hi <- 5.2                       # 找峰頂用的粗略區間
caf_idx <- which(tvec >= caf_lo & tvec <= caf_hi)
apex_i  <- function(v) caf_idx[which.max(v[caf_idx])]
apex_rt <- function(v) tvec[apex_i(v)]

# std-0 是空白，該區間沒有咖啡因峰，峰頂只是雜訊 -> 不納入診斷統計
has_peak <- all_files != "std-0.txt"

ref_apex <- apex_rt(REF)
diag <- data.frame(
  file      = all_files,
  cor_before = vapply(all_files, function(fn) similarity(REF, prep[[fn]]$signal), 0),
  cor_after  = vapply(all_files, function(fn) similarity(REF, aligned[[fn]]), 0),
  apex_before = vapply(all_files, function(fn) apex_rt(prep[[fn]]$signal), 0),
  apex_after  = vapply(all_files, function(fn) apex_rt(aligned[[fn]]), 0)
)
diag$dev_before <- ifelse(has_peak, diag$apex_before - ref_apex, NA)
diag$dev_after  <- ifelse(has_peak, diag$apex_after  - ref_apex, NA)
diag[!has_peak, c("apex_before", "apex_after")] <- NA   # 空白無峰

# ---- 固定積分視窗（取所有對齊後圖譜峰邊界的聯集）----------------------------
# 樣品峰比標準品寬（A/H 17-19 vs 11），若只用參考圖決定視窗會把樣品峰切掉，
# 因此改用聯集；並限制最多離參考峰頂 ±1.0 min，避免碰到鄰近的基質峰。
CAP <- round(1.0 * (N - 1) / T_END)
ref_apex_i <- apex_i(REF)
peak_bounds <- function(v, frac = 0.01) {
  ap <- apex_i(v)
  thr <- frac * v[ap]
  a <- ap; while (a > max(1, ref_apex_i - CAP) && v[a] > thr) a <- a - 1
  b <- ap; while (b < min(N, ref_apex_i + CAP) && v[b] > thr) b <- b + 1
  c(a, b)
}
bounds <- t(vapply(all_files[has_peak],
                   function(fn) peak_bounds(aligned[[fn]]), c(0, 0)))
cat("\n對齊後各檔案的咖啡因峰邊界（訊號降到峰高 1% 處）:\n")
for (k in seq_len(nrow(bounds)))
  cat(sprintf("  %-13s %.3f - %.3f min\n", rownames(bounds)[k],
              tvec[bounds[k, 1]], tvec[bounds[k, 2]]))

i0 <- max(1, min(bounds[, 1]) - 2)                  # 各留 1 秒緩衝
i1 <- min(N, max(bounds[, 2]) + 2)
cat(sprintf("=> 固定積分視窗: %.3f - %.3f min  (index %d-%d, %d 點)\n",
            tvec[i0], tvec[i1], i0, i1, i1 - i0 + 1))

# 基線已扣除，直接等距梯形積分 (uV*sec)
integrate_fixed <- function(v) {
  seg <- v[i0:i1]
  DT_SEC * (sum(seg) - (seg[1] + seg[length(seg)]) / 2)
}

area_aligned  <- vapply(all_files, function(fn) integrate_fixed(aligned[[fn]]), 0)
area_unaligned <- vapply(all_files, function(fn) integrate_fixed(prep[[fn]]$signal), 0)

# ---- 檢量線與定量 -----------------------------------------------------------
quantify <- function(area_vec, tag) {
  cal <- data.frame(file = names(standards),
                    conc = as.numeric(standards),
                    area = area_vec[names(standards)])
  fit <- lm(area ~ conc, data = cal)
  sm  <- summary(fit)
  b0 <- unname(coef(fit)[1]); b1 <- unname(coef(fit)[2])
  n <- nrow(cal); dof <- n - 2
  x <- cal$conc; xbar <- mean(x); Sxx <- sum((x - xbar)^2)
  s <- sm$sigma; tcrit <- qt(0.975, dof)

  cal$back_calc <- (cal$area - b0) / b1
  cal$recovery  <- ifelse(cal$conc > 0, cal$back_calc / cal$conc * 100, NA)

  sp <- data.frame(file = samples, area = area_vec[samples])
  sp$conc   <- (sp$area - b0) / b1
  sp$se     <- (s / b1) * sqrt(1 + 1/n + (sp$conc - xbar)^2 / Sxx)
  sp$ci_low <- sp$conc - tcrit * sp$se
  sp$ci_high<- sp$conc + tcrit * sp$se

  list(tag = tag, cal = cal, samp = sp, b0 = b0, b1 = b1,
       r2 = sm$r.squared, s = s, lod = 3.3*s/b1, loq = 10*s/b1,
       xbar = xbar, Sxx = Sxx, n = n)
}

q_before <- quantify(area_unaligned, "對齊前（固定視窗）")
q_after  <- quantify(area_aligned,  "對齊後（固定視窗）")

# ---- 視窗寬度敏感度：對齊到底買到了什麼 -------------------------------------
# 以參考峰頂為中心、不同半寬的固定視窗重跑定量。
# 穩健的方法應該在各種視窗寬度下都給出接近的濃度。
half_widths <- c(0.20, 0.25, 0.30, 0.40, 0.50)
sens <- lapply(half_widths, function(hw) {
  k <- round(hw * (N - 1) / T_END)
  a <- max(1, ref_apex_i - k); b <- min(N, ref_apex_i + k)
  integ <- function(v) { s <- v[a:b]; DT_SEC * (sum(s) - (s[1] + s[length(s)]) / 2) }
  list(hw = hw,
       before = quantify(vapply(all_files, function(fn) integ(prep[[fn]]$signal), 0), "")$samp$conc,
       after  = quantify(vapply(all_files, function(fn) integ(aligned[[fn]]), 0), "")$samp$conc)
})

# ---- 輸出 -------------------------------------------------------------------
write.csv(diag, file.path(out_dir, "alignment_diagnostics.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

sig_out <- data.frame(time_min = tvec)
for (fn in all_files) {
  id <- sub("\\.txt$", "", fn)
  sig_out[[paste0(id, "_bc")]]      <- round(prep[[fn]]$signal, 2)
  sig_out[[paste0(id, "_aligned")]] <- round(aligned[[fn]], 2)
}
write.csv(sig_out, file.path(out_dir, "aligned_signals.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cmp <- data.frame(
  file        = samples,
  conc_orig   = c(92.78, 83.97, 86.72, 87.33),     # 原始（儀器積分邊界，未對齊）
  conc_before = q_before$samp$conc,
  conc_after  = q_after$samp$conc
)
write.csv(cmp, file.path(out_dir, "alignment_comparison.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# 供互動教學網頁使用
write.csv(do.call(rbind, lapply(sens, function(r) data.frame(
            half_width = r$hw, sample = sub("\\.txt$", "", samples),
            before = r$before, after = r$after))),
          file.path(out_dir, "alignment_sensitivity.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(data.frame(
            key = c("reference", "ref_apex_min", "win_lo_min", "win_hi_min",
                    "slack", "cluster_gap", "de_np", "de_itermax", "alignDE_version"),
            value = c(sub("\\.txt$", "", REFERENCE), sprintf("%.4f", ref_apex),
                      sprintf("%.4f", tvec[i0]), sprintf("%.4f", tvec[i1]),
                      SLACK, CLUSTER_GAP, DE_NP, DE_ITERMAX,
                      as.character(packageVersion("alignDE")))),
          file.path(out_dir, "alignment_meta.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---- 報告 -------------------------------------------------------------------
con <- file(file.path(out_dir, "report_alignment.txt"), open = "wt", encoding = "UTF-8")
emit <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); writeLines(s, con) }
fmt <- function(v, d = 2) formatC(v, format = "f", big.mark = ",", digits = d)

emit(strrep("=", 78))
emit("HPLC 層析圖峰對齊 (alignDE) 與重新定量")
emit(strrep("=", 78))
emit(sprintf("套件: alignDE %s  (Zhang, Chen & Liang 2011, Talanta 83:1108)",
             as.character(packageVersion("alignDE"))))
emit(sprintf("參考層析圖: %s，咖啡因峰頂 %.3f min", REFERENCE, ref_apex))
emit(sprintf("參數: SNR.Th=%g, ridgeLength=%d, lambda=%g, threshold=%g, ",
             SNR_TH, RIDGE_LEN, BL_LAMBDA, BL_THRESH))
emit(sprintf("      clusterGap=%d, slack=%d, n=%d, NP=%d, itermax=%d, seed=20161017",
             CLUSTER_GAP, SLACK, GROUP_N, DE_NP, DE_ITERMAX))
emit("")
emit("[1] 對齊診斷")
emit(strrep("-", 78))
emit(sprintf("%-13s %11s %11s %12s %12s %11s", "檔案", "對齊前相關", "對齊後相關",
             "峰頂前(min)", "峰頂後(min)", "與參考差"))
for (i in seq_len(nrow(diag))) {
  na2 <- function(v, f) if (is.na(v)) "       —" else sprintf(f, v)
  emit(sprintf("%-13s %11.4f %11.4f %12s %12s %11s",
               diag$file[i], diag$cor_before[i], diag$cor_after[i],
               na2(diag$apex_before[i], "%.3f"), na2(diag$apex_after[i], "%.3f"),
               na2(diag$dev_after[i], "%+.3f")))
}
emit("  （std-0 為空白，該區間無咖啡因峰，峰頂欄位不適用）")
emit("")
emit(sprintf("  峰頂與參考的偏差 (絕對值平均):  對齊前 %.4f min -> 對齊後 %.4f min",
             mean(abs(diag$dev_before), na.rm = TRUE),
             mean(abs(diag$dev_after), na.rm = TRUE)))
emit(sprintf("  相關係數平均 (不含空白):        對齊前 %.4f  -> 對齊後 %.4f",
             mean(diag$cor_before[has_peak]), mean(diag$cor_after[has_peak])))
emit("")
emit(sprintf("[2] 固定積分視窗 %.3f - %.3f min（對齊後全部檔案共用）",
             tvec[i0], tvec[i1]))
emit(strrep("-", 78))
emit(sprintf("%-13s %8s %16s %16s %10s", "檔案", "濃度", "對齊前面積", "對齊後面積", "變化%"))
for (fn in all_files) {
  cc <- if (fn %in% names(standards)) sprintf("%.0f", standards[[fn]]) else "樣品"
  emit(sprintf("%-13s %8s %16s %16s %10.2f", fn, cc,
               fmt(area_unaligned[[fn]], 0), fmt(area_aligned[[fn]], 0),
               (area_aligned[[fn]] - area_unaligned[[fn]]) / area_unaligned[[fn]] * 100))
}
emit("")
for (q in list(q_before, q_after)) {
  emit(sprintf("[3] 檢量線 — %s", q$tag))
  emit(strrep("-", 78))
  emit(sprintf("  Area = %s x C %s %s", fmt(q$b1), ifelse(q$b0 >= 0, "+", "-"), fmt(abs(q$b0))))
  emit(sprintf("  R^2 = %.6f,  s = %s,  LOD = %.2f ppm,  LOQ = %.2f ppm",
               q$r2, fmt(q$s), q$lod, q$loq))
  emit(sprintf("  標準品回收率: %s",
               paste(sprintf("%.0fppm:%.1f%%", q$cal$conc[-1], q$cal$recovery[-1]),
                     collapse = "  ")))
  emit("")
}
emit("[4] 樣品濃度比較 (ppm)")
emit(strrep("-", 78))
emit(sprintf("%-13s %20s %20s %20s", "樣品", "原始(儀器積分邊界)",
             "固定視窗-對齊前", "固定視窗-對齊後"))
for (i in seq_along(samples)) {
  emit(sprintf("%-13s %20.2f %20.2f %20.2f", samples[i],
               cmp$conc_orig[i], cmp$conc_before[i], cmp$conc_after[i]))
}
emit("")
emit("[5] 視窗寬度敏感度：對齊買到了什麼")
emit(strrep("-", 78))
emit("  以參考峰頂 4.658 min 為中心、不同半寬的固定視窗重跑定量 (ppm)。")
emit("  穩健的方法在任何合理視窗下都該給出接近的濃度。")
emit("")
emit(sprintf("%-10s %-30s %-30s", "視窗半寬", "未對齊  s1 / s2 / s3 / s4",
             "對齊後  s1 / s2 / s3 / s4"))
for (r in sens) {
  emit(sprintf("±%.2f min  %-30s %-30s", r$hw,
               paste(sprintf("%.1f", r$before), collapse = " / "),
               paste(sprintf("%.1f", r$after),  collapse = " / ")))
}
spread <- function(g) {
  m <- do.call(rbind, lapply(sens, function(r) r[[g]]))
  mean(apply(m, 2, function(v) max(v) - min(v)))
}
emit("")
emit(sprintf("  各樣品在不同視窗下的濃度全距平均: 未對齊 %.2f ppm -> 對齊後 %.2f ppm",
             spread("before"), spread("after")))
emit("")
emit("[6] 備註")
emit(strrep("-", 78))
emit("  * alignDE 只平移峰的位置、峰間以線性內插伸縮，峰形與峰面積依設計保留，")
emit("    因此對齊本身不會改變面積；它帶來的好處是讓所有檔案能共用同一個固定")
emit("    積分視窗，免除逐檔設定積分邊界造成的不一致。")
emit("  * 此處面積以 Whittaker 懲罰最小平方基線扣除後積分，與儀器的線性基線積分")
emit("    不同，故絕對數值與原始報告不同；但標準品與樣品同樣處理，檢量線會吸收")
emit("    這個系統差異，最終濃度仍可直接比較。")
emit(strrep("=", 78))
close(con)

# ---- 繪圖 -------------------------------------------------------------------
png(file.path(out_dir, "alignment_r.png"), width = 1200, height = 1000, res = 120)
par(mfrow = c(2, 1), mar = c(4, 5, 3, 1))
zoom <- which(tvec >= 4.0 & tvec <= 5.3)
cols <- setNames(c("#94a3b8","#38bdf8","#34d399","#fbbf24","#fb7185","#a78bfa",
                   "#2563eb","#0891b2","#059669","#d97706"), all_files)
for (tag in c("before", "after")) {
  ys <- lapply(all_files, function(fn)
    (if (tag == "before") prep[[fn]]$signal else aligned[[fn]])[zoom] / 1000)
  plot(NA, xlim = range(tvec[zoom]), ylim = range(unlist(ys)),
       xlab = "Retention time (min)", ylab = "Signal (mV, baseline corrected)",
       main = if (tag == "before") "Before alignment (baseline corrected)"
              else "After alignDE alignment")
  abline(v = ref_apex, col = "grey60", lty = 2)
  for (i in seq_along(all_files))
    lines(tvec[zoom], ys[[i]], col = cols[all_files[i]], lwd = 1.6)
  legend("topright", legend = sub("\\.txt$", "", all_files),
         col = cols[all_files], lty = 1, lwd = 1.6, cex = .62, ncol = 2, bty = "n")
}
dev.off()

cat("\n輸出檔案位於: ", out_dir, "\n", sep = "")
