# -----------------------------------------------------------------------------
# HPLC 咖啡因檢量線與樣品濃度計算 (R 版)
#
# 輸入：Shimadzu LabSolutions ASCII 匯出檔 (std-*.txt, sample-*.txt)
# 輸出：outputs/ 下的 CSV 與 PNG
#
# 參考文獻:
#   檢量線建立與線性評估
#     Danzer K, Currie LA (1998) Pure Appl Chem 70(4):993-1014.
#     doi:10.1351/pac199870040993
#   LOD = 3.3 s/b, LOQ = 10 s/b
#     ICH Q2(R2) Validation of Analytical Procedures (2023);
#     Currie LA (1995) Pure Appl Chem 67(10):1699-1723.
#     doi:10.1351/pac199567101699
#   反推濃度的標準誤與信賴區間
#     Miller JN, Miller JC (2010) Statistics and Chemometrics for Analytical
#     Chemistry, 6th ed., Ch. 5. Pearson.
#   咖啡因 HPLC 定量
#     DiNunzio JE (1985) J Chem Educ 62(5):446. doi:10.1021/ed062p446
#     Naik JP, Nagalakshmi S (1997) J Agric Food Chem 45(10):3973-3975.
#     doi:10.1021/jf970147i
#
# 用法:  Rscript caffeine_calibration.R
# 僅使用 base R，不需額外套件。
# -----------------------------------------------------------------------------

options(stringsAsFactors = FALSE)

base_dir <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
out_dir <- file.path(base_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE)

# ---- 設定 -------------------------------------------------------------------
standards <- c("std-0.txt" = 0, "std-25.txt" = 25, "std-50.txt" = 50,
               "std-100.txt" = 100, "std-150.txt" = 150, "std-200.txt" = 200)
samples <- c("sample-1.txt", "sample-2.txt", "sample-3.txt", "sample-4.txt")

RT_LO <- 4.20   # 咖啡因判峰視窗 (min)
RT_HI <- 4.90
AREA_NOISE_FRAC <- 0.005

# ---- 解析 Shimadzu 匯出檔（共用模組）---------------------------------------
source(file.path(base_dir, "hplc_io.R"))

load_run <- function(fn) load_run_file(file.path(base_dir, fn))

# ---- 判峰與重新積分 ---------------------------------------------------------
pick_caffeine <- function(pk) {
  # 不採用儀器自動 ID（多筆誤判），改取滯留時間視窗內面積最大的峰
  if (is.null(pk) || !nrow(pk)) return(NULL)
  sel <- pk[!is.na(pk$R.Time) & pk$R.Time >= RT_LO & pk$R.Time <= RT_HI, , drop = FALSE]
  if (!nrow(sel)) return(NULL)
  sel[which.max(sel$Area), , drop = FALSE]
}

reintegrate <- trapz_baseline   # 見 hplc_io.R

# ---- 讀入全部檔案 -----------------------------------------------------------
all_files <- c(names(standards), samples)
runs <- lapply(all_files, load_run)
names(runs) <- all_files

peak_tab <- do.call(rbind, lapply(all_files, function(fn) {
  p <- pick_caffeine(runs[[fn]]$peaks)
  if (is.null(p)) {
    data.frame(file = fn, rt = NA_real_, area = 0, height = 0,
               t_start = NA_real_, t_end = NA_real_, area_recalc = 0)
  } else {
    data.frame(file = fn, rt = p$R.Time, area = p$Area, height = p$Height,
               t_start = p$I.Time, t_end = p$F.Time,
               area_recalc = reintegrate(runs[[fn]]$chrom, p$I.Time, p$F.Time))
  }
}))
rownames(peak_tab) <- peak_tab$file
peak_tab$detected <- peak_tab$area > AREA_NOISE_FRAC * max(peak_tab$area)

# ---- 檢量線 -----------------------------------------------------------------
cal <- peak_tab[names(standards), ]
cal$conc <- as.numeric(standards[cal$file])

fit  <- lm(area ~ conc, data = cal)
cf   <- coef(fit)
b0   <- unname(cf[1]); b1 <- unname(cf[2])
sm   <- summary(fit)
se0  <- sm$coefficients[1, 2]; se1 <- sm$coefficients[2, 2]
r2   <- sm$r.squared
s_res <- sm$sigma
n     <- nrow(cal); dof <- n - 2
x     <- cal$conc; xbar <- mean(x); Sxx <- sum((x - xbar)^2)
tcrit <- qt(0.975, dof)

lod <- 3.3 * s_res / b1
loq <- 10  * s_res / b1

cal$area_fit     <- fitted(fit)
cal$residual     <- resid(fit)
cal$residual_pct <- ifelse(cal$area != 0, cal$residual / cal$area * 100, NA)
cal$back_calc    <- (cal$area - b0) / b1
cal$recovery_pct <- ifelse(cal$conc > 0, cal$back_calc / cal$conc * 100, NA)

# 二次式對照（檢查線性偏離）
fit2 <- lm(area ~ conc + I(conc^2), data = cal)
r2_quad <- summary(fit2)$r.squared

# ---- 樣品濃度（反推 + 95% 信賴區間，Miller & Miller）------------------------
samp <- peak_tab[samples, c("file", "rt", "area", "height", "area_recalc")]
samp$conc_ppm <- (samp$area - b0) / b1
samp$se       <- (s_res / b1) * sqrt(1 + 1 / n + (samp$conc_ppm - xbar)^2 / Sxx)
samp$ci_low   <- samp$conc_ppm - tcrit * samp$se
samp$ci_high  <- samp$conc_ppm + tcrit * samp$se
samp$in_range <- ifelse(samp$conc_ppm >= min(x) & samp$conc_ppm <= max(x), "是", "否(外插)")
samp$above_loq <- ifelse(samp$conc_ppm >= loq, "是", "否")

# ---- 輸出 CSV ---------------------------------------------------------------
write.csv(cal[, c("file", "conc", "rt", "area", "area_recalc", "area_fit",
                  "residual", "residual_pct", "back_calc", "recovery_pct")],
          file.path(out_dir, "calibration_r.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")
write.csv(samp, file.path(out_dir, "samples_r.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

# ---- 文字報告 ---------------------------------------------------------------
con <- file(file.path(out_dir, "report_r.txt"), open = "wt", encoding = "UTF-8")
emit <- function(...) {
  s <- paste0(...)
  cat(s, "\n", sep = "")
  writeLines(s, con)
}
fmt <- function(v, d = 2) formatC(v, format = "f", big.mark = ",", digits = d)

# 中日韓字元在等寬終端佔 2 欄，sprintf 的寬度以字元數計，需自行補齊
dwidth <- function(s) {
  ch <- strsplit(as.character(s), "")[[1]]
  if (!length(ch)) return(0L)
  sum(ifelse(grepl("[ᄀ-ᅟ⺀-꓏가-힣豈-﫿︰-﹏＀-｠￠-￦]", ch), 2L, 1L))
}
padw <- function(s, w, left = FALSE) {
  s <- as.character(s)
  p <- strrep(" ", max(w - dwidth(s), 0))
  if (left) paste0(s, p) else paste0(p, s)
}

emit(strrep("=", 74))
emit("HPLC 咖啡因定量分析 (R)")
emit(strrep("=", 74))
emit("偵測器 Detector A-Ch1, 280 nm；注射量 10 uL；稀釋倍數 1")
emit(sprintf("咖啡因判峰視窗: %.2f - %.2f min（取視窗內最大峰）", RT_LO, RT_HI))
emit("")
emit("[1] 檢量線資料點")
emit(strrep("-", 74))
emit(paste0(padw("檔案", 13, TRUE), padw("濃度(ppm)", 11), padw("RT(min)", 9),
            padw("峰面積", 14), padw("重積分面積", 14), padw("殘差%", 9),
            padw("回算濃度", 11)))
for (i in seq_len(nrow(cal))) {
  emit(paste0(
    padw(cal$file[i], 13, TRUE), padw(sprintf("%.1f", cal$conc[i]), 11),
    padw(ifelse(is.na(cal$rt[i]), "-", sprintf("%.3f", cal$rt[i])), 9),
    padw(fmt(cal$area[i], 0), 14), padw(fmt(cal$area_recalc[i], 0), 14),
    padw(ifelse(is.na(cal$residual_pct[i]), "-", sprintf("%.2f", cal$residual_pct[i])), 9),
    padw(sprintf("%.2f", cal$back_calc[i]), 11)))
}
emit("")
emit("[2] 線性迴歸結果  Area = a + b x Conc")
emit(strrep("-", 74))
emit(sprintf("  斜率 b        = %s  ± %s  (area / ppm)", fmt(b1), fmt(se1)))
emit(sprintf("  截距 a        = %s  ± %s", fmt(b0), fmt(se0)))
emit(sprintf("  相關係數 r    = %.6f", sqrt(r2) * sign(b1)))
emit(sprintf("  判定係數 R^2  = %.6f", r2))
emit(sprintf("  迴歸標準誤    = %s   (自由度 %d)", fmt(s_res), dof))
emit(sprintf("  檢量線方程式  : Area = %s x C %s %s",
             fmt(b1), ifelse(b0 >= 0, "+", "-"), fmt(abs(b0))))
emit(sprintf("  反算式        : C(ppm) = (Area %s %s) / %s",
             ifelse(b0 >= 0, "-", "+"), fmt(abs(b0)), fmt(b1)))
emit("")
emit(sprintf("  截距是否顯著異於 0: t = %.3f, t(0.975,%d) = %.3f -> %s",
             b0 / se0, dof, tcrit,
             ifelse(abs(b0 / se0) > tcrit, "顯著", "不顯著（可視為過原點）")))
emit(sprintf("  LOD (3.3 s/b) = %.2f ppm", lod))
emit(sprintf("  LOQ (10  s/b) = %.2f ppm", loq))
emit(sprintf("  二次式對照 R^2 = %.6f (與線性 R^2 差 %.6f)", r2_quad, r2_quad - r2))
emit("")
emit("[3] 樣品濃度")
emit(strrep("-", 74))
emit(paste0(padw("檔案", 13, TRUE), padw("RT(min)", 9), padw("峰面積", 14),
            padw("濃度(ppm)", 12), padw("95%CI下限", 13), padw("95%CI上限", 13),
            padw("在範圍內", 11)))
for (i in seq_len(nrow(samp))) {
  emit(paste0(
    padw(samp$file[i], 13, TRUE), padw(sprintf("%.3f", samp$rt[i]), 9),
    padw(fmt(samp$area[i], 0), 14), padw(sprintf("%.2f", samp$conc_ppm[i]), 12),
    padw(sprintf("%.2f", samp$ci_low[i]), 13), padw(sprintf("%.2f", samp$ci_high[i]), 13),
    padw(samp$in_range[i], 11)))
}
emit("")
emit("[4] ANOVA — 檢量線線性顯著性")
emit(strrep("-", 74))
av <- anova(fit)
emit(sprintf("  F = %.1f, p = %.3e", av$`F value`[1], av$`Pr(>F)`[1]))
emit("")
emit("[5] 備註")
emit(strrep("-", 74))
emit("  * 儀器自動指認的化合物峰不可靠（std-25 指到 4.177 min、sample-3 指到")
emit("    3.242 min 的微小雜峰），故本程式改以滯留時間視窗內最大峰判定咖啡因。")
emit("  * 「重積分面積」為以峰起訖點連線作基線、對原始層析圖梯形積分之結果。")
emit("  * Dilution Factor = 1；若實際有稀釋，濃度須再乘上稀釋倍數。")
emit(strrep("=", 74))
close(con)

# ---- 繪圖 -------------------------------------------------------------------
png(file.path(out_dir, "calibration_r.png"), width = 1200, height = 520, res = 120)
par(mfrow = c(1, 2), mar = c(4.5, 5, 3, 1))
plot(cal$conc, cal$area, pch = 19, col = "#2c3e50", cex = 1.2,
     xlab = "Caffeine concentration (ppm)", ylab = "Peak area",
     main = sprintf("Calibration curve (R2 = %.5f)", r2),
     xlim = c(-5, 215), ylim = c(0, max(cal$area) * 1.05))
abline(fit, col = "#c0392b", lwd = 2)
points(samp$conc_ppm, samp$area, pch = 17, col = "#2980b9", cex = 1.3)
# 樣品濃度相近，標籤拉到右下方錯開避免重疊
lab_y <- 2.1e6 + (seq_len(nrow(samp)) - 1) * 0.55e6
segments(samp$conc_ppm, samp$area, 128, lab_y, col = "#2980b955", lwd = .8)
text(130, lab_y, sprintf("%s: %.1f ppm", sub("\\.txt$", "", samp$file), samp$conc_ppm),
     pos = 4, cex = .65, col = "#2980b9")
legend("topleft", bty = "n", cex = .8,
       legend = c(sprintf("y = %.0f x %s %.0f", b1, ifelse(b0 >= 0, "+", "-"), abs(b0)),
                  "standards", "samples"),
       col = c("#c0392b", "#2c3e50", "#2980b9"),
       lty = c(1, NA, NA), pch = c(NA, 19, 17))
grid()

plot(cal$conc, cal$residual, type = "h", lwd = 3, col = "#2980b9",
     xlab = "Caffeine concentration (ppm)", ylab = "Residual (area)",
     main = "Residual plot")
points(cal$conc, cal$residual, pch = 19, col = "#2980b9")
abline(h = 0, col = "grey40")
grid()
dev.off()

png(file.path(out_dir, "chromatograms_r.png"), width = 1200, height = 900, res = 120)
par(mfrow = c(2, 1), mar = c(4, 5, 3, 1))
for (grp in list(list(f = names(standards), t = "Standards"),
                 list(f = samples, t = "Samples"))) {
  ch <- lapply(grp$f, function(fn) runs[[fn]]$chrom)
  ok <- !vapply(ch, is.null, TRUE)
  ch <- ch[ok]; labs <- sub("\\.txt$", "", grp$f[ok])
  ylim <- range(unlist(lapply(ch, function(d) d$intensity / 1000)))
  plot(NA, xlim = c(0, 6), ylim = ylim, xlab = "Retention time (min)",
       ylab = "Intensity (mV)", main = paste(grp$t, "- Detector A-Ch1, 280 nm"))
  rect(RT_LO, ylim[1], RT_HI, ylim[2], col = "#f1c40f33", border = NA)
  cols <- rainbow(length(ch), v = .85)
  for (i in seq_along(ch)) lines(ch[[i]]$time_min, ch[[i]]$intensity / 1000, col = cols[i])
  legend("topright", legend = labs, col = cols, lty = 1, cex = .7, ncol = 2, bty = "n")
  grid()
}
dev.off()

cat("\n輸出檔案位於: ", out_dir, "\n", sep = "")
