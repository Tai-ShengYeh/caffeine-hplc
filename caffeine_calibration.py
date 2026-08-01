# -*- coding: utf-8 -*-
"""
HPLC 咖啡因檢量線與樣品濃度計算 (Python 版)

輸入：Shimadzu LabSolutions ASCII 匯出檔 (std-*.txt, sample-*.txt)
輸出：outputs/ 下的 CSV 與圖檔

參考文獻:
    檢量線建立與線性評估
        Danzer K, Currie LA (1998) Pure Appl Chem 70(4):993-1014.
        doi:10.1351/pac199870040993
    LOD = 3.3 s/b, LOQ = 10 s/b
        ICH Q2(R2) Validation of Analytical Procedures (2023);
        Currie LA (1995) Pure Appl Chem 67(10):1699-1723.
        doi:10.1351/pac199567101699
    反推濃度的標準誤與信賴區間
        Miller JN, Miller JC (2010) Statistics and Chemometrics for
        Analytical Chemistry, 6th ed., Ch. 5. Pearson.
    咖啡因 HPLC 定量
        DiNunzio JE (1985) J Chem Educ 62(5):446. doi:10.1021/ed062p446
        Naik JP, Nagalakshmi S (1997) J Agric Food Chem 45(10):3973-3975.
        doi:10.1021/jf970147i

用法:
    python caffeine_calibration.py
"""

import os
import re
import sys
import math
import unicodedata

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- 設定
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "outputs")

# 標準品檔名 -> 標稱濃度 (ppm, 即 ug/mL)
STANDARDS = {
    "std-0.txt": 0.0,
    "std-25.txt": 25.0,
    "std-50.txt": 50.0,
    "std-100.txt": 100.0,
    "std-150.txt": 150.0,
    "std-200.txt": 200.0,
}
SAMPLES = ["sample-1.txt", "sample-2.txt", "sample-3.txt", "sample-4.txt"]

# 咖啡因滯留時間視窗 (min)。標準品約 4.65、樣品因基質稍早約 4.43~4.60
RT_WINDOW = (4.20, 4.90)
# 判定為「未檢出」的最低面積門檻 (相對於最高標準品面積)
AREA_NOISE_FRAC = 0.005


# ---------------------------------------------------------------- 排版工具
def dwidth(s):
    """中日韓字元在等寬終端佔 2 欄，計算實際顯示寬度。"""
    return sum(2 if unicodedata.east_asian_width(c) in ("W", "F") else 1 for c in str(s))


def pad(s, w, left=False):
    s = str(s)
    sp = " " * max(w - dwidth(s), 0)
    return s + sp if left else sp + s


# ---------------------------------------------------------------- 解析
def read_text(path):
    """Shimadzu 匯出檔為 Big5 (繁中 Windows)，失敗時退回 latin-1。"""
    for enc in ("big5", "cp950", "utf-8"):
        try:
            with open(path, "r", encoding=enc) as f:
                return f.read()
        except UnicodeDecodeError:
            continue
    with open(path, "r", encoding="latin-1") as f:
        return f.read()


def split_sections(text):
    """把 [Section Name] ... 切成 dict: {section: [lines]}"""
    sections, current = {}, None
    for line in text.splitlines():
        m = re.match(r"^\[(.+)\]\s*$", line)
        if m:
            current = m.group(1)
            sections[current] = []
        elif current is not None:
            sections[current].append(line)
    return sections


def parse_kv(lines):
    d = {}
    for line in lines:
        if "\t" in line:
            k, _, v = line.partition("\t")
            d[k.strip()] = v.strip()
    return d


def parse_peak_table(lines):
    """回傳 peak table DataFrame（可能為空）。"""
    if not lines:
        return pd.DataFrame()
    n_peaks = 0
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("# of Peaks"):
            n_peaks = int(line.split("\t")[1])
        elif line.startswith("Peak#"):
            header_idx = i
            break
    if header_idx is None or n_peaks == 0:
        return pd.DataFrame()
    cols = lines[header_idx].split("\t")
    rows = [l.split("\t") for l in lines[header_idx + 1: header_idx + 1 + n_peaks]]
    df = pd.DataFrame(rows, columns=cols)
    for c in ("R.Time", "I.Time", "F.Time", "Area", "Height", "Conc."):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def parse_chromatogram(lines):
    """回傳 (DataFrame[time, intensity_uV], meta)。"""
    meta, data_idx = {}, None
    for i, line in enumerate(lines):
        if line.startswith("R.Time"):
            data_idx = i + 1
            break
        if "\t" in line:
            k, _, v = line.partition("\t")
            meta[k.strip()] = v.strip()
    if data_idx is None:
        return pd.DataFrame(), meta
    rows = []
    for line in lines[data_idx:]:
        if not line.strip():
            break
        parts = line.split("\t")
        if len(parts) < 2:
            break
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except ValueError:
            break
    return pd.DataFrame(rows, columns=["time_min", "intensity"]), meta


def load_run(path):
    text = read_text(path)
    sec = split_sections(text)
    chrom, chrom_meta = parse_chromatogram(sec.get("LC Chromatogram(Detector A-Ch1)", []))
    return {
        "file": os.path.basename(path),
        "info": parse_kv(sec.get("Sample Information", [])),
        "peaks": parse_peak_table(sec.get("Peak Table(Detector A-Ch1)", [])),
        "chrom": chrom,
        "chrom_meta": chrom_meta,
    }


# ---------------------------------------------------------------- 判峰 / 積分
def pick_caffeine_peak(peaks, rt_lo, rt_hi):
    """在滯留時間視窗內取面積最大的峰（不信任儀器自動 ID）。"""
    if peaks.empty:
        return None
    sel = peaks[(peaks["R.Time"] >= rt_lo) & (peaks["R.Time"] <= rt_hi)]
    if sel.empty:
        return None
    return sel.loc[sel["Area"].idxmax()]


def reintegrate(chrom, t_start, t_end):
    """以起訖點連線為基線，梯形法重新積分。回傳面積 (uV*sec)。"""
    if chrom.empty:
        return np.nan
    m = (chrom["time_min"] >= t_start) & (chrom["time_min"] <= t_end)
    seg = chrom.loc[m]
    if len(seg) < 3:
        return np.nan
    t = seg["time_min"].to_numpy() * 60.0      # min -> sec
    y = seg["intensity"].to_numpy()            # uV
    baseline = np.linspace(y[0], y[-1], len(y))
    return float(np.trapezoid(y - baseline, t)) if hasattr(np, "trapezoid") \
        else float(np.trapz(y - baseline, t))


# ---------------------------------------------------------------- 主流程
def main():
    # Windows 主控台預設 cp950，強制以 UTF-8 輸出中文
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    os.makedirs(OUT_DIR, exist_ok=True)

    runs = {}
    for fn in list(STANDARDS) + SAMPLES:
        p = os.path.join(BASE_DIR, fn)
        if not os.path.exists(p):
            sys.exit(f"找不到檔案: {p}")
        runs[fn] = load_run(p)

    # ---- 擷取每個檔案的咖啡因峰
    rows = []
    for fn, run in runs.items():
        pk = pick_caffeine_peak(run["peaks"], *RT_WINDOW)
        if pk is None:
            rows.append(dict(file=fn, rt=np.nan, area=0.0, height=0.0,
                             t_start=np.nan, t_end=np.nan, area_recalc=0.0))
        else:
            rows.append(dict(
                file=fn, rt=pk["R.Time"], area=pk["Area"], height=pk["Height"],
                t_start=pk["I.Time"], t_end=pk["F.Time"],
                area_recalc=reintegrate(run["chrom"], pk["I.Time"], pk["F.Time"]),
            ))
    peak_df = pd.DataFrame(rows).set_index("file")

    max_area = peak_df["area"].max()
    peak_df["detected"] = peak_df["area"] > AREA_NOISE_FRAC * max_area

    # ---- 檢量線 (最小平方法, 面積 vs 濃度)
    cal = peak_df.loc[list(STANDARDS)].copy()
    cal["conc"] = [STANDARDS[f] for f in cal.index]
    x = cal["conc"].to_numpy(float)
    y = cal["area"].to_numpy(float)
    n = len(x)

    fit = stats.linregress(x, y)
    slope, intercept, r = fit.slope, fit.intercept, fit.rvalue
    y_hat = intercept + slope * x
    resid = y - y_hat
    dof = n - 2
    s_res = math.sqrt(np.sum(resid ** 2) / dof)        # 迴歸標準誤 s_y/x
    x_bar = x.mean()
    Sxx = float(np.sum((x - x_bar) ** 2))
    se_slope, se_intercept = fit.stderr, fit.intercept_stderr
    t_crit = stats.t.ppf(0.975, dof)

    lod = 3.3 * s_res / slope      # 偵測極限 (ICH, 以 s_y/x 估計)
    loq = 10.0 * s_res / slope     # 定量極限

    cal["area_fit"] = y_hat
    cal["residual"] = resid
    with np.errstate(divide="ignore", invalid="ignore"):
        cal["residual_pct"] = np.where(y != 0, resid / np.where(y == 0, 1, y) * 100, np.nan)
    cal["back_calc_conc"] = (y - intercept) / slope
    cal["recovery_pct"] = np.where(x > 0, cal["back_calc_conc"] / x * 100, np.nan)

    # ---- 樣品濃度 (反推 + 95% 信賴區間)
    srows = []
    for fn in SAMPLES:
        pk = peak_df.loc[fn]
        y0 = float(pk["area"])
        x0 = (y0 - intercept) / slope
        se_x0 = (s_res / slope) * math.sqrt(1.0 / 1 + 1.0 / n + (x0 - x_bar) ** 2 / Sxx)
        srows.append(dict(
            file=fn, rt_min=pk["rt"], area=y0, height=pk["height"],
            area_recalc=pk["area_recalc"],
            conc_ppm=x0, se=se_x0,
            ci_low=x0 - t_crit * se_x0, ci_high=x0 + t_crit * se_x0,
            in_range="是" if (x.min() <= x0 <= x.max()) else "否(外插)",
            above_loq="是" if x0 >= loq else "否",
        ))
    samp = pd.DataFrame(srows).set_index("file")

    # ---- 二次式對照（檢查線性偏離）
    quad = np.polyfit(x, y, 2)
    y_q = np.polyval(quad, x)
    r2_quad = 1 - np.sum((y - y_q) ** 2) / np.sum((y - y.mean()) ** 2)

    # ---- 輸出
    cal_out = cal[["conc", "rt", "area", "area_recalc", "area_fit",
                   "residual", "residual_pct", "back_calc_conc", "recovery_pct"]]
    cal_out.to_csv(os.path.join(OUT_DIR, "calibration_py.csv"),
                   encoding="utf-8-sig", float_format="%.4f")
    samp.to_csv(os.path.join(OUT_DIR, "samples_py.csv"),
                encoding="utf-8-sig", float_format="%.4f")

    with open(os.path.join(OUT_DIR, "report_py.txt"), "w", encoding="utf-8") as fh:
        def emit(s=""):
            print(s)
            fh.write(s + "\n")

        emit("=" * 74)
        emit("HPLC 咖啡因定量分析 (Python)")
        emit("=" * 74)
        emit(f"偵測器 Detector A-Ch1, 280 nm；注射量 10 uL；稀釋倍數 1")
        emit(f"咖啡因判峰視窗: {RT_WINDOW[0]:.2f} - {RT_WINDOW[1]:.2f} min（取視窗內最大峰）")
        emit()
        emit("[1] 檢量線資料點")
        emit("-" * 74)
        emit(pad("檔案", 13, left=True) + pad("濃度(ppm)", 11) + pad("RT(min)", 9)
             + pad("峰面積", 14) + pad("重積分面積", 14) + pad("殘差%", 9)
             + pad("回算濃度", 11))
        for f_, row in cal.iterrows():
            rt = "-" if np.isnan(row["rt"]) else f"{row['rt']:.3f}"
            rp = "-" if pd.isna(row["residual_pct"]) else f"{row['residual_pct']:.2f}"
            emit(pad(f_, 13, left=True) + pad(f"{row['conc']:.1f}", 11) + pad(rt, 9)
                 + pad(f"{row['area']:,.0f}", 14) + pad(f"{row['area_recalc']:,.0f}", 14)
                 + pad(rp, 9) + pad(f"{row['back_calc_conc']:.2f}", 11))
        emit()
        emit("[2] 線性迴歸結果  Area = a + b x Conc")
        emit("-" * 74)
        emit(f"  斜率 b        = {slope:,.2f}  ± {se_slope:,.2f}  (area / ppm)")
        emit(f"  截距 a        = {intercept:,.2f}  ± {se_intercept:,.2f}")
        emit(f"  相關係數 r    = {r:.6f}")
        emit(f"  判定係數 R^2  = {r ** 2:.6f}")
        emit(f"  迴歸標準誤    = {s_res:,.2f}   (自由度 {dof})")
        emit(f"  檢量線方程式  : Area = {slope:,.2f} x C {'+' if intercept >= 0 else '-'} "
             f"{abs(intercept):,.2f}")
        emit(f"  反算式        : C(ppm) = (Area {'-' if intercept >= 0 else '+'} "
             f"{abs(intercept):,.2f}) / {slope:,.2f}")
        emit()
        emit(f"  截距是否顯著異於 0: t = {intercept / se_intercept:.3f}, "
             f"t(0.975,{dof}) = {t_crit:.3f} -> "
             f"{'顯著' if abs(intercept / se_intercept) > t_crit else '不顯著（可視為過原點）'}")
        emit(f"  LOD (3.3 s/b) = {lod:.2f} ppm")
        emit(f"  LOQ (10  s/b) = {loq:.2f} ppm")
        emit(f"  二次式對照 R^2 = {r2_quad:.6f} "
             f"(二次項係數 {quad[0]:,.2f}；與線性 R^2 差 {r2_quad - r ** 2:.6f})")
        emit()
        emit("[3] 樣品濃度")
        emit("-" * 74)
        emit(pad("檔案", 13, left=True) + pad("RT(min)", 9) + pad("峰面積", 14)
             + pad("濃度(ppm)", 12) + pad("95%CI下限", 13) + pad("95%CI上限", 13)
             + pad("在範圍內", 11))
        for f_, row in samp.iterrows():
            emit(pad(f_, 13, left=True) + pad(f"{row['rt_min']:.3f}", 9)
                 + pad(f"{row['area']:,.0f}", 14) + pad(f"{row['conc_ppm']:.2f}", 12)
                 + pad(f"{row['ci_low']:.2f}", 13) + pad(f"{row['ci_high']:.2f}", 13)
                 + pad(row["in_range"], 11))
        emit()
        emit("[4] 備註")
        emit("-" * 74)
        emit("  * 儀器自動指認的化合物峰不可靠（std-25 指到 4.177 min、sample-3 指到")
        emit("    3.242 min 的微小雜峰），故本程式改以滯留時間視窗內最大峰判定咖啡因。")
        emit("  * 「重積分面積」為以峰起訖點連線作基線、對原始層析圖梯形積分之結果，")
        emit("    用於交叉驗證儀器報告之面積。")
        emit("  * 樣品濃度單位同標準品配製單位（此處以 ppm = ug/mL 表示）；")
        emit("    Dilution Factor = 1，若實際有稀釋須再乘上稀釋倍數。")
        emit("=" * 74)

    # ---- 繪圖
    plot_calibration(x, y, slope, intercept, r ** 2, samp, cal)
    plot_chromatograms(runs)
    print(f"\n輸出檔案位於: {OUT_DIR}")

    return cal, samp


def plot_calibration(x, y, slope, intercept, r2, samp, cal):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    xs = np.linspace(-5, 215, 200)
    ax1.plot(xs, intercept + slope * xs, "-", color="#c0392b", lw=1.6,
             label=f"y = {slope:,.0f}x {'+' if intercept >= 0 else '-'} {abs(intercept):,.0f}")
    ax1.plot(x, y, "o", ms=8, color="#2c3e50", label="Standards")
    for xi, yi in zip(x, y):
        ax1.annotate(f"{xi:.0f}", (xi, yi), textcoords="offset points",
                     xytext=(8, -10), fontsize=8, color="#2c3e50")
    # 樣品濃度相近，標籤錯開避免重疊
    for i, (f_, row) in enumerate(samp.iterrows()):
        ax1.plot(row["conc_ppm"], row["area"], "^", ms=9, color="#2980b9",
                 label="Samples" if i == 0 else None)
        ax1.annotate(f"{f_.replace('.txt', '')}: {row['conc_ppm']:.1f} ppm",
                     (row["conc_ppm"], row["area"]),
                     xytext=(130, 2.1e6 + i * 0.55e6), textcoords="data",
                     fontsize=8, color="#2980b9",
                     arrowprops=dict(arrowstyle="-", color="#2980b9",
                                     lw=.6, alpha=.5))
    ax1.set_xlabel("Caffeine concentration (ppm)")
    ax1.set_ylabel("Peak area")
    ax1.set_title(f"Calibration curve (R$^2$ = {r2:.5f})")
    ax1.legend(loc="upper left", fontsize=9)
    ax1.grid(alpha=.3)

    ax2.axhline(0, color="#888", lw=1)
    ax2.stem(x, cal["residual"].to_numpy(), linefmt="C0-", markerfmt="C0o", basefmt=" ")
    ax2.set_xlabel("Caffeine concentration (ppm)")
    ax2.set_ylabel("Residual (area)")
    ax2.set_title("Residual plot")
    ax2.grid(alpha=.3)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "calibration_py.png"), dpi=150)
    plt.close(fig)


def plot_chromatograms(runs):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    for fn, run in runs.items():
        if run["chrom"].empty:
            continue
        ax = axes[0] if fn.startswith("std") else axes[1]
        ax.plot(run["chrom"]["time_min"], run["chrom"]["intensity"] / 1000.0,
                lw=1, label=fn.replace(".txt", ""))
    for ax, title in zip(axes, ["Standards", "Samples"]):
        ax.axvspan(RT_WINDOW[0], RT_WINDOW[1], color="#f1c40f", alpha=.15)
        ax.set_ylabel("Intensity (mV)")
        ax.set_title(f"{title} — Detector A-Ch1, 280 nm")
        ax.legend(fontsize=8, ncol=3)
        ax.grid(alpha=.3)
    axes[1].set_xlabel("Retention time (min)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "chromatograms_py.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
