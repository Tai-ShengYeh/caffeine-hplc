# -*- coding: utf-8 -*-
"""
把原始 Shimadzu 匯出檔轉成互動教學網頁 hplc_teaching.html

會把層析圖原始資料、峰表、檢量線結果全部內嵌進單一 HTML（不依賴外部檔案/CDN）。

用法:
    python build_teaching_html.py
"""

import os
import sys
import json

import numpy as np
import pandas as pd
from scipy import stats

from caffeine_calibration import (
    BASE_DIR, STANDARDS, SAMPLES, RT_WINDOW,
    load_run, pick_caffeine_peak, reintegrate,
)

TEMPLATE = os.path.join(BASE_DIR, "hplc_teaching_template.html")
OUTPUT = os.path.join(BASE_DIR, "hplc_teaching.html")
BASIC_TEMPLATE = os.path.join(BASE_DIR, "hplc_basic_template.html")
BASIC_OUTPUT = os.path.join(BASE_DIR, "hplc_basic.html")
OUT_DIR = os.path.join(BASE_DIR, "outputs")
MARKER = "/*__HPLC_DATA__*/{}"

# 各檔案在圖上的顏色
COLORS = {
    "std-0": "#94a3b8", "std-25": "#38bdf8", "std-50": "#34d399",
    "std-100": "#fbbf24", "std-150": "#fb7185", "std-200": "#a78bfa",
    "sample-1": "#2563eb", "sample-2": "#0891b2",
    "sample-3": "#059669", "sample-4": "#d97706",
}


def load_alignment():
    """讀入 align_chromatograms.R 的輸出（若尚未執行則回傳 None）。"""
    need = ["aligned_signals.csv", "alignment_diagnostics.csv",
            "alignment_sensitivity.csv", "alignment_meta.csv"]
    paths = {n: os.path.join(OUT_DIR, n) for n in need}
    missing = [n for n, p in paths.items() if not os.path.exists(p)]
    if missing:
        print("！ 尚未找到 " + ", ".join(missing) +
              "；請先執行 Rscript align_chromatograms.R。教學頁將略過對齊章節。")
        return None

    sig = pd.read_csv(paths["aligned_signals.csv"])
    diag = pd.read_csv(paths["alignment_diagnostics.csv"])
    sens = pd.read_csv(paths["alignment_sensitivity.csv"])
    meta = pd.read_csv(paths["alignment_meta.csv"]).set_index("key")["value"]

    traces = {}
    for col in sig.columns:
        if col == "time_min":
            continue
        rid, kind = col.rsplit("_", 1)          # 例: "std-100_bc" / "std-100_aligned"
        key = {"bc": "bc", "aligned": "al"}[kind]
        traces.setdefault(rid, {})[key] = [int(round(v)) for v in sig[col]]

    diag["file"] = diag["file"].str.replace(".txt", "", regex=False)
    # NaN 在 JSON 中是非法字面值，轉成 null（空白樣品沒有峰頂）
    diag = diag.astype(object).where(pd.notna(diag), None)
    hw = sorted(sens["half_width"].unique())
    return dict(
        traces=traces,
        diag=diag.where(pd.notna(diag), None).to_dict("records"),
        sens=[dict(hw=float(h),
                   before=[float(v) for v in sens[sens.half_width == h]["before"]],
                   after=[float(v) for v in sens[sens.half_width == h]["after"]],
                   samples=list(sens[sens.half_width == h]["sample"]))
              for h in hw],
        refApex=float(meta["ref_apex_min"]),
        winLo=float(meta["win_lo_min"]),
        winHi=float(meta["win_hi_min"]),
        reference=str(meta["reference"]),
        params=dict(slack=int(meta["slack"]), clusterGap=int(meta["cluster_gap"]),
                    np=int(meta["de_np"]), itermax=int(meta["de_itermax"]),
                    version=str(meta["alignDE_version"])),
    )


def load_airpls():
    """讀入 baseline_airpls.R 的輸出（若尚未執行則回傳 None）。"""
    need = ["airpls_baselines.csv", "airpls_demo.csv", "airpls_drift.csv",
            "airpls_scan.csv", "airpls_meta.csv"]
    paths = {n: os.path.join(OUT_DIR, n) for n in need}
    missing = [n for n, p in paths.items() if not os.path.exists(p)]
    if missing:
        print("！ 尚未找到 " + ", ".join(missing) +
              "；請先執行 Rscript baseline_airpls.R。教學頁將略過基線章節。")
        return None

    bl = pd.read_csv(paths["airpls_baselines.csv"])
    demo = pd.read_csv(paths["airpls_demo.csv"])
    drift = pd.read_csv(paths["airpls_drift.csv"])
    scan = pd.read_csv(paths["airpls_scan.csv"])
    meta = pd.read_csv(paths["airpls_meta.csv"]).set_index("key")["value"]

    # 欄名格式為 "<run>|<lambda 序號>"
    baselines = {}
    for col in bl.columns:
        if col == "time_min":
            continue
        rid, k = col.split("|")
        baselines.setdefault(rid, {})[int(k)] = [int(round(v)) for v in bl[col]]
    # 轉成依 lambda 序號排序的 list
    baselines = {rid: [d[k] for k in sorted(d)] for rid, d in baselines.items()}

    dr = []
    for (run, name), g in drift.groupby(["run", "drift"], sort=False):
        g = g.sort_values("lambda")
        dr.append(dict(run=run, drift=name,
                       lambdas=[float(v) for v in g["lambda"]],
                       err=[float(v) for v in g["err_pct"]]))

    return dict(
        baselines=baselines,
        demo=demo.to_dict("records"),
        drift=dr,
        scan=scan.to_dict("records"),
        bestLambda=float(meta["best_lambda"]),
        bestDiff=int(meta["best_differences"]),
        halfWidth=float(meta["half_width_min"]),
        ref=dict(r2=float(meta["ref_r2"]), maxRecErr=float(meta["ref_max_rec_err"]),
                 conc=[float(meta[f"ref_s{i}"]) for i in range(1, 5)]),
        version=str(meta["airPLS_version"]),
    )


def build_payload():
    runs_out = []
    areas = {}

    for fn in list(STANDARDS) + SAMPLES:
        run = load_run(os.path.join(BASE_DIR, fn))
        rid = fn.replace(".txt", "")
        chrom = run["chrom"]
        if chrom.empty:
            sys.exit(f"{fn} 沒有層析圖資料")

        # 時間軸為等間隔 (0.5 s)，只需存強度陣列即可還原
        y = [int(round(v)) for v in chrom["intensity"].tolist()]

        peaks = []
        pk_tbl = run["peaks"]
        if not pk_tbl.empty:
            for _, r in pk_tbl.iterrows():
                peaks.append(dict(
                    n=int(r["Peak#"]), rt=float(r["R.Time"]),
                    t0=float(r["I.Time"]), t1=float(r["F.Time"]),
                    area=float(r["Area"]), h=float(r["Height"]),
                    named=bool(str(r.get("Name", "")).strip()),
                    conc=float(r["Conc."]),
                ))

        pk = pick_caffeine_peak(pk_tbl, *RT_WINDOW)
        caf = None
        if pk is not None:
            caf = dict(rt=float(pk["R.Time"]), t0=float(pk["I.Time"]),
                       t1=float(pk["F.Time"]), area=float(pk["Area"]),
                       h=float(pk["Height"]),
                       recalc=float(reintegrate(chrom, pk["I.Time"], pk["F.Time"])))
        areas[rid] = caf["area"] if caf else 0.0

        runs_out.append(dict(
            id=rid, kind="std" if fn in STANDARDS else "sample",
            conc=STANDARDS.get(fn), color=COLORS[rid],
            acquired=run["info"].get("Acquired", ""),
            y=y, peaks=peaks, caf=caf,
        ))

    # 檢量線（與 caffeine_calibration.py 相同的算法，供頁面顯示預設值）
    x = np.array([STANDARDS[f] for f in STANDARDS], float)
    yv = np.array([areas[f.replace(".txt", "")] for f in STANDARDS], float)
    fit = stats.linregress(x, yv)
    resid = yv - (fit.intercept + fit.slope * x)
    s_res = float(np.sqrt(np.sum(resid ** 2) / (len(x) - 2)))

    return dict(
        meta=dict(
            nPoints=len(runs_out[0]["y"]), dtSec=0.5, tEnd=6.0,
            wavelength=280, injection=10,
            rtLo=RT_WINDOW[0], rtHi=RT_WINDOW[1],
        ),
        runs=runs_out,
        calib=dict(slope=float(fit.slope), intercept=float(fit.intercept),
                   r2=float(fit.rvalue ** 2), s=s_res,
                   lod=3.3 * s_res / float(fit.slope),
                   loq=10.0 * s_res / float(fit.slope)),
        align=load_alignment(),
        airpls=load_airpls(),
    )


def render(template, output, payload, label):
    if not os.path.exists(template):
        sys.exit(f"找不到模板: {template}")
    with open(template, "r", encoding="utf-8") as f:
        html = f.read()
    if MARKER not in html:
        sys.exit(f"{template} 中找不到資料佔位符 {MARKER}")
    blob = json.dumps(payload, separators=(",", ":"), ensure_ascii=False)
    with open(output, "w", encoding="utf-8") as f:
        f.write(html.replace(MARKER, blob))
    print(f"已產生 {output}  ({os.path.getsize(output)/1024:.0f} KB) — {label}")


def main():
    payload = build_payload()
    render(TEMPLATE, OUTPUT, payload,
           f"進階版，內嵌 {len(payload['runs'])} 個層析圖 + 對齊/基線資料")

    # 入門版只需要層析圖與檢量線，不含 alignDE / airPLS 的大型陣列
    basic = {k: v for k, v in payload.items() if k in ("meta", "runs", "calib")}
    render(BASIC_TEMPLATE, BASIC_OUTPUT, basic, "入門版")

    # --site DIR：另外產生一份給 GitHub Pages 用的版本
    # 網站上入門版是 index.html、進階版是 advanced.html，兩頁互連的檔名要跟著改
    if "--site" in sys.argv:
        site_dir = sys.argv[sys.argv.index("--site") + 1]
        os.makedirs(site_dir, exist_ok=True)
        rewrites = {
            os.path.join(site_dir, "index.html"):
                (BASIC_TEMPLATE, basic, {"hplc_teaching.html": "advanced.html"}, "網站入門版"),
            os.path.join(site_dir, "advanced.html"):
                (TEMPLATE, payload, {"hplc_basic.html": "index.html"}, "網站進階版"),
        }
        for out, (tpl, data, subs, label) in rewrites.items():
            with open(tpl, "r", encoding="utf-8") as f:
                html = f.read()
            for a, b in subs.items():
                if a not in html:
                    sys.exit(f"{tpl} 中找不到待改寫的連結 {a}")
                html = html.replace(f'href="{a}"', f'href="{b}"')
            blob = json.dumps(data, separators=(",", ":"), ensure_ascii=False)
            with open(out, "w", encoding="utf-8") as f:
                f.write(html.replace(MARKER, blob))
            print(f"已產生 {out}  ({os.path.getsize(out)/1024:.0f} KB) — {label}")


if __name__ == "__main__":
    main()
