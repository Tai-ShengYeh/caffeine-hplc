# HPLC 咖啡因檢量線與樣品定量

以 R 與 Python 兩套程式，從 Shimadzu LabSolutions ASCII 匯出檔直接建立檢量線並計算樣品咖啡因濃度。兩者結果完全一致。

## 執行方式

主定量（R 與 Python 各一套，結果一致）：

```bash
python caffeine_calibration.py
```

```bash
Rscript caffeine_calibration.R
```

> 以上指令假設 `Rscript` 已在系統 PATH 中。Windows 使用者若尚未把 R 加入 PATH，
> 可改用自己 R 安裝路徑下 `Rscript.exe` 的完整路徑來取代 `Rscript`（例如 R 預設安裝目錄下的
> `bin\Rscript.exe`），依實際安裝位置調整即可。

峰對齊（alignDE）與對齊後重新定量：

```bash
Rscript align_chromatograms.R
```

airPLS 基線校正可行性測試：

```bash
Rscript baseline_airpls.R
```

產生互動教學網頁（入門版 `hplc_basic.html` + 進階版 `hplc_teaching.html`）：

```bash
python build_teaching_html.py
```

Python 需 `numpy / pandas / scipy / matplotlib`；`caffeine_calibration.R` 只用 base R，
`align_chromatograms.R` 需 `alignDE`（相依 `DEoptim`、`Matrix`），
`baseline_airpls.R` 另需 `airPLS`。

安裝 alignDE：

```r
install.packages(c("DEoptim", "Matrix", "remotes"))
remotes::install_github("Tai-ShengYeh/alignDE")
packageVersion("alignDE")   # 需要 >= 3.0.1
```

> **不要安裝 `zmzhang/alignDE`。** 那是原作者的上游 repo，停留在 2012 年的
> 2.0.1：`widthEstimationCWT()` 在訊號有平台時會丟出
> `number of items to replace is not a multiple of replacement length`，
> 內建的純 R 版 DEoptim 也會讓對齊從幾秒變成幾分鐘。
> 請安裝 `Tai-ShengYeh/alignDE`（v3.0.1，LGPL >= 2）。

| 檔案 | 角色 |
|---|---|
| `hplc_io.R` | Shimadzu 匯出檔解析器（R 各程式共用） |
| `caffeine_calibration.py` / `.R` | 主定量流程 |
| `align_chromatograms.R` | alignDE 峰對齊 + 對齊後定量 |
| `baseline_airpls.R` | airPLS 基線校正參數掃描與可行性測試 |
| `build_teaching_html.py` | 產生兩個互動教學頁 |
| `hplc_basic_template.html` | 入門版模板 → `hplc_basic.html` |
| `hplc_teaching_template.html` | 進階版模板 → `hplc_teaching.html` |

## 資料

| 檔案 | 用途 | 標稱濃度 (ppm) |
|---|---|---|
| `std-0.txt` … `std-200.txt` | 標準品 | 0, 25, 50, 100, 150, 200 |
| `sample-1.txt` … `sample-4.txt` | 未知樣品 | — |

條件：Detector A-Ch1、280 nm、注射量 10 µL、Sample Amount = 1、Dilution Factor = 1、跑程 6 min。

## 分析流程

1. **解析**：檔案為 Big5 編碼，切出 `[Sample Information]`、`[Peak Table(Detector A-Ch1)]`、`[LC Chromatogram(Detector A-Ch1)]` 三個區段。
2. **判峰**：**不採用儀器的自動化合物指認**，改取滯留時間視窗 **4.20–4.90 min** 內面積最大的峰。
3. **重新積分**：以峰起訖點（I.Time / F.Time）連線為基線，對原始層析圖做梯形積分，交叉驗證儀器報告的面積。
4. **檢量線**：面積對濃度做最小平方線性迴歸（含 0 ppm 空白點，n = 6）。
5. **樣品**：以檢量線反推濃度，並依 Miller & Miller 公式給 95% 信賴區間。

### 為什麼要自己判峰

儀器的自動指認在多個檔案是錯的：

| 檔案 | 儀器指認 | 實際咖啡因峰 |
|---|---|---|
| `std-25` | 4.177 min，面積 174,985 | 4.653 min，面積 1,080,560 |
| `sample-1` | 未指認 | 4.446 min，面積 4,378,116 |
| `sample-3` | 3.242 min，面積 21,355 | 4.594 min，面積 4,085,132 |
| `sample-4` | 4.028 min，面積 17,636 | 4.578 min，面積 4,114,384 |

若照用儀器的 `Conc.` 欄位，sample-3、sample-4 會得到負濃度。層析圖（`outputs/chromatograms_*.png`）可見標準品與樣品在 4.2–4.9 min 都只有單一主峰，空白 `std-0` 在該區間無峰，佐證此判峰視窗正確。

## 結果

### 檢量線

```
Area = 48,393.42 × C − 111,613.06        R² = 0.999435   (n = 6)
C(ppm) = (Area + 111,613.06) / 48,393.42
```

| 項目 | 值 |
|---|---|
| 斜率 b | 48,393.42 ± 575.50 area/ppm |
| 截距 a | −111,613.06 ± 64,610.38（t = −1.73，未達顯著，可視為過原點）|
| 迴歸標準誤 s | 99,158.97（df = 4）|
| LOD (3.3 s/b) | 6.76 ppm |
| LOQ (10 s/b) | 20.49 ppm |
| ANOVA | F = 7071，p = 1.2 × 10⁻⁷ |

各標準品回算濃度與標稱值相差在 −1.7% ~ +3.1%（0 ppm 空白回算為 2.31 ppm，在 LOD 以下）。

### 樣品濃度

| 樣品 | RT (min) | 峰面積 | 濃度 (ppm) | 95% CI |
|---|---|---|---|---|
| sample-1 | 4.446 | 4,378,116 | **92.78** | 86.63 – 98.92 |
| sample-2 | 4.428 | 3,952,115 | **83.97** | 77.83 – 90.12 |
| sample-3 | 4.594 | 4,085,132 | **86.72** | 80.58 – 92.87 |
| sample-4 | 4.578 | 4,114,384 | **87.33** | 81.18 – 93.47 |

四個樣品都落在檢量線範圍內且遠高於 LOQ。單位同標準品配製單位（ppm = µg/mL）；因 Dilution Factor = 1，若前處理實際有稀釋，需再乘上稀釋倍數。

## 峰對齊（alignDE）

用 [alignDE](https://github.com/Tai-ShengYeh/alignDE) 3.0.1（Zhang, Chen & Liang 2011, *Talanta* 83:1108）
對 10 個層析圖做滯留時間校準：CWT（Mexican Hat）找峰 → Whittaker 懲罰最小平方扣基線 →
峰群合併 → 差分演化搜尋位移使與參考圖的相關係數最大。參考圖為 `std-100`（峰頂 4.658 min），
參數 `slack = 45`（點）、`n = 4`、`NP = 90`、`itermax = 200`，固定亂數種子 20161017。

| 檔案 | 對齊前相關 | 對齊後相關 | 峰頂前 (min) | 峰頂後 (min) |
|---|---|---|---|---|
| std-25 | 0.9987 | 0.9989 | 4.658 | 4.667 |
| std-50 | 0.9968 | 0.9999 | 4.650 | 4.658 |
| std-150 | 0.9736 | 0.9995 | 4.683 | 4.658 |
| std-200 | 0.9925 | 0.9993 | 4.675 | 4.667 |
| sample-1 | 0.1491 | 0.8887 | 4.450 | 4.683 |
| sample-2 | 0.1009 | 0.8865 | 4.433 | 4.675 |
| sample-3 | 0.6961 | 0.8704 | 4.600 | 4.692 |
| sample-4 | 0.6496 | 0.8830 | 4.583 | 4.683 |

峰頂與參考的偏差（絕對值平均）由 0.0685 min 降到 0.0130 min；相關係數平均 0.7286 → 0.9474。
`std-0` 是空白，該區間本來就沒有咖啡因峰，不納入統計。

### 對齊買到了什麼

alignDE 依設計會保留峰形與峰面積（只平移峰、峰間線性內插），所以**對齊本身不改變面積**。
它的價值是讓所有檔案能共用同一個固定積分視窗。以參考峰頂為中心、不同半寬的固定視窗重跑定量：

| 視窗半寬 | 未對齊 s1 / s2 / s3 / s4 | 對齊後 s1 / s2 / s3 / s4 |
|---|---|---|
| ±0.20 min | 47.1 / 39.9 / 75.9 / 79.3 | 87.7 / 81.0 / 80.9 / 84.1 |
| ±0.25 min | 61.5 / 54.5 / 81.8 / 84.8 | 90.9 / 82.7 / 84.3 / 86.3 |
| ±0.30 min | 74.4 / 68.0 / 85.0 / 87.2 | 92.3 / 83.4 / 85.7 / 87.1 |
| ±0.40 min | 88.9 / 82.0 / 86.8 / 87.9 | 93.3 / 83.9 / 86.6 / 87.7 |
| ±0.50 min | 92.6 / 83.9 / 86.9 / 87.9 | 93.6 / 84.0 / 86.9 / 87.9 |

各樣品在五種視窗下的濃度**全距平均：未對齊 27.3 ppm → 對齊後 4.7 ppm**。

### 對齊後的定量結果

用聯集決定的固定視窗 4.250–5.117 min（樣品峰比標準品寬，若只照參考圖決定視窗會把樣品峰切掉）：

| 樣品 | 原始（儀器積分邊界） | 固定視窗・未對齊 | 固定視窗・對齊後 |
|---|---|---|---|
| sample-1 | 92.78 | 89.48 | **93.42** |
| sample-2 | 83.97 | 82.38 | **83.98** |
| sample-3 | 86.72 | 86.83 | **86.78** |
| sample-4 | 87.33 | 87.89 | **87.84** |

對齊後的檢量線 `Area = 47,914.16 × C − 129,973.59`（R² = 0.999396）。
斜率與原始的 48,393 差約 1%，來自基線模型不同（Whittaker vs 儀器的線性基線）；
標準品與樣品同樣處理，檢量線會吸收這個系統差異，最終濃度與原始流程在 0.7% 內一致。

## 互動教學網頁

已上線：

| 版本 | 網址 |
|---|---|
| 入門版 | https://tai-shengyeh.github.io/hplc-caffeine/ |
| 進階版 | https://tai-shengyeh.github.io/hplc-caffeine/advanced.html |

原始碼與資料：https://github.com/Tai-ShengYeh/Tai-ShengYeh.github.io/tree/main/hplc-caffeine

部署指令（`--site` 會把兩頁另存為 `index.html` / `advanced.html` 並自動改寫互連檔名）：

```bash
python build_teaching_html.py --site <repo>/hplc-caffeine
```

兩個版本皆為單一 HTML 檔、無外部相依，可直接用瀏覽器開或分享；兩頁互相連結。

### 入門版 `hplc_basic.html`（73 KB）

給第一次接觸 HPLC 定量的人，只保留核心觀念鏈：**峰 → 面積 → 檢量線 → 濃度**。
不出現基線校正、峰對齊、統計檢定等進階內容。五個步驟：

1. **層析圖** — 拖滑桿換六個標準品濃度，看峰長高；可疊圖比較
2. **峰面積** — 一鍵把峰「壓扁變寬」，親眼看峰高掉 43% 但面積只變 0.07%，
   解釋為什麼定量一律用面積
3. **檢量線** — 按鈕逐點加入標準品，線與 R² 即時更新
4. **測樣品** — 面積 → 檢量線 → 濃度的三段動畫，附「換你算算看」計算練習
5. **小測驗** — 三題選擇題（為何用面積、空白的意義、超出檢量線範圍怎麼辦），答錯有解釋

結尾附 ppm 單位換算與完整流程回顧，並導向進階版。

### 進階版 `hplc_teaching.html`（310 KB）

九個步驟：

1. 定量的邏輯（流程圖）
2. 層析圖 — 可勾選疊圖、縮放、游標讀值
3. **基線校正** — λ 滑桿即時看 airPLS 基線爬進峰底下，含漂移回收測試圖
4. **峰對齊** — 對齊前/後切換疊圖、診斷表、視窗寬度敏感度圖
5. 峰積分 — 可拖曳積分邊界、調整梯形數，即時比對儀器面積
6. 檢量線 — 可切換是否納入空白點／強制過原點，含殘差圖
7. 樣品定量 — 反推動畫與 95% 預測區間帶
8. 判峰陷阱 — 儀器誤判峰對照
9. 結果總表 + 參考文獻（14 篇，行內引註可點跳轉）

兩頁的圖形都由內嵌的原始 721 點資料即時繪製，計算與 R／Python 程式一致
（例如入門版最後算出的 sample-2 = 84.0 ppm，與主流程的 83.97 ppm 相同）。

`build_teaching_html.py` 一次產生兩頁：進階版吃完整資料，入門版只吃層析圖與檢量線。
建置順序：先跑 `caffeine_calibration.py`、`align_chromatograms.R`、`baseline_airpls.R`，
再跑 `build_teaching_html.py`（缺哪一步，進階版會自動略過對應章節並提示；入門版不受影響）。

## 基線校正：airPLS 可行性測試

[baseline_airpls.R](baseline_airpls.R) 實測 airPLS 3.0.0 是否適用於本資料
（airPLS 與 alignDE 是同一組作者；`airPLS()` 回傳的是**基線**，校正訊號 = x − airPLS(x)）。

```bash
Rscript baseline_airpls.R
```

三段測試：[A] λ 與 differences 參數掃描、[B] 疊加已知漂移的回收測試、
[C] 用 airPLS 取代 alignDE 內建 `baselineCorrectionCWT` 跑完整流程。

### 結論：可以用，但 λ 要挑對

`λ = 10⁴–10⁵`、`differences = 2` 表現良好——峰高保留 99.7–100.2%（沒把峰吃掉），
檢量線 R² 與線性基線相當，最終樣品濃度與現有兩套流程差 ≤ 1 ppm。

| λ (differences = 2) | 檢量線 R² | 回收率最大偏差 | 峰高保留 | 斜率比 | 樣品濃度 s1/s2/s3/s4 |
|---|---|---|---|---|---|
| 10³ | 0.999300 | 4.3% | 98.8% | 0.973 | 84.9 / 81.5 / 83.6 / 83.3 |
| 10⁴ | 0.999351 | 4.9% | 99.7% | 0.990 | 91.0 / 83.0 / 85.8 / 86.1 |
| **10⁵** | **0.999389** | **2.9%** | **100.2%** | **0.996** | **92.4 / 83.3 / 85.9 / 87.1** |
| 10⁶ | 0.999372 | 3.0% | 100.4% | 1.000 | 92.7 / 84.3 / 86.3 / 87.2 |
| 10⁷ | 0.999453 | 3.0% | 100.5% | 1.001 | 92.7 / 84.5 / **92.1** / 87.1 |

基準（線性基線，等同儀器做法）：R² = 0.999404、回收率最大偏差 4.0%、
樣品濃度 92.03 / 83.40 / 85.86 / 87.21。

兩端都會出事：

- **λ 太小**（10³）：基線爬進峰底下，峰高只剩 98.8%、斜率低 2.7%，sample-1 從 92.0 掉到 84.9 ppm。
- **λ 太大**（≥10⁶ 且 d=2，或 ≥10⁴ 且 d=1）：基線變成近乎水平線，sample-3 被高估到 92–118 ppm。
- **differences = 1** 在本資料整體偏不穩，建議一律用 `differences = 2`。

### 漂移回收測試

在原始訊號疊上已知漂移，比較有／無疊加兩次校正後的面積變化（越接近 0 越好）：

| 漂移型態 | λ=10³ | λ=10⁴ | λ=10⁵ | λ=10⁶ | λ=10⁷ |
|---|---|---|---|---|---|
| 線性斜坡 (std-100) | −11.8% | −1.7% | −1.1% | −0.3% | −0.4% |
| 指數衰減 (std-100) | −5.4% | −1.7% | −0.1% | +16.7% | +17.4% |
| 緩弧半週期 (std-100) | −11.8% | −6.3% | −1.1% | −0.3% | +86.4% |
| 起伏 1.5 週期 (sample-1) | −16.2% | −1.9% | **+51.7%** | +57.4% | +58.3% |

`λ = 10⁴` 對各種漂移形狀最穩健（誤差都 ≤ 2%）；`λ = 10⁵` 雖然在本資料上與基準最接近，
但遇到週期約 4 分鐘的起伏型基線會失效。**基線只是平緩漂移用 10⁵；不確定就用 10⁴。**

### 與現有做法的關係

alignDE 內建的 `baselineCorrectionCWT` 也是 Whittaker 懲罰最小平方（同樣源自 Eilers 2003），
差別在它用 CWT 偵測到的峰區決定哪裡該擬合，airPLS 則靠自適應加權自動判斷、不需先偵測峰。
完整流程（airPLS + alignDE 對齊 + 固定視窗）的結果：

| 樣品 | 原始流程 | 現行（CWT 基線） | airPLS 基線 |
|---|---|---|---|
| sample-1 | 92.78 | 93.42 | 92.34 |
| sample-2 | 83.97 | 83.98 | 83.28 |
| sample-3 | 86.72 | 86.78 | 85.83 |
| sample-4 | 87.33 | 87.84 | 87.13 |

三者差異 ≤ 1 ppm，可互相替換。

## 參考文獻

書目資料已透過 Crossref 以 DOI 逐筆核對（第 11、12 項為法規指引與教科書，無 DOI）。

### A · 峰對齊與訊號前處理

1. Zhang, Z.-M., Chen, S., & Liang, Y.-Z. (2011). Peak alignment using wavelet pattern matching and differential evolution. *Talanta*, 83(4), 1108–1117. https://doi.org/10.1016/j.talanta.2010.08.008
   — alignDE 方法原始論文，`align_chromatograms.R` 的依據。
2. Du, P., Kibbe, W. A., & Lin, S. M. (2006). Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching. *Bioinformatics*, 22(17), 2059–2065. https://doi.org/10.1093/bioinformatics/btl355
   — CWT ridge-line 峰偵測；alignDE 的 `cwt()` / `getRidge()` / `identifyMajorPeaks()` 改寫自此。
3. Eilers, P. H. C. (2003). A perfect smoother. *Analytical Chemistry*, 75(14), 3631–3636. https://doi.org/10.1021/ac034173t
   — Whittaker 懲罰最小平方；`WhittakerSmooth()`、`baselineCorrectionCWT()` 與 airPLS 的共同理論基礎。
4. Zhang, Z.-M., Chen, S., & Liang, Y.-Z. (2010). Baseline correction using adaptive iteratively reweighted penalized least squares. *The Analyst*, 135(5), 1138–1146. https://doi.org/10.1039/b922045c
   — airPLS 原始論文，`baseline_airpls.R` 的依據。
5. Storn, R., & Price, K. (1997). Differential evolution — A simple and efficient heuristic for global optimization over continuous spaces. *Journal of Global Optimization*, 11(4), 341–359. https://doi.org/10.1023/A:1008202821328
6. Mullen, K., Ardia, D., Gil, D., Windover, D., & Cline, J. (2011). DEoptim: An R package for global optimization by differential evolution. *Journal of Statistical Software*, 40(6), 1–26. https://doi.org/10.18637/jss.v040.i06
   — alignDE 3.0.0 實際呼叫的最佳化引擎。
7. Nielsen, N.-P. V., Carstensen, J. M., & Smedsgaard, J. (1998). Aligning of single and multiple wavelength chromatographic profiles for chemometric data analysis using correlation optimised warping. *Journal of Chromatography A*, 805(1–2), 17–35. https://doi.org/10.1016/S0021-9673(98)00021-1
   — COW 原始論文，層析圖對齊的基準方法。
8. Tomasi, G., van den Berg, F., & Andersson, C. (2004). Correlation optimized warping and dynamic time warping as preprocessing methods for chromatographic data. *Journal of Chemometrics*, 18(5), 231–241. https://doi.org/10.1002/cem.859
   — COW 與 DTW 的比較；說明伸縮型方法為何會改變峰面積，而 alignDE 不會。

### B · 檢量線、偵測極限與方法確效

9. Danzer, K., & Currie, L. A. (1998). Guidelines for calibration in analytical chemistry. Part I. Fundamentals and single component calibration (IUPAC Recommendations 1998). *Pure and Applied Chemistry*, 70(4), 993–1014. https://doi.org/10.1351/pac199870040993
10. Currie, L. A. (1995). Nomenclature in evaluation of analytical methods including detection and quantification capabilities (IUPAC Recommendations 1995). *Pure and Applied Chemistry*, 67(10), 1699–1723. https://doi.org/10.1351/pac199567101699
11. ICH (2023). *ICH Harmonised Guideline Q2(R2): Validation of Analytical Procedures.* International Council for Harmonisation. https://database.ich.org/sites/default/files/ICH_Q2%28R2%29_Guideline_2023_1130.pdf
    — 本專案 LOD = 3.3 σ/S、LOQ = 10 σ/S 的來源（σ 以迴歸標準誤估計、S 為斜率）。
12. Miller, J. N., & Miller, J. C. (2010). *Statistics and Chemometrics for Analytical Chemistry* (6th ed., Ch. 5). Pearson Education, Harlow.
    — 反推濃度標準誤 sC = (s/b)√(1/m + 1/n + (C−C̄)²/Sxx) 與信賴區間的推導。

### C · 咖啡因的 HPLC 定量

13. DiNunzio, J. E. (1985). Determination of caffeine in beverages by high performance liquid chromatography. *Journal of Chemical Education*, 62(5), 446. https://doi.org/10.1021/ed062p446
14. Naik, J. P., & Nagalakshmi, S. (1997). Determination of caffeine in tea products by an improved high-performance liquid chromatography method. *Journal of Agricultural and Food Chemistry*, 45(10), 3973–3975. https://doi.org/10.1021/jf970147i

### 軟體

| 軟體 | 版本 |
|---|---|
| R | 4.6.1 |
| alignDE | 3.0.0 |
| DEoptim | 2.2.8 |
| Python | 3.12 |
| NumPy / SciPy / pandas / Matplotlib | — |

## 注意事項

- **std-0 與 sample-1、sample-2 是 2016/10/14 跑的，其餘為 10/17。** 本分析把 10/17 的五點標準品加上 10/14 的空白合併成單一檢量線。sample-2 儀器原本報 93.26 ppm（用它自己 10/14 的舊曲線，斜率 46,424），本程式以合併曲線算得 83.97 ppm。若要嚴格分批定量，應各自建線——但 10/14 只有空白一點，無法單獨建線。
- 檢量線有輕微正曲率（響應因子由 43,222 升到 48,154 area/ppm）。二次式擬合 R² = 0.999861，僅比線性高 0.0004，故仍採線性。
- 樣品峰比標準品峰寬（A/H ≈ 17–19 vs 11，理論板數 1400–1800 vs 約 4000）且滯留時間略早，屬基質效應；面積積分不受影響。

## 輸出

`outputs/` 下：

| 檔案 | 內容 |
|---|---|
| `report_py.txt` / `report_r.txt` | 完整文字報告 |
| `calibration_py.csv` / `calibration_r.csv` | 檢量線各點、殘差、回算濃度 |
| `samples_py.csv` / `samples_r.csv` | 樣品濃度與信賴區間 |
| `calibration_py.png` / `calibration_r.png` | 檢量線 + 殘差圖 |
| `chromatograms_py.png` / `chromatograms_r.png` | 標準品與樣品層析圖疊圖 |

## 資料來源

層析資料由 Shimadzu LCsolution 匯出為 ASCII。依檔案表頭的 `Acquired` 欄位，
採集時間為 **2016/10/14 與 10/17**（`Output Date` 顯示的 2019/1/7 是 ASCII
匯出日期，不是量測日期）。六個標準品（0、25、50、100、150、200 ppm）
與四個樣品，樣品以編號標示、未指涉任何商品或品牌。

資料照原樣提供，用途是教學示範——說明檢量線建立與樣品定量的流程——
不是可引用的分析方法確效資料。

## 授權

| 範圍 | 授權 |
|---|---|
| 實驗資料與教材內容（`std-*.txt`、`sample-*.txt`、`hplc_*.html`、`outputs/`、本 README） | [CC BY 4.0](LICENSE) |
| 程式碼（`*.R`、`*.py`） | [MIT](LICENSE-CODE) |

alignDE 套件本身不在此授權範圍內：方法與原始實作為 Zhi-Min Zhang 等人所有，
本專案使用的維護分支以 LGPL (>= 2) 散布於
<https://github.com/Tai-ShengYeh/alignDE>。airPLS 亦有其自身條款。
