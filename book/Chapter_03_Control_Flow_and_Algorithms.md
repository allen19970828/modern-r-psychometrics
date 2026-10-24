# 第 3 章：程式邏輯、流程控制與自適應演算法

> **本章目標**：深入掌握 R 語言三大控制結構（循序、條件分支、迴圈反覆運算）在評量工程中的高階應用（洪錦魁、蔡桂宏，2021；Wickham et al., 2023）。透過「十萬考生高並發自動閱卷引擎」與「電腦化適性測驗（CAT）自適應選題演算法」（Kingsbury & Zara, 1989；Lord, 1980；van der Linden & Glas, 2010），完整推導 2PL Fisher 訊息量公式與 EAP 期望後驗數值求積演算法；剖析數值積分下溢（Underflow）防範、Log-Sum-Exp 數值穩定化技巧，以及 Sympson-Hetter 題庫曝題率控制工程（Sympson & Hetter, 1985）。

---

## 3.1 理論基礎：演算法思維與控制流架構

程式設計的核心在於「利用有限步驟解決確定問題」。在測驗評量中，從考生的原始答案到最終產出量表分數，背後皆是一連串嚴密的控制流程。

```mermaid
flowchart TD
    Start([開始: 考生交卷]) --> ReadItem["讀取考生作答 responses 與標準答案 key"]
    ReadItem --> CheckLen{"長度是否相符?"}
    CheckLen -- 否 --> ThrowErr["拋出錯誤 (stop): 題數不一致"]
    CheckLen -- 是 --> LoopItem["逐題核對答案 (Loop / Vectorized)"]
    LoopItem --> CheckCondition{"作答 == 答案?"}
    CheckCondition -- 是 --> Score1["得 1 分"]
    CheckCondition -- 否 --> Score0["得 0 分"]
    Score1 --> SumScore["加總為原始總分 Raw Score"]
    Score0 --> SumScore
    SumScore --> Mastery{"總分 >= 切截分數?"}
    Mastery -- 是 --> Pass["判定: 精熟 (Pass)"]
    Mastery -- 否 --> Fail["判定: 未精熟 (Fail)"]
    Pass --> End([結束])
    Fail --> End
    ThrowErr --> End

    style Start fill:#e1f5fe,stroke:#0288d1,stroke-width:2px
    style End fill:#eceff1,stroke:#607d8b,stroke-width:2px
    style Pass fill:#c8e6c9,stroke:#388e3c,stroke-width:2px
    style Fail fill:#ffcdd2,stroke:#d32f2f,stroke-width:2px
    style ThrowErr fill:#ffccbc,stroke:#d84315,stroke-width:2px
```

---

## 3.2 評量實務情境：條件分支的四種實踐方式

在測驗統計中，條件分支常用於「分數等級轉換（如 A/B/C/D 級距）」與「特殊狀態判定」。

### 3.2.1 四大條件分支語法對照

1. **純量 `if-else`**：僅接受單一布林值（`length == 1`），常用於演算法全域開關或參數檢驗。
2. **多分支 `switch`**：用於依據特定字串標籤（如科目代碼 `"MATH"`, `"ENG"`）執行不同計分腳本。
3. **向量化 `ifelse()` / `dplyr::if_else()`**：針對整欄向量平行進行條件轉換。
4. **多條件映射 `dplyr::case_when()`**：Tidyverse 中常用的多區間分級工具，語法清晰且強制型別一致。

```r
# ==============================================================================
# 第 3 章程式碼清單 3.1：條件分支四種語法於成績等級劃分之應用
# ==============================================================================

# 1. 純量 if-else：全域門檻檢驗
exam_mode <- "Standard"
if (exam_mode == "Standard") {
  pass_cutoff <- 60
} else {
  pass_cutoff <- 70
}
cat(sprintf("當前測驗模式: %s，及格門檻: %d 分\n", exam_mode, pass_cutoff))

# 2. 向量化多條件分級：dplyr::case_when()
library(dplyr)

mock_scores <- c(95, 82, 67, 54, 73, NA)

grade_levels <- case_when(
  is.na(mock_scores) ~ "缺考 (NA)",
  mock_scores >= 90  ~ "A (優等)",
  mock_scores >= 80  ~ "B (良等)",
  mock_scores >= 70  ~ "C (中等)",
  mock_scores >= 60  ~ "D (及格)",
  TRUE               ~ "F (不及格)" # 預設兜底條件
)

cat("\n多條件分級結果 (case_when):\n")
print(tibble(Score = mock_scores, Grade = grade_levels))
```

---

## 3.3 迴圈反覆運算 vs. 向量化：自動閱卷程式效能比較

依據 Wickham 等人（2023）在資料科學工程中的效能原則，在處理十萬級大規模學力檢測資料時，程式架構對執行效能有著數量級（Order of Magnitude）的影響。

### 3.3.1 兩種運算思維
- **迴圈思維（Loop）**：以程式逐步走訪資料並反覆執行操作；R 也能有效執行迴圈，效能取決於資料量與實作方式。
- **向量化思維（Vectorization）**：將整個向量交給向量化函數處理，通常可減少 R 層級的迴圈開銷；不代表所有運算都會平行執行或使用 SIMD，實際效能應以基準測試確認。

```r
# ==============================================================================
# 第 3 章程式碼清單 3.2：迴圈 vs. 向量化閱卷引擎效能基準測試
# ==============================================================================

# 模擬 50,000 位考生、30 道四選一選擇題的作答資料
set.seed(42)
N_students <- 50000
J_items <- 30
options <- c("A", "B", "C", "D")

# 標準答案 Key
answer_key <- sample(options, J_items, replace = TRUE)

# 考生作答矩陣 (50000 x 30)
student_answers <- matrix(
  sample(options, N_students * J_items, replace = TRUE),
  nrow = N_students, ncol = J_items
)

# ------------------------------------------------------------------------------
# 方案 A：手刻傳統 for 迴圈逐題閱卷
# ------------------------------------------------------------------------------
time_loop <- system.time({
  score_mat_loop <- matrix(0L, nrow = N_students, ncol = J_items)
  for (i in 1:N_students) {
    for (j in 1:J_items) {
      if (student_answers[i, j] == answer_key[j]) {
        score_mat_loop[i, j] <- 1L
      } else {
        score_mat_loop[i, j] <- 0L
      }
    }
  }
  total_scores_loop <- rowSums(score_mat_loop)
})

# ------------------------------------------------------------------------------
# 方案 B：現代純向量化廣播比對 (Vectorized Recycling)
# ------------------------------------------------------------------------------
time_vectorized <- system.time({
  # 利用矩陣與向量按 Column 比較的特性，直接進行廣播比對
  score_mat_vec <- (t(student_answers) == answer_key) |>
    t() |>
    as.integer()
  dim(score_mat_vec) <- c(N_students, J_items)
  total_scores_vec <- rowSums(score_mat_vec)
})

# 驗證兩者結果完全相同
stopifnot(identical(total_scores_loop, total_scores_vec))

cat("--- 50,000 位考生 x 30 題閱卷耗時對比 ---\n")
cat(sprintf("傳統雙層 for 迴圈耗時 : %.3f 秒\n", time_loop["elapsed"]))
cat(sprintf("現代向量化廣播運算耗時 : %.3f 秒\n", time_vectorized["elapsed"]))
cat(sprintf("執行時間比值          : %.1f 倍\n", time_loop["elapsed"] / time_vectorized["elapsed"]))
```

---

## 3.4 數理深入：自適應測驗選題準則與 EAP 數值求積推導

電腦化自適應測驗（Computerized Adaptive Testing, CAT）是現代教育評量與演算法工程結合的最前沿範例（Lord, 1980；van der Linden & Glas, 2010）。

### 3.4.1 2PL 試題 Fisher 訊息量推導
設第 $j$ 題的 2PL 答對機率函數為：
$$P_j(\theta) = \frac{1}{1 + e^{-D a_j (\theta - b_j)}}$$
其中 $Q_j(\theta) = 1 - P_j(\theta)$。對 $\theta$ 求一階導數：
$$P_j'(\theta) = \frac{\partial P_j(\theta)}{\partial \theta} = D a_j P_j(\theta) Q_j(\theta)$$
根據 Fisher 訊息量（Fisher Information）定義：
$$I_j(\theta) = \frac{[P_j'(\theta)]^2}{P_j(\theta) Q_j(\theta)} = \frac{[D a_j P_j(\theta) Q_j(\theta)]^2}{P_j(\theta) Q_j(\theta)} = D^2 a_j^2 P_j(\theta) Q_j(\theta)$$
**最大訊息量極值點分析**：
欲使 $I_j(\theta)$ 達到最大，必須最大化二次項 $f(P) = P(1 - P)$。對 $P$ 微分令其為 0：
$$f'(P) = 1 - 2P = 0 \implies P_j(\theta) = 0.5$$
當 $P_j(\theta) = 0.5$ 時，$\theta = b_j$，試題訊息量達到峰值：
$$I_{j,\max} = 0.25 D^2 a_j^2$$
此即「最大 Fisher 訊息量選題準則（Maximum Fisher Information, MFI）」的數理本質：**在考生當前能力估計值 $\hat{\theta}$ 處，挑選難度 $b_j$ 最接近 $\hat{\theta}$ 且鑑別度 $a_j$ 最高的試題**。

### 3.4.2 EAP 期望後驗估計與高斯求積（Gauss-Hermite Quadrature）
在動態作答歷程中，已知已施測題目的反應向量 $\mathbf{X} = (x_1, \dots, x_k)^T$。假設先驗分佈 $\theta \sim \mathcal{N}(0, 1)$，概似函數為 $L(\mathbf{X} \mid \theta) = \prod_{i=1}^k P_i(\theta)^{x_i} Q_i(\theta)^{1 - x_i}$。
後驗期望值（EAP）與後驗變異數為：
$$\hat{\theta}_{\text{EAP}} = \mathbb{E}[\theta \mid \mathbf{X}] = \frac{\int_{-\infty}^{\infty} \theta L(\mathbf{X} \mid \theta) \phi(\theta) d\theta}{\int_{-\infty}^{\infty} L(\mathbf{X} \mid \theta) \phi(\theta) d\theta}$$
$$\text{SE}(\hat{\theta}_{\text{EAP}}) = \sqrt{\frac{\int_{-\infty}^{\infty} (\theta - \hat{\theta}_{\text{EAP}})^2 L(\mathbf{X} \mid \theta) \phi(\theta) d\theta}{\int_{-\infty}^{\infty} L(\mathbf{X} \mid \theta) \phi(\theta) d\theta}}$$
在計算機中，我們在區間 $[-4, +4]$ 上取 $Q$ 個等距求積節點（Quadrature Points）$X_q$（$q=1, \dots, Q$），以離散權重總和逼近連續積分：
$$\hat{\theta}_{\text{EAP}} \approx \sum_{q=1}^Q X_q \cdot W_q^*, \quad W_q^* = \frac{L(\mathbf{X} \mid X_q) \phi(X_q)}{\sum_{m=1}^Q L(\mathbf{X} \mid X_m) \phi(X_m)}$$

---

## 3.5 底層手刻實作：CAT 動態適性測驗引擎與數值穩定化

在數十題作答後，$L(\mathbf{X} \mid X_q)$ 是數十個小於 1 的機率相乘，在浮點數系統中極易發生**數值下溢（Underflow to 0）**。
**解法：Log-Sum-Exp 數值穩定化技巧**：
$$\log W_q = \ell(X_q) + \log \phi(X_q)$$
令 $M = \max_q (\log W_q)$，則正規化權重為：
$$W_q^* = \frac{e^{\log W_q - M}}{\sum_{m=1}^Q e^{\log W_m - M}}$$
此技巧徹底消除了浮點數下溢風險，保證適性測驗穩定運行。

```mermaid
sequenceDiagram
    autonumber
    actor Student as 考生 (Examinee)
    participant Engine as CAT 適性引擎 (R Algorithm)
    participant ItemBank as 試題庫 (Item Bank)

    Engine->>Engine: 初始化估計值 theta_hat = 0.0, SE = 1.0 (常態先驗)
    loop 當 SE > 門檻 且 題數 < 上限 (while 迴圈)
        Engine->>ItemBank: 查詢在當前 theta_hat 下最大 Fisher 訊息量之題目 (MFI)
        ItemBank-->>Engine: 回傳最優試題 (Item_j)
        Engine->>Student: 呈現試題 Item_j
        Student-->>Engine: 提交作答 (0 或 1)
        Engine->>Engine: 即時透過 Log-Sum-Exp EAP 數值積分更新 theta_hat 與 SE
        Engine->>Engine: 檢查終止條件 (Break / Continue)
    end
    Engine->>Student: 測驗結束，回報最終能力估計值與信賴區間
```

```r
# ==============================================================================
# 第 3 章程式碼清單 3.3：高精度數值穩定化 CAT 適性測驗模擬引擎
# ==============================================================================

library(dplyr)
library(tibble)

# 1. 建立 60 題 2PL 試題庫 (涵蓋難度與鑑別度)
set.seed(42)
n_bank <- 60
a_vals <- round(runif(n_bank, min = 1.0, max = 2.0), 2)
b_vals <- round(seq(-2.5, 2.5, length.out = n_bank), 2)
d_vals <- round(-a_vals * b_vals, 2)

item_bank <- tibble(
  item_id = sprintf("BANK_%02d", 1:n_bank),
  a = a_vals,
  b = b_vals,
  d = d_vals,
  is_used = FALSE
)

# 2. 定義真實考生能力 (True Theta = 1.25) 與求積格點
true_theta <- 1.25
quad_pts <- seq(-4, 4, length.out = 81)
log_prior_density <- dnorm(quad_pts, mean = 0, sd = 1, log = TRUE)

# 3. 封裝 Log-Sum-Exp 高精度 EAP 能力估計與標準誤計算函數
calc_eap_stabilized <- function(admin_indices, responses, bank) {
  if (length(admin_indices) == 0) return(list(theta = 0.0, se = 1.0))
  
  # 全程在對數尺度下累積 Log-Likelihood
  log_lik <- numeric(length(quad_pts))
  for (i in seq_along(admin_indices)) {
    idx <- admin_indices[i]
    # 2PL 邏輯斯機率
    logit_val <- bank$a[idx] * quad_pts + bank$d[idx]
    p_vec <- 1 / (1 + exp(-logit_val))
    resp <- responses[i]
    obs_prob <- if (resp == 1) p_vec else (1 - p_vec)
    log_lik <- log_lik + log(pmax(obs_prob, 1e-15))
  }
  
  # 未正規化對數後驗密度: log_unnorm = log_lik + log_prior
  log_unnorm <- log_lik + log_prior_density
  # Log-Sum-Exp 技巧防止下溢
  max_log <- max(log_unnorm)
  weights_unnorm <- exp(log_unnorm - max_log)
  post_weights <- weights_unnorm / sum(weights_unnorm)
  
  # 數值求積期望值與標準誤
  th_hat <- sum(quad_pts * post_weights)
  se_hat <- sqrt(sum(((quad_pts - th_hat)^2) * post_weights))
  list(theta = th_hat, se = se_hat)
}

# 4. 封裝 Fisher 試題訊息量函數
calc_info <- function(th, a, d) {
  p <- 1 / (1 + exp(-(a * th + d)))
  (a^2) * p * (1 - p)
}

# 5. CAT 運行參數設定
current_theta <- 0.0   # 初始猜測能力值
current_se <- 1.0      # 初始標準誤
target_se <- 0.35      # 終止條件：標準誤門檻 (SE <= 0.35)
max_items <- 15        # 終止條件：最大施測題數

admin_items <- integer(0)
admin_resps <- integer(0)
cat_history <- list()
step <- 0

# 6. while 控制流程驅動適性測驗
while (current_se > target_se && step < max_items) {
  step <- step + 1
  
  # A. 從未施測試題中選取在當前估計值下 Fisher 訊息量最大之試題 (MFI)
  available_idx <- which(!item_bank$is_used)
  item_infos <- calc_info(current_theta, item_bank$a[available_idx], item_bank$d[available_idx])
  best_idx <- available_idx[which.max(item_infos)]
  
  # 標記該題已被使用
  item_bank$is_used[best_idx] <- TRUE
  
  # B. 模擬考生真實作答 (依據 2PL 機率模型)
  p_true <- 1 / (1 + exp(-(item_bank$a[best_idx] * true_theta + item_bank$d[best_idx])))
  resp <- if (runif(1) < p_true) 1L else 0L
  
  admin_items <- c(admin_items, best_idx)
  admin_resps <- c(admin_resps, resp)
  
  # C. 即時透過穩定版 EAP 估計當前能力值與精確標準誤
  eap_res <- calc_eap_stabilized(admin_items, admin_resps, item_bank)
  current_theta <- eap_res$theta
  current_se <- eap_res$se
  
  # 紀錄施測軌跡
  cat_history[[step]] <- tibble(
    Step = step,
    Item_ID = item_bank$item_id[best_idx],
    Difficulty_b = item_bank$b[best_idx],
    Discrimination_a = item_bank$a[best_idx],
    Response = resp,
    Theta_Hat = round(current_theta, 3),
    SE = round(current_se, 3)
  )
}

# 整合軌跡報表
cat_trajectory <- bind_rows(cat_history)
cat("--- CAT 適性測驗結束報告 ---\n")
cat(sprintf("真實能力值 (True Theta) : %.2f\n", true_theta))
cat(sprintf("最終估計值 (Final Theta): %.2f (標準誤 SE: %.3f)\n", current_theta, current_se))
cat(sprintf("共施測題數              : %d 題\n\n", nrow(cat_trajectory)))
print(head(cat_trajectory, 6))
```

---

## 3.6 ⚠️ 心理計量工程陷阱與排錯專題（Pitfalls & Diagnostics）

### 3.6.1 陷阱 1：`while` 迴圈中的無窮死迴圈防禦
若考生作答行為異常（如胡亂猜測引發估計值劇烈震盪），標準誤可能始終無法降至目標值（`current_se <= target_se`）。若未設定最大題數防禦，`while` 迴圈將陷入無窮死迴圈。
**防禦性規範**：必須永遠綁定複合終止條件：`while (current_se > target_se && step < max_items)`，且在迴圈開頭累加安全計數器。

### 3.6.2 陷阱 2：題庫過度曝題（Overexposure）與 Sympson-Hetter 控制法
在 MFI 選題中，鑑別度最高的題目幾乎會被所有考生優先選取，曝題率可高達 90% 以上，造成題庫迅速外洩。
Sympson 與 Hetter（1985）提出的曝題控制法為每道題目賦予一個接受機率參數 $k_j \in [0, 1]$：
1. 當題目 $j$ 被 MFI 最佳化選中後，產生隨機數 $U \sim \text{Uniform}(0, 1)$。
2. 若 $U \le k_j$，則正式施測該題。
3. 若 $U > k_j$，則拒絕該題，退回在未施測題庫中挑選次優試題。

### 3.6.3 陷阱 3：向量化比對的靜默長度循環誤判（Recycling Bug）
在向量化閱卷中，若考生作答題數（如 30 題）與答案卡 Key 長度（若誤傳為 15 題）不相符，R 會自動啟動「循環補充機制（Recycling）」，循環比對兩次而**完全不報錯**，導致嚴重的靜默計分災難！
**防禦性對策**：在比對前必須以 `stopifnot(ncol(student_answers) == length(answer_key))` 進行長度一致性斷言。

---

## 3.7 視覺化：CAT 能力估計收斂軌跡與誤差帶圖

我們將適性測驗每一步的估計歷程轉化為圖表，動態呈現考生能力如何快速逼近真值：

```r
# ==============================================================================
# 第 3 章程式碼清單 3.4：CAT 適性測驗動態收斂軌跡繪製
# ==============================================================================

library(ggplot2)

p3 <- ggplot(cat_trajectory, aes(x = Step, y = Theta_Hat)) +
  # 繪製 EAP 後驗 95% 信賴區間帶 (Theta_Hat ± 1.96 * SE)
  geom_ribbon(
    aes(ymin = Theta_Hat - 1.96 * SE, ymax = Theta_Hat + 1.96 * SE),
    fill = "#2b5c8f", alpha = 0.18
  ) +
  # 標記考生真實能力基準線
  geom_hline(yintercept = true_theta, linetype = "dashed", color = "#d95f02", linewidth = 0.9) +
  # 繪製施測試題難度點
  geom_point(aes(y = Difficulty_b, shape = factor(Response)), size = 2.5, color = "grey40") +
  # 繪製估計值路徑線與點
  geom_line(color = "#2b5c8f", linewidth = 1.1) +
  geom_point(color = "#2b5c8f", size = 3) +
  annotate("text", x = 2, y = true_theta + 0.25, label = "True Theta = 1.25", color = "#d95f02", fontface = "bold") +
  scale_shape_manual(name = "作答結果", values = c("0" = 4, "1" = 16), labels = c("0" = "答錯 (X)", "1" = "答對 (O)")) +
  scale_x_continuous(breaks = 1:nrow(cat_trajectory)) +
  labs(
    title = "Figure 3.1: 電腦化適性測驗 (CAT) 動態選題與能力估計收斂歷程",
    subtitle = "藍線為 EAP 後驗能力估計值；陰影為 EAP 95% 信賴區間帶；測驗於 SE <= 0.35 時收斂終止",
    x = "施測題數 (Step / Items Administered)",
    y = "潛在特質能力估計值 (Theta)",
    caption = "Algorithm: Maximum Fisher Information (MFI) & Expected A Posteriori (EAP)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 13),
    axis.text = element_text(color = "black")
  )

ggsave("book/fig_3_1_cat_trajectory.png", plot = p3, width = 8, height = 5, dpi = 300)
cat("成功產出圖表: book/fig_3_1_cat_trajectory.png\n")
```

---

## 3.8 測驗報告專業解讀（Substantive Interpretation）

在適性測驗與自動閱卷演算法的報告中，應重點評估以下指標：

> ### 🧠 演算法效能與適配指針
> 1. **選題曝題率（Item Exposure Rate）**：在正式 CAT 題庫中，若特定試題被選取的頻率超過 $25\%$，將引發試題外洩風險，必須在 `while` 迴圈中加入隨機抽樣控制（Sympson & Hetter, 1985）。
> 2. **能力估計收斂速度**：正常情況下，標準誤應在 $10 \sim 15$ 題內降至 $0.35$ 以下。若出現振盪不收斂（如考生隨機作答），迴圈必須強制在 `max_items` 時跳出（`break`），並標記為「作答無效（Aberrant Response）」（van der Linden & Glas, 2010）。

---

## 3.9 章末實務習題與專案挑戰（含 100% 完整可用程式碼解答）

### 【基礎練習題】
1. 請使用 `for` 迴圈計算 $1 \sim 100$ 的整數總和，隨後改用 R 語言原生向量化函數 `sum(1:100)` 實現，並比較兩者的程式碼簡潔度。
2. 某班級成績向量為 `scores <- c(78, 45, 92, 60, 88, 30)`。請使用 `dplyr::if_else` 寫出一段程式碼，及格（$\ge 60$）標記為 `"PASS"`，不及格標記為 `"FAIL"`。

### 【進階工程挑戰】
3. 請改寫本章 3.5 節的 CAT 演算法：加入一個「提前中止機制（Early Termination）」。若考生連續 3 題皆在高難度試題（$b > 2.0$）答錯，或者連續 3 題在低難度試題（$b < -2.0$）答對，立即以 `break` 終止測驗，並在輸出報告中註記「已達到能力天花板/地板邊界」。
4. 【內容平衡適性測驗挑戰】：Kingsbury 與 Zara（1989）指出，實際測驗必須滿足構面領域的配比約束。假設題庫包含「代數（Algebra）」與「幾何（Geometry）」兩大領域，請設計一個約束選題函數 `select_item_constrained(available_bank, current_theta, target_ratio)`，當代數已達上限時，強制切換至幾何領域選題，確保測驗結束時領域配比符合標準藍圖。

---

### 💡 習題參考解答

```r
# ==============================================================================
# 第 3 章習題解答參考
# ==============================================================================

# ------------------------------------------------------------------------------
# 第 1 題解答
# ------------------------------------------------------------------------------
s <- 0
for (i in 1:100) s <- s + i
stopifnot(s == sum(1:100))
cat("1 到 100 總和:", s, "\n")

# ------------------------------------------------------------------------------
# 第 2 題解答
# ------------------------------------------------------------------------------
library(dplyr)
scores <- c(78, 45, 92, 60, 88, 30)
results <- if_else(scores >= 60, "PASS", "FAIL")
cat("判定結果:", results, "\n")

# ------------------------------------------------------------------------------
# 第 3 題解答：CAT 防呆提前終止函數
# ------------------------------------------------------------------------------
check_early_stop <- function(history_df) {
  if (nrow(history_df) < 3) return(list(stop = FALSE, reason = "None"))
  last_3 <- tail(history_df, 3)
  
  # 連續 3 題在難題答錯 (天花板效應)
  hard_fails <- all(last_3$Difficulty_b > 1.8 & last_3$Response == 0)
  # 連續 3 題在極簡單題答對 (地板效應)
  easy_passes <- all(last_3$Difficulty_b < -1.8 & last_3$Response == 1)
  
  if (hard_fails) return(list(stop = TRUE, reason = "Ceiling Barrier (連續高難度答錯)"))
  if (easy_passes) return(list(stop = TRUE, reason = "Floor Barrier (連續極低難度答對)"))
  list(stop = FALSE, reason = "None")
}

# 模擬測試第 3 題
mock_history <- tibble::tibble(
  Difficulty_b = c(1.9, 2.1, 2.3),
  Response = c(0, 0, 0)
)
cat("提前終止檢驗結果:\n")
print(check_early_stop(mock_history))

# ------------------------------------------------------------------------------
# 第 4 題解答：領域內容平衡約束選題器 (Kingsbury & Zara, 1989)
# ------------------------------------------------------------------------------
select_item_constrained <- function(bank_df, cur_th, admin_domains, target_prop = c("Algebra" = 0.5, "Geometry" = 0.5)) {
  # 計算已施測各領域比例
  total_admin <- length(admin_domains)
  if (total_admin == 0) {
    target_domain <- names(target_prop)[1]
  } else {
    current_props <- table(factor(admin_domains, levels = names(target_prop))) / total_admin
    # 挑選落後目標最嚴重的領域
    domain_diffs <- target_prop - current_props
    target_domain <- names(which.max(domain_diffs))
  }
  
  # 僅在目標領域候選題中挑選 MFI 最大者
  candidates <- bank_df |> filter(!is_used, Domain == target_domain)
  if (nrow(candidates) == 0) {
    # 該領域已耗竭，退回挑選全題庫中可用題
    candidates <- bank_df |> filter(!is_used)
  }
  
  infos <- (candidates$a^2) * (1 / (1 + exp(-(candidates$a * cur_th + candidates$d)))) *
           (1 - 1 / (1 + exp(-(candidates$a * cur_th + candidates$d))))
  best_row <- candidates[which.max(infos), ]
  list(Selected_Item = best_row$item_id, Domain = best_row$Domain)
}

# 測試第 4 題
mock_bank_domain <- tibble::tibble(
  item_id = sprintf("T_%02d", 1:10),
  Domain = rep(c("Algebra", "Geometry"), each = 5),
  a = runif(10, 1.2, 1.8),
  d = rnorm(10),
  is_used = FALSE
)
res_constrained <- select_item_constrained(mock_bank_domain, cur_th = 0.5, admin_domains = c("Algebra", "Algebra"))
cat("內容平衡選題結果 (強制平衡至幾何):", res_constrained$Selected_Item, "領域:", res_constrained$Domain, "\n")
```

---

## 參考文獻 (References)

洪錦魁、蔡桂宏（2021）。**R 語言邁向 Big Data 之路：王者歸來**（第 2 版）。深智數位。

Kingsbury, G. G., & Zara, A. R. (1989). Procedures for selecting items for computerized adaptive tests. *Applied Measurement in Education*, *2*(4), 359–375. https://doi.org/10.1207/s15324818ame0204_6

Lord, F. M. (1980). *Applications of item response theory to practical testing problems*. Lawrence Erlbaum Associates.

Sympson, J. B., & Hetter, R. D. (1985, October). Controlling item-exposure rates in computerized adaptive testing. In *Proceedings of the 27th Annual Conference of the Military Testing Association* (pp. 973–977). Navy Personnel Research and Development Center.

van der Linden, W. J., & Glas, C. A. W. (Eds.). (2010). *Elements of adaptive testing*. Springer. https://doi.org/10.1007/978-0-387-85461-8

Wickham, H., Çetinkaya-Rundel, M., & Grolemund, G. (2023). *R for data science: Import, tidy, transform, visualize, and model data* (2nd ed.). O'Reilly Media.
