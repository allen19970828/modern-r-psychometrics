# ==============================================================================
# 測驗統計程式設計：現代 R 語言導向
# 全書貫穿式模擬資料集生成腳本：AssessmentData
# ==============================================================================

set.seed(20261023)
suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

n_students <- 1000

# 1. 產生學生背景與潛在特質 (Latent Traits)
# 潛在能力 theta (認知) 與潛在學習動機 eta (情意)，設定相關係數約 0.45
mu <- c(0, 0)
sigma_mat <- matrix(c(1.0, 0.45, 0.45, 1.0), nrow = 2)
chol_mat <- chol(sigma_mat)
raw_traits <- matrix(rnorm(n_students * 2), ncol = 2) %*% chol_mat

theta <- as.vector(raw_traits[, 1])
eta <- as.vector(raw_traits[, 2])

# 性別、學校類型與社經地位
gender <- sample(c("F", "M"), size = n_students, replace = TRUE, prob = c(0.51, 0.49))
school_type <- sample(c("Public", "Private"), size = n_students, replace = TRUE, prob = c(0.75, 0.25))
ses_level <- sample(c("Low", "Medium", "High"), size = n_students, replace = TRUE, prob = c(0.25, 0.55, 0.20))

# 2. 產生 15 題二分選擇題 (Cognitive Items 0/1) - 基於 2PL IRT 模型
# P(X = 1 | theta) = 1 / (1 + exp(-1.7 * a * (theta - b)))
a_params <- c(0.8, 1.1, 1.4, 0.9, 1.3, 1.6, 0.7, 1.2, 1.5, 1.0, 1.2, 1.7, 0.9, 1.4, 1.2)
b_params <- c(-1.8, -1.3, -0.9, -0.6, -0.3, 0.0, 0.2, 0.5, 0.8, 1.1, 1.3, 1.6, -0.5, 0.3, 0.0)

cognitive_mat <- matrix(0, nrow = n_students, ncol = 15)
for (j in 1:15) {
  # 第 15 題人為注入性別差異試題功能 (DIF): 男性的難度低 0.75 (對男性有利的均勻 DIF)
  b_eff <- if (j == 15) {
    ifelse(gender == "M", b_params[j] - 0.75, b_params[j])
  } else {
    rep(b_params[j], n_students)
  }
  prob <- 1 / (1 + exp(-1.7 * a_params[j] * (theta - b_eff)))
  cognitive_mat[, j] <- ifelse(runif(n_students) < prob, 1, 0)
}
colnames(cognitive_mat) <- sprintf("item_%02d", 1:15)

# 3. 產生 10 題李克特 5 點量表 (Affective Scale: 1 ~ 5)
# scale_01 ~ scale_05: 內在動機 (Intrinsic Motivation)
# scale_06 ~ scale_10: 自我效能 (Self-Efficacy)
# scale_03 與 scale_08 為反向題 (Negative items)
affective_mat <- matrix(0, nrow = n_students, ncol = 10)
for (j in 1:10) {
  loading <- runif(1, 0.65, 0.85)
  latent_score <- if (j <= 5) eta else (0.6 * eta + 0.4 * theta)
  cont_val <- loading * latent_score + rnorm(n_students, mean = 0, sd = sqrt(1 - loading^2))
  
  # 切分門檻轉為 1~5 點
  cuts <- quantile(cont_val, probs = c(0.15, 0.35, 0.65, 0.85))
  score <- findInterval(cont_val, cuts) + 1
  
  # 若為反向題 (第 3 與第 8 題)，進行反向轉換 (1->5, 2->4, 3->3, 4->2, 5->1)
  if (j %in% c(3, 8)) {
    score <- 6 - score
  }
  affective_mat[, j] <- score
}
colnames(affective_mat) <- sprintf("scale_%02d", 1:10)

# 4. 產生作答歷程變項：反應時間 (Response Times, 秒)
# 對數常態分佈，能力較高者作答時間平均稍短
rt_mat <- matrix(0, nrow = n_students, ncol = 15)
for (j in 1:15) {
  base_time <- 25 + abs(b_params[j]) * 10
  log_mu <- log(base_time) - 0.15 * theta
  rt_mat[, j] <- round(rlnorm(n_students, meanlog = log_mu, sdlog = 0.35), 1)
}
colnames(rt_mat) <- sprintf("rt_%02d", 1:15)
total_time <- rowSums(rt_mat)

# 5. 學業落後預警標籤 (At Risk): 總答對題數 <= 5 題標記為 1 (高風險)
total_score <- rowSums(cognitive_mat)
at_risk <- ifelse(total_score <= 5, 1, 0)

# 6. 整合為現代標準 tibble
assessment_df <- tibble(
  student_id = sprintf("STU_%04d", 1:n_students),
  gender = gender,
  school_type = school_type,
  ses_level = ses_level,
  latent_ability = round(theta, 3),
  latent_motivation = round(eta, 3),
  total_score = total_score,
  total_time = round(total_time, 1),
  at_risk = at_risk
) |>
  bind_cols(as_tibble(cognitive_mat)) |>
  bind_cols(as_tibble(affective_mat)) |>
  bind_cols(as_tibble(rt_mat))

# 7. 隨機引入 1.5% 缺失值 (NA) 以模擬真實測驗中的跳題漏答情況
# 僅在 item_11 ~ item_14 注入少量 NA
for (col_idx in sprintf("item_%02d", 11:14)) {
  na_mask <- runif(n_students) < 0.015
  assessment_df[[col_idx]][na_mask] <- NA_real_
}

# 儲存為 CSV
write.csv(assessment_df, "book/data/assessment_data.csv", row.names = FALSE)
cat(sprintf("成功生成模擬資料集！樣本數：%d，欄位數：%d\n", nrow(assessment_df), ncol(assessment_df)))
cat("檔案儲存至: book/data/assessment_data.csv\n")
