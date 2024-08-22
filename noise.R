noise = residential_greenness[,c(1,28:32)]
summary(noise1)

colnames(noise) = c('eid', '16h', '24h', 'daytime', 'evening', 'nighttime')
apply(noise,2, pMiss)
table(noise1$education)
noise = drop_na(noise, 4)
sum(noise$followup_time_pso < 0)
noise = merge(noise, Covariates1, by = 'eid')
count()
noise1 = subset(noise, followup_time_pso >=0)
colnames(noise1)
###########persd######
noise1 = noise1[,-c(19:23)]
z_score <- function(col) {
  (col - mean(col, na.rm = TRUE)) / sd(col, na.rm = TRUE)
}

# 批量应用此函数计算z得分
perSD_columns <- as.data.frame(lapply(noise1 %>% select(25), z_score))

# 将新列名设置为原列名后加上_perSD
colnames(perSD_columns) <- paste0(colnames(perSD_columns), "_perSD")

# 将新列添加到原始数据框中
noise1 <- cbind(noise1, perSD_columns)
colnames(noise1)
library(gtsummary)
?tbl_summary
table1_noise <- noise1 %>%
  select(c(18,17,10,13,14,12,16,11,15,9,7,3,6,34,35,36)) %>%
  tbl_summary(
    by = source_of_pso,
    statistic = list(
      all_continuous() ~ '{median} ({p25}, {p75})',  # 默认统计方法
      Age_at_recruitment ~ '{mean} ({sd})'  # 特定于第18列的统计方法
    ),
    digits = list(
      all_continuous() ~ 2,
      Age_at_recruitment ~ 2  # 确保第18列数字显示两位小数
    )
  ) %>%
  add_p() %>%
  add_overall() %>%
  add_stat_label() %>%
  modify_header(label = "**Variable**") %>%
  bold_labels()

# 输出表格到Word文档
table1_noise %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = './table/table1_baseline1.docx')



#########四分位##########
colnames(noise1)[c(2:3)] = c('X16h', 'X24h')
noise1 = noise1 %>%
  mutate(
    X24h_C = case_when(
      X24h <= 55 ~ 'low',
      X24h > 55 & X24h <= 60 ~ 'Low_medium',
      X24h > 60 & X24h <=65 ~ 'medium_high',
      X24h > 65 ~ 'high'
    )
  )
table(noise1$X24h_C)
noise1$X24h_C = factor(noise1$X24h_C, levels = c('low', 'Low_medium', 'medium_high', 'high'))
noise1 = noise1 %>%
  mutate(
    nighttime_C = case_when(
      nighttime <= 45 ~ 'low',
      nighttime > 45 & nighttime <= 50 ~ 'Low_medium',
      nighttime > 50 & nighttime <=55 ~ 'medium_high',
      nighttime > 55 ~ 'high'
    )
  )

table(noise1$nighttime_C)
noise1$nighttime_C = factor(noise1$nighttime_C, levels = c('low', 'Low_medium', 'medium_high', 'high'))
###########P_trend########
colnames(noise1)
summary(noise1)
noise1 <- noise1 %>%
  mutate(across(all_of(27:28), 
                ~recode(.,'low' = 1, 'Low_medium' = 2, 'medium_high' = 3, 'high' = 4), 
                .names = "numeric_{.col}"))


library(dplyr)

results_list <- list()

# 迭代第34到52列进行多变量Cox回归
for (i in 32:33) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1)
  
  # 提取第i列的系数、标准误差、置信区间
  coef <- coef(cox_model)[colnames(noise1)[i]]
  se <- sqrt(diag(vcov(cox_model)))[colnames(noise1)[i]]
  lower <- coef - 1.96 * se
  upper <- coef + 1.96 * se
  
  # 提取p值
  p_value <- summary(cox_model)$coef[colnames(noise1)[i], "Pr(>|z|)"]
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR = exp(coef),
    Lower = exp(lower),
    Upper = exp(upper),
    PValue = p_value
  )
}
# 将结果汇总到一个数据框中
results_df_p_trend <- do.call(rbind, results_list)
write.csv(results_df_p_trend, './unicox/p_trend_mulcox_snp.csv')
#############unicox#############
results_list <- list()
noise1 = as.data.frame(noise1)
noise1$time = noise1$followup_time_pso/365
# 迭代第34到52列进行Cox回归
for (i in c(19:23,27:28)) {
  # 单因素Cox回归模型
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i])
  cox_model <- coxph(as.formula(formula_str), data = noise1)
  
  # 提取系数、标准误差、置信区间
  coef <- coef(cox_model)
  se <- sqrt(diag(vcov(cox_model)))
  lower <- coef - 1.96 * se
  upper <- coef + 1.96 * se
  # 提取p值
  p_value <- summary(cox_model)$coef[,'Pr(>|z|)']
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR = exp(coef),
    Lower = exp(lower),
    Upper = exp(upper),
    PValue = p_value
  )
}

# 将结果汇总到一个数据框中
results_df_unicox_noise <- do.call(rbind, results_list)
results_df_unicox_noise$adjusted_pvalue <- p.adjust(results_df_unicox_noise$PValue, method = "BH")
results_df_unicox_noise <- results_df_unicox_noise %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", Lower),
    Upper = sprintf("%0.2f", Upper),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_unicox_noise <- results_df_unicox_noise %>%
  mutate(
    PValue = sprintf("%0.3f", PValue)
  )
write.csv(results_df_unicox_noise, './unicox/results_df_unicox_noise.csv')
#########model1##########
results_list <- list()
summary(noise1)
str(noise1)
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(17,18)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(19:23,27:28)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1)
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}

# 将结果汇总到一个数据框中
results_df_multi_little_noise <- do.call(rbind, results_list)
results_df_multi_little_noise$Variable = row.names(results_df_multi_little_noise)
results_df_multi_little_noise$Variable <- gsub(".*\\.", "", results_df_multi_little_noise$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_little_noise <- results_df_multi_little_noise %>%
  filter(grepl(pattern, Variable))

results_df_multi_little_noise$adjusted_pvalue <- p.adjust(results_df_multi_little_noise$pvalue, method = "BH")
write.csv(results_df_multi_little_noise, './noise/mulcox/multicox_little.csv')
results_df_multi_little_noise <- results_df_multi_little_noise %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_little_noise <- results_df_multi_little_noise %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_multi_little_noise, './mulcox/results_df_multi_little_noise.csv')
############model2###### ====
colnames(noise1)
results_list <- list()
setwd('./noise/')
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:18)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(19:23,27:28)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1, subset = c(time >= 2))
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_multi_noise <- do.call(rbind, results_list)
results_df_multi_noise$Variable = row.names(results_df_multi_noise)
results_df_multi_noise$Variable <- gsub(".*\\.", "", results_df_multi_noise$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_noise <- results_df_multi_noise %>%
  filter(grepl(pattern, Variable))
results_df_multi_noise$adjusted_pvalue <- p.adjust(results_df_multi_noise$pvalue, method = "BH")
results_df_multi_noise <- results_df_multi_noise %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_noise <- results_df_multi_noise %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_multi_noise, './mulcox/results_df_multi_noise.csv')

############model3###### ====
colnames(noise1)
air_cov = residential_greenness1[,c(1,11,15)]
noise1 = merge(noise1, air_cov, by = 'eid')
results_list <- list()
setwd('./noise/')
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:18,34)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(19:23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1)
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_model3 <- do.call(rbind, results_list)
results_model3$Variable = row.names(results_model3)
results_model3$Variable <- gsub(".*\\.", "", results_model3$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_model3 <- results_model3 %>%
  filter(grepl(pattern, Variable))
results_model3$adjusted_pvalue <- p.adjust(results_model3$pvalue, method = "BH")
results_model3 <- results_model3 %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_model3 <- results_model3 %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_model3, './mulcox/results_model3.csv')
############model4###### ====
colnames(noise1)
results_list <- list()
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:18,35)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(19:23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1)
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_model4 <- do.call(rbind, results_list)
results_model4$Variable = row.names(results_model4)
results_model4$Variable <- gsub(".*\\.", "", results_model4$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_model4 <- results_model4 %>%
  filter(grepl(pattern, Variable))
results_model4$adjusted_pvalue <- p.adjust(results_model4$pvalue, method = "BH")
results_model4 <- results_model4 %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_model4 <- results_model4 %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_model4, './mulcox/results_model4.csv')

############model5###### ====
colnames(noise1)
livetime1 = fread('live_time_participant.csv')
str(livetime1)
head(livetime1)
livetime = livetime[,1:3]
colnames(livetime) = c('eid', 'date_location', 'date_attend')
str(livetime)
livetime$date_location = as.IDate(livetime$date_location)
livetime$live_time = livetime$date_attend - livetime$date_location
livetime$live_time = livetime$live_time/365 
livetime = livetime[,c(1,4)]
noise1 = merge(noise1, livetime, by = 'eid')
summary(noise1)
# Replace negative values in the 'live_time' column with 0
noise1$live_time <- ifelse(noise1$live_time < 0, 0, noise1$live_time)
library(mice)
miss_data = noise1[,c(18,35:36)]
colnames(miss_data)
str(miss_data)
imputed_data <- mice(miss_data, method = c('pmm','pmm','pmm'), m = 1)
completed_data <- complete(imputed_data, 1)  # 获取第一个插补数据集
pMiss = function(x){sum(is.na(x))/length(x)*100}
apply(miss_data,2, pMiss)
noise1 = noise1[,-c(35:36)]
noise1 = cbind(noise1, completed_data[,2:3])



results_list <- list()
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:18,36)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(19:23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1)
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_model5 <- do.call(rbind, results_list)
results_model5$Variable = row.names(results_model5)
results_model5$Variable <- gsub(".*\\.", "", results_model5$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_model5 <- results_model5 %>%
  filter(grepl(pattern, Variable))
results_model5$adjusted_pvalue <- p.adjust(results_model5$pvalue, method = "BH")
results_model5 <- results_model5 %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_model5 <- results_model5 %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_model5, './mulcox/results_model5.csv')

#########sensitive#########
colnames(noise1)
results_list <- list()
setwd('./noise/')
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:18)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(19:23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1, subset = c(time >= 2))
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_multi_noise_2y <- do.call(rbind, results_list)
results_df_multi_noise_2y$Variable = row.names(results_df_multi_noise_2y)
results_df_multi_noise_2y$Variable <- gsub(".*\\.", "", results_df_multi_noise_2y$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_noise_2y <- results_df_multi_noise_2y %>%
  filter(grepl(pattern, Variable))
results_df_multi_noise_2y$adjusted_pvalue <- p.adjust(results_df_multi_noise_2y$pvalue, method = "BH")
results_df_multi_noise_2y <- results_df_multi_noise_2y %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_noise_2y <- results_df_multi_noise_2y %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_multi_noise_2y, './mulcox/results_df_multi_noise_2y.csv')



results_list <- list()
setwd('./noise/')
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:18)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(19:23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1, subset = c(time >= 2.5))
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_multi_noise_5y <- do.call(rbind, results_list)
results_df_multi_noise_5y$Variable = row.names(results_df_multi_noise_5y)
results_df_multi_noise_5y$Variable <- gsub(".*\\.", "", results_df_multi_noise_5y$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_noise_5y <- results_df_multi_noise_5y %>%
  filter(grepl(pattern, Variable))
results_df_multi_noise_5y$adjusted_pvalue <- p.adjust(results_df_multi_noise_5y$pvalue, method = "BH")
results_df_multi_noise_5y <- results_df_multi_noise_5y %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_noise_5y <- results_df_multi_noise_5y %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_multi_noise_5y, './mulcox/results_df_multi_noise_5y.csv')
############sensitive_more#############
hearing = fread('noise/hearing_loss_participant.csv')
head(hearing)
colnames(hearing) = c('eid', 'date_hl', 'source_hl', 'date_ohl', 'source_ohl', 'ukb_centre')
table(hearing$source_hl)
table(hearing$source_ohl)
table(hearing$ukb_centre)
hearing = hearing[,c(1,3,5,6)]
hearing <- hearing %>%
  mutate(
    across(.cols = 2:3, .fns = ~if_else(. == "" | is.na(.), "no", "yes"))
  )
hearing <- hearing %>%
  mutate(
    hearing_loss = if_else(source_hl == "yes" | source_ohl == "yes", "yes", "no")
  )

hearing$hearing_loss = factor(hearing$hearing_loss, levels = c('no', 'yes'))
table(hearing$hearing_loss)
hearing1 = hearing[,c(1,4,5)]

summary(hearing1)
noise1 = merge(noise1, hearing1, by = 'eid')
library(survival)
results_list <- list()
setwd('./noise/')
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:18,38)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(19:23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1)
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_multi_noise_hearing_loss <- do.call(rbind, results_list)
results_df_multi_noise_hearing_loss$Variable = row.names(results_df_multi_noise_hearing_loss)
results_df_multi_noise_hearing_loss$Variable <- gsub(".*\\.", "", results_df_multi_noise_hearing_loss$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_noise_hearing_loss <- results_df_multi_noise_hearing_loss %>%
  filter(grepl(pattern, Variable))
results_df_multi_noise_hearing_loss$adjusted_pvalue <- p.adjust(results_df_multi_noise_hearing_loss$pvalue, method = "BH")
results_df_multi_noise_hearing_loss <- results_df_multi_noise_hearing_loss %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_noise_hearing_loss <- results_df_multi_noise_hearing_loss %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_multi_noise_hearing_loss, './mulcox/results_df_multi_noise_hearing_loss.csv')



results_list <- list()
setwd('./noise/')
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:18,37)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(19:23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1)
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_multi_noise_center <- do.call(rbind, results_list)
results_df_multi_noise_center$Variable = row.names(results_df_multi_noise_center)
results_df_multi_noise_center$Variable <- gsub(".*\\.", "", results_df_multi_noise_center$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_noise_center <- results_df_multi_noise_center %>%
  filter(grepl(pattern, Variable))
results_df_multi_noise_center$adjusted_pvalue <- p.adjust(results_df_multi_noise_center$pvalue, method = "BH")
results_df_multi_noise_center <- results_df_multi_noise_center %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_noise_center <- results_df_multi_noise_center %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_multi_noise_center, './mulcox/results_df_multi_noise_center.csv')


###########interaction########
colnames(noise1)
library(survival)
library(broom)
# 定义自变量和协变量的列索引
independent_vars <- c(20,23,27,28)
covariate_vars <- c(9:18,26,34,35,36)

# 初始化数据框存储所有模型结果
all_results_df <- data.frame()

# 生成包含所有协变量的字符串（除了当前正在分析的那一个）
for (ind in independent_vars) {
  for (cov in covariate_vars) {
    # 除了当前分析的协变量外，获取其他所有协变量
    other_covariates <- covariate_vars[covariate_vars != cov]
    other_covs_formula <- paste(colnames(noise1)[other_covariates], collapse=" + ")
    
    # 构建包含交互作用和其他协变量的模型公式
    formula_str <- as.formula(paste("Surv(time, source_of_pso) ~", 
                                    paste(colnames(noise1)[ind], "*", colnames(noise1)[cov], sep=""),
                                    "+", other_covs_formula))
    
    # 拟合Cox回归模型
    cox_model <- coxph(formula_str, data = noise1)
    
    # 使用broom包的tidy()函数提取模型摘要
    model_summary <- tidy(cox_model)
    
    # 为结果数据框添加标识列
    model_summary$Model <- paste(colnames(noise1)[ind], "x", colnames(noise1)[cov], sep="")
    
    # 将当前模型结果追加到总结果数据框中
    all_results_df <- bind_rows(all_results_df, model_summary)
  }
}


all_results_df = all_results_df %>%
  filter(grepl(":", term))
write.csv(all_results_df, './table/interaction_cov.csv')
############PRS#############
pso_PRS = fread('PRS/merged_dataset.best')
head(pso_PRS)
noise1 = merge(noise1, pso_PRS, by.x = 'eid', by.y = 'FID')
summary(noise1)
noise1 = noise1[,-c(25:26)]
colnames(noise1)
noise1 = noise1 %>%
  mutate(PRS_C = 
    case_when(
      PRS <= median(PRS) ~ 'low',
      PRS > median(PRS) ~ 'high'
    )
  )
table(noise1$PRS_C)
noise1$PRS_C = factor(noise1$PRS_C, levels = c('low', 'high'))

results_list <- list()

# 迭代第34到52列进行多变量Cox回归
for (i in c(26,29)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i])
  cox_model <- coxph(as.formula(formula_str), data = noise1)
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_PRS_unicox <- do.call(rbind, results_list)
results_df_PRS_unicox$Variable = row.names(results_df_PRS_unicox)
results_df_PRS_unicox$Variable <- gsub(".*\\.", "", results_df_PRS_unicox$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(26,29)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_PRS_unicox <- results_df_PRS_unicox %>%
  filter(grepl(pattern, Variable))
results_df_PRS_unicox$adjusted_pvalue <- p.adjust(results_df_PRS_unicox$pvalue, method = "BH")
results_df_PRS_unicox <- results_df_PRS_unicox %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_PRS_unicox <- results_df_PRS_unicox %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_PRS_unicox, './mulcox/results_df_PRS_unicox.csv')

results_list <- list()

# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:18)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(26,29)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1)
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_PRS <- do.call(rbind, results_list)
results_df_PRS$Variable = row.names(results_df_PRS)
results_df_PRS$Variable <- gsub(".*\\.", "", results_df_PRS$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(26,29)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_PRS <- results_df_PRS %>%
  filter(grepl(pattern, Variable))
results_df_PRS$adjusted_pvalue <- p.adjust(results_df_PRS$pvalue, method = "BH")
results_df_PRS <- results_df_PRS %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_PRS <- results_df_PRS %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_PRS, './mulcox/results_df_PRS.csv')
#############JOINT##########
noise1 <- noise1 %>%
  mutate(
    X24h_C_PRS_C = paste(X24h_C, PRS_C, sep = "_"),
    nighttime_C_PRS_C = paste(nighttime_C, PRS_C, sep = "_")
  )

levels_combinations <- c("low_low",  "Low_medium_low", 
                         "medium_high_low",  "high_low","low_high","Low_medium_high","medium_high_high", "high_high")
library(dplyr)

noise1 <- noise1 %>%
  mutate(
    X24h_C_PRS_C = factor(X24h_C_PRS_C, levels = levels_combinations),
    nighttime_C_PRS_C = factor(nighttime_C_PRS_C, levels = levels_combinations)
  )
colnames(noise1)
results_list <- list()

# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:18)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(30:31)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1)
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_injoint <- do.call(rbind, results_list)
results_df_injoint$Variable = row.names(results_df_injoint)
results_df_injoint$Variable <- gsub(".*\\.", "", results_df_injoint$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(30:31)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_injoint <- results_df_injoint %>%
  filter(grepl(pattern, Variable))

write.csv(results_df_injoint, 'results_df_injoint.csv')


#########sex#########
colnames(noise1)
results_list <- list()
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:16,18)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(19:23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1, subset = c(sex == 'male'))
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_multi_noise_male <- do.call(rbind, results_list)
results_df_multi_noise_male$Variable = row.names(results_df_multi_noise_male)
results_df_multi_noise_male$Variable <- gsub(".*\\.", "", results_df_multi_noise_male$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_noise_male <- results_df_multi_noise_male %>%
  filter(grepl(pattern, Variable))
results_df_multi_noise_male$adjusted_pvalue <- p.adjust(results_df_multi_noise_male$pvalue, method = "BH")
results_df_multi_noise_male <- results_df_multi_noise_male %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_noise_male <- results_df_multi_noise_male %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_multi_noise_male, './mulcox/results_df_multi_noise_male.csv')



results_list <- list()
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:16,18)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(19:23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1, subset = c(sex=='female'))
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_multi_noise_female <- do.call(rbind, results_list)
results_df_multi_noise_female$Variable = row.names(results_df_multi_noise_female)
results_df_multi_noise_female$Variable <- gsub(".*\\.", "", results_df_multi_noise_female$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_noise_female <- results_df_multi_noise_female %>%
  filter(grepl(pattern, Variable))
results_df_multi_noise_female$adjusted_pvalue <- p.adjust(results_df_multi_noise_female$pvalue, method = "BH")
results_df_multi_noise_female <- results_df_multi_noise_female %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_noise_female <- results_df_multi_noise_female %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_multi_noise_female, './mulcox/results_df_multi_noise_female.csv')

#########age#########
colnames(noise1)
results_list <- list()
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:17)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(19:23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1, subset = c(Age_at_recruitment < 60))
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_multi_noise_young <- do.call(rbind, results_list)
results_df_multi_noise_young$Variable = row.names(results_df_multi_noise_young)
results_df_multi_noise_young$Variable <- gsub(".*\\.", "", results_df_multi_noise_young$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_noise_young <- results_df_multi_noise_young %>%
  filter(grepl(pattern, Variable))
results_df_multi_noise_young$adjusted_pvalue <- p.adjust(results_df_multi_noise_young$pvalue, method = "BH")
results_df_multi_noise_young <- results_df_multi_noise_young %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_noise_young <- results_df_multi_noise_young %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_multi_noise_young, './mulcox/results_df_multi_noise_young.csv')



colnames(noise1)
results_list <- list()
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:17)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(19:23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1, subset = c(Age_at_recruitment >= 60))
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_multi_noise_old <- do.call(rbind, results_list)
results_df_multi_noise_old$Variable = row.names(results_df_multi_noise_old)
results_df_multi_noise_old$Variable <- gsub(".*\\.", "", results_df_multi_noise_old$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_noise_old <- results_df_multi_noise_old %>%
  filter(grepl(pattern, Variable))
results_df_multi_noise_old$adjusted_pvalue <- p.adjust(results_df_multi_noise_old$pvalue, method = "BH")
results_df_multi_noise_old <- results_df_multi_noise_old %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_noise_old <- results_df_multi_noise_old %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_multi_noise_old, './mulcox/results_df_multi_noise_old.csv')

#########PRS#########
colnames(noise1)
table(noise1$PRS_C)
library(survival)
results_list <- list()
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:18)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(20,23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1, subset = c(PRS_C == 'low'))
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_multi_noise_low_prs <- do.call(rbind, results_list)
results_df_multi_noise_low_prs$Variable = row.names(results_df_multi_noise_low_prs)
results_df_multi_noise_low_prs$Variable <- gsub(".*\\.", "", results_df_multi_noise_low_prs$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(20,23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_noise_low_prs <- results_df_multi_noise_low_prs %>%
  filter(grepl(pattern, Variable))
results_df_multi_noise_low_prs$adjusted_pvalue <- p.adjust(results_df_multi_noise_low_prs$pvalue, method = "BH")
results_df_multi_noise_low_prs <- results_df_multi_noise_low_prs %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_noise_low_prs <- results_df_multi_noise_low_prs %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
setwd('./noise/')
write.csv(results_df_multi_noise_low_prs, './mulcox/results_df_multi_noise_low_prs.csv')



colnames(noise1)
results_list <- list()
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:18)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(20,23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1, subset = c(PRS_C == 'high'))
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_multi_noise_high_prs <- do.call(rbind, results_list)
results_df_multi_noise_high_prs$Variable = row.names(results_df_multi_noise_high_prs)
results_df_multi_noise_high_prs$Variable <- gsub(".*\\.", "", results_df_multi_noise_high_prs$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_noise_high_prs <- results_df_multi_noise_high_prs %>%
  filter(grepl(pattern, Variable))
results_df_multi_noise_high_prs$adjusted_pvalue <- p.adjust(results_df_multi_noise_high_prs$pvalue, method = "BH")
results_df_multi_noise_high_prs <- results_df_multi_noise_high_prs %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_noise_high_prs <- results_df_multi_noise_high_prs %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_multi_noise_high_prs, './mulcox/results_df_multi_noise_high_prs.csv')

#########townsend#########
colnames(noise1)
median(noise1$Townsend)
results_list <- list()
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:13,15:18)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(20,23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1, subset = c(Townsend < -2.17))
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_multi_noise_low_town <- do.call(rbind, results_list)
results_df_multi_noise_low_town$Variable = row.names(results_df_multi_noise_low_town)
results_df_multi_noise_low_town$Variable <- gsub(".*\\.", "", results_df_multi_noise_low_town$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(20,23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_noise_low_town <- results_df_multi_noise_low_town %>%
  filter(grepl(pattern, Variable))
results_df_multi_noise_low_town$adjusted_pvalue <- p.adjust(results_df_multi_noise_low_town$pvalue, method = "BH")
results_df_multi_noise_low_town <- results_df_multi_noise_low_town %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_noise_low_town <- results_df_multi_noise_low_town %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_multi_noise_low_town, './mulcox/results_df_multi_noise_low_town.csv')



colnames(noise1)
results_list <- list()
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:13,15:18)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(20,23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1, subset = c(Townsend >= -2.17))
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_multi_noise_high_town <- do.call(rbind, results_list)
results_df_multi_noise_high_town$Variable = row.names(results_df_multi_noise_high_town)
results_df_multi_noise_high_town$Variable <- gsub(".*\\.", "", results_df_multi_noise_high_town$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(20,23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_noise_high_town <- results_df_multi_noise_high_town %>%
  filter(grepl(pattern, Variable))
results_df_multi_noise_high_town$adjusted_pvalue <- p.adjust(results_df_multi_noise_high_town$pvalue, method = "BH")
results_df_multi_noise_high_town <- results_df_multi_noise_high_town %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_noise_high_town <- results_df_multi_noise_high_town %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_multi_noise_high_town, './mulcox/results_df_multi_noise_high_town.csv')

#########income#########
colnames(noise1)
results_list <- list()
# 创建包含所有协变量的字符串
covariates <- paste(colnames(noise1)[c(9:10,12:18)], collapse=" + ")

# 迭代第34到52列进行多变量Cox回归
for (i in c(19:23,27:28,32:33)) {
  # 多变量Cox回归模型，包含协变量
  formula_str <- paste("Surv(time, source_of_pso) ~", colnames(noise1)[i], "+", covariates)
  cox_model <- coxph(as.formula(formula_str), data = noise1, subset = c(income == 'Greater_than_100,000'))
  
  # 提取第i列的系数、标准误差、置信区间
  a = summary(cox_model)
  
  # 存储结果
  results_list[[colnames(noise1)[i]]] <- data.frame(
    Variable = colnames(noise1)[i],
    HR=a$conf.int[,"exp(coef)"],
    L95CI=a$conf.int[,"lower .95"],
    H95CI=a$conf.int[,"upper .95"],
    pvalue=a$coefficients[,"Pr(>|z|)"])
}


# 将结果汇总到一个数据框中
results_df_multi_noise_greater_100000 <- do.call(rbind, results_list)
results_df_multi_noise_greater_100000$Variable = row.names(results_df_multi_noise_greater_100000)
results_df_multi_noise_greater_100000$Variable <- gsub(".*\\.", "", results_df_multi_noise_greater_100000$Variable)
# 定义要筛选的字符
characters_to_include <- colnames(noise1)[c(19:23,27:28,32:33)]

# 创建正则表达式，匹配任何包含上述字符的字符串
pattern <- paste(characters_to_include, collapse = "|")

# 使用dplyr的filter和grepl函数筛选包含特定字符的行
results_df_multi_noise_greater_100000 <- results_df_multi_noise_greater_100000 %>%
  filter(grepl(pattern, Variable))
results_df_multi_noise_greater_100000$adjusted_pvalue <- p.adjust(results_df_multi_noise_greater_100000$pvalue, method = "BH")
results_df_multi_noise_greater_100000 <- results_df_multi_noise_greater_100000 %>%
  mutate(
    HR = sprintf("%0.2f", HR),
    Lower = sprintf("%0.2f", L95CI),
    Upper = sprintf("%0.2f", H95CI),
    ci = paste(HR, " [", Lower, ";", Upper, "]", sep = "")
  )
results_df_multi_noise_greater_100000 <- results_df_multi_noise_greater_100000 %>%
  mutate(
    PValue = sprintf("%0.3f", pvalue)
  )
write.csv(results_df_multi_noise_greater_100000, './mulcox/results_df_multi_noise_greater_100000.csv')



