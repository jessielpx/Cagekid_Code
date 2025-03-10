#### Oncoplot + CNV ####
# 加载必要包 ----------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(readxl)

# 数据预处理 ----------------------------------------------------------------
# 1. 读取临床数据并筛选 VHL = 0 的患者 -----------------------------------
cagekid_clin <- read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin') %>% filter(VHL == 0)

vhl0_patients <- unique(cagekid_clin$PatientID)

# 2. 读取 CNA 数据并处理 -----------------------------------------------
cna_data <- read_excel('/Users/jessie/Desktop/Cagekid/merged_seg.xlsx') %>%
  filter(PatientID %in% vhl0_patients) %>%
  mutate(
    chrom = paste0("chr", chrom),
    Corrected_Copy_Number = ifelse(is.infinite(Corrected_Copy_Number), 2, Corrected_Copy_Number) # 处理 Inf 值
  )

# 3. 读取突变数据 ------------------------------------------------------
mutation_data <- read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_mut') %>%
  filter(PatientID %in% vhl0_patients)

# 4. 创建染色体着丝粒位置数据框（hg19，使用中间值） -------------------
centromeres <- tibble(
  chrom = paste0("chr", 1:22),
  centromere_mid = c(132000000, 93900000, 91000000, 50000000, 48750000, 60550000,
                     61150000, 45500000, 49450000, 39800000, 53400000, 35500000,
                     17550000, 17150000, 19000000, 36950000, 24000000, 17200000,
                     26150000, 28050000, 11950000, 15550000)
)

# 5. 处理 CNA 数据（按 p/q 臂分类） ------------------------------------
split_arm_segments <- function(data) {
  data %>%
    mutate(
      centromere_mid = centromeres$centromere_mid[match(chrom, centromeres$chrom)],
      p_length = pmax(0, pmin(end, centromere_mid) - pmax(start, 0)),
      q_length = pmax(0, pmin(end, 1e9) - pmax(start, centromere_mid))
    ) %>%
    pivot_longer(cols = c(p_length, q_length), names_to = "arm_type", values_to = "arm_length") %>%
    filter(arm_length > 0) %>%
    mutate(
      arm = ifelse(arm_type == "p_length", "p", "q"),
      chrom_arm = paste0(chrom, arm),
      weighted_logR = Corrected_Copy_Number * arm_length
    )
}

# 6. 计算 CNA 加权平均值 ------------------------------------------------
cna_data <- cna_data %>%
  mutate(
    Corrected_Copy_Number = as.numeric(Corrected_Copy_Number),  # Ensure numeric type
    num.mark = as.numeric(num.mark)  # Ensure num.mark is numeric
  )

cna_agg <- cna_data %>%
  split_arm_segments() %>%
  group_by(chrom_arm, PatientID) %>%
  summarise(
    arm_logR = {
      total_weight = sum(num.mark * arm_length, na.rm = TRUE)  # Use num.mark * arm_length as weight
      if(total_weight == 0) {
        2  # When total weight is zero, assume diploid
      } else {
        sum(weighted_logR * num.mark, na.rm = TRUE) / total_weight  # Weighted by num.mark & arm_length
      }
    },
    .groups = "drop"
  )


# 7. 构建 CNA 矩阵 ------------------------------------------------------
chrom_order <- paste0("chr", rep(1:22, each = 2), c("p", "q"))

cna_matrix <- cna_agg %>%
  pivot_wider(
    names_from = PatientID, 
    values_from = arm_logR, 
    values_fill = 2  # 确保缺失值填充为2（二倍体）
  ) %>%
  column_to_rownames("chrom_arm") %>%
  as.matrix() %>%
  .[intersect(chrom_order, rownames(.)), ] %>%
  replace(., is.na(.) | is.infinite(.), 2)  # 二次替换异常值

if ("chr3p" %in% rownames(cna_matrix)) {
  sample_order <- names(sort(cna_matrix["chr3p", ], decreasing = TRUE))  # 降序排列
  cna_matrix <- cna_matrix[, sample_order, drop = FALSE]  # 重新排序CNA矩阵
}


cnv_breaks <- c(0, 1, 2, 3, 5, 7)
cnv_labels <- c("Homo DEL(0)", "Hemi DEL(1)", "Diploid(2)", "Low AMP(3)", "High AMP(5)", "Focal AMP(7)")
cnv_col_fun <- colorRamp2(
  breaks = cnv_breaks,
  colors = c("#084594", "#4292c6", "white", "#fb6a4a", "#cb181d", "#99000d")
)

# 8. 处理突变矩阵 ------------------------------------------------------
mutation_colors <- c(Missense_Mutation = "#FF7F00", Frame_Shift_Del = "#4DAF4A")

# 步骤1：创建包含所有患者和基因的全组合框架
all_patients_genes <- expand_grid(
  Hugo_Symbol = unique(mutation_data$Hugo_Symbol),
  PatientID = vhl0_patients  # 确保包含所有VHL=0患者
)

# 步骤2：构建突变矩阵（保留所有患者）
# 修复后的突变矩阵构建代码
mutation_matrix <- read_xlsx(str_glue("/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx"), "Cagekid_mut") %>%
  filter(PatientID %in% colnames(cna_matrix)) %>%  
  # 创建全组合框架时显式指定列名
  right_join(
    expand_grid(
      Hugo_Symbol = unique(.$Hugo_Symbol),
      PatientID = colnames(cna_matrix)
    ),
    by = c("Hugo_Symbol", "PatientID")  # 显式指定连接键
  ) %>%  
  # 分组处理时保留必要列
  group_by(Hugo_Symbol, PatientID) %>%
  summarise(
    Mutation = if_else(
      all(is.na(Variant_Classification)), 
      "", 
      str_c(unique(na.omit(Variant_Classification)), collapse = ";")
    ),
    .groups = "drop"
  ) %>% 
  # 关键修复：明确指定数值列
  pivot_wider(
    names_from = PatientID, 
    values_from = Mutation,  # 显式指定数值来源列
    values_fill = ""
  ) %>%
  column_to_rownames("Hugo_Symbol") %>% 
  as.matrix()
mutation_matrix[mutation_matrix == "Silent"] <- ""  # Convert "Silent" to an empty string


# 9. 确保样本一致 ------------------------------------------------------
# Ensure both matrices have the same column order
sample_order <- names(sort(cna_matrix["chr3p", ], decreasing = FALSE))
common_samples <- intersect(colnames(cna_matrix), colnames(mutation_matrix))
cna_matrix <- cna_matrix[, common_samples, drop = FALSE]
mutation_matrix <- mutation_matrix[, common_samples, drop = FALSE]
cna_matrix <- cna_matrix[, sample_order, drop = FALSE]  # Apply sorting to CNA
mutation_matrix <- mutation_matrix[, sample_order, drop = FALSE]

ht_cnv <- Heatmap(
  cna_matrix, name = "CNA", col = cnv_col_fun,
  cluster_rows = FALSE,  # Allow row clustering but not column
  cluster_columns = FALSE,  # Disable column clustering (fixed 3p order)
  row_names_side = "left", column_names_side = "bottom",
  column_names_rot = 45, show_column_names = TRUE,
  heatmap_legend_param = list(
    title = "Log2 Ratio", 
    at = c(0, 1, 2, 3, 5, 7),
    labels = c("Homo DEL(0)", "Hemi DEL(1)", "Diploid(2)", 
               "Low AMP(3)", "High AMP(5)", "Focal AMP(7)")
  )
)



alter_fun <- list(
  background = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = "#F5F5F5", col = NA)),
  Missense_Mutation = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = "#FF7F00", col = NA)),
  Frame_Shift_Del = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = "#4DAF4A", col = NA))
)


ht_mut <- oncoPrint(
  mutation_matrix, 
  alter_fun = alter_fun, 
  alter_fun_is_vectorized = TRUE, 
  col = mutation_colors,
  remove_empty_columns = FALSE, 
  remove_empty_rows = FALSE,
  row_names_side = "left", 
  pct_side = "right",  
  top_annotation = NULL,
  column_order = colnames(mutation_matrix),  # Ensure sorting follows chr3p
  heatmap_legend_param = list(
    title = "Mutation Types", 
    at = names(mutation_colors),
    labels = c("Missense", "Frameshift Deletion")
  )
)

# 12. 组合并绘制图形 ----------------------------------------------------
# 组合并绘制图形
combined_ht <- ht_mut %v% ht_cnv
draw(
  combined_ht, 
  merge_legends = TRUE,
  column_title = "Integrated Genomic Profile (Sorted by chr3p CNA)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  padding = unit(c(2, 2, 12, 2), "mm")
)