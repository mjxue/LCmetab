library(ggrepel)
library(ggplot2)
library(dplyr)

prepare_plot_data <- function(df) {
  df %>%
    # 1. 计算 Beta 值 (log HR)
    mutate(
      Beta_train = log(HR_train),
      Beta_valid = log(HR_valid),
      
      # 2. 根据你的三大条件进行分组
      Significance = case_when(
        fdr_train <= 0.05 & P_Value_valid <= 0.05 & ((HR_train>1&HR_valid>1)|(HR_train<1&HR_valid<1))~ "Validated",
        fdr_train <= 0.05 ~ "Discovery only",
        TRUE ~ "Non-significant"
      )
    ) %>%
    # 3. 将 Significance 转化为因子，并设定水平顺序
    # 这一步极其重要！它确保在画图时，灰色的不显著点画在最底层，最显著的验证点画在最顶层，防止被遮挡。
    mutate(Significance = factor(Significance, levels = c(
      "Non-significant",
      "Discovery only",
      "Validated"
    ))) %>%
    arrange(Significance)
}

# =====================================================================
# 2. 核心绘图函数
# =====================================================================

create_beta_plot <- function(df_plot, title_text, label_col = "Index") {
  
  # 计算 X 和 Y 轴的统一范围
  lims <- range(c(df_plot$Beta_train, df_plot$Beta_valid), na.rm = TRUE)
  
  p <- ggplot(df_plot, aes(x = Beta_train, y = Beta_valid, color = Significance)) +
    geom_point(alpha = 0.7, size = 1.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
    
    # 颜色映射
    scale_color_manual(values = c(
      "Non-significant" = "#E0E0E0",
      "Discovery only" = "#4DBBD5",
      "Validated" = "#DC0000"
    )) +
    
    # === 新增：添加 ggrepel 文本图层 ===
    # data 参数中只筛选出 Significance == "Validated" 的行进行标记
    # aes 中动态使用传入的 label_col 列名
    # geom_text_repel(
    #   data = df_plot %>% filter(Significance == "Validated"),
    #   aes(label = .data[[label_col]]), 
    #   size = 3.5,                 # 字体大小
    #   color = "black",            # 标签统一用黑色，或者改为 "#DC0000" 和点同色
    #   box.padding = 1,          # 标签与标签/点之间的避让空间
    #   point.padding = 0.3,        # 标签与点之间的距离
    #   segment.color = "grey50",   # 连接线的颜色
    #   max.overlaps = 20,          # 允许的最大重叠计算次数（如果点太多，设大一点以显示所有标签）
    #   show.legend = FALSE         # 不在图例中显示文本图层
    # ) +
    # ==================================
  
  labs(
    title = title_text,
    x = expression(paste("Beta value in discovery set")),
    y = expression(paste("Beta value in validation set")),
    color = NULL
  ) +
    scale_x_continuous(limits = lims) +
    scale_y_continuous(limits = lims) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = c(0.5,0.15),
      legend.direction = "horizontal",
      legend.background = element_rect(fill = NULL),
      legend.title = element_text(size = 9),
      legend.text =  element_text(size = 9),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8),   # x 轴标题字体大小
      axis.text.y = element_text(size = 8),   # y 轴标题字体大小
      axis.title.x = element_text(size = 10),   # x 轴标题字体大小
      axis.title.y = element_text(size = 10),   # y 轴标题字体大小
    )
  
  return(p)
}