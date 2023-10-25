###Author Marco Sanna (git:marco-sanna) 25/10/2023


library("HTqPCR")
library(dplyr)
library(tidyr)
library(stringi)
library("tibble")
library(ggplot2)
library(cowplot)
library(stringr)
library(ctrlGene)
library(readxl)
library(factoextra)
library(ggstar)
library(ggforce)
library(ggpubr)
library(ggsignif)
library(reshape2)
library(gridExtra)
library(grid)
library(EnhancedVolcano)
library("ggsci")
library(patchwork)
library(car)
library(ggpubr)
library(reshape2)
library(reshape)
library(pheatmap)
library(EnhancedVolcano)
library(pathfindR)
library(enrichplot)
library(doBy)

################################################################################# Functions
nth_element <- function(vector, starting_position, n) { 
  vector[seq(starting_position, length(vector), n)] }

test_group_difference <- function(df, group_col, selected_col) {
  # Load necessary packages
  library(tidyverse)
  library(car)
  df = filter(df, !is.na(!!sym(selected_col)))
  df[[group_col]] = factor(df[[group_col]])
  
  # Check if the selected column and group column exist in the data frame
  if (!(selected_col %in% names(df))) {
    stop(paste("Error: Selected column", selected_col, "does not exist in data frame."))
  }
  if (!(group_col %in% names(df))) {
    stop(paste("Error: Group column", group_col, "does not exist in data frame."))
  }
  
  # Split the data frame by the group column
  groups <- df %>% 
    split(.[[group_col]])
  
  # Check if there are only two groups
  if (length(groups) != 2) {
    stop("Error: This function is designed to test differences between two groups only.")
  }
  
  # Test for normality and equal variance of the two groups
  normality_test <- tapply(df[[selected_col]], df[[group_col]], shapiro.test)
  equal_var_test <- leveneTest(df[[selected_col]] ~ df[[group_col]])
  
  # If normality and equal variance assumptions are met, perform Welch's t-test
  if (all(sapply(normality_test, function(x) x$p.value >= 0.05)) & equal_var_test$'Pr(>F)'[1] >= 0.05) {
    t_test <- t.test(df[[selected_col]] ~ df[[group_col]], var.equal = TRUE, na.rm = T)
    return(t_test)
  }
  # If normality assumptions are not met, perform Wilcoxon rank sum test
  wilcox_test <- wilcox.test(df[[selected_col]] ~ df[[group_col]], na.rm = T, exact = FALSE)
  return(wilcox_test)
}

##########################################################################  LOAD INPUT DATA  ###############################################################################################

men_sba = read.csv("MenB_hSBA_data.csv")
df_frequencies = read.csv("MenB_frequency_data.csv")
inp_men = read.delim("Combined_Chip_Run_MeN.csv", 
                     skip = 11, 
                     sep = ",", 
                     colClasses= "character")

men_meta = read.csv("MenB_Fluidigm_meta.csv", row.names = 1)

#LOAD FLUIDIGM Ct DATA
men_ct = readCtData(files = "Combined_Chip_Run_MeN.csv", 
                    n.features = 96,
                    n.data = length(samplenames_men),
                    samples = samplenames_men,
                    format = "BioMark")

############################################################################# FIGURE 2 ########################################################################################################

############################################################     Figure 2A: COMBINED HC-HIV hSBA CUTOFFS CURVES
df = as.data.frame(men_sba)

df = df %>% 
  rename(
    T0 = SBA_T0,
    T2 = SBA_T2
  )

# Filter out rows with missing values for T0
df_filtered_t0 <- df[, c("StudyPID", "Group", "T0")] %>% 
  filter(!is.na(T0), !is.na(Group))

# Filter out rows with missing values for T2
df_filtered_t2 <- df[, c("StudyPID", "Group", "T2")] %>% 
  filter(!is.na(T2), !is.na(Group))

# Combine the filtered data frames
df_combined <- merge(df_filtered_t0, df_filtered_t2, all.x = T)
df_combined <- df_combined[complete.cases(df_combined),]

cutoffs <- c(2, 4, 8, 16, 32, 64, 128, 256, 512)

# Reshape the data into the desired format

df_long <- df_combined %>%
  pivot_longer(cols = c(T0, T2), names_to = "Value", values_to = "TempValue") %>%
  mutate(Group = paste(Group, Value, sep = "_"),
         Value = TempValue) %>%
  select(-TempValue)

freq_df <- as.data.frame(aggregate(Value ~ Group, data = df_long, FUN = function(x) {
  sapply(cutoffs, function(cutoff) {
    sum(x >= cutoff) / length(x) * 100
  })
}))

treats = freq_df$Group

freq_df = data.frame(freq_df$Value)
colnames(freq_df) = cutoffs
freq_df$Group = treats


# Reshape the data into the desired format
freq_df_long <- freq_df %>%
  gather(key = "Cutoff", value = "Value", -Group)

# Order the levels of the "Cutoff" variable
freq_df_long$Cutoff <- factor(freq_df_long$Cutoff, levels = as.character(cutoffs))

# Plot the data
plot_1 = ggplot(freq_df_long, aes(x = Cutoff, y = Value, color = Group)) +
  geom_vline(xintercept = "4", 
             linetype = "dashed", 
             color = "darkgrey", 
             size = 1) +  # Add horizontal line
  geom_line(size = 1.5, 
            aes(group = Group, linetype = Group)) +
  geom_point(size = 4, 
             aes(shape = Group)) +
  geom_rect(aes(xmin = 0, xmax = 2, ymin = -Inf, ymax = Inf), 
            color = NA, 
            fill='gainsboro', 
            alpha= 0.02)+
  labs(x = "hSBA", 
       y = "Patients (%)", 
       color = "Group") +
  scale_color_manual(values = c("HC_T0" = "blue", 
                                "PHIV_T0" = "red", 
                                "HC_T2" = "blue", 
                                "PHIV_T2" = "red")) +
  scale_linetype_manual(values = c("HC_T0" = "dashed", 
                                   "PHIV_T0" = "dashed", 
                                   "HC_T2" = "solid", 
                                   "PHIV_T2" = "solid")) +
  ggtitle("Patients (%) per hSBA titers") +
  theme_classic(base_size = 30) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

############################################################       Figure 2D: COMBINED HIV
df = as.data.frame(men_sba)
df$Treatment[df$Treatment == "early treated"] = "Early"
df$Treatment[df$Treatment == "late treated"] = "Late"

df = df %>% 
  rename(
    T0 = SBA_T0,
    T2 = SBA_T2
  )

# Filter out rows with missing values for T0
df_filtered_t0 <- df[, c("StudyPID", "Treatment", "T0")] %>% 
  filter(!is.na(T0), !is.na(Treatment))

# Filter out rows with missing values for T2
df_filtered_t2 <- df[, c("StudyPID", "Treatment", "T2")] %>% 
  filter(!is.na(T2), !is.na(Treatment))

# Combine the filtered data frames
df_combined <- merge(df_filtered_t0, df_filtered_t2, all.x = T)
df_combined <- df_combined[complete.cases(df_combined),]


cutoffs <- c(2, 4, 8, 16, 32, 64, 128, 256, 512)

# Reshape the data into the desired format

df_long <- df_combined %>%
  pivot_longer(cols = c(T0, T2), names_to = "Value", values_to = "TempValue") %>%
  mutate(Treatment = paste(Treatment, Value, sep = "_"),
         Value = TempValue) %>%
  select(-TempValue)

freq_df <- as.data.frame(aggregate(Value ~ Treatment, data = df_long, FUN = function(x) {
  sapply(cutoffs, function(cutoff) {
    sum(x >= cutoff) / length(x) * 100
  })
}))

treats = freq_df$Treatment

freq_df = data.frame(freq_df$Value)
colnames(freq_df) = cutoffs
freq_df$Treatment = treats


# Reshape the data into the desired format
freq_df_long <- freq_df %>%
  gather(key = "Cutoff", value = "Value", -Treatment)

# Order the levels of the "Cutoff" variable
freq_df_long$Cutoff <- factor(freq_df_long$Cutoff, levels = as.character(cutoffs))
 
colnames(freq_df_long)[1] = "ART start"

# Plot the data
plot_2 = ggplot(freq_df_long, aes(x = Cutoff, y = Value, color = `ART start`, linetype = `ART start`)) +
  geom_vline(xintercept = "4", 
             linetype = "dashed", 
             color = "darkgrey", 
             size = 1) +  # Add horizontal line
  geom_line(size = 1.5, 
            aes(group = `ART start`, linetype = `ART start`)) +
  geom_point(size = 4, 
             aes(shape = `ART start`)) +
  geom_rect(aes(xmin=0, xmax=2,ymin=-Inf,ymax=Inf), 
            color = NA, 
            fill='gainsboro', 
            alpha= 0.02)+
  labs(x = "hSBA", 
       y = "Patients (%)", 
       color = "ART start", 
       linetype = "ART start") +
  scale_color_manual(values = c("Early_T0" = "darkcyan", 
                                "Late_T0" = "darkorange", 
                                "Early_T2" = "darkcyan", 
                                "Late_T2" = "darkorange")) +
  scale_linetype_manual(values = c("Early_T0" = 2, 
                                   "Late_T0" = 2, 
                                   "Early_T2" = 1, 
                                   "Late_T2" = 1)) +
  ggtitle("PHIV Patients (%) per hSBA titers") +
  theme_classic(base_size = 30) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


############################################################          Figure 2B: hSBA per TP and Group
# create a new dataframe for plotting
df = as.data.frame(men_sba) 

df_plot <- df %>%
  pivot_longer(cols = c(SBA_T0, SBA_T2, SBA_T3),
               names_to = "time",
               values_to = "value") %>%
  mutate(time = factor(time, levels = c("SBA_T0", "SBA_T2", "SBA_T3")))

df_plot$Group <- as.factor(df_plot$Group)

t0 = test_group_difference(df_plot[df_plot$time=="SBA_T0",], "Group", "value")
t2 = test_group_difference(df_plot[df_plot$time=="SBA_T2",], "Group", "value")
t3 = test_group_difference(df_plot[df_plot$time=="SBA_T3",], "Group", "value")

# add comparisons of means for each group separately
cmp_t0_t2_HC <- wilcox.test(df_plot$value[df_plot$Group == "HC" & df_plot$time == "SBA_T0"], 
                            df_plot$value[df_plot$Group == "HC" & df_plot$time =="SBA_T2"], 
                            exact = F)
cmp_t2_t3_HC <- wilcox.test(df_plot$value[df_plot$Group == "HC" & df_plot$time =="SBA_T2"], 
                            df_plot$value[df_plot$Group == "HC" & df_plot$time =="SBA_T3"],  
                            exact = F)
cmp_t0_t2_HIV <- wilcox.test(df_plot$value[df_plot$Group == "PHIV" & df_plot$time =="SBA_T0"], 
                             df_plot$value[df_plot$Group == "PHIV" & df_plot$time =="SBA_T2"], 
                             exact = F)
cmp_t2_t3_HIV <- wilcox.test(df_plot$value[df_plot$Group == "PHIV" & df_plot$time =="SBA_T2"], 
                             df_plot$value[df_plot$Group == "PHIV" & df_plot$time =="SBA_T3"], 
                             exact = F)

# plot the values with Group as hue
fc_tp_1 = ggplot(df_plot, aes(x = time, y = value, color = Group, fill = Group)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(values = c("HC" = "blue", "PHIV" = "red")) +
  scale_color_manual(values = c("HC" = "blue", "PHIV" = "red")) +
  labs(x = "Time Points", y = "hSBA", color = "Group", title = "hSBA Values by Group") +
  stat_summary(fun = "median", 
               aes(colour = Group), 
               geom = "crossbar", 
               width = 0.6,
               position = position_dodge(width = 0.75)) +
  stat_summary(fun = median, 
               geom = "line", 
               aes(group = Group, color = Group), 
               linewidth = 1.5, 
               position = position_dodge(width = 0.75)) +
  geom_signif(y_position = c(500, 525, 500, 525, 555,555, 555), 
              xmin = c(0.8, 1.2, 1.85, 2.2, 0.8, 1.8, 2.8), 
              xmax = c(1.8, 2.15, 2.8, 3.2, 1.2, 2.2, 3.2), 
              aes(size = 20),
              annotation = c("< 0.001",
                             "< 0.001",
                             "0.004",
                             "0.015",
                             round(t0$p.value,2),
                             round(t2$p.value, 2),
                             round(t3$p.value, 2)),
              tip_length = 0.01, color = "black") +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.75), alpha = 0.6, size = 3) +  # Add dodged points to the plot()
  scale_x_discrete(labels = c("T0", "T2", "T3")) +
  theme_classic(base_size = 30)


############################################################       Figure 2E : hSBA per TP and Group TREATMENT

# create a new dataframe for plotting
df = as.data.frame(men_sba)

df$Treatment[df$Treatment == "early treated"] = "Early"
df$Treatment[df$Treatment == "late treated"] = "Late"

df <- df[!is.na(df$Treatment), ]

df_plot <- df %>%
  pivot_longer(cols = c(SBA_T0, SBA_T2, SBA_T3),
               names_to = "time",
               values_to = "value") %>%
  mutate(time = factor(time, levels = c("SBA_T0", "SBA_T2", "SBA_T3")))

df_plot$Treatment <- as.factor(df_plot$Treatment)

t0 = test_group_difference(df_plot[df_plot$time=="SBA_T0",], "Treatment", "value")
t2 = test_group_difference(df_plot[df_plot$time=="SBA_T2",], "Treatment", "value")
t3 = test_group_difference(df_plot[df_plot$time=="SBA_T3",], "Treatment", "value")

# add comparisons of means for each group separately
cmp_t0_t2_HC <- wilcox.test(df_plot$value[df_plot$Treatment == "Early" & df_plot$time == "SBA_T0"], 
                            df_plot$value[df_plot$Treatment == "Early" & df_plot$time =="SBA_T2"], 
                            exact = F)
cmp_t2_t3_HC <- wilcox.test(df_plot$value[df_plot$Treatment == "Early" & df_plot$time =="SBA_T2"], 
                            df_plot$value[df_plot$Treatment == "Early" & df_plot$time =="SBA_T3"], 
                            exact = F)
cmp_t0_t2_HIV <- wilcox.test(df_plot$value[df_plot$Treatment == "Late" & df_plot$time =="SBA_T0"], 
                             df_plot$value[df_plot$Treatment == "Late" & df_plot$time =="SBA_T2"], 
                             exact = F)
cmp_t2_t3_HIV <- wilcox.test(df_plot$value[df_plot$Treatment == "Late" & df_plot$time =="SBA_T2"], 
                             df_plot$value[df_plot$Treatment == "Late" & df_plot$time =="SBA_T3"], 
                             exact = F)

# plot the values with Group as hue
colnames(df_plot)[4] = "ART start"

fc_tp_2 = ggplot(df_plot, aes(x = time, y = value, color = `ART start`, fill = `ART start`)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(values = c("Early" = "darkcyan", "Late" = "darkorange")) +
  scale_color_manual(values = c("Early" = "darkcyan", "Late" = "darkorange")) +
  labs(x = "Time Points", y = "hSBA", color = "ART start", title = "PHIV hSBA Values by ART start") +
  stat_summary(fun = "median", 
               aes(colour = `ART start`), 
               geom = "crossbar", 
               width = 0.6,
               position = position_dodge(width = 0.75)) +
  stat_summary(fun = median, 
               geom = "line", 
               aes(group = `ART start`, color = `ART start`), 
               linewidth = 1.5, 
               position = position_dodge(width = 0.75)) +
  geom_signif(y_position = c(500, 525, 500, 525, 555,555, 555), 
              xmin = c(0.8, 1.2, 1.85, 2.2, 0.8, 1.8, 2.8), 
              xmax = c(1.8, 2.15, 2.8, 3.2, 1.2, 2.2, 3.2), 
              aes(size = 20),
              annotation = c("< 0.001",
                             "< 0.001",
                             "0.02",
                             "0.12",
                             round(t0$p.value,2),
                             round(t2$p.value, 2),
                             round(t3$p.value, 2)),
              tip_length = 0.01, color = "black") +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), 
              alpha = 0.6, 
              size = 3) +  # Add dodged points to the plot() 
  scale_x_discrete(labels = c("T0", "T2", "T3")) +
  theme_classic(base_size = 30)

############################################################          Figure 2C: FC T0-T2

# Create plot
df = as.data.frame(men_sba) 

df$SBA_FC = df$SBA_T2 / df$SBA_T0
df$Seroprotection <- ifelse(df$SBA_T0 < 4, "NSP", "SP")

df = df[complete.cases(df$SBA_FC), ]

fc_t0t2_1 = ggplot(df, aes(x = Group, y = SBA_FC, color = Group)) +
  geom_jitter(data = df, 
              aes(fill = Group, 
                  size = SBA_T0), 
              position = position_jitter(width = 0.2, height = .05), 
              alpha = .35) +  # Use a different shape for specific points
  geom_hline(yintercept = 4, linetype = "dashed", color = "red", size = 1) +  # Add horizontal line
  geom_text(x = 0.65, 
            y = 4.3, 
            label = paste0("SC=", sum(df$SBA_FC[df$Group == "HC"] >= 4, na.rm = T), 
                            " (", sum(df$SBA_FC[df$Group == "HC"] >= 4, na.rm = T)/sum(df$Group == "HC", na.rm = TRUE) * 100, "%)"), 
            size = 5,  
            hjust = 0.5, 
            vjust = -0.75, 
            color = "black") +
  geom_text(x = 1.65, 
            y = 4.3, 
            label = paste0("SC=", sum(df$SBA_FC[df$Group == "PHIV"] >= 4,  na.rm = T), 
                            " (", sum(df$SBA_FC[df$Group == "PHIV"] >= 4, na.rm = T)/sum(df$Group == "PHIV", na.rm = TRUE) * 100, "%)"), 
            size = 5,  hjust = 0.5, vjust = -0.75, color = "black") +
  geom_text(x = 0.65, 
            y = -2, 
            label = paste0("NSC=", sum(df$SBA_FC[df$Group == "HC"] < 4, na.rm = T), 
                            " (", sum(df$SBA_FC[df$Group == "HC"] < 4, na.rm = T)/sum(df$Group == "HC", na.rm = TRUE) * 100, "%)"), 
            size = 5,  
            hjust = 0.5, 
            vjust = -0.75, 
            color = "black") +
  geom_text(x = 1.65, 
            y = -2, 
            label = paste0("NSC=", sum(df$SBA_FC[df$Group == "PHIV"] < 4, na.rm = T), 
                            " (", sum(df$SBA_FC[df$Group == "PHIV"] < 4, na.rm = T)/sum(df$Group == "PHIV", na.rm = TRUE) * 100, "%)"), 
            size = 5,  
            hjust = 0.5, 
            vjust = -0.75, 
            color = "black") +
  stat_compare_means(comparisons = list(c("HC", "PHIV")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     exact = FALSE,
                     size = 6) +
  stat_summary(fun.data = "mean_cl_boot", 
               colour = "black", 
               geom = "crossbar", 
               width = 0.2) +
  labs(x = "Group", 
       y = "hSBA FC T0-T2", 
       title = "hSBA Fold Change T0-T2") +
  scale_color_manual(values = c("HC" = "blue", "PHIV" = "red")) +
  scale_fill_manual(values = c("HC" = "blue", "PHIV" = "red")) +
  scale_size_continuous(range = c(3, 26), breaks = c(2, 4, 8, 16, 32, 64, 128, 256, 512)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  guides(shape = guide_legend(override.aes = list(size = 6))) +
  theme_classic(base_size = 30)

############################################################            Figure 2F: FC T0-T2 early-late

df = as.data.frame(men_sba)
df$Treatment[df$Treatment == "early treated"] = "Early"
df$Treatment[df$Treatment == "late treated"] = "Late"

df$SBA_FC = df$SBA_T2 / df$SBA_T0
df$Seroprotection <- ifelse(df$SBA_T0 < 4, "NSP", "SP")


df <- df[!is.na(df$Treatment), ]
df =  df[!is.na(df$SBA_FC), ]

df = df[complete.cases(df$SBA_FC), ]

fc_t0t2_2 = ggplot(df, aes(x = `ART start`, y = SBA_FC, color = `ART start`)) +
  geom_jitter(data = df, 
              aes(fill = `ART start`, 
                  size = SBA_T0), 
              position = position_jitter(width = 0.2, 
                                         height = .05),
              alpha = .35) +  # Use a different shape for specific points
  geom_hline(yintercept = 4, 
             linetype = "dashed", 
             color = "red", 
             size = 1) +  # Add horizontal line
  geom_text(x = 0.65, 
            y = 4.3, 
            label = paste0("SC=", sum(df$SBA_FC[df$`ART start` == "Early"] >= 4, na.rm = T), 
                    " (", sum(df$SBA_FC[df$`ART start` == "Early"] >= 4, na.rm = T) / sum(df$`ART start` == "Early", na.rm = TRUE) * 100, "%)"), 
            size = 5,  
            hjust = 0.5, 
            vjust = -0.75, 
            color = "black") +
  geom_text(x = 1.65, 
            y = 4.3, 
            label = paste0("SC=", sum(df$SBA_FC[df$`ART start` == "Late"] >= 4, na.rm = T), 
                            " (", round(sum(df$SBA_FC[df$`ART start` == "Late"] >= 4, na.rm = T)/sum(df$`ART start` == "Late", na.rm = TRUE) * 100, digits = 0), "%)"), 
            size = 5,  
            hjust = 0.5, 
            vjust = -0.75, 
            color = "black") +
  geom_text(x = 0.65, 
            y = -2, 
            label = paste0("NSC=", sum(df$SBA_FC[df$`ART start` == "Early"] < 4, na.rm = T), 
                            " (", sum(df$SBA_FC[df$`ART start` == "Early"] < 4, na.rm = T)/sum(df$`ART start` == "Early", na.rm = TRUE) * 100, "%)"), 
            size = 5,  
            hjust = 0.5, 
            vjust = -0.75, 
            color = "black") +
  geom_text(x = 1.65, 
            y = -2, 
            label = paste0("NSC=", sum(df$SBA_FC[df$`ART start` == "Late"] < 4, na.rm = T), 
                            " (", round(sum(df$SBA_FC[df$`ART start` == "Late"] < 4, na.rm = T)/sum(df$`ART start` == "Late", na.rm = TRUE) * 100, digits = 0), "%)"), 
            size = 5,  
            hjust = 0.5, 
            vjust = -0.75, 
            color = "black") +
  stat_compare_means(comparisons = list(c("Early", "Late")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     exact = FALSE, 
                     size = 6) +
  stat_summary(fun.data = "mean_cl_boot", 
               colour = "black", 
               geom = "crossbar", 
               width = 0.2) +
  labs(x = "ART start", 
       y = "hSBA FC T0-T2", 
       title = "PHIV hSBA Fold Change T0-T2", 
       color = "ART start") +
  scale_color_manual(values = c("Early" = "darkcyan", 
                                "Late" = "darkorange")) +
  scale_fill_manual(values = c("Early" = "darkcyan", 
                               "Late" = "darkorange")) +
  scale_size_continuous(range = c(3, 26), 
                        breaks = c(2, 4, 8, 16, 32, 64, 128, 256, 512)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  guides(shape = guide_legend(override.aes = list(size = 6))) +
  theme_classic(base_size = 30)


combined_plot_fig2a <- grid.arrange(
  grobs = list(plot_1, fc_tp_1, fc_t0t2_1, plot_2, fc_tp_2, fc_t0t2_2),
  widths = c(1, 1, 1, 0.1),
  layout_matrix = rbind(c(1, 2, 3, 3),
                        c(4, 5, 6, 6)))
ggsave("Figure_2.png",  combined_plot_fig2a , height = 20, width = 33, dpi = 600)

############################################################################# FIGURE 3 ########################################################################################################
############################################################     Figure 3A

df = df_frequencies
selection = c("CD19+", 
      "CD27-IgD+ (Naive)",
      "CD27+IgD+ (USM)",
      "CD27+IgD- (SM)",
      "CD27-IgD- (DN)",
      "CD27+CD21- (AM)",
      "CD27-CD21- (TLM)")
df = df_frequencies[df_frequencies$cell_populations %in% selection,]
# Reorder the levels of the cell_populations variable
order = c("CD19+", 
  "CD27-IgD+ (Naive)",
  "CD27+IgD+ (USM)",
  "CD27+IgD- (SM)",
  "CD27-IgD- (DN)",
  "CD27+CD21- (AM)",
  "CD27-CD21- (TLM)")

df_long$cell_populations <- factor(df_long$cell_populations, levels = order)

df_long$Group = as.factor(df_long$Group)

# Plot for CD19+ variable
plot_cd19 <- ggplot(df_long[df_long$cell_populations == "CD19+", ], 
                    aes(x = cell_populations, y = frequencies, fill = Group)) +
  geom_violin(position = position_dodge(width = 1), 
              scale = "width", 
              alpha = 0.7) +
  geom_boxplot(width = 0.2, 
               position = position_dodge(width = 1), 
               alpha = 0) +
  scale_fill_manual(values = c("HC" = "blue", "PHIV" = "red")) +
  geom_jitter(shape = 16, 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1)) +
  stat_compare_means(aes(group = Group), 
                     size = 6, 
                     label = "p.format", 
                     method = "wilcox.test", 
                     vjust = -.5) +
  geom_segment(aes(x = 0.8, xend = 1.2, y = 90, yend = 90)) +
  xlab("") +
  ylab("% of Lymphocytes") +
  labs(fill = "Group    ") +
  theme_classic(base_size = 30) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1, size = 20),
        axis.title.y = element_text(size = 25),
        plot.margin = margin(10, 10, 61.5, 10, unit = "pt"),
        panel.spacing = unit(10, "pt"))

# Plot for remaining variables
plot_remaining <- ggplot(df_long[df_long$cell_populations != "CD19+", ], 
                         aes(x = cell_populations, y = frequencies, fill = Group)) +
  geom_violin(position = position_dodge(width = 1), 
              scale = "width", 
              alpha = 0.7) +
  geom_boxplot(width = 0.2, 
               position = position_dodge(width = 1), 
               alpha = 0) +
  scale_fill_manual(values = c("HC" = "blue", "PHIV" = "red")) +
  geom_jitter(shape = 16, 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1)) +
  stat_compare_means(aes(group = Group), 
                     size = 6, 
                     label = "p.format", 
                     method = "wilcox.test", 
                     vjust = -.5) +
  geom_segment(aes(x=0.8, xend=1.2, y=90, yend=90)) +
  geom_segment(aes(x=1.8, xend=2.2, y=90, yend=90)) +
  geom_segment(aes(x=2.8, xend=3.2, y=90, yend=90)) +
  geom_segment(aes(x=3.8, xend=4.2, y=90, yend=90)) +
  geom_segment(aes(x=4.8, xend=5.2, y=90, yend=90)) +
  geom_segment(aes(x=5.8, xend=6.2, y=90, yend=90)) +
  xlab("") +
  ylab("% of B cells") +
  labs(fill = "Group    ") +
  theme_classic(base_size = 30) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1, size = 20),
        axis.title.y = element_text(size = 25),
        plot.margin = margin(10, 10, 0, 10, unit = "pt"),
        panel.spacing = unit(10, "pt")) +
  scale_y_continuous(position = "right") +
  theme(legend.position = "none", )

# Combine the plots

merged_plot <- ggarrange(plot_cd19, plot_remaining, 
                         nrow = 1, 
                         ncol = 2, 
                         widths = c(.3, 1.1), 
                         common.legend = TRUE, 
                         legend = "right")

############################################################    FIGURE 3 B

df = df_frequencies
df$Treatment[df$Treatment == "early treated"] = "Early"
df$Treatment[df$Treatment == "late treated"] = "Late"
selection = c("CD19+", 
      "CD27-IgD+ (Naive)",
      "CD27+IgD+ (USM)",
      "CD27+IgD- (SM)",
      "CD27-IgD- (DN)",
      "CD27+CD21- (AM)",
      "CD27-CD21- (TLM)")
df = df_frequencies[df_frequencies$cell_populations %in% selection,]
# Reorder the levels of the cell_populations variable
order = c("CD19+", 
  "CD27-IgD+ (Naive)",
  "CD27+IgD+ (USM)",
  "CD27+IgD- (SM)",
  "CD27-IgD- (DN)",
  "CD27+CD21- (AM)",
  "CD27-CD21- (TLM)")

# Reorder the levels of the cell_populations variable
order = c("CD19+", 
          "CD27-IgD+ (Naive)",
          "CD27+IgD+ (USM)",
          "CD27+IgD- (SM)",
          "CD27-IgD- (DN)",
          "CD27+CD21- (AM)",
          "CD27-CD21- (TLM)")
          
df_long$cell_populations <- factor(df_long$cell_populations, levels = order)

df_long$Treatment = as.factor(df_long$Treatment)

# Plot for CD19+ variable
plot_cd19_2 <- ggplot(df_long[df_long$cell_populations == "CD19+", ], 
                    aes(x = cell_populations, y = frequencies, fill = Treatment)) +
  geom_violin(position = position_dodge(width = 1), 
              scale = "width", 
              alpha = 0.7) +
  geom_boxplot(width = 0.2, 
               position = position_dodge(width = 1), 
               alpha = 0) +
  scale_fill_manual(values = c("Early" = "darkcyan", "Late" = "darkorange")) +
  geom_jitter(shape = 16, 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1)) +
  stat_compare_means(aes(group = Treatment), 
                     size = 6, 
                     label = "p.format", 
                     method = "wilcox.test", 
                     vjust = -.5) +
  geom_segment(aes(x = 0.8, xend = 1.2, y = 90, yend = 90)) +
  xlab("") +
  ylab("% of Lymphocytes") +
  labs(fill = "ART start") +
  #ggtitle("CD19+ B-cell Subset by Groups") +
  theme_classic(base_size = 30) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(hjust = 1, size = 20),
        axis.title.y = element_text(size = 25),
        plot.margin = margin(10, 10, 61.5, 10, unit = "pt"),
        panel.spacing = unit(10, "pt"))

# Plot for remaining variables
plot_remaining_2 <- ggplot(df_long[df_long$cell_populations != "CD19+", ], 
                         aes(x = cell_populations, y = frequencies, fill = Treatment)) +
  geom_violin(position = position_dodge(width = 1), 
              scale = "width", 
              alpha = 0.7) +
  geom_boxplot(width = 0.2, 
               position = position_dodge(width = 1), 
               alpha = 0) +
  scale_fill_manual(values = c("Early" = "darkcyan", "Late" = "darkorange")) +
  geom_jitter(shape = 16, 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1)) +
  stat_compare_means(aes(group = Treatment), 
                     size = 6, 
                     label = "p.format", 
                     method = "wilcox.test", 
                     vjust = -.5) +
  geom_segment(aes(x=0.8, xend=1.2, y=90, yend=90)) +
  geom_segment(aes(x=1.8, xend=2.2, y=90, yend=90)) +
  geom_segment(aes(x=2.8, xend=3.2, y=90, yend=90)) +
  geom_segment(aes(x=3.8, xend=4.2, y=90, yend=90)) +
  geom_segment(aes(x=4.8, xend=5.2, y=90, yend=90)) +
  geom_segment(aes(x=5.8, xend=6.2, y=90, yend=90)) +
  xlab("") +
  ylab("% of B cells") +
  labs(fill = "ART start") +
  theme_classic(base_size = 30) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(hjust = 1, size = 20),
        axis.title.y = element_text(size = 25),
        plot.margin = margin(10, 10, 0, 10, unit = "pt"),
        panel.spacing = unit(10, "pt")) +
  scale_y_continuous(position = "right") +
  theme(legend.position = "none")

# Combine the plots

merged_plot_2 <- ggarrange(plot_cd19_2, plot_remaining_2, 
                         nrow = 1, 
                         ncol = 2, 
                         widths = c(.3, 1.1), 
                         common.legend = TRUE, 
                         legend = "right")

combined_vio = grid.arrange(merged_plot, merged_plot_2, ncol = 1) #v2
ggsave("Figure_3.png", combined_vio, height = 16, width = 16, dpi = 600)

############################################################################# FIGURE 4B ########################################################################################################

orig_locale <- Sys.getlocale("LC_NUMERIC")
Sys.setlocale("LC_NUMERIC", "C")

df = df_frequencies
df$Treatment[df$Treatment == "early treated"] = "Early"
df$Treatment[df$Treatment == "late treated"] = "Late"

df$cell_populations<- gsub("fHbp\\+ CD27\\+IgD-IgG\\+ \\(Switched Memory\\)", "fHbp+ SM", df$cell_populations, ignore.case = TRUE)
df$cell_populations <- gsub("fHbp\\+ CD27\\+ CD21\\- \\(Activated Memory\\)", "fHbp+ AM", df$cell_populations, ignore.case = TRUE)
df$cell_populations <- gsub("fHbp\\+ CD27\\- CD21\\- \\(Tissue like Memory\\)", "fHbp+ TLM", df$cell_populations, ignore.case = TRUE)

selection = c("fHbp+ SM", 
      "fHbp+ AM",
      "fHbp+ TLM")
df = df_frequencies[df_frequencies$cell_populations %in% selection,]

# Reorder the levels of the cell_populations variable
order = c("fHbp+ SM",
          "fHbp+ AM",
          "fHbp+ TLM")

df_long$cell_populations <- factor(df_long$cell_populations, levels = order)

df_long$Group = as.factor(df_long$Group)

el_plot_1 = ggplot(df_long, aes(x = cell_populations, y = frequencies, fill = Group)) +
  geom_violin(position = position_dodge(width = 1), 
              scale = "width", 
              alpha = 0.7) +
  geom_boxplot(width = 0.2, 
               position = position_dodge(width = 1), 
               alpha = 0) +
  scale_fill_manual(values = c("HC" = "blue", "PHIV" = "red")) +
  geom_jitter(shape = 16, 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1)) +
  stat_compare_means(aes(group = Group), 
                     label = "p.format", 
                     size = 5, 
                     method = "wilcox.test") +
  xlab("") +
  ylab("Frequency (%)") +
  theme_classic(base_size = 30) +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~ cell_populations, scales = "free", ncol = 3) +
  theme(strip.text = element_text(size = 25))

ggsave("Figure_4B.png", el_plot_1, height = 6, width = 18, dpi = 600)

############################################################################# FIGURE 4C & 4D ########################################################################################################
 
df_wide <- df_frequencies %>%
  pivot_wider(
    id_cols = c(StudyPID, Group, Gender, Treatment),
    names_from = cell_populations,
    values_from = frequencies
  )

# Pivot the dataframe to the wide format
igg_wide <- pivot_wider(men_igg, 
                       id_cols = StudyPID, 
                       names_from = Ag, 
                       values_from = c(T0, T2), 
                       names_sep = "_"
)

# Rename the columns to match the desired format
colnames(igg_wide) <- c("StudyPID", "fHbp_T0", "NadA_T0", "NHBA_T0", "fHbp_T2", "NadA_T2", "NHBA_T2")

df <- merge(df_wide, igg_wide, by= c("StudyPID"), all.x = T)

df$SBA_res = ifelse(df$SBA_T2 >= 4*df$SBA_T0, "R", "UR")
Grp_levels <- levels(factor(df$Group.x))
resp_levels = levels(factor(df$SBA_res))
df$Grp <- factor(df$Group.x, levels = Grp_levels)
df$SBA_res = factor(df$SBA_res, levels = resp_levels)
df$SBA_FC = df$SBA_T2 / df$SBA_T0
df$fHbp_FC = df$fHbp_T2 / df$fHbp_T0
df$NHBA_FC = df$NHBA_T2 / df$NHBA_T0
df$NadA_FC = df$NadA_T2 / df$NadA_T0


# Define the variables to loop over
variables_to_plot <- c('fHbp_T2',
                       'fHbp_FC', 
                       'NadA_T2',
                       'NadA_FC',
                       'NHBA_T2',
                       'NHBA_FC',
                       'SBA_T2',
                       'SBA_FC',
                       'fHbp+ CD27+IgD-IgG+ (Switched Memory)',
                       'fHbp+ CD27+IgD-IgG+ (Switched Memory)',
                       'fHbp+ CD27- CD21- (Tissue like Memory)')
                       )

####HC
df_hc = df[df$Group.x == "HC",][variables_to_plot]

for (i in colnames(df_hc)[9:length(colnames(df_hc))]) {
  df_hc[, i] = as.numeric(gsub(",", ".", df_hc[, i]))
}

df_hc[, variables_to_plot] <- lapply(df_hc[, variables_to_plot], as.numeric)

colnames(df_hc) <- gsub("fHbp\\+ CD27\\+IgD-IgG\\+ \\(Switched Memory\\)", "fHbp+ SM", colnames(df_hc), ignore.case = TRUE)
colnames(df_hc) <- gsub("fHbp\\+ CD27\\+ CD21\\- \\(Activated Memory\\)", "fHbp+ AM", colnames(df_hc), ignore.case = TRUE)
colnames(df_hc) <- gsub("fHbp\\+ CD27\\- CD21\\- \\(Tissue like Memory\\)", "fHbp+ TLM", colnames(df_hc), ignore.case = TRUE)
colnames(df_hc) <- gsub("SBA", "hSBA", colnames(df_hc))


# Compute correlation matrix only HIV
cor_hc_mat <- cor(df_hc, use = "pairwise.complete.obs")
p_values_hc = cor.mtest(df_hc, use = "pairwise.complete.obs")$p

####HIV

df_hiv = df[df$Group.x == "HIV",][variables_to_plot]

for (i in colnames(df_hiv)[9:length(colnames(df_hiv))]) {
  df_hiv[, i] = as.numeric(gsub(",", ".", df_hiv[, i]))
}

df_hiv[, variables_to_plot] <- lapply(df_hiv[, variables_to_plot], as.numeric)

colnames(df_hiv) <- gsub("fHbp\\+ CD27\\+IgD-IgG\\+ \\(Switched Memory\\)", "fHbp+ SM", colnames(df_hiv), ignore.case = TRUE)
colnames(df_hiv) <- gsub("fHbp\\+ CD27\\+ CD21\\- \\(Activated Memory\\)", "fHbp+ AM", colnames(df_hiv), ignore.case = TRUE)
colnames(df_hiv) <- gsub("fHbp\\+ CD27\\- CD21\\- \\(Tissue like Memory\\)", "fHbp+ TLM", colnames(df_hiv), ignore.case = TRUE)
colnames(df_hiv) <- gsub("SBA", "hSBA", colnames(df_hiv))


# Compute correlation matrix only HIV
cor_hiv_mat <- cor(df_hiv, use = "pairwise.complete.obs")
p_values_hiv = cor.mtest(df_hiv, use = "pairwise.complete.obs")$p

####EARLY

df_early = df[df$Treatment == "early treated",][variables_to_plot]

for (i in colnames(df_early)[9:length(colnames(df_early))]) {
  df_early[, i] = as.numeric(gsub(",", ".", df_early[, i]))
}

df_early[, variables_to_plot] <- lapply(df_early[, variables_to_plot], as.numeric)

colnames(df_early) <- gsub("fHbp\\+ CD27\\+IgD-IgG\\+ \\(Switched Memory\\)", "fHbp+ SM", colnames(df_early), ignore.case = TRUE)
colnames(df_early) <- gsub("fHbp\\+ CD27\\+ CD21\\- \\(Activated Memory\\)", "fHbp+ AM", colnames(df_early), ignore.case = TRUE)
colnames(df_early) <- gsub("fHbp\\+ CD27\\- CD21\\- \\(Tissue like Memory\\)", "fHbp+ TLM", colnames(df_early), ignore.case = TRUE)
colnames(df_early) <- gsub("SBA", "hSBA", colnames(df_early))


# Compute correlation matrix only HIV
cor_early_mat <- cor(df_early, use = "pairwise.complete.obs")
p_values_early = cor.mtest(df_early, use = "pairwise.complete.obs")$p

####LATE

df_late = df[df$Treatment == "late treated",][variables_to_plot]

for (i in colnames(df_late)[9:length(colnames(df_late))]) {
  df_late[, i] = as.numeric(gsub(",", ".", df_late[, i]))
}

df_late[, variables_to_plot] <- lapply(df_late[, variables_to_plot], as.numeric)

colnames(df_late) <- gsub("fHbp\\+ CD27\\+IgD-IgG\\+ \\(Switched Memory\\)", "fHbp+ SM", colnames(df_late), ignore.case = TRUE)
colnames(df_late) <- gsub("fHbp\\+ CD27\\+ CD21\\- \\(Activated Memory\\)", "fHbp+ AM", colnames(df_late), ignore.case = TRUE)
colnames(df_late) <- gsub("fHbp\\+ CD27\\- CD21\\- \\(Tissue like Memory\\)", "fHbp+ TLM", colnames(df_late), ignore.case = TRUE)
colnames(df_late) <- gsub("SBA", "hSBA", colnames(df_late))

# Compute correlation matrix only HIV
cor_late_mat <- cor(df_late, use = "pairwise.complete.obs")
p_values_late = cor.mtest(df_late, use = "pairwise.complete.obs")$p

############################################################ HIV vs HC : FIGURE 4 C

# List of variable names to process
variable_names <- c("fHbp_T2", "NHBA_T2", "NadA_T2", "hSBA_T2", "hSBA_FC")
process_variable_data <- function(cor_data1, cor_data2, pval_data1, pval_data2, variable_name) {
  cat("Processing variable:", variable_name, "\n")
  
  cor_data1$Group <- "PHIV"
  cor_data2$Group <- "HC"
  cor_data1$variables <- rownames(cor_data1)
  cor_data2$variables <- rownames(cor_data2)
  
  pval_data1$Group <- "PHIV"
  pval_data2$Group <- "HC"
  pval_data1$variables <- rownames(pval_data1)
  pval_data2$variables <- rownames(pval_data2)
  
  cor_data1_wide <- reshape(cor_data1[, c(variable_name, "variables", "Group"), drop = FALSE], idvar = "Group", timevar = "variables", direction = "wide")
  cor_data2_wide <- reshape(cor_data2[, c(variable_name, "variables", "Group"), drop = FALSE], idvar = "Group", timevar = "variables", direction = "wide")
  pval_data1_wide <- reshape(pval_data1[, c(variable_name, "variables", "Group"), drop = FALSE], idvar = "Group", timevar = "variables", direction = "wide")
  pval_data2_wide <- reshape(pval_data2[, c(variable_name, "variables", "Group"), drop = FALSE], idvar = "Group", timevar = "variables", direction = "wide")
  
  combined_data <- rbind(cor_data1_wide, cor_data2_wide)
  combined_data$Ag <- gsub("_T2", "", variable_name)
  rownames(combined_data) <- combined_data$Group
  
  
  combined_pval <- rbind(pval_data1_wide, pval_data2_wide)
  combined_pval$Ag <- gsub("_T2", "", variable_name)
  rownames(combined_pval) <- combined_pval$Group
  
  
  return(list(combined_data, combined_pval))
}
# List to store the processed data for each variable
processed_data_list <- lapply(variable_names, function(var) {
  cor_data1 <- data.frame(cor_hiv_mat)
  cor_data2 <- data.frame(cor_hc_mat)
  pval_data1 <- data.frame(p_values_hiv)
  pval_data2 <- data.frame(p_values_hc)
  
  process_variable_data(cor_data1, cor_data2, pval_data1, pval_data2, var)
})

# Extract the combined data frames and p-value data frames from the processed_data_list
combined_data_list <- lapply(processed_data_list, function(res) as.matrix(res[[1]]))
combined_pval_list <- lapply(processed_data_list, function(res) as.matrix(res[[2]]))

# Merge the combined data frames and p-value data frames into one single data frame
combined_data <- as.data.frame(do.call(rbind, combined_data_list))
combined_pval <- as.data.frame(do.call(rbind, combined_pval_list))


# Convert all non-numeric columns in combined_data to numeric
combined_data[, 2:(dim(combined_data)[2] - 1)] <- lapply(
  combined_data[, 2:(dim(combined_data)[2] - 1)], 
  function(x) as.numeric(as.character(x))
)

# Convert all non-numeric columns in combined_pval to numeric
combined_pval[, 2:(dim(combined_pval)[2] - 1)] <- lapply(
  combined_pval[, 2:(dim(combined_pval)[2] - 1)], 
  function(x) as.numeric(as.character(x))
)

grp = combined_pval$Group

# Recode p-values into different categories
combined_pval[combined_pval <= 0.001] <- "***"
combined_pval[0.001 < combined_pval & combined_pval <= 0.014] <- "**"
combined_pval[combined_pval <= 0.054 & combined_pval > 0.014] <- "*"
combined_pval[combined_pval > 0.055] <- " "

# Remove the common prefix "fHbp_T2." from all column names
colnames(combined_data) <- sub("^fHbp_T2\\.", "", colnames(combined_data))
colnames(combined_pval) <- sub("^fHbp_T2\\.", "", colnames(combined_pval))

combined_pval$Group = grp

#Versione Solo Frequenze
combined_data = combined_data[combined_data$Group == "PHIV", c(1, (ncol(combined_data) - 3):ncol(combined_data))]
combined_pval = combined_pval[combined_pval$Group == "PHIV", c(1, (ncol(combined_pval) - 3):ncol(combined_pval))]

# Apply the same color scheme for annotations as before
annot <- data.frame(Group = combined_data$Group, row.names = rownames(combined_data))
colnames(annot)[1] <- "Group    "
annoCol <- list("Group    " = c(HC = "blue", PHIV = "red"))#Correlates = c(fHbp = "yellow", NHBA = "turquoise", NadA = "#CC99FF", hSBA = "lightgreen", hSBA_FC = "darkgreen"), 


# Create the pheatmap plot
cor_hm = pheatmap(t(combined_data[, 2:(dim(combined_data)[2] - 1)]), 
                  annotation_col = annot, 
                  annotation_colors = annoCol,
                  cellheight = 25,
                  cluster_cols = FALSE,
                  cluster_rows = FALSE, 
                  scale = "none",
                  show_colnames = F,
                  fontsize = 12,
                  fontsize_number = 25, 
                  #main = "Correlation Between Variables",
                  display_numbers = t(combined_pval[, 2:(dim(combined_pval)[2] - 1)]))

############################################################ Early vs Late : FIGURE 4D

# List of variable names to process
variable_names <- c("fHbp_T2", "NHBA_T2", "NadA_T2", "hSBA_T2", "hSBA_FC")
process_variable_data <- function(cor_data1, cor_data2, pval_data1, pval_data2, variable_name) {
  cat("Processing variable:", variable_name, "\n")
  
  cor_data1$Group <- "Early"
  cor_data2$Group <- "Late"
  cor_data1$variables <- rownames(cor_data1)
  cor_data2$variables <- rownames(cor_data2)
  
  pval_data1$Group <- "Early"
  pval_data2$Group <- "Late"
  pval_data1$variables <- rownames(pval_data1)
  pval_data2$variables <- rownames(pval_data2)
  
  cor_data1_wide <- reshape(cor_data1[, c(variable_name, "variables", "Group"), drop = FALSE], idvar = "Group", timevar = "variables", direction = "wide")
  cor_data2_wide <- reshape(cor_data2[, c(variable_name, "variables", "Group"), drop = FALSE], idvar = "Group", timevar = "variables", direction = "wide")
  pval_data1_wide <- reshape(pval_data1[, c(variable_name, "variables", "Group"), drop = FALSE], idvar = "Group", timevar = "variables", direction = "wide")
  pval_data2_wide <- reshape(pval_data2[, c(variable_name, "variables", "Group"), drop = FALSE], idvar = "Group", timevar = "variables", direction = "wide")
  
  combined_data <- rbind(cor_data1_wide, cor_data2_wide)
  combined_data$Ag <- gsub("_T2", "", variable_name)
  rownames(combined_data) <- combined_data$Group
  
  
  combined_pval <- rbind(pval_data1_wide, pval_data2_wide)
  combined_pval$Ag <- gsub("_T2", "", variable_name)
  rownames(combined_pval) <- combined_pval$Group
  
  
  return(list(combined_data, combined_pval))
}
# List to store the processed data for each variable
processed_data_list <- lapply(variable_names, function(var) {
  cor_data1 <- data.frame(cor_early_mat)
  cor_data2 <- data.frame(cor_late_mat)
  pval_data1 <- data.frame(p_values_early)
  pval_data2 <- data.frame(p_values_late)
  
  process_variable_data(cor_data1, cor_data2, pval_data1, pval_data2, var)
})

# Extract the combined data frames and p-value data frames from the processed_data_list
combined_data_list <- lapply(processed_data_list, function(res) as.matrix(res[[1]]))
combined_pval_list <- lapply(processed_data_list, function(res) as.matrix(res[[2]]))

# Merge the combined data frames and p-value data frames into one single data frame
combined_data <- as.data.frame(do.call(rbind, combined_data_list))
combined_pval <- as.data.frame(do.call(rbind, combined_pval_list))


# Convert all non-numeric columns in combined_data to numeric
combined_data[, 2:(dim(combined_data)[2] - 1)] <- lapply(
  combined_data[, 2:(dim(combined_data)[2] - 1)], 
  function(x) as.numeric(as.character(x))
)

# Convert all non-numeric columns in combined_pval to numeric
combined_pval[, 2:(dim(combined_pval)[2] - 1)] <- lapply(
  combined_pval[, 2:(dim(combined_pval)[2] - 1)], 
  function(x) as.numeric(as.character(x))
)



# Recode p-values into different categories
combined_pval[combined_pval <= 0.001] <- "***"
combined_pval[0.001 < combined_pval & combined_pval <= 0.014] <- "**"
combined_pval[combined_pval <= 0.054 & combined_pval > 0.014] <- "*"
combined_pval[combined_pval > 0.055] <- " "

# Remove the common prefix "fHbp_T2." from all column names
colnames(combined_data) <- sub("^fHbp_T2\\.", "", colnames(combined_data))
colnames(combined_pval) <- sub("^fHbp_T2\\.", "", colnames(combined_pval))

combined_data = combined_data[, c(1, (ncol(combined_data) - 3):ncol(combined_data))]
combined_pval = combined_pval[, c(1, (ncol(combined_pval) - 3):ncol(combined_pval))]


annot <- data.frame(Treatment = combined_data$Group, row.names = rownames(combined_data))

colnames(annot)[1] <- "ART start"

annoCol <- list("ART start" = c(Early = "darkcyan", Late = "darkorange"))

# Create the pheatmap plot
cor_el = pheatmap(t(combined_data[, 2:(dim(combined_data)[2] - 1)]), 
                  annotation_col = annot, 
                  annotation_colors = annoCol,
                  #cellwidth = .5, 
                  cellheight = 25,
                  cluster_cols = FALSE,
                  cluster_rows = FALSE, 
                  scale = "none",
                  show_colnames = F,
                  fontsize = 12,
                  fontsize_number = 25, 
                  #main = "Correlation Between Variables",
                  display_numbers = t(combined_pval[, 2:(dim(combined_pval)[2] - 1)]),
                  margin = c(10,10))

comb_hm = grid.arrange(cor_hm[[4]], cor_el[[4]], ncol = 3, widths = c(0.2, 1, 0.2), layout_matrix = rbind(c(NA,1,NA), c(NA, 2, NA)))
ggsave("Figure_4_CD.png", comb_hm, height = 6, width = 12, dpi = 600)

############################################################################# FIGURE 5 ########################################################################################################
############################################################   FIGURE 5 A-B-C

#sample names
samplenames_men = nth_element(inp_men$Name, 1, 96)

#METADATA CLEANING
samplenames_men <- gsub("MenB", "MenB_0", samplenames_men)
samplenames_men <- gsub("MenC", "MenC_", samplenames_men)
samplenames_men <- gsub("HHIV_040_Ag", "HIV_H040_Ag", samplenames_men)
samplenames_men <- gsub("HIV_011_CD19", "HIV_H011_CD19", samplenames_men)

pData(men_ct) = men_meta

#############     set categories
men_cat = setCategory(men_ct, groups = NULL, Ct.min = 13, Ct.max = 40 - 0.0001,flag = TRUE, flag.out = "Flag")

#############     FLUIDIGM: filter undetermined categories
filter_und = filterCategory(men_cat, na.categories = c("Undetermined"))

filt_samplenames = colnames(filter_und)
expr_matr = exprs(filter_und)
colnames(expr_matr) = filt_samplenames

# Calculate the proportion of non-NA values in each row
prop_non_na <- rowSums(!is.na(expr_matr)) / ncol(expr_matr)

#Filter genes with less than 30% non-NA values 
df_subset <- expr_matr[round(prop_non_na, 1) >= 0.3,]

#Compute % of NA in each column (samples)
non_na_prop <- colMeans(!is.na(df_subset))

#Filter samples with less than 60% non-NA values
M_subset <- df_subset[, round(non_na_prop, 1) >= 0.6]

# remove bad samples manually
mat_clean <- M_subset[, -grep("Ale", colnames(M_subset))]
mat_clean = mat_clean[, -grep("HC_MenB_013_CD19", colnames(mat_clean))]

#### FLUIDIGM DATA PREPROCESSING:

# Set NA values to 40
mat_clean[is.na(mat_clean)] <- 40

# Create vector of weights for cyclicloess normalization
weigt_men = c(rownames(mat_clean))

# Set all genes to same weight
weight_men[weight_men != 50] = 1
weight_men = as.numeric(weight_men)

# Convert Ct values in Et (expression tresholds: 40 - ct)
et_men = apply(mat_clean,c(1,2), function(x) 40 - x)

# Assign new values
mat_clean = et_men

#### set Et values < 13 to 13
mu = et_men
mu[mu < 13] = 13
mat_clean = mu

# NORMALIZATION
bwta_men = normalizeBetweenArrays(mat_clean, method="cyclicloess", weights = weigt_men)

raw_m <- new("qPCRset", exprs = bwta_men, featureCategory =
               featureCategory(men_cat)[rownames(mat_clean), colnames(mat_clean)])

featureNames(raw_m) <- rownames(mat_clean)
pData(raw_m) = pData(men_cat)[colnames(mat_clean),]

##### LIMMA DEGs

base = exprs(raw_m)
colnames(base) = rownames(pData(raw_m))

#### KEEP REQUIRED METADATA
samples <- pData(raw_m)[, c(1,3,4,5, 12)]
samples = na.omit(samples)

samples$pairs = str_c(samples$Grp, "_", samples$Cond)
samples$full = str_c(samples$pairs, "_", samples$SBA_res)

pids = samples$StudyPID
gruppi = samples$full
response= samples$SBA_res
gender = samples$Gender

# Define the factor levels for each variable
levels(pids) <- levels(factor(pids))
levels(gruppi)<- levels(factor(gruppi))
levels(response) <- levels(factor(response))
levels(gender) = levels(factor(gender))


for (levels in levels(gruppi)) {
  print(levels)
  print(length(grep(levels, gruppi, value = FALSE)))
}

### design matrix for LIMMA
design <- model.matrix(~ 0 + gruppi + pids, 
                       data = as.data.frame(base))

####contrast matrix
contrasts_men <- makeContrasts(
  #PHIV vs HC
  PHIV_CD19_RvsPHIV_CD19_UR = gruppiPHIV_CD19_R-gruppiPHIV_CD19_UR,
  HC_CD19_RvsHC_CD19_UR = gruppiHC_CD19_R-gruppiHC_CD19_UR,
  PHIV_CD19_RvsHC_CD19_R = gruppiPHIV_CD19_R-gruppiHC_CD19_R,
  PHIV_CD19_URvsHC_CD19_UR = gruppiPHIV_CD19_UR-gruppiHC_CD19_UR,
  PHIV_Naive_RvsPHIV_Naive_UR = gruppiPHIV_Naive_R-gruppiPHIV_Naive_UR,
  HC_Naive_RvsHC_Naive_UR = gruppiHC_Naive_R-gruppiHC_Naive_UR,
  PHIV_Naive_RvsHC_Naive_R = gruppiPHIV_Naive_R-gruppiHC_Naive_R,
  PHIV_Naive_URvsHC_Naive_UR = gruppiPHIV_Naive_UR-gruppiHC_Naive_UR,
  PHIV_NotAg_RvsPHIV_NotAg_UR = gruppiPHIV_NotAg_R-gruppiPHIV_NotAg_UR,
  HC_NotAg_RvsHC_NotAg_UR = gruppiHC_NotAg_R-gruppiHC_NotAg_UR,
  PHIV_NotAg_RvsHC_NotAg_R = gruppiPHIV_NotAg_R-gruppiHC_NotAg_R,
  PHIV_NotAg_URvsHC_NotAg_UR = gruppiPHIV_NotAg_UR-gruppiHC_NotAg_UR,
  PHIV_Ag_RvsPHIV_Ag_UR = gruppiPHIV_Ag_R-gruppiPHIV_Ag_UR,
  HC_Ag_RvsHC_Ag_UR = gruppiHC_Ag_R-gruppiHC_Ag_UR,
  PHIV_Ag_RvsHC_Ag_R = gruppiPHIV_Ag_R-gruppiHC_Ag_R,
  PHIV_Ag_URvsHC_Ag_UR = gruppiPHIV_Ag_UR-gruppiHC_Ag_UR,
  levels = colnames(design))

fit <- lmFit(select(as.data.frame(base), c(rownames(samples))),
             design)
fit2 <- contrasts.fit(fit, contrasts_men)
fit2 <- eBayes(fit2,robust = TRUE, trend = TRUE)

# Prepare output
out <- list()

# Brief summary across all contrasts/tests
res <- decideTests(fit2)
res
rownames(res) <- rownames(topTable(fit2, sort="none", n=nrow(fit2)))
out[["Summary"]] <- res
colnames(out$Summary[, colSums(out$Summary != 0) > 0])

volcano_list <- list()
diff_genes = list ()
for (coef in colnames(out$Summary)) {
  print(coef)
  if (coef == "HC_Ag_RvsHC_Ag_UR"){
    titlenam = "fHbp+ Sw.Memory DEGs profile in HC SC vs NSC"
    ydx = 4
    
  } else if (coef == "PHIV_Ag_RvsPHIV_Ag_UR") {
    titlenam = "fHbp+ Sw.Memory DEGs profile in PHIV SC vs NSC"
    ydx = 3.5
    }
  else {
    titlenam = "fHbp+ Sw.Memory DEGs profile in PHIV NSC vs HC NSC"
    ydx = 4.5
  }
  # Extract a table of differentially expressed genes
  tab <- topTable(fit2, coef=coef, number=Inf, adjust = "BH")
  
  # Extract differentially expressed genes
  de_genes <- tab %>%
    filter(adj.P.Val < 0.054)
  
  # Extract differentially expressed genes
  diff_genes[[coef]] <- de_genes
  keyvals.colour <- ifelse(
    tab$logFC <= -1.5 & tab$adj.P.Val <= 0.054, 'orchid4',
    ifelse(tab$logFC >= 1.5 & tab$adj.P.Val <= 0.054, 'yellowgreen',
           ifelse(abs(tab$logFC) > 1.5 & tab$P.Val <= 0.054, 'gray67',
                  'gray67')))
                  keyvals.colour[is.na(keyvals.colour)] <- 'gray67'
                    names(keyvals.colour)[keyvals.colour == 'yellowgreen'] <- 'Upregulated'
                    names(keyvals.colour)[keyvals.colour == 'gray67'] <- 'NS'
                    names(keyvals.colour)[keyvals.colour == 'orchid4'] <- 'Downregulated'
                    
                    volcano_plot = EnhancedVolcano(tab,
                                                   lab = rownames(tab),
                                                   x = "logFC",
                                                   y = "adj.P.Val",
                                                   ylim = c(0, ydx),#11 o 17   #5 per le Ag HIv e HC
                                                   title = "",
                                                   titleLabSize = 15,
                                                   subtitle = "",
                                                   caption = bquote("Cutoffs: "~Log[2]~ "FC: 1.5 ; p-value" [adj] ~ ": 0.05"),
                                                   captionLabSize = 10,
                                                   colCustom = keyvals.colour,
                                                   boxedLabels = TRUE,
                                                   drawConnectors = TRUE,
                                                   legendPosition = 'bottom',
                                                   widthConnectors = 0.75,
                                                   colConnectors = 'black',
                                                   gridlines.major = FALSE,
                                                   gridlines.minor = FALSE,
                                                   cutoffLineWidth = 0.8,
                                                   pCutoff = 0.054,
                                                   FCcutoff = 1.5,
                                                   pointSize = 4,
                                                   labSize = 7,
                                                   colAlpha = 1
                    )
                    volcano_plot <- volcano_plot + theme_classic(base_size = 30)
                    
                    volcano_list[[coef]] <- volcano_plot
}

vol_1 = volcano_list[14]$HC_Ag_RvsHC_Ag_UR
vol_2 = volcano_list[13]$PHIV_Ag_RvsPHIV_Ag_UR
vol_3 = volcano_list[16]$PHIV_Ag_URvsHC_Ag_UR

# Hide legends in vol_2 and vol_3 plots
vol_1 <- vol_1 + theme(legend.position = "none")
vol_2 <- vol_2 + theme(legend.position = "bottom") + guides(color = guide_legend(title = "", override.aes = list(size = 5)))
vol_3 <- vol_3 + theme(legend.position = "none")


plprov_h = vol_1 + plot_spacer() + vol_2 + plot_spacer() + vol_3 + plot_layout(widths = c(1, 0.1, 1, 0.1, 1))

ggsave("Figure_5_top.png", plprov_h, height = 12, width = 25, dpi = 600)

############################################################     FIGURE 5 D-E

genes_symbols = c("FOXP1", "CD74", "CXCR4", "CD83", "PAX5",
                  "POLH", "CD40", "TNFSF13", "IL10RA", "PRDM1",
                  "MKI67", "BTK", "MX1", "CD22", "CD86",
                  "BTLA", "CD79B", "DOCK8", "IL6R", "MTOR",
                  "SELPLG", "CYBB", "TLR9", "FYN", "BAX",
                  "TNFSF13B", "NFKB1", "CAV1", "GAPDH", "NOTCH2",
                  "SOCS1", "BATF", "STAT5A", "TRIM5", "MAPK3",
                  "CD80", "IL10", "STAT3", "IGHM", "PLCG1",
                  "TIRAP", "IFIT2", "RUNX3", "ETS1", "IGHD",
                  "IL21R", "CAMK4", "MZB1", "STAT1", "TNFRSF13B",
                  "FOXO1", "BCL2", "LILRB1", "IL2RA", "SIVA1",
                  "PIK3C2B", "BCL6", "IFNAR2", "TNSF14", "IRAK4",
                  "IL4R", "PPP1R13B", "CD27", "ZBTB32", "SLAMF1",
                  "CXCR3", "BST2", "IRF4", "IKBKG", "SAMHD1",
                  "PPP3CA", "NOD2", "MYD88", "ABCB1", "NKRF",
                  "CXCR5", "FOXO3", "SYK", "APOBEC3G", "STAT4",
                  "PTEN", "ITCH", "FAS", "TCF3", "CD38")

genes_symbols_hiv = c("TNFRSF13B", "IL4R", "IKBKG", "NOD2", "PPP1R13B")

de = mapIds(org.Hs.eg.db, genes_symbols_pv, 'ENTREZID', 'SYMBOL')
#Per HIV
de_hiv = mapIds(org.Hs.eg.db, genes_symbols_hiv, 'ENTREZID', 'SYMBOL')

edo <- enrichGO(de, 'org.Hs.eg.db', ont = "BP",  pAdjustMethod = "fdr", minGSSize = 2)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
edox2 <- pairwise_termsim(edox)

edo_hiv <- enrichGO(de_hiv, 'org.Hs.eg.db', ont = "BP",  pAdjustMethod = "fdr", minGSSize = 2)
edox_hiv <- setReadable(edo_hiv, 'org.Hs.eg.db', 'ENTREZID')
edox2_hiv <- pairwise_termsim(edox_hiv)

go_hc = mutate(edox2, GeneRatio) %>% 
  dotplot(x="GeneRatio", showCategory = 10) +
  scale_color_continuous(type = "gradient", low = "blue", high = "blue") + 
  scale_size_continuous(range = c(15, 35), breaks = c(1, 2, 3, 4, 5), labels = c("1", "2", "3", "4", "5")) +
  theme_classic(base_size = 30) + theme(legend.position = "bottom")

go_hc <- go_hc + scale_y_discrete(labels = function(labels) str_wrap(labels, width = 50))

# Remove the legends from individual plots
go_hc_no_legend <- go_hc + theme_classic(base_size = 30) + theme(legend.position = "right")

go_hiv = mutate(edox2_hiv, GeneRatio) %>% 
  dotplot(x="GeneRatio", showCategory = 10) +
  scale_color_continuous(type = "gradient", low = "red", high = "red") + 
  scale_size_continuous(range = c(15, 25), breaks = c(1, 2, 3), labels = c("1", "2", "3")) +
  theme_classic(base_size = 30)

# Apply the function to the y-axis labels
go_hiv <- go_hiv + scale_y_discrete(labels = function(labels) str_wrap(labels, width = 50))
# Customize the appearance of ticks
go_hiv <- go_hiv + theme(axis.text = element_text(face = "bold"))
go_hiv_no_legend <- go_hiv + theme_classic(base_size = 30) + theme(legend.position = "none")

# Combine the plots and specify the position of the legend
go_plots <- go_hc_no_legend + plot_spacer() + go_hiv_no_legend +
  plot_layout(guides = 'collect', widths = c(1, 0.1, 1)) 

#versione orizzontale
ggsave("Figure_5_bot.png", go_plots, height = 15, width = 43, dpi = 600)

