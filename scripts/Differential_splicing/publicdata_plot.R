setwd("/work/users/l/i/limjw123/External/cluster_OA")
library(limma)
library(magrittr)
library(data.table)
library(dplyr)
library(leafviz)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(plotgardener)
library(grid)
library(ggrepel)
library(leafcutter)
library(colorspace)
library(ggtext)

#leafviz("../cluster_next500_OA/norm_vs_oa.Rdata")
#calculate deltapis for all of the pbs vs fnf
ratios_fnf <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/clu_fnf/ratio_fnf.txt")

ratios_fnf$MeanDiff <- apply(ratios_fnf, 1, function(row) {
  ctl_cols <- grep("CTL", colnames(ratios_fnf))
  fnf_cols <- grep("FNF", colnames(ratios_fnf))
  
  ctl_mean <- mean(as.numeric(row[ctl_cols]), na.rm = TRUE)
  fnf_mean <- mean(as.numeric(row[fnf_cols]), na.rm = TRUE)
  
  if (is.na(ctl_mean) || is.na(fnf_mean)) {
    NA
  } else {
    fnf_mean - ctl_mean
  }
})

ratios_fnf <- ratios_fnf %>%
  separate(Junction, into = c("part1", "part2", "part3","part4"), sep = ":", fill = "right", extra = "drop") %>%
  mutate(loc = paste(part1, part2, part3, sep = ":")) %>%
           dplyr::select(-part1, -part2, -part3, -part4)

norm_oa <- load("../cluster_next500_OA/norm_vs_oa.Rdata")

for (name in norm_oa) {
  # Construct the new name by prefixing with "fnf_"
  new_name <- paste0(name,"_public")
  
  # Assign the object to the new name in the global environment
  assign(new_name, get(name))
  remove(name)
}
normalize_column <- function(x) {
  x / sum(x, na.rm = TRUE)
}
ratios_public <- counts_public  %>%
  mutate(clu = str_split_fixed(rownames(counts_public ), ":", 4)[,4]) %>%
  group_by(clu) %>%
  mutate_all(normalize_column) %>%
  ungroup() %>%
  as.data.frame() %>%
  set_rownames(rownames(counts_public )) %>%
  dplyr::select(-clu)

ratios_public_qc = ratios_public[rowMeans(is.na(ratios_public)) <= 0.4,,drop=F ] #Try to remove 40% or less NA values in each rows
#10304  to 9677 --> 627 rows dropped
row_means_public = rowMeans(ratios_public_qc, na.rm = T) #calculate the mean without NAs in the rows. 
row_means_outer_public = outer(row_means_public, rep(1,ncol(ratios_public_qc)))# making outlier to the NAs
ratios_public_qc[is.na(ratios_public_qc)] = row_means_outer_public[is.na(ratios_public_qc)] # instead of Na, add the rowmean
colnames(ratios_public_qc) <- colnames(ratios_public_qc)
ratios_public_qc <- cbind(rownames(ratios_public_qc), ratios_public_qc)
colnames(ratios_public_qc)[1] <- c('Junction')



ratios_public_qc$MeanDiff <- apply(ratios_public_qc, 1, function(row) {
  norm_cols <- grep("_1", colnames(ratios_public_qc))
  oa_cols <- grep("_4", colnames(ratios_public_qc))
  
  norm_mean <- mean(as.numeric(row[norm_cols]), na.rm = TRUE)
  oa_mean <- mean(as.numeric(row[oa_cols]), na.rm = TRUE)
  
  if (is.na(norm_mean) || is.na(oa_mean)) {
    NA
  } else {
    oa_mean - norm_mean
  }
})

ratios_public_qc <- ratios_public_qc %>%
  separate(Junction, into = c("part1", "part2", "part3","part4"), sep = ":", fill = "right", extra = "drop") %>%
  mutate(loc = paste(part1, part2, part3, sep = ":")) %>%
  dplyr::select(-part1, -part2, -part3, -part4)

dim(ratios_public_qc)
introns_public_subset_sig <- introns_public[abs(introns_public$deltapsi) >= 0.2,]
introns_public_subset_sig$loc <- paste0(introns_public_subset_sig$chr,":",introns_public_subset_sig$start,":",introns_public_subset_sig$end)

introns_public_subset_sig_up <- introns_public_subset_sig %>% dplyr::filter(deltapsi > 0) %>% 
  select(loc,deltapsi)
introns_public_subset_sig_down <- introns_public_subset_sig %>% dplyr::filter(deltapsi < 0) %>% 
  select(loc,deltapsi)

fnf_all_public_subset <- bind_rows(introns_public_subset_sig_up |> mutate(group = "Up in GSE114007"),
                                   introns_public_subset_sig_down |> mutate(group = "Down in GSE114007")) |>
  mutate(group = factor(group, levels = c("Up in GSE114007", "Down in GSE114007")))

sig_public_subset_psi20 <- fnf_all_public_subset %>%
  inner_join(ratios_fnf %>% select(loc, MeanDiff), by = "loc")

public_from_all_fnf_up <- ratios_fnf %>% filter(loc %in% introns_public_subset_sig_up$loc) %>%
  dplyr::select(loc,MeanDiff)
public_from_all_fnf_down <- ratios_fnf %>% filter(loc %in% introns_public_subset_sig_down$loc) %>%
  dplyr::select(loc,MeanDiff)

# wilcox test

up_test_public <- wilcox.test(x = as.numeric(introns_public_subset_sig_up$deltapsi),
                          mu=0,
                          alternative = "greater")

down_test_public <- wilcox.test(x = as.numeric(introns_public_subset_sig_down$deltapsi),
                                mu=0,
                            alternative = "less")

cat("p-value:", format.pval(up_test_public$p.value, digits = 3), "\n")
cat("p-value:", format.pval(down_test_public$p.value, digits = 3), "\n")

gse11407_boxplot_nextseq <- ggplot(sig_public_subset_psi20 ,aes(x = group, y = MeanDiff, fill = group))  +
  geom_hline(yintercept = 0, lty = 2, color = "grey25", linewidth = 0.25) +
  geom_jitter(width = 0.2, color = "grey40", size = 0.25) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5, alpha = 0.4,width = 0.5,color=c("#FFB81C",'#005587')) +
  #facet_wrap(vars(group), nrow = 2, strip.position = "top", scales = "free_x") 
  stat_boxplot(geom = "errorbar", width = 0.5,color=c("#FFB81C",'#005587'))+
  scale_color_manual(values = c(darken("#FFB81C", 0.3),darken('#005587', 0.3))) +
  scale_fill_manual(values = c("#FFB81C",'#005587')) +
  scale_y_continuous(name = "deltaPSI in response to FN-f",
                     limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3,0.2)) +
  coord_cartesian(clip = "off") +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(size = 6, family = "Helvetica",
                                        margin = margin(r = -15)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 8, margin = margin(b = -1)),
        strip.background = element_blank(),
        strip.text.x.bottom = element_markdown(size = 8, margin = margin(t = 1)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, family = "Helvetica",
                                  size = 10, margin = margin(b=-3)))

ggsave(filename = "/work/users/s/e/seyoun/CQTL_sQTL/output/results_plots/Figure1_GSE114007_boxplot_nextseq.pdf",
       plot = gse11407_boxplot_nextseq, width = 5, height = 5, units = "in")
save(gse11407_boxplot_nextseq, file = "/work/users/s/e/seyoun/CQTL_sQTL/output/results_plots/Figure1_differntial_splicing/gse11407_nextseq_boxplot.rda")


#-------------------------------------------------------------------------------
#katsoula

setwd("/work/users/s/e/seyoun/CQTL_sQTL/external_data")
katsoula <- fread("./katsoula.suptable8.txt")

katsoula$loc <- paste0(katsoula$chr,":",katsoula$start,":",katsoula$end)
katsoula_subset_sig <- katsoula[abs(katsoula$deltapsi) >= 0.05,]

katsoula_up <- katsoula_subset_sig %>% dplyr::filter(deltapsi > 0) %>% 
  select(loc,deltapsi)
katsoula_down <- katsoula_subset_sig %>% dplyr::filter(deltapsi < 0) %>% 
  select(loc,deltapsi)

fnf_all_public_subset <- bind_rows(katsoula_up |> mutate(group = "Up in Katsoula et.al"),
                                   katsoula_down |> mutate(group = "Down Katsoula et.al")) |>
  mutate(group = factor(group, levels = c("Up in Katsoula et.al", "Down Katsoula et.al")))

sig_public_subset <- fnf_all_public_subset %>%
  inner_join(ratios_fnf %>% select(loc, MeanDiff), by = "loc")

public_from_all_fnf_up <- ratios_fnf %>% filter(loc %in% katsoula_up$loc) %>%
  dplyr::select(loc,MeanDiff)
public_from_all_fnf_down <- ratios_fnf %>% filter(loc %in% katsoula_down$loc) %>%
  dplyr::select(loc,MeanDiff)


fnf_subset_katsoula_douwn <- sig_public_subset %>%
  filter(group == "Down Katsoula et.al") %>%
  select(MeanDiff) %>%
  unlist() %>%
  as.numeric()

# wilcox test

up_test_public <- wilcox.test(x = as.numeric(katsoula_up$deltapsi),
                              y = as.numeric(ratios_fnf$MeanDiff),
                              alternative = "greater")

down_test_public <- wilcox.test(x = as.numeric(katsoula_down$deltapsi),
                                y = sig_public_subset$MeanDiff,
                                alternative = "less")
katsoula_boxplot <- ggplot(sig_public_subset ,aes(x = group, y = MeanDiff, fill = group))  +
  geom_hline(yintercept = 0, lty = 2, color = "grey25", linewidth = 0.25) +
  geom_jitter(width = 0.2, color = "grey40", size = 0.25) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5, alpha = 0.4,width = 0.5,color=c("#FFB81C",'#005587')) +
  #facet_wrap(vars(group), nrow = 2, strip.position = "top", scales = "free_x") 
  stat_boxplot(geom = "errorbar", width = 0.5,color=c("#FFB81C",'#005587'))+
  scale_color_manual(values = c(darken("#FFB81C", 0.3),darken('#005587', 0.3))) +
  scale_fill_manual(values = c("#FFB81C",'#005587')) +
  scale_y_continuous(name = "deltaPSI in response to FN-f",
                     limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4,0.1)) +
  coord_cartesian(clip = "off") +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(size = 6, family = "Helvetica",
                                        margin = margin(r = -15)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 8, margin = margin(b = -1)),
        strip.background = element_blank(),
        strip.text.x.bottom = element_markdown(size = 8, margin = margin(t = 1)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, family = "Helvetica",
                                  size = 10, margin = margin(b=-3)))


ggsave(filename = "/work/users/s/e/seyoun/CQTL_sQTL/output/results_plots/Figure1_Katsoula_boxplot.pdf",
       plot = katsoula_boxplot, width = 5, height = 5, units = "in")
save(katsoula_boxplot, file = "/work/users/s/e/seyoun/CQTL_sQTL/output/results_plots/Figure1_differntial_splicing/katsoula_boxplot.rda")


