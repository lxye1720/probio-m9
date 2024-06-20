#(E) The heatmap illustrates Pearson’s correlation between significant differential abundant KOs and species-level genome bins (SGBs), while the network shows the correlations between cytokine levels, including interferon (IFN)-γ, interleukin (IL)-10, IL-1β, IL-2, and IL-6, and these KOs and SGBs, evaluated by the partial Mantel test.

##安装包
# install.packages("devtools")
devtools::install_github("Hy4m/linkET", force = TRUE)
packageVersion("linkET")
##作图
library(linkET)
library(dplyr)
library(ggplot2)

mantel <- mantel_test(varespec, varechem,
                      spec_select = list(Spec01 = 1:7,
                                         Spec02 = 8:18,
                                         Spec03 = 19:37,
                                         Spec04 = 38:44)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


set_corrplot_style()     ##colours = c("red", "white", "blue")设置红白蓝配色，scale = ggplot2::scale_fill_viridis_c()scale函数设置其他配色；
set_default_style()    ##还原样式
qcorrplot(correlate(varechem), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel, 
              curvature = nice_curvature()) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))