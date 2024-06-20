#(E) The heatmap illustrates Pearson’s correlation between significant differential abundant KOs and species-level genome bins (SGBs), while the network shows the correlations between cytokine levels, including interferon (IFN)-γ, interleukin (IL)-10, IL-1β, IL-2, and IL-6, and these KOs and SGBs, evaluated by the partial Mantel test.

##安装包
# install.packages("devtools")
devtools::install_github("Hy4m/linkET", force = TRUE)
packageVersion("linkET")
##作图
library(linkET)
library(dplyr)
library(ggplot2)
setwd("path/to/your/directory")  #设置工作路径
varespec <- read.delim("varespec.txt", header = TRUE, sep = "\t")
varechem <- read.delim("varechem.txt", header = TRUE, sep = "\t")
mantel <- mantel_test(varespec, varechem,
                      spec_select = list(Spec01 = 1,
                                         Spec02 = 2,
                                         Spec03 = 3,
                                         Spec04 = 4,
                                         Spec05 = 5)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


set_corrplot_style(scale = ggplot2::scale_fill_viridis_c())
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
