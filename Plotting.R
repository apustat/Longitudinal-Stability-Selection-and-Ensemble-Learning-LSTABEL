library(ggplot2)
library(dplyr)

d1=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/indep-n=600-p=200-sig2u=0.8.csv")
d2=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/indep-n=600-p=800-sig2u=0.8.csv")
d3=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/indep-n=600-p=200-sig2u=1.2.csv")
d4=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/indep-n=600-p=800-sig2u=1.2.csv")

apply(d1, 2, mean); apply(d1, 2, sd)
apply(d2, 2, mean); apply(d2, 2, sd)
apply(d3, 2, mean); apply(d3, 2, sd)
apply(d4, 2, mean); apply(d4, 2, sd)

data1 <- data.frame(
  Methods = rep(c("LSTABEL", "GLMMLASSO", "GLMMLASSO-EL", "LSTABEL-GLMM", "LSTABEL-RF"), each = 400),
  Sigma_squared = rep(rep(c(0.8, 1.2), each = 200), times = 5),
  p = rep(rep(c(200, 800), each = 100), times = 5 * 2),
  MSE =c(d1$pmse.el, d2$pmse.el, d3$pmse.el, d4$pmse.el, d1$pmse.lasso, d2$pmse.lasso, d3$pmse.lasso, d4$pmse.lasso, d1$pmse.lasso.el, 
         d2$pmse.lasso.el, d3$pmse.lasso.el, d4$pmse.lasso.el, d1$pmse.glm, d2$pmse.glm, d3$pmse.glm, d4$pmse.glm,
         d1$pmse.rf, d2$pmse.rf, d3$pmse.rf, d4$pmse.rf)
)

data1$Method <- factor(data1$Method, levels = c("LSTABEL", "GLMMLASSO", "GLMMLASSO-EL", "LSTABEL-GLMM", "LSTABEL-RF"))  #
custom_colors <- c("0.8" = "steelblue1", "1.2" = "palegreen")  # Replace with your desired colors

ggplot(data1, aes(x = Method, y = MSE, fill = factor(Sigma_squared))) +
  geom_boxplot() +
  facet_wrap(~ p, scales = "free_x", nrow = 1, labeller = labeller(p = function(x) paste0("p=", x))) +
  labs( y = "MSPE", fill = expression(sigma[b]^2) ) +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 13),  # Increase x-axis label size
    axis.title.y = element_text(size = 13),  # Increase y-axis label size
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),  # Rotate x-axis labels
    strip.text = element_text(size = 13),  # Adjust facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),  # Increase legend title size
    #axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    plot.title = element_blank(),           # Remove plot title
  )



######################################################################################
library(ggplot2)
library(dplyr)

d5=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/indep-n=1000-p=500-sig2u=0.8.csv")
d6=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/indep-n=1000-p=1500-sig2u=0.8.csv")
d7=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/indep-n=1000-p=500-sig2u=1.2.csv")
d8=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/indep-n=1000-p=1500-sig2u=1.2.csv")

apply(d5, 2, mean); apply(d5, 2, sd)
apply(d6, 2, mean); apply(d6, 2, sd)
apply(d7, 2, mean); apply(d7, 2, sd)
apply(d8, 2, mean); apply(d8, 2, sd)

data2 <- data.frame(
  Methods = rep(c("LSTABEL", "GLMMLASSO", "GLMMLASSO-EL", "LSTABEL-GLMM", "LSTABEL-RF"), each = 400),
  Sigma_squared = rep(rep(c(0.8, 1.2), each = 200), times = 5),
  p = rep(rep(c(500, 1500), each = 100), times = 5 * 2),
  MSE =c(d5$pmse.el, d6$pmse.el, d7$pmse.el, d8$pmse.el, d5$pmse.lasso, d6$pmse.lasso, d7$pmse.lasso, d8$pmse.lasso, 
         d5$pmse.lasso.el, d6$pmse.lasso.el, d7$pmse.lasso.el, d8$pmse.lasso.el, d5$pmse.glm, d6$pmse.glm, d7$pmse.glm, d8$pmse.glm,
         d5$pmse.rf, d6$pmse.rf, d7$pmse.rf, d8$pmse.rf)
)

data2$Method <- factor(data2$Method, levels = c("LSTABEL", "GLMMLASSO", "GLMMLASSO-EL", "LSTABEL-GLMM", "LSTABEL-RF"))  #
custom_colors <- c("0.8" = "steelblue1", "1.2" = "palegreen")  # Replace with your desired colors

ggplot(data2, aes(x = Method, y = MSE, fill = factor(Sigma_squared))) +
  geom_boxplot() +
  facet_wrap(~ p, scales = "free_x", nrow = 1, labeller = labeller(p = function(x) paste0("p=", x))) +
  labs( y = "MSPE", fill = expression(sigma[b]^2) ) +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 13),  # Increase x-axis label size
    axis.title.y = element_text(size = 13),  # Increase y-axis label size
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),  # Rotate x-axis labels
    strip.text = element_text(size = 13),  # Adjust facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),  # Increase legend title size
    #axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    plot.title = element_blank(),           # Remove plot title
  )

###########################################################################
library(ggplot2)
library(dplyr)

d9=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/corr-n=600-p=200-sig2u=0.8.csv")
d10=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/corr-n=600-p=800-sig2u=0.8.csv")
d11=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/corr-n=600-p=200-sig2u=1.2.csv")
d12=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/corr-n=600-p=800-sig2u=1.2.csv")

apply(d9, 2, mean); apply(d9, 2, sd)
apply(d10, 2, mean); apply(d10, 2, sd)
apply(d11, 2, mean); apply(d11, 2, sd)
apply(d12, 2, mean); apply(d12, 2, sd)

data3 <- data.frame(
  Methods = rep(c("LSTABEL", "GLMMLASSO", "GLMMLASSO-EL", "LSTABEL-GLMM", "LSTABEL-RF"), each = 400),
  Sigma_squared = rep(rep(c(0.8, 1.2), each = 200), times = 5),
  p = rep(rep(c(200, 800), each = 100), times = 5 * 2),
  MSE =c( d9$pmse.el, d10$pmse.el, d11$pmse.el, d12$pmse.el,  d9$pmse.lasso, d10$pmse.lasso, d11$pmse.lasso, d12$pmse.lasso, 
          d9$pmse.lasso.el, d10$pmse.lasso.el, d11$pmse.lasso.el, d12$pmse.lasso.el, d9$pmse.glm, d10$pmse.glm, d11$pmse.glm, d12$pmse.glm,
         d9$pmse.rf, d10$pmse.rf, d11$pmse.rf, d12$pmse.rf)
)

data3$Method <- factor(data3$Method, levels = c("LSTABEL", "GLMMLASSO", "GLMMLASSO-EL", "LSTABEL-GLMM", "LSTABEL-RF"))  #
custom_colors <- c("0.8" = "steelblue1", "1.2" = "palegreen")  # Replace with your desired colors

ggplot(data3, aes(x = Method, y = MSE, fill = factor(Sigma_squared))) +
  geom_boxplot() +
  facet_wrap(~ p, scales = "free_x", nrow = 1, labeller = labeller(p = function(x) paste0("p=", x))) +
  labs( y = "MSPE", fill = expression(sigma[b]^2) ) +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 13),  # Increase x-axis label size
    axis.title.y = element_text(size = 13),  # Increase y-axis label size
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),  # Rotate x-axis labels
    strip.text = element_text(size = 13),  # Adjust facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),  # Increase legend title size
    #axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    plot.title = element_blank(),           # Remove plot title
  )


######################################################################################
library(ggplot2)
library(dplyr)

d13=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/corr-n=1000-p=500-sig2u=0.8.csv")
d14=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/corr-n=1000-p=1500-sig2u=0.8.csv")
d15=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/corr-n=1000-p=500-sig2u=1.2.csv")
d16=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Longitudinal Project/R code/corr-n=1000-p=1500-sig2u=1.2.csv")

apply(d13, 2, mean); apply(d13, 2, sd)
apply(d14, 2, mean); apply(d14, 2, sd)
apply(d15, 2, mean); apply(d15, 2, sd)
apply(d16, 2, mean); apply(d16, 2, sd)

data4 <- data.frame(
  Methods = rep(c("LSTABEL", "GLMMLASSO", "GLMMLASSO-EL", "LSTABEL-GLMM", "LSTABEL-RF"), each = 400),
  Sigma_squared = rep(rep(c(0.8, 1.2), each = 200), times = 5),
  p = rep(rep(c(500, 1500), each = 100), times = 5 * 2),
  MSE =c(d13$pmse.el, d14$pmse.el, d15$pmse.el, d16$pmse.el, d13$pmse.lasso, d14$pmse.lasso, d15$pmse.lasso, d16$pmse.lasso, 
         d13$pmse.lasso.el, d14$pmse.lasso.el, d15$pmse.lasso.el, d16$pmse.lasso.el, d13$pmse.glm, d14$pmse.glm, d15$pmse.glm, d16$pmse.glm,
         d13$pmse.rf, d14$pmse.rf, d15$pmse.rf, d16$pmse.rf)
)

data4$Method <- factor(data4$Method, levels = c("LSTABEL", "GLMMLASSO", "GLMMLASSO-EL", "LSTABEL-GLMM", "LSTABEL-RF"))  #
custom_colors <- c("0.8" = "steelblue1", "1.2" = "palegreen")  # Replace with your desired colors

ggplot(data4, aes(x = Method, y = MSE, fill = factor(Sigma_squared))) +
  geom_boxplot() +
  facet_wrap(~ p, scales = "free_x", nrow = 1, labeller = labeller(p = function(x) paste0("p=", x))) +
  labs( y = "MSPE", fill = expression(sigma[b]^2) ) +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 13),  # Increase x-axis label size
    axis.title.y = element_text(size = 13),  # Increase y-axis label size
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),  # Rotate x-axis labels
    strip.text = element_text(size = 13),  # Adjust facet label size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),  # Increase legend title size
    #axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    plot.title = element_blank(),           # Remove plot title
  )
