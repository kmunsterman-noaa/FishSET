
library(tidyverse)

exp_outputs <- read_csv("~/Documents/GitHub/FishSET/data/confidential/rds/exp_outputs.csv")

exp_outputs$port <- factor(exp_outputs$port, levels = c("brookings", "eureka", "fort bragg", "morro", "crescent city"))
exp_outputs$temp.window <- factor(exp_outputs$temp.window, levels = c(7, 14, 30, 90, 180))

ggplot(exp_outputs, aes(x = temp.window, y = exp.1)) +
  geom_point(color = "black", size = 3, alpha = 0.8) +
  facet_wrap(~ port,  nrow = 1) + # Facet in one row
  labs(
    x = "temp.window",
    y = "exp. coefficient"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_rect(fill = "grey", color = "grey"),
    strip.text = element_text(face = "bold", color = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1) # Add black border to facets
  )
