
library(tidyverse)

exp_outputs <- read_csv("~/Documents/GitHub/FishSET/data/confidential/rds/exp_outputs_rtmb.csv")

EUR_MOR <- exp_outputs %>%
  filter(port == "Morro" | port == "Eureka")

EUR_MOR$port <- factor(EUR_MOR$port, levels = c("Eureka", "Morro"))

EUR_MOR$temp.window <- factor(EUR_MOR$temp.window, levels = c(7, 14, 30, 90))

ggplot(EUR_MOR, aes(x = temp.window, y = catch)) +
  geom_point(color = "black", size = 3, alpha = 0.8) +
  facet_wrap(~port,  nrow = 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1) # Add black border to facets
  )
