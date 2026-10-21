library(ggplot2)

p <- readRDS("C:/Users/fvsk236/Desktop/Results_Paper_13/Figures/N=240/Model_2_nT_12_nB_5_aP_1/Model_2_nT_12_nB_5_aP_1_a_1")

p <- p +
  labs(
    x = "Randomisation time to target",
    fill = "Strategy"
  ) +
  scale_fill_discrete(
    breaks = c("Wal25", "SPP2", "Non.rand", "IBD", "Flexible", "Flex.SPP2"),
    labels = c("E1", "E2", "E3", "N1", "N2", "N3")
  )

p