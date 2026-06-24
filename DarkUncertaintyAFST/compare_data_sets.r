
##---Plots--------------------------------------
# plot real BACON2 data
ggplot(ratiodf, aes(x=as.character(date), y=offset)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=offset-statistical_unc, ymax=offset+statistical_unc), width=0) +
  ylim(c(-110,-90)) +
  theme_bw() +
  ggtitle("BACON2 data")

# plot simulated data set
# plot includes mu_hat with uncertainty bars showing k*u, where k=1 and true mu (red)
p1 <- ggplot(sim_dat, aes(x=day, y=x)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=x-u, ymax=x+u), width=0) +
  geom_hline(aes(yintercept = mu, color = "true mu", linetype = "true mu")) + 
  geom_hline(aes(yintercept = mu_hat_b, color = "weighted mean", linetype = "weighted mean")) + 
  geom_hline(aes(yintercept = mu_hat_b + u_mu_hat_b, color = "unadjusted Birge bounds", linetype = "unadjusted Birge bounds")) + 
  geom_hline(aes(yintercept = mu_hat_b - u_mu_hat_b, color = "unadjusted Birge bounds", linetype = "unadjusted Birge bounds")) +
  geom_hline(aes(yintercept = mu_hat_b + u_mu_hat_b_corrected, color = "adjusted Birge bounds", linetype = "adjusted Birge bounds")) + 
  geom_hline(aes(yintercept = mu_hat_b - u_mu_hat_b_corrected, color = "adjusted Birge bounds", linetype = "adjusted Birge bounds")) + 
  # ylim(c(-6,6)) +
  # ylim(c(-110,-90)) + # if use mu=mean(BACON2 data)
  scale_color_manual(name = "Legend",
              breaks = c("true mu", 
               "weighted mean", 
               "unadjusted Birge bounds", 
               "adjusted Birge bounds"),
              values = c("true mu" = "red", 
                "weighted mean" = "purple", 
                "unadjusted Birge bounds" = "purple", 
                "adjusted Birge bounds" = "purple"),
              labels = c("true mu" = expression(mu), 
               "weighted mean" = expression(hat(mu)[WM]),
               "unadjusted Birge bounds" = "unadjusted Birge bounds",
               "adjusted Birge bounds" = "adjusted Birge bounds")) +
  scale_linetype_manual(name = "Legend",
              breaks = c("true mu", 
               "weighted mean", 
               "unadjusted Birge bounds", 
               "adjusted Birge bounds"),
              values = c("true mu" = 1, 
                "weighted mean" = 1, 
                "unadjusted Birge bounds" = 2, 
                "adjusted Birge bounds" = 3),
              labels = c("true mu" = expression(mu), 
               "weighted mean" = expression(hat(mu)[WM]),
               "unadjusted Birge bounds" = "unadjusted Birge bounds",
               "adjusted Birge bounds" = "adjusted Birge bounds")) +
  theme_bw() +
  ylab(expression(x[i] %+-% u(x[i]))) + xlab("Day") + 
  ggtitle(paste("Simulated data, N =", N, ", c =", birg_constant))

#ggsave("DarkUncertaintyAFST/mul_birge_plot.png", plot = p1, width = 6, height = 4, units = "in", dpi = 300)


