library(tidyverse)
library(gridExtra)

### function to generate x_i and sigma_i for each day i
plot_data <- function(df){
# df is a data frame with columns "value" and "day"
    summary_info = df %>% 
                    group_by(day) %>% 
                    summarise(mean_val = mean(value), sd_val = sd(value))
    day = 1:5
    point_est = summary_info$mean_val
    sd_est = summary_info$sd_val
    LB = summary_info$mean_val - summary_info$sd_val
    UB = summary_info$mean_val + summary_info$sd_val
    return(data.frame(cbind(day, LB, point_est, UB, sd_est)))
}

### Inconsistent data
set.seed(409)
day1_inc <- rnorm(11, 0, 0.85)
day2_inc <- rnorm(15, -1, 1)
day3_inc <- rnorm(3, 5, 0.7)
day4_inc <- rnorm(10, 4, 1)
day5_inc <- rnorm(15, -2.5, 0.9)

days_inc <- c(rep(1, length(day1_inc)),
        rep(2, length(day2_inc)),
        rep(3, length(day3_inc)),
        rep(4, length(day4_inc)),
        rep(5, length(day5_inc))
        )

df_inc <- cbind(c(day1_inc, day2_inc, day3_inc, day4_inc, day5_inc), as.factor(days_inc))
df_inc <- as.data.frame(df_inc)
colnames(df_inc) <- c("value", "day")
data_inc <- plot_data(df_inc)

# plot1 <- ggplot(df_inc) + 
#     geom_boxplot(aes(x=as.factor(day), y = value), outlier.shape = NA) + 
#     ylim(-5, 7) +
#     xlab(" ") + ylab(" ") +
#     theme_classic() +
#     theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank()
#     )


### Consistent data
day1_con <- day1_inc
day2_con <- day2_inc + 0.5
day3_con <- day3_inc - 4.5
day4_con <- day4_inc - 3.5
day5_con <- day5_inc + 2.3

days_con <- c(rep(1, length(day1_con)),
        rep(2, length(day2_con)),
        rep(3, length(day3_con)),
        rep(4, length(day4_con)),
        rep(5, length(day5_con))
        )

df_con <- cbind(c(day1_con, day2_con, day3_con, day4_con, day5_con), as.factor(days_con))
df_con <- as.data.frame(df_con)
colnames(df_con) <- c("value", "day")
data_con <- plot_data(df_con)

# plot2 <- ggplot(df_con) +
#     geom_boxplot(aes(x=as.factor(day), y = value), outlier.shape = NA) + 
#     ylim(-5, 7) +
#     xlab(" ") + ylab(" ") +
#     theme_classic() +
#     theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank()
#     )

# grid.arrange(plot1, plot2, nrow = 1)


p1 <- ggplot(data_inc) +
    geom_segment(aes(x = day, y = LB, xend = day, yend = UB)) + 
    geom_point(aes(x = day, y = point_est)) +
    ylim(-5, 7) +
    xlab("Day") + ylab("Inconsistent Measurements") +
    theme_classic() +
    theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
    )

p2 <- ggplot(data_con) +
    geom_segment(aes(x = day, y = LB, xend = day, yend = UB)) + 
    geom_point(aes(x = day, y = point_est)) +
    ylim(-5, 7) +
    xlab("Day") + ylab("Consistent Measurements") +
    theme_classic() +
    theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
    )


path = "/Users/smt3/Documents/GitHub/atomic-clock/"
figfolder = "DarkUncertaintyAFST/figures/" 

pdf(paste0(path, figfolder, "intro_fig.pdf"), width=8, height=6)
grid.arrange(p1, p2, nrow = 1)
dev.off()

##Reduced chi-square (Methods section of BACON1)
## \chi^{2}_{red} = (1/(N-1))\sum_{i=1}^{N}((x_{i} - \bar{x})^{2}/\sigma_{i}^{2})
chi_sq_red_fn <- function(df){
    N <- length(df[,1])
    res <- ((df$point_est - mean(df$point_est))^2)/df$sd_est
    res <- (1/(N-1))*sum(res)
    return(res)
}

chi_sq_red_fn(data_inc)
chi_sq_red_fn(data_con)
