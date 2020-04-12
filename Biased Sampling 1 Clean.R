### DDI as Function of Sampling Scheme

library(ggplot2)

## Parameters

#set.seed(2020)
N <- 10^3 # ease burden
mu <- 1
sigma <- 1
alpha_dist <- 0.1
dev_alpha <- 3
n_min <- 100
n_dist <- 100
n_rep <- 50 # rep for smoothing (was 100)
n_max_var <- 50 # new coefficient

## Storage Setup

alpha_vec <- seq(-dev_alpha, dev_alpha, alpha_dist)
n_vec <- seq(n_min, N - n_min, n_dist)

result_df1 <- expand.grid(alpha_vec, n_vec)
res_rows <- nrow(result_df1)
result_df <- cbind(result_df1, as.data.frame(matrix(0, ncol=4, nrow=res_rows)))
colnames(result_df) <- c("alpha", "n", "ddi", "error", "max_var", "di_do")

## Simulation

for(i in 1:res_rows){
  
  error_inter <- numeric(n_rep)
  corr_inter <- numeric(n_rep)
  max_var_inter <- numeric(n_rep) # new thing added
  
  for(j in 1:n_rep){
    
    # Data Generation
    
    pop_data <- rexp(N, mu)
    prob_vec <- pop_data^(result_df$alpha[i])
    
    samp_data <- sample(pop_data, result_df$n[i], replace=F, prob=prob_vec)
    samp_ind <- pop_data %in% samp_data
    
    # Compute Error/Corr
    
    y_star <- max(pop_data)
    w_tot = sum(prob_vec)
    
    corr_inter[j] <- cor(samp_ind, pop_data)
    error_inter[j] <- y_star - max(samp_data) # was max
    
    y_rand <- sample(pop_data, n_max_var, replace=F)
    p_t <- numeric(n_max_var)
    for(k in n_max_var){
      lt_vec <- y_star - pop_data > y_rand[k]
      p_t[k] <- sum(prob_vec[lt_vec]/w_tot)
    }
    
    max_var_inter[j] <- mean(p_t)
  }
  
  result_df$error[i] <- mean(error_inter^2) # MSE for straight error for easier comparison
  result_df$ddi[i] <- mean(corr_inter^2)
  result_df$max_var[i] <- mean(result_df$n[i] * max_var_inter) # based on heuristic calc
}

result_df$di_do <- (N - result_df$n)/result_df$n * result_df$ddi

ggplot(data=result_df, mapping=aes(max_var, error)) + geom_point()

# ggplot(data=result_df, mapping=aes(di_do, error)) + geom_point() +
#   labs(x="D_I * D_O", y="MSE of Max") +
#   ggtitle("MSE of Mean versus MSE of Max") +
#   theme(plot.title = element_text(hjust = 0.5))

# ggplot(data=result_df[result_df$alpha==2,], mapping=aes(n, error)) + geom_point() +
#   labs(x="D_I * D_O", y="MSE of Max") +
#   ggtitle("Sample Size vs. MSE of Max") +
#   theme(plot.title = element_text(hjust = 0.5))

# ggplot(data=result_df, mapping=aes(alpha, n/N, fill=ddi)) + geom_tile() + 
#   labs(x="Alpha", y = "Sampling Frequency") +
#   ggtitle("Heatplot of DDI vs. Sampling Scheme") + 
#   theme(plot.title = element_text(hjust = 0.5))

# ggplot(data=result_df[result_df$n == 200,], mapping=aes(alpha, ddi)) + geom_point()
# ggplot(data=result_df, mapping=aes(ddi, error)) + geom_point()