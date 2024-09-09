rm(list = ls())

library(ggplot2)
library(MASS)
library(viridis)
library(gifski)
library(ggExtra)

# Define the target distribution (Mixture of two bivariate Gaussians)
target_distribution <- function(x, y) {
  mu1 <- c(2.5, 1)
  mu2 <- c(-2.5, -1)
  sigma <- matrix(c(2, 0, 0, 3), nrow = 2)  
  
  p1 <- 0.5 * dmvnorm(cbind(x, y), mean = mu1, sigma = sigma)
  p2 <- 0.5 * dmvnorm(cbind(x, y), mean = mu2, sigma = sigma)
  
  return(p1 + p2)
}

# Calculate theoretical marginal densities
# marginal_density_x <- function(x) {
#   mu1 <- 1
#   mu2 <- -1
#   sigma <- 1
  
#   p1 <- dnorm(x, mean = mu1, sd = sigma)
#   p2 <- dnorm(x, mean = mu2, sd = sigma)
  
#   return(0.5 * p1 + 0.5 * p2)
# }

# marginal_density_y <- function(y) {
#   mu1 <- 1
#   mu2 <- -1
#   sigma <- 1
  
#   p1 <- dnorm(y, mean = mu1, sd = sigma)
#   p2 <- dnorm(y, mean = mu2, sd = sigma)
  
#   return(0.5 * p1 + 0.5 * p2)
# }

# Create a grid over x and y for the contour plot
x_seq <- seq(-6, 6, length.out = 100)
y_seq <- seq(-6, 6, length.out = 100)
grid <- expand.grid(x = x_seq, y = y_seq)
grid$z <- target_distribution(grid$x, grid$y)

# Metropolis-Hastings algorithm 
metropolis_hastings_bivariate <- function(n_iter, init_val, proposal_sd) {
  samples <- matrix(NA, nrow = n_iter, ncol = 2)  # Store bivariate samples
  samples[1, ] <- init_val  # Initial value
  
  for (i in 2:n_iter) {
    # Propose new values for both dimensions
    proposed_val <- mvrnorm(1, mu = samples[i - 1, ], Sigma = proposal_sd * diag(2))
    
    # Compute acceptance ratio
    numerator <- target_distribution(proposed_val[1], proposed_val[2])
    denominator <- target_distribution(samples[i - 1, 1], samples[i - 1, 2])
    acceptance_ratio <- numerator / denominator
    
    # Accept or reject the new value
    if (runif(1) < acceptance_ratio) {
      samples[i, ] <- proposed_val
    } else {
      samples[i, ] <- samples[i - 1, ]
    }
  }
  
  return(samples)
}

# Set MCMC parameters
n_iter <- 1000  # Number of iterations
initial_value <- c(0, 0)
proposal_sd <- 0.5  # Proposal distribution's standard deviation

# Run the Metropolis-Hastings algorithm
mcmc_samples <- metropolis_hastings_bivariate(n_iter, initial_value, proposal_sd)

# Create data frame for the samples
samples_df <- data.frame(
  x = mcmc_samples[, 1],
  y = mcmc_samples[, 2],
  iteration = 1:n_iter
)

# Calculate limits for the marginal histograms
x_limits <- range(grid$x)
y_limits <- range(grid$y)

# Create an empty list to store frame file names
frame_files <- character(n_iter)

# Loop over each iteration to create frames
for (i in 1:n_iter) {
  # Samples up to iteration i
  samples_i <- samples_df[1:i, ]
  
  # Main plot with contour, sample points, and path
  main_plot <- ggplot(samples_i, aes(x = x, y = y)) +
    geom_raster(data = grid, aes(x = x, y = y, fill = z)) +
    geom_contour(data = grid, aes(x = x, y = y, z = z), bins = 10, color = '#851e3e') +
    scale_fill_gradient(low = "#eee3e7", high = "#f6abb6") +  # Gradient with midpoint
    geom_path(color = "#de70ec", size = 0.8) +
    geom_point(color = "#4c005a", size = 1.2) +  # Add path of consecutive points
    labs(title = "Markov Chain Monte Carlo Visualization",
         subtitle = paste("Iteration:", i), x = "X", y = "Y") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  
  # Add marginal histograms
  combined_plot <- ggMarginal(main_plot, type = "histogram", fill = "#63ace5", alpha = 0.6, size = 5)
  
  # Save the frame to a temporary file
  frame_file <- sprintf("frame_%03d.png", i)
  ggsave(frame_file, combined_plot, width = 6, height = 6, dpi = 300, bg = "white")
  
  # Store the file name
  frame_files[i] <- frame_file
  print(i)
}

# Combine frames into a GIF
fps <- 20  # Frames per second
gifski(png_files = frame_files, gif_file = "mcmc_animation.gif", delay = 1/fps, width = 600, height = 600)

# Clean up temporary files
file.remove(frame_files)
