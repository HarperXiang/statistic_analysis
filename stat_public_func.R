###################################
# Title   : "Statistical Functions"
# Author  : Harper
# Release : 20190915
###################################

# FUNCTION: Sample Statistics
stat_Sample <- function(sample_data) {
  stat_ <- list(
    "size" = length(sample_data), 
    "mean" = mean(sample_data), 
    "var"  = var(sample_data)
  )
  # return the statistics of the Sample
  return (stat_)
}

# FUNCTION: create the Sampling Distribution of the Mean
dis_Mean <- function(sample) {
  samp_stat <- stat_Sample(sample)
  # return the parameters of the sampling distribution
  u   <- samp_stat$mean
  df  <- samp_stat$size - 1
  se  <- sqrt(samp_stat$var / samp_stat$size)
  dis_mean <- list(
    "u"   = u, 
    "df"  = df, 
    "se"  = se
  )
  return (dis_mean)
}

# FUNCTION: create the Sampling Distribution of the Difference between Means
dis_Diff_means <- function(sample_1, sample_2) {
  samp1_stat <- stat_Sample(sample_1)
  samp2_stat <- stat_Sample(sample_2)
  
  df <- samp1_stat$size-1 + samp2_stat$size-1
  if (samp1_stat$size != samp2_stat$size) {
    sse <- sum((sample_1 - samp1_stat$mean)^2) + sum((sample_2 - samp2_stat$mean)^2)
    mse <- sse / df
    n   <- 2 / (1/samp1_stat$size + 1/samp2_stat$size)
  } else {
    mse <- (samp1_stat$var + samp2_stat$var) / 2
    n   <- samp1_stat$size
  }
  se <- sqrt(2*mse/n)
  u  <- samp1_stat$mean - samp2_stat$mean
  # return parameters of the Sampling Distribution
  dis_diff <- list(
    "u"   = u, 
    "df"  = df, 
    "mse" = mse, 
    "n"   = n, 
    "se"  = se
  )
  return (dis_diff)
}


# FUNCTION: t-test on Sampling Distribution of Mean
# Input :
# sample (vector), 
# Null_hypothesis value, one-/two-tailed, direction of one-tailed
# Output:
# t_ {
#   t_test   { # t-test result }
#   dis_mean { # parameters of the sampling distribution }
# }
t_Mean_run <- function(sample, h0_u=0, tailed=2, direction=1) {
  if ((tailed!=1)&&(tailed!=2) || 
      ((direction!=1)&&(direction!=-1)&&(direction!=0))) 
  { return (NULL) }
  
  dis_mean <- dis_Mean(sample)
  t_t <- (dis_mean$u - h0_u) / dis_mean$se
  if (tailed==2) { # Two-tailed
    t_p <- 2 * pt(abs(t_t), dis_mean$df, lower.tail=FALSE)
  } else {         # One-tailed
    if (direction==1) {         # positive direction
      t_p <- pt(t_t, dis_mean$df, lower.tail=FALSE)
    } else if (direction==-1) { # negative direction
      t_p <- pt(t_t, dis_mean$df, lower.tail=TRUE)
    } else {                    # either positive or negative direction
      t_p <- pt(abs(t_t), dis_mean$df, lower.tail=FALSE)
    }
  }
  # return t-test result and the sampling distribution parameters
  t_test <- list(
    "t" = t_t, 
    "p" = t_p
  )
  t_ <- list(
    "t_test"   = t_test, 
    "dis_mean" = dis_mean
  )
  return (t_)
}


# FUNCTION: t-test on Sampling Distribution of the Difference between Means
# Input :
# sample1 (vector), sampe2 (vector), 
# Null_hypothesis value, one-/two-tailed, direction of one-tailed
# Output:
# t_ {
#   t_test   { # t-test result }
#   dis_diff { # parameters of the sampling distribution }
# }
t_Diff_means_run <- function(sample_1, sample_2, h0_diff_u=0, tailed=2, direction=1) {
  if ((tailed!=1)&&(tailed!=2) || 
      ((direction!=1)&&(direction!=-1)&&(direction!=0))) 
  { return (NULL) }
  
  dis_diff <- dis_Diff_means(sample_1, sample_2)
  t_t <- (dis_diff$u - h0_diff_u) / dis_diff$se
  if (tailed==2) { # Two-tailed
    t_p <- 2 * pt(abs(t_t), dis_diff$df, lower.tail=FALSE)
  } else {         # One-tailed
    if (direction==1) {         # positive direction
      t_p <- pt(t_t, dis_diff$df, lower.tail=FALSE)
    } else if (direction==-1) { # negative direction
      t_p <- pt(t_t, dis_diff$df, lower.tail=TRUE)
    } else {                    # either positive or negative direction
      t_p <- pt(abs(t_t), dis_diff$df, lower.tail=FALSE)
    }
  }
  # return t-test result and the sampling distribution parameters
  t_test <- list(
    "t" = t_t, 
    "p" = t_p
  )
  t_ <- list(
    "t_test"   = t_test, 
    "dis_diff" = dis_diff
  )
  return (t_)
}


# FUNCTION: Confidence Interval for Mean (or Means)
# Assuming “homogeneity of variance”; both gender means are normal distributed; and each score is independent
# Input :
# mean, standard error, degree of freedom, Confidence Interval level
# Output:
# CI_ { # the Lower and Upper Critical Values of mean }
CI_Mean <- function(u, se, df, ci=0.95, tailed=2, direction=1) {
  if (((ci<0)||(ci>1)) || 
      (tailed!=1)&&(tailed!=2) || 
      ((direction!=1)&&(direction!=-1)&&(direction!=0))) 
  { return (NULL) }
  
  if (tailed==2) { # Two-tailed
    t_crit <- qt((1-ci)/2, df, lower.tail=FALSE)
    ci_L <- u - se * t_crit
    ci_U <- u + se * t_crit
  } else {         # One-tailed
    if (direction==1) {         # positive direction
      t_crit <- qt(1-ci, df, lower.tail=FALSE)
      ci_L <- u + se * t_crit
      ci_U <- NULL
    } else if (direction==-1) { # negative direction
      t_crit <- qt(ci-1, df, lower.tail=TRUE)
      ci_L <- NULL
      ci_U <- u - se * t_crit
    } else {                    # either positive or negative direction
      t_crit <- qt(1-ci, df, lower.tail=FALSE)
      ci_L <- u - se * t_crit
      ci_U <- u + se * t_crit
    }
  }
  # return the Confidence Interval
  CI_ <- list(
    "t_crit" = t_crit, 
    "ci_L"    = ci_L, 
    "ci_U"    = ci_U 
  )
  return (CI_)
}


# FUNCTION: Confidence Interval calculation for Pearson Correaltion
# Input :
# sample1 (vector), sampe2 (vector), Confidence Interval level
# Output:
# CI_r { # the Lower and Upper limits of r }
CI_Correlation <- function(sample_1, sample_2, ci=0.95) {
  n <- length(sample_1)
  if ((length(sample_2)!=n) || ((ci<0)||(ci>1))) { return (NULL) }
  
  # Correlation r
  r      <- cor(sample_1, sample_2)
  # Fisher’s Transformation to z'
  z_     <- log((1.0+r)/(1.0-r)) / 2.0
  z_se   <- 1.0 / sqrt(n-3.0)
  # z' Critical Value
  z_crit <- qnorm((1.0-ci)/2, mean=0, sd=1, lower.tail=FALSE)
  # CI on z'
  z_L <- z_ - z_crit * z_se
  z_U <- z_ + z_crit * z_se
  # Fisher's Transformation back to r
  r_L <- exp(2*z_L)
  r_L <- (r_L-1)/(r_L+1)
  r_U <- exp(2*z_U)
  r_U <- (r_U-1)/(r_U+1)
  
  CI_r <- list(
    "r_L" = r_L, 
    "r_U" = r_U
  )
  return (CI_r)
}

