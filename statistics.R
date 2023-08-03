#!/usr/bin/R


# Bootstrap test statistics for 2-sample comparison
# Algorithm 1 in the paper
critval_twosamplecomparison = function(f1, f2, alpha=0.05, B=999) {

  g_mean = mean(c(f1,f2))
  fuzz1_ = f1 - mean(f1) + g_mean
  fuzz2_ = f2 - mean(f2) + g_mean
  
  Ts = c()
  for (i in 1:B) {
    fuzz1_b = sample(fuzz1_,length(fuzz1_),replace = TRUE)
    fuzz2_b = sample(fuzz2_,length(fuzz2_),replace = TRUE)
    
    var_f1 = var(fuzz1_b)
    var_f2 = var(fuzz2_b)
    s_p = sqrt(var_f1/length(fuzz1_) + var_f2/length(fuzz2_))
    
    T_b = abs(mean(fuzz1_b)-mean(fuzz2_b))/s_p
    Ts[i] = T_b
  }
  Ts_sort = sort(Ts)
  c = Ts_sort[round((B+1)*(1-alpha))]
  c
}


# Test decision for 2-sample comparison
# Algorithm 2 in the paper
dec_twosamplecomparison=function(f1, f2, alpha=0.05, B=999) {
  c = critval_twosamplecomparison (f1,f2,alpha,B)
  var_f1 = var(f1)
  var_f2 = var(f2)
  s_p = sqrt(var_f1/length(f1) + var_f2/length(f2))
  
  T_b = abs(mean(f1)-mean(f2))/s_p

  res = c(T_b,c,(T_b>=c))  
  res
# Output is test statistics, critical value and test decision.
}


# Bootstrap test statistics for >=3-sample comparison ANOVA and an associated posthoc test
# Algorithm 3 in the paper
critval_anova = function(fs, alpha=0.05, B=999) {
  
  fs.an = stack(fs)
  g_mean_n = mean(fs.an$values)
  d <- length(fs)
  fs.an = fs
  for (j in 1:d) {
    fs.an[[j]] = fs[[j]] - mean(fs[[j]]) + g_mean_n
  }
#  fuzz1_3 = fuzz1 - mean(fuzz1) + g_mean_3
#  fuzz2_3 = fuzz2 - mean(fuzz2) + g_mean_3
#  fuzz3_3 = fuzz3 - mean(fuzz3) + g_mean_3
  
  Ts_anova = c()
  Ts_posthoc = c()
  Ts_posthoc_diff = c()
  
  fb = fs
  for (i in 1:B) {
    for (j in 1:d) {
      fb[[j]] = sample(fs.an[[j]],length(fs.an[[j]]),replace = TRUE)
    }
    
    fb_long = stack(fb)
    res = oneway.test(fb_long$values ~ fb_long$ind)
    Ts_anova[i] = res$statistic
    
    diffmax = 0.
    diff_mean_max = 0.
    n = 0
    for (j1 in 1:(d-1)) {
      for (j2 in (j1+1):d) {
        n = n + 1
        res = t.test(fb[[j1]],fb[[j2]])
        diffmax = max(c(diffmax,res$statistic))
        diff_mean_max = max(c(diff_mean_max,abs(mean(fb[[j1]])-mean(fb[[j2]]))))
      }
    }
    Ts_posthoc[i] = diffmax
    
    Ts_posthoc_diff[i] = diff_mean_max
  }
  
  Ts_anova_sort = sort(Ts_anova)
  c_anova = Ts_anova_sort[round((B+1)*(1-alpha))]
  
  Ts_posthoc_sort = sort(Ts_posthoc)
  c_posthoc = Ts_posthoc_sort[round((B+1)*(1-alpha))]
  
  Ts_posthoc_diff_sort = sort(Ts_posthoc_diff)
  c_posthoc_diff = Ts_posthoc_diff_sort[round((B+1)*(1-alpha))]
  res = c(c_anova,c_posthoc,c_posthoc_diff)
  res
  
}


# Test decision for 3-sample comparison ANOVA and an associated posthoc test
# Algorihtms 4 and 5 in the paper
dec_anova=function(fs, alpha=0.05, B=999) {
  c = critval_anova(fs,alpha,B)
  d = length(fs)
  
  fs_long = stack(fs)
  res_data = oneway.test(fs_long$values ~ fs_long$ind)

  test_statistics = res_data$statistic  
  
  decision_anova = (res_data$statistic >= c[1])
  
  decisions_posthoc = data.frame(decisions = rep(0,d*(d-1)/2), test_statistics = rep(0,d*(d-1)/2), sample1 = rep(0,d*(d-1)/2), sample2 = rep(0,d*(d-1)/2))
  decisions_posthoc_mean = data.frame(decisions = rep(0,d*(d-1)/2), test_statistics = rep(0,d*(d-1)/2), sample1 = rep(0,d*(d-1)/2), sample2 = rep(0,d*(d-1)/2))
  index = 0
  for (j1 in 1:(d-1)) {
    for (j2 in (j1+1):d) {
      index = index + 1
      decisions_posthoc[index,] = c((abs(t.test(fs[[j1]],fs[[j2]])$statistic) >= c[2]),t.test(fs[[j1]],fs[[j2]])$statistic,j1,j2)
      decisions_posthoc_mean[index,] = c((abs(mean(fs[[j1]])-mean(fs[[j2]])) >= c[3]),abs(mean(fs[[j1]])-mean(fs[[j2]])),j1,j2)
    }
  }
  res = list(decision_anova = decision_anova, test_statistics_anova = res_data$statistic, critvals = c, decisions_posthoc = decisions_posthoc, decisions_posthoc_mean = decisions_posthoc_mean)
  res
  
  # Output is test decision for the anova, the test statistics, a vector with critical values for the ANOVA and the posthoc test based on the t-test and the 
  # difference of means, and 
  # two data frames which give test decision in first column and the test statistics in the second column for the comparison
  # of samples with numbers given in columns 3 and 4 for the posthoc test based on the t-test and on
  # comparison of means, respectively.
}


# Examples
# ========

# Data input by the user. Here 3 fuzzers in comparison.

fuzz1 = c(19223, 18156, 21213, 18217, 21224, 17204, 19232, 18782, 19881, 20327)
fuzz2 = c(14177, 14731, 14755, 14295, 11202, 14708, 12291, 14578, 17352, 14147)
fuzz3 = c(17127, 18923, 18782, 19221, 15521, 18358, 17845, 15252, 18417, 14484)

# Two-sample comparison

dec_twosamplecomparison(fuzz1,fuzz2)

# Anova

fs = list(fuzz1=fuzz1,fuzz2=fuzz2,fuzz3=fuzz3)
res = dec_anova(fs)
