### Need to source the function.R file###
#source("C:/Users/alber/Dropbox/5. år UNI/PUK Rolf/funktioner.R")
#Make the grid of values
set.seed(44440)
N  <- 10000
S0 <- 1

r_vec     <- runif(N,-0.02,0.05)
delta_vec <- runif(N,-0.01,0.03)
capT_vec  <- runif(N,1/52,2)
K_vec     <- runif(N,0.8,1.2)
sigma_vec <- runif(N,0.05,0.55)

eps         <- .Machine$double.eps
max_digits  <- -log10(eps)

uniroot_sigma <- numeric(N)
NR_sigma      <- numeric(N)
IG_sigma      <- numeric(N)
uniroot_status <- character(N)
NR_status      <- character(N)
IG_status      <- character(N)

obsprice_vec = numeric(N)
start_uniroot_vec = logical(N)
start_IG_vec = logical(N)

for (i in seq_len(N)){
    r <- r_vec[i]
    delta <- delta_vec[i]
    capT <- capT_vec[i]
    K <- K_vec[i]
    sigma <- sigma_vec[i]
    obsprice_vec[i] <- BS_call_sim(S0,0,K,capT,r,delta,sigma)
    
    ## Check if Bisection is able to start. Demands that the value is of opposite signs ##
    f_lo <- IV_func(0.001, obsprice_vec[i], S0,0,K,capT,r,delta,"call")
    f_hi <- IV_func(10,    obsprice_vec[i], S0,0,K,capT,r,delta,"call")
    start_uniroot_vec[i] <- (f_lo*f_hi < 0)
    
    
    Forward <- exp((r-delta)*capT)*S0
    D       <- exp(-r*capT)
    k       <-log(K/Forward)
    m       <- pmin(1, exp(k))
    c= obsprice_vec[i]/(D*Forward)
    
    ### IG demands that 1-c/m  has to be between 0 and 1 
    start_IG_vec[i] <- ((1-c)/m) > 0 & ((1-c)/m)<1
}
    
time_uniroot <- system.time({
  for (i in seq_len(N)){
    r <- r_vec[i]; delta <- delta_vec[i]; capT <- capT_vec[i]
    K <- K_vec[i]; obsprice <- obsprice_vec[i]
    if (start_uniroot_vec[i]){
      res_UR <- tryCatch(list(sig=IV_uniroot(obsprice,S0,0,K,capT,r,delta,"call"), status="attempted"),
                         error=function(e) list(sig=NA_real_, status="Crash"))
    } else res_UR <- list(sig=NA_real_, status="Non-start")
    uniroot_sigma[i] <- res_UR$sig
    uniroot_status[i] <- res_UR$status
  }
})[3] 
    
time_NR <- system.time({
  for (i in seq_len(N)){
    r <- r_vec[i]; delta <- delta_vec[i]; capT <- capT_vec[i]
    K <- K_vec[i]; obsprice <- obsprice_vec[i]
    res_NR <- tryCatch(list(sig=IV_NR(obsprice,S0,0,K,capT,r,delta,"call"), status="attempted"),
                      error=function(e) list(sig=NA_real_, status="Crash"))
    NR_sigma[i] <- res_NR$sig
    NR_status[i] <- res_NR$status}
})[3]

time_IG <- system.time({
  for (i in seq_len(N)){
    r <- r_vec[i]; delta <- delta_vec[i]; capT <- capT_vec[i]
    K <- K_vec[i]; obsprice <- obsprice_vec[i]
    if (start_IG_vec[i]){
      res_IG <- tryCatch(list(sig=IV_IG(obsprice,S0,0,K,capT,r,delta,"call"), status="attempted"),
                      error=function(e) list(sig=NA_real_, status="Crash"))
    } else res_IG <- list(sig=NA_real_, status="Non-start")
    IG_sigma[i] <- res_IG$sig
    IG_status[i] <- res_IG$status
  }
})[3]

time_table = data.frame(uniroot=time_uniroot,NR=time_NR,IG =time_IG, 
                           uniroot_ms_pr_it = time_uniroot/N * 10^6,
                           NR_ms_pr_it=time_NR/N *10^6,
                           IG_ms_pr_it = time_IG/N*10^6)

# Looking at the final status of each of the three methods # 
final_status <- function(status, sigma_hat, sigma_true){
  ifelse(status != "attempted", status,
         ifelse(!is.finite(sigma_hat) | abs(sigma_hat-sigma_true)>1, "Diverged", "Accurate"))
}
uniroot_status <- final_status(uniroot_status, uniroot_sigma, sigma_vec)
NR_status      <- final_status(NR_status,      NR_sigma,      sigma_vec)
IG_status      <- final_status(IG_status,      IG_sigma,      sigma_vec)


# Looking at number of correct digits
digits_correct <- function(sigma_hat, sigma_true){
  rel_err <- abs(sigma_hat - sigma_true) / sigma_true
  ifelse(rel_err <= 0, (max_digits), pmin((-log10(rel_err)),(max_digits)))
}
uniroot_digits <- ifelse(uniroot_status=="Accurate", digits_correct(uniroot_sigma, sigma_vec), NA)
NR_digits      <- ifelse(NR_status=="Accurate",      digits_correct(NR_sigma, sigma_vec),      NA)
IG_digits      <- ifelse(IG_status=="Accurate",      digits_correct(IG_sigma, sigma_vec),      NA)


long_df <- rbind(
  data.frame(method="uniroot", sigma_true=sigma_vec, sigma_hat=uniroot_sigma,
             status=uniroot_status, digits=uniroot_digits),
  data.frame(method="NR",      sigma_true=sigma_vec, sigma_hat=NR_sigma,
             status=NR_status,      digits=NR_digits),
  data.frame(method="IG",      sigma_true=sigma_vec, sigma_hat=IG_sigma,
             status=IG_status,      digits=IG_digits)
)
long_df$status <- factor(long_df$status, levels=c("Accurate","Crash","Non-start","Diverged"))

d1_vec = (log(S0/K_vec)+(r_vec-delta_vec+sigma_vec^2/2) *(capT_vec))/(sigma_vec*sqrt(capT_vec))
Delta_vec = pnorm(d1_vec)
long_df$Delta = rep(Delta_vec,3)
status_summary <- as.data.frame(table(method=long_df$method, status=long_df$status))
status_summary$pct <- 100*status_summary$Freq/N
#   
status_cols <- c(Accurate="grey60",Diverged="blue", Crash="red", `Non-start`="purple")


# 3) Alle punkter, farvet efter status, faceted pr. metode
plot_y <- ifelse(long_df$status=="Accurate", long_df$digits, -5)
long_df$plot_y <- plot_y

p <- ggplot(long_df, aes(x=Delta, y=plot_y, colour=status)) +
  geom_point(alpha=0.5) +
  facet_wrap(~method) +
  geom_hline(yintercept=max_digits, colour="green4", linewidth=0.8) +
  scale_colour_manual(values=status_cols) +
  labs(title="10000 calculations of implied volatility for the 3 methods ",
       x="Delta", y=expression(-log[10](epsilon))) +
  theme_minimal()


p2 <- ggplot(long_df[long_df$status=="Accurate",], aes(x=Delta, y=plot_y)) +
  geom_point(alpha=0.5, color ="grey") +
  facet_wrap(~method) +
  geom_hline(yintercept=max_digits, colour="green4", linewidth=0.8) +
  labs(title="10000 calculations of implied volatility for the 3 methods ",
       x="Delta", y=expression(-log[10](epsilon))) +
  theme_minimal()

p2

table_status_overleaf<- status_summary %>% 
  mutate(
    pct = sprintf("%.2f%%", pct),
    method = recode(method, uniroot = "Uniroot")
  ) %>%
  select(status, method, pct) %>%
  pivot_wider(
    names_from = method,
    values_from = pct
  )
row_speed = data.frame(status ="Iter./microsecond", 
                       IG=time_table$IG_ms_pr_it,
                       NR=time_table$NR_ms_pr_it,
                       Uniroot=time_table$uniroot_ms_pr_it)
table_status_overleaf = rbind(table_status_overleaf,row_speed)
kable(table_status_overleaf,format="latex")


