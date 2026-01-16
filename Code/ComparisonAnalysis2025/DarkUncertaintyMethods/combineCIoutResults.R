rm(list=ls())

library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
########### make CIout

collectRuns=function(filepath,runDate,howDatSim){
  allRunDateFiles=list.files(filepath, pattern=runDate)
  
  CIoutFiles=grep(paste(howDatSim,"CIout",runDate,sep=""),allRunDateFiles)
  
  CIout=read.csv(paste(filepath,allRunDateFiles[CIoutFiles[1]],sep=""))#,row.names = 1)
  
  for(i in CIoutFiles[-1]){
    tmp=read.csv(paste(filepath,allRunDateFiles[i],sep=""))#,row.names = 1)
    dim(tmp)
    CIout=bind_rows(CIout,tmp)
  }

  print(table(CIout$Iteration))
  
  return(CIout)
  
}

myfilepath="/home/aak3/NIST/atomic-clock/Code/ComparisonAnalysis2025/DarkUncertaintyMethods/CIouts/"

###################
CIout=collectRuns(filepath = myfilepath,runDate = "251208",howDatSim = "normal") 

################################
## How are runs
################################

checkBadRuns=function(CIout){
  bad1=dim(filter(CIout,parameter=="mu",min.bfmi<.3))[1] 
  bad2=dim(filter(CIout,parameter=="mu",min.bfmi<.3 | min.neff <300))[1]
  bad3=dim(filter(CIout,parameter=="mu",min.bfmi<.3 | min.neff <300 | bad.rhat > 0))[1]
  bad4=dim(filter(CIout,parameter=="mu",min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0))[1]
  bad5=dim(filter(CIout,parameter=="mu",min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0))[1]
  
  print(paste("BFMI",bad1))
  print(paste("neff",bad2-bad1))
  print(paste("rhat",bad3-bad2))
  print(paste("div",bad4-bad3))
  print(paste("tree",bad5-bad4))
  print(paste("total dropped",bad5))
}

checkBadRuns(CIout)

CIout=CIout %>% group_by(method,parameter,tau,N) %>% 
  mutate(Repeat=1:n())

CIoutRepVal=max(CIout$Repeat)

bad.CIout=filter(CIout,min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0)
final.CIout=filter(CIout,is.na(min.bfmi) | !(min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0))

print((dim(final.CIout)[1]+dim(bad.CIout)[1]) == dim(CIout)[1])

keptPercent=final.CIout %>% group_by(method,parameter,tau,N) %>% 
  summarise(kept=n()/CIoutRepVal)
View(keptPercent)

final.CIout$method = factor(final.CIout$method, 
                            levels=c("DL","SJ","REML","PM",
                                     "GaussGauss","Birge Ratio (Multiplicative)", "Dark Uncertainty (Additive)"))
outfilepath=paste(myfilepath,"summaryPlots/",sep="")


the.N.vals=sort(unique(final.CIout$N))

the.tau.vals=sort(unique(final.CIout$tau))

dets="YbSr"
#######proportional bias

pdf(paste(outfilepath,"proportionalBiasFor",final.CIout$simulation[1],
          "SimDataWith",max(final.CIout$Repeat),"Repeats",dets,".pdf",sep=""),width = 7,height = 5.5)

for (aTau in the.tau.vals){
  ######### For tau
  propBias=final.CIout %>% group_by(method,parameter,tau,N) %>%
    filter(parameter=="tau")%>%
    mutate(proportionalBias=(estimate-tau)/tau*100)
  
  AvgPropBias=propBias %>% 
    summarize(AvgProportionalBias=mean(proportionalBias))
  
  print(ggplot(filter(AvgPropBias,tau == aTau),aes(N,AvgProportionalBias,col=method))+
          geom_point()+
          geom_line(aes(group=method))+
          ggtitle(paste("tau, ","true tau=",aTau,sep=""))+
          geom_hline(yintercept = 0)
        # scale_x_continuous(trans='log')
  )
  
  propBias$N = as.factor(propBias$N)
  
  print(ggplot(filter(propBias,tau == aTau),aes(N,proportionalBias))+
          geom_violin()+
          stat_summary(fun=mean, geom="point",color="cyan")+
          facet_wrap(~method,ncol=4)+
          ggtitle(paste("tau, ","true tau=",aTau,sep=""))+
          geom_hline(yintercept = 0)
  )
  
  ######### For mu
  propBiasMu=final.CIout %>% group_by(method,parameter,tau,N) %>%
    filter(parameter=="mu")%>%
    mutate(proportionalBias=(estimate-mu)/mu*100)
  
  AvgPropBiasMu=propBiasMu %>% 
    summarize(AvgProportionalBias=mean(proportionalBias))
  
  print(ggplot(filter(AvgPropBiasMu,tau == aTau),aes(N,AvgProportionalBias,col=method))+
          geom_point()+
          geom_line(aes(group=method))+
          ggtitle(paste("mu, ","true tau=",aTau,sep=""))+
          geom_hline(yintercept = 0)
  )
  
  propBiasMu$N = as.factor(propBiasMu$N)
  
  print(ggplot(filter(propBiasMu,tau == aTau),aes(N,proportionalBias))+
          geom_violin()+
          stat_summary(fun=mean, geom="point",color="cyan")+
          facet_wrap(~method,ncol=4)+
          ggtitle(paste("mu, ","true tau=",aTau,sep=""))+
          geom_hline(yintercept = 0)
  )
  
}

dev.off()




propBiasViolinPlotsForSelectN=function(CIout,Nvals,Nlabs,tauVals,tauLabs,dets,plotht,plotwd){
  
  CIout=CIout %>% group_by(method,parameter,tau,N) %>% 
    mutate(Repeat=1:n())
  
  CIoutRepVal=max(CIout$Repeat)
  
  bad.CIout=filter(CIout,min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0)
  final.CIout=filter(CIout,is.na(min.bfmi) | !(min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0))
  
  keptPercent=final.CIout %>% group_by(method,parameter,tau,N) %>% 
    summarise(kept=n()/CIoutRepVal)
  
  final.CIout$method = factor(final.CIout$method, 
                              levels=c("DL","SJ","REML","PM",
                                       "GaussGauss","Birge Ratio (Multiplicative)", "Dark Uncertainty (Additive)"))
  outfilepath=paste(myfilepath,"summaryPlots/",sep="")
  
  the.tau.vals=sort(unique(final.CIout$tau))
  
  #######focus on small N
  
  pdf(paste(outfilepath,"propBiasViolinPlotsForNis",paste(Nvals,collapse = "_"),"and",final.CIout$simulation[1],
            "SimDataWith",max(final.CIout$Repeat),"Repeats",dets,".pdf",sep=""),width = plotwd,height = plotht)
  
  
  ######### For mu
  propBiasMu=final.CIout %>% group_by(method,parameter,tau,N) %>%
    filter(parameter=="mu")%>%
    mutate(proportionalBias=(estimate-mu)/mu*100)
  propBiasMu$N = as.factor(propBiasMu$N)
  
  propBiasMu$tau <- factor(propBiasMu$tau, levels = tauVals,
                           ordered = TRUE, labels=tauLabs)
  propBiasMu$Nlab <- factor(propBiasMu$N, levels = Nvals,
                            ordered = TRUE, labels=Nlabs)
  
  
  plot1=ggplot(filter(propBiasMu, N %in% Nvals),aes(method,proportionalBias))+
    geom_violin()+
    stat_summary(fun=mean, geom="point",color="cyan")+
    ggtitle(expression(mu))+
    geom_hline(yintercept = 0)+
    facet_wrap(~Nlab+tau,labeller = labeller(.cols = label_parsed, .multi_line = FALSE),ncol=2)+
    theme(axis.text.x=element_text(angle = 45,hjust=1))
  ######### For tau
  propBias=final.CIout %>% group_by(method,parameter,tau,N) %>%
    filter(parameter=="tau")%>%
    mutate(proportionalBias=(estimate-tau)/tau*100)
  propBias$N = as.factor(propBias$N)
  
  propBias$tau <- factor(propBias$tau, levels = tauVals,
                         ordered = TRUE, labels=tauLabs)
  propBias$Nlab <- factor(propBias$N, levels = Nvals,
                          ordered = TRUE, labels=Nlabs)
  
  
  plot2=ggplot(filter(propBias, N %in% Nvals),aes(method,proportionalBias))+
    geom_violin()+
    stat_summary(fun=mean, geom="point",color="cyan")+
    ggtitle(expression(tau))+
    geom_hline(yintercept = 0)+
    facet_wrap(~Nlab+tau,labeller = labeller(.cols = label_parsed, .multi_line = FALSE),ncol=2)+
    theme(axis.text.x=element_text(angle = 45,hjust=1))
  print(grid.arrange(plot1, plot2, ncol=2))
  
  # }
  
  dev.off()
  
}



# N_set=c(13,13+20,100)
# tau_set=c(3,10)



tauValz=c(3,10)
tauLabz=c(expression(paste(tau,"=3")),expression(paste(tau,"=10")))

Nvalz=c(13,33,100)
Nlabz=c(expression(paste("N=13")), expression(paste("N=33")),expression(paste("N=100")))

propBiasViolinPlotsForSelectN(CIout,Nvals = Nvalz,
                              Nlabs = Nlabz,tauVals = tauValz,tauLabs = tauLabz,"TEST",
                              plotht=5,plotwd=8)









##################################



GeomSplitViolin <- ggplot2::ggproto("GeomSplitViolin", ggplot2::GeomViolin, 
                                    draw_group = function(self, data, ..., draw_quantiles = NULL) {
                                      data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                                      grp <- data[1, "group"]
                                      newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                                      newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                                      newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                                      
                                      if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                                        quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                                        aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                                        aesthetics$alpha <- rep(1, nrow(quantiles))
                                        both <- cbind(quantiles, aesthetics)
                                        quantile_grob <- ggplot2::GeomPath$draw_panel(both, ...)
                                        ggplot2:::ggname("geom_split_violin", grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                                      }
                                      else {
                                        ggplot2:::ggname("geom_split_violin", ggplot2::GeomPolygon$draw_panel(newdata, ...))
                                      }
                                    }
)

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
                 position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
                 params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
# ----------------------------------------------------------------------


library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

# --- Define the GeomSplitViolin and geom_split_violin (MUST BE RUN) ---
# Assuming you have these definitions from the previous response:
# GeomSplitViolin <- ggplot2::ggproto(...)
# geom_split_violin <- function(...)
# ... (The definitions for GeomSplitViolin and geom_split_violin go here) ...
# I will omit them here for brevity but ensure they are run in your session.

## 🎨 Modified Function: Color by Method
propBiasSplitViolinForTwoMethods_ColorByMethod <- function(CIout, Nvals, Nlabs, tauVals, tauLabs, methodsToCompare, dets, plotht, plotwd, myfilepath) {
  
  if (length(methodsToCompare) != 2) {
    stop("Please provide exactly two methods to compare in 'methodsToCompare'.")
  }
  
  CIout_filtered <- CIout %>%
    filter(method %in% methodsToCompare) %>%
    group_by(method, parameter, tau, N) %>%
    mutate(Repeat = 1:n())
  
  # Filter for good runs
  final.CIout <- CIout_filtered %>%
    filter(is.na(min.bfmi) | !(min.bfmi < .3 | min.neff < 300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex > 0))
  
  # Ensure the methods are factors
  final.CIout$method <- factor(final.CIout$method, 
                               levels = methodsToCompare)
  
  # Calculate proportional bias for both 'mu' and 'tau'
  propBias_data <- final.CIout %>%
    group_by(method, parameter, tau, N) %>%
    mutate(
      proportionalBias = case_when(
        parameter == "mu" ~ (estimate - mu) / mu * 100,
        parameter == "tau" ~ (estimate - tau) / tau * 100,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(proportionalBias))
  
  # Factor N and tau for plotting
  propBias_data$Nlab <- factor(propBias_data$N, levels = Nvals, ordered = TRUE, labels = Nlabs)
  propBias_data$tau <- factor(propBias_data$tau, levels = tauVals, ordered = TRUE, labels = tauLabs)
  
  # Set the order of parameters for the X-axis/split categories
  propBias_data$parameter_lab <- factor(propBias_data$parameter, levels = c("mu", "tau"), 
                                        labels = c(expression(mu), expression(tau)))
  
  plot_data <- filter(propBias_data, N %in% Nvals)
  
  pdf(paste(myfilepath, "SplitViolinBiasPlot_ColorByMethod_Methods", paste(methodsToCompare, collapse = "_vs_"), ".pdf", sep = ""), 
      width = plotwd, height = plotht)
  
  # 5. Create the plot with swapped aesthetics 
  plot_split <- ggplot(plot_data, 
                       aes(x = parameter_lab,        # X-axis is now the Parameter (mu/tau)
                           y = proportionalBias, 
                           fill = method,            # Fill/Color is now the Method
                           color = method)) +
    geom_split_violin() +
    # Use position_dodge to ensure means align correctly on the split violin
    stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.5), size = 2,color="black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ Nlab + tau+parameter,scales="free", 
               labeller = labeller(.cols = label_parsed, .multi_line = FALSE), 
               ncol = 2) +
    labs(
      x = "Parameter",
      y = "Proportional Bias (%)",
      fill = "Method",
      color = "Method",
      title = paste("Proportional Bias Comparison by Method:", methodsToCompare[1], "vs.", methodsToCompare[2])
    ) +
    # Optional: Customize colors
    # scale_fill_manual(values = c("DL" = "blue", "PM" = "red")) +
    # scale_color_manual(values = c("DL" = "darkblue", "PM" = "darkred")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0), # Simplified X-axis labels (mu/tau)
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  print(plot_split)
  
  dev.off()
}

propBiasSplitViolinForTwoMethods_ColorByMethod(
  CIout = CIout,
  Nvals = Nvalz,
  Nlabs = Nlabz,
  tauVals = tauValz,
  tauLabs = tauLabz,
  methodsToCompare = c("DL", "Dark Uncertainty (Additive)"), # <--- Specify the two methods here
  dets = NA,
  plotht = 8,
  plotwd = 10,
  myfilepath = outfilepath
)

rawMuSplitViolinPlot(
  CIout = CIout,
  Nvals = Nvalz,
  Nlabs = Nlabz,
  tauVals = tauValz,
  tauLabs = tauLabz,
  methodsToCompare = c("DL", "Dark Uncertainty (Additive)"), # <--- Specify the two methods here
  dets = NA,
  plotht = 8,
  plotwd = 10,
  myfilepath = outfilepath
)
  


# --- Ensure GeomSplitViolin and geom_split_violin are defined in your session ---

## 🛠️ Function for Split Violin Plot of Raw Mu Estimates
rawMuSplitViolinPlot <- function(CIout, Nvals, Nlabs, tauVals, tauLabs, methodsToCompare, dets, plotht, plotwd, myfilepath) {
  
  if (length(methodsToCompare) != 2) {
    stop("Please provide exactly two methods to compare in 'methodsToCompare'.")
  }
  
  # 1. Filter and Clean Data for mu
  plot_data <- CIout %>%
    filter(method %in% methodsToCompare, parameter == "mu") %>%
    group_by(method, parameter, tau, N) %>%
    mutate(Repeat = 1:n()) %>%
    filter(is.na(min.bfmi) | !(min.bfmi < .3 | min.neff < 300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex > 0))
  
  # Ensure the methods are factors
  plot_data$method <- factor(plot_data$method, levels = methodsToCompare)
  
  # 2. Factor N and tau for plotting
  plot_data$Nlab <- factor(plot_data$N, levels = Nvals, ordered = TRUE, labels = Nlabs)
  plot_data$tau <- factor(plot_data$tau, levels = tauVals, ordered = TRUE, labels = tauLabs)
  
  # 3. Create a common x-axis variable for splitting (we just use a placeholder)
  plot_data$x_group <- factor("Mu Estimate")
  
  plot_data_filtered <- filter(plot_data, N %in% Nvals)
  
  pdf(paste(myfilepath, "SplitViolin_RawMuEstimates_Methods", paste(methodsToCompare, collapse = "_vs_"), ".pdf", sep = ""), 
      width = plotwd, height = plotht)
  
  # 4. Create the Plot
  plot_mu <- ggplot(plot_data_filtered, 
                    aes(x = x_group, 
                        y = estimate, 
                        fill = method, 
                        color = method)) +
    geom_split_violin(alpha = 0.6) +
    # Add mean and median
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", 
                 position = position_dodge(width = 0.5), show.legend = FALSE) +
    stat_summary(fun = median, geom = "point", shape = 18, size = 3, 
                 position = position_dodge(width = 0.5), show.legend = FALSE) +
    
    # Add true value line, mapping 'mu' to the y-intercept
    geom_hline(aes(yintercept = mu), linetype = "dashed", color = "red", data = distinct(plot_data_filtered, N, tau, mu)) +
    
    facet_wrap(~ Nlab + tau, 
               labeller = labeller(.cols = label_parsed, .multi_line = FALSE), 
               ncol = 2, scales = "free_y") +
    labs(
      x = "",
      y = expression(paste("Estimate of Mean (", hat(mu), ")")),
      fill = "Method",
      color = "Method",
      title = paste("Distribution of Raw", expression(hat(mu)), "Estimates:", methodsToCompare[1], "vs.", methodsToCompare[2])
    ) +
    theme_bw() +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  print(plot_mu)
  
  dev.off()
}

## 🛠️ Function for Split Violin Plot of Raw Tau Estimates (Final Fix)
rawTauSplitViolinPlot <- function(CIout, Nvals, Nlabs, tauVals, tauLabs, methodsToCompare, dets, plotht, plotwd, myfilepath) {
  
  if (length(methodsToCompare) != 2) {
    stop("Please provide exactly two methods to compare in 'methodsToCompare'.")
  }
  
  # 1. Prepare and Filter Data
  plot_data <- CIout %>%
    filter(method %in% methodsToCompare, parameter == "tau") %>%
    group_by(method, parameter, tau, N) %>%
    mutate(Repeat = 1:n()) %>%
    filter(is.na(min.bfmi) | !(min.bfmi < .3 | min.neff < 300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex > 0))
  
  # Ensure the methods are factors
  plot_data$method <- factor(plot_data$method, levels = methodsToCompare)
  
  # 2. Factor N and tau for plotting
  plot_data$Nlab <- factor(plot_data$N, levels = Nvals, ordered = TRUE, labels = Nlabs)
  
  # Store the numeric tau value before making the factor used for facetting
  plot_data$tau_numeric <- plot_data$tau 
  
  # Factor tau for facetting labels
  plot_data$tau <- factor(plot_data$tau, levels = tauVals, ordered = TRUE, labels = tauLabs)
  
  # 3. Create a common x-axis variable for splitting
  plot_data$x_group <- factor("Tau Estimate")
  
  plot_data_filtered <- filter(plot_data, N %in% Nvals)
  
  # 4. --- Data for geom_hline (Crucial Fix) ---
  # We extract the true numeric tau value and the factor labels used for facetting.
  # This ensures the numeric intercept value is matched correctly to the facet label.
  true_tau_data <- plot_data_filtered %>%
    distinct(Nlab, tau, tau_numeric) %>%
    rename(yintercept = tau_numeric)
  # -------------------------------------------
  
  pdf(paste(myfilepath, "SplitViolin_RawTauEstimates_Methods", paste(methodsToCompare, collapse = "_vs_"), ".pdf", sep = ""), 
      width = plotwd, height = plotwd)
  
  # 5. Create the Plot
  plot_tau <- ggplot(plot_data_filtered, 
                     aes(x = x_group, 
                         y = estimate, 
                         fill = method, 
                         color = method)) +
    geom_split_violin(alpha = 0.6) +
    
    # Add mean and median
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", 
                 position = position_dodge(width = 0.5), show.legend = T) +
    stat_summary(fun = median, geom = "point", shape = 18, size = 3, 
                 position = position_dodge(width = 0.5), show.legend = T) +
    
    # Use the new true_tau_data. ggplot automatically matches the Nlab and tau factors
    # to the facet panels and uses the yintercept value.
    geom_hline(aes(yintercept = yintercept), linetype = "dashed", color = "red", data = true_tau_data) +
    
    # Facet by the factor variables Nlab and tau
    facet_wrap(~ Nlab + tau, 
               labeller = labeller(.cols = label_parsed, .multi_line = FALSE), 
               ncol = 2, scales = "free_y") +
    
    labs(
      x = "",
      y = expression(paste("Estimate of Heterogeneity (", hat(tau), ")")),
      fill = "Method",
      color = "Method",
      title = paste("Distribution of Raw", expression(hat(tau)), "Estimates:", methodsToCompare[1], "vs.", methodsToCompare[2])
    ) +
    theme_bw() +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  print(plot_tau)
  
  dev.off()
}
rawMuSplitViolinPlot(
  CIout = CIout,
  Nvals = Nvalz,
  Nlabs = Nlabz,
  tauVals = tauValz,
  tauLabs = tauLabz,
  methodsToCompare = c("DL", "Dark Uncertainty (Additive)"), # <--- Specify the two methods here
  dets = NA,
  plotht = 8,
  plotwd = 10,
  myfilepath = outfilepath
)

rawTauSplitViolinPlot(
  CIout = CIout,
  Nvals = Nvalz,
  Nlabs = Nlabz,
  tauVals = tauValz,
  tauLabs = tauLabz,
  methodsToCompare = c("DL", "Dark Uncertainty (Additive)"), # <--- Specify the two methods here
  dets = NA,
  plotht = 8,
  plotwd = 10,
  myfilepath = outfilepath
)











# 
# 
# 
# 
# 
# processResults=function(CIout,dets){
#   
#   checkBadRuns(CIout)
#   
#   CIout=CIout %>% group_by(method,parameter,tau,N) %>% 
#     mutate(Repeat=1:n())
# 
#   CIoutRepVal=max(CIout$Repeat)
#     
#   bad.CIout=filter(CIout,min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0)
#   # View(bad.CIout)
#   # test=filter(CIout,min.bfmi >=.3 & min.neff < 300 & bad.rhat == 0 & num.Divergent == 0 & tree.Depth.Ex == 0)
#   # View(test)
#   # good.CIout=filter(CIout,min.bfmi >=.3 & min.neff >= 300 & bad.rhat == 0 & num.Divergent == 0 & tree.Depth.Ex == 0)
#   final.CIout=filter(CIout,is.na(min.bfmi) | !(min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0))
#   
#   print((dim(final.CIout)[1]+dim(bad.CIout)[1]) == dim(CIout)[1])
#   
#   keptPercent=final.CIout %>% group_by(method,parameter,tau,N) %>% 
#     summarise(kept=n()/CIoutRepVal)
#   View(keptPercent)
#   
#   final.CIout$method = factor(final.CIout$method, 
#                                     levels=c("DL","SJ","REML","PM",
#                                              "LaplaceLaplace","GaussGauss","SkewtGaussian","LaplaceGaussian"))
#   outfilepath=myfilepath
#   
#   the.N.vals=sort(unique(final.CIout$N))
#   the.tau.vals=sort(unique(final.CIout$tau))
#   
#   pdf(paste(outfilepath,"intervalsFor",final.CIout$simulation[1],
#             "SimDataWith",max(final.CIout$Repeat),"Repeats",dets,".pdf",sep=""),width = 8.5,height = 11)
#   
#   for(aN in the.N.vals){
#     for (aTau in the.tau.vals){
#       print(aN)
#       print(aTau)
#       print(ggplot(filter(final.CIout,parameter=="mu" & N==aN & tau == aTau),aes(Repeat,estimate,ymin=ci.lb,ymax=ci.ub,col=inInt))+
#               geom_point()+
#               geom_errorbar()+
#               facet_wrap(~method,ncol = 2)+
#               geom_hline(aes(yintercept=mu))+
#               ggtitle(paste("mu, N=",aN,", tau=",aTau,sep="")))
#       if(dim(filter(final.CIout,parameter=="mu" & N==aN & tau==aTau & !inInt))[1]!=0){
#         print(ggplot(filter(final.CIout,parameter=="mu" & N==aN & tau == aTau & !inInt),aes(Repeat,estimate,ymin=ci.lb,ymax=ci.ub,col=inInt))+
#                 geom_point()+
#                 geom_errorbar()+
#                 facet_wrap(~method,ncol = 2)+
#                 geom_hline(aes(yintercept=mu))+
#                 ggtitle(paste("mu, N=",aN,", tau=",aTau,sep="")))
#         
#       }
#       
#       print(ggplot(filter(final.CIout,parameter=="tau" & N==aN & tau==aTau),aes(Repeat,estimate,ymin=ci.lb,ymax=ci.ub,col=inInt))+
#               geom_point()+
#               geom_errorbar()+
#               facet_wrap(~method,ncol = 2)+
#               geom_hline(aes(yintercept=tau))+
#               ggtitle(paste("tau, N=",aN,", tau=",aTau,sep="")))
#       if(dim(filter(final.CIout,parameter=="tau" & N==aN & tau==aTau & !inInt))[1]!=0){
#         print(ggplot(filter(final.CIout,parameter=="tau" & N==aN & tau==aTau & !inInt),aes(Repeat,estimate,ymin=ci.lb,ymax=ci.ub,col=inInt))+
#                 geom_point()+
#                 geom_errorbar()+
#                 facet_wrap(~method,ncol = 2)+
#                 geom_hline(aes(yintercept=tau))+
#                 ggtitle(paste("tau, N=",aN,", tau=",aTau,sep="")))
#       }
#       
#     }
#   }
#   
#   dev.off()
#   
#   ########### coverage probs
#   
#   CPsum=final.CIout %>% group_by(method,parameter,tau,N) %>% 
#     summarise(numberOfSims=n(),coverage=sum(inInt)/numberOfSims)
#   
#   pdf(paste(outfilepath,"coverageProbsFor",final.CIout$simulation[1],
#             "SimDataWith",max(final.CIout$Repeat),"Repeats",dets,".pdf",sep=""),width = 7,height = 5.5)
#   
#   for (aTau in the.tau.vals){
#     print(
#       ggplot(filter(CPsum,parameter=="mu" & tau == aTau),aes(N,coverage))+
#         geom_point()+
#         facet_wrap(~method,ncol=4)+
#         geom_hline(aes(yintercept=.95))+
#         ggtitle(paste("mu, tau=",aTau,sep=""))
#     )
#     
#     print(
#       ggplot(filter(CPsum,parameter=="tau" & tau == aTau),aes(N,coverage))+
#         geom_point()+
#         facet_wrap(~method,ncol=4)+
#         geom_hline(aes(yintercept=.95))+
#         ggtitle(paste("tau, tau=",aTau,sep=""))
#     )
#     
#     #######all plotted together, 4 rma results will overlap
#     print(
#       ggplot(filter(CPsum,parameter=="mu" & tau == aTau),aes(N,coverage,col=method))+
#         geom_point()+
#         geom_line(aes(group=method))+
#         geom_hline(aes(yintercept=.95))+
#         ggtitle(paste("mu, tau=",aTau,sep=""))
#     )
#     
#     print(
#       ggplot(filter(CPsum,parameter=="tau" & tau == aTau),aes(N,coverage,col=method))+
#         geom_point()+
#         geom_line(aes(group=method))+
#         geom_hline(aes(yintercept=.95))+
#         ggtitle(paste("tau, tau=",aTau,sep=""))
#     )
#     
#   }
#   
#   dev.off()
#   
#   #######proportional bias
#   
#   pdf(paste(outfilepath,"proportionalBiasFor",final.CIout$simulation[1],
#             "SimDataWith",max(final.CIout$Repeat),"Repeats",dets,".pdf",sep=""),width = 7,height = 5.5)
#   
#   for (aTau in the.tau.vals){
#     ######### For tau
#     propBias=final.CIout %>% group_by(method,parameter,tau,N) %>%
#       filter(parameter=="tau")%>%
#       mutate(proportionalBias=(estimate-tau)/tau*100)
#     
#     AvgPropBias=propBias %>% 
#       summarize(AvgProportionalBias=mean(proportionalBias))
#     
#     print(ggplot(filter(AvgPropBias,tau == aTau),aes(N,AvgProportionalBias,col=method))+
#             geom_point()+
#             geom_line(aes(group=method))+
#             ggtitle(paste("tau, ","true tau=",aTau,sep=""))+
#             geom_hline(yintercept = 0)
#             # scale_x_continuous(trans='log')
#     )
#     
#     propBias$N = as.factor(propBias$N)
#     
#     print(ggplot(filter(propBias,tau == aTau),aes(N,proportionalBias))+
#       geom_violin()+
#       stat_summary(fun=mean, geom="point",color="cyan")+
#       facet_wrap(~method,ncol=4)+
#       ggtitle(paste("tau, ","true tau=",aTau,sep=""))+
#       geom_hline(yintercept = 0)
#     )
#     
#     ######### For mu
#     propBiasMu=final.CIout %>% group_by(method,parameter,tau,N) %>%
#       filter(parameter=="mu")%>%
#       mutate(proportionalBias=(estimate-mu)/mu*100)
#     
#     AvgPropBiasMu=propBiasMu %>% 
#       summarize(AvgProportionalBias=mean(proportionalBias))
#     
#     print(ggplot(filter(AvgPropBiasMu,tau == aTau),aes(N,AvgProportionalBias,col=method))+
#             geom_point()+
#             geom_line(aes(group=method))+
#             ggtitle(paste("mu, ","true tau=",aTau,sep=""))+
#             geom_hline(yintercept = 0)
#     )
# 
#     propBiasMu$N = as.factor(propBiasMu$N)
#     
#     print(ggplot(filter(propBiasMu,tau == aTau),aes(N,proportionalBias))+
#       geom_violin()+
#       stat_summary(fun=mean, geom="point",color="cyan")+
#       facet_wrap(~method,ncol=4)+
#         ggtitle(paste("mu, ","true tau=",aTau,sep=""))+
#         geom_hline(yintercept = 0)
#     )
#     
#   }
#   
#   dev.off()
#   
#   return(final.CIout)
# }
# 
# 
# propBiasViolinPlotsForSelectN=function(CIout,Nvals,Nlabs,tauVals,tauLabs,dets,plotht,plotwd){
#   # CIout=filter(CIout,method != "LaplaceLaplace")
#                
#   CIout=CIout %>% group_by(method,parameter,tau,N) %>% 
#     mutate(Repeat=1:n())
#   
#   CIoutRepVal=max(CIout$Repeat)
#   
#   bad.CIout=filter(CIout,min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0)
#   # View(bad.CIout)
#   # test=filter(CIout,min.bfmi >=.3 & min.neff < 300 & bad.rhat == 0 & num.Divergent == 0 & tree.Depth.Ex == 0)
#   # View(test)
#   # good.CIout=filter(CIout,min.bfmi >=.3 & min.neff >= 300 & bad.rhat == 0 & num.Divergent == 0 & tree.Depth.Ex == 0)
#   final.CIout=filter(CIout,is.na(min.bfmi) | !(min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0))
#   
#   keptPercent=final.CIout %>% group_by(method,parameter,tau,N) %>% 
#     summarise(kept=n()/CIoutRepVal)
#   
#   final.CIout$method = factor(final.CIout$method, 
#                               levels=c("DL","SJ","REML","PM",
#                                        "LaplaceLaplace","GaussGauss","SkewtGaussian","LaplaceGaussian"))
#   outfilepath=myfilepath
#   
#   the.tau.vals=sort(unique(final.CIout$tau))
#   
#   #######focus on small N
#   
#   pdf(paste(outfilepath,"propBiasViolinPlotsForNis",paste(Nvals,collapse = "_"),"and",final.CIout$simulation[1],
#             "SimDataWith",max(final.CIout$Repeat),"Repeats",dets,".pdf",sep=""),width = plotwd,height = plotht)
#   
#   # for (aTau in the.tau.vals){
#     # propBias=final.CIout %>% group_by(method,parameter,tau,N) %>%
#     #   mutate(proportionalBiasMu=(estimate-mu)/mu*100, proportionalBiasTau=(estimate-tau)/tau*100)
#     # propBias$N = as.factor(propBias$N)
#     # 
#     # print(ggplot(filter(propBias,N==Nval),aes(method,proportionalBiasMu))+
#     #         geom_violin()+
#     #         stat_summary(fun=mean, geom="point",color="cyan")+
#     #         facet_wrap(~method+tau,ncol=4)+
#     #         ggtitle(paste("mu, ","true tau=",aTau,sep=""))+
#     #         geom_hline(yintercept = 0)
#     # )
#     # 
#     # 
#     
#     ######### For mu
#     propBiasMu=final.CIout %>% group_by(method,parameter,tau,N) %>%
#       filter(parameter=="mu")%>%
#       mutate(proportionalBias=(estimate-mu)/mu*100)
#     propBiasMu$N = as.factor(propBiasMu$N)
#     
#     propBiasMu$tau <- factor(propBiasMu$tau, levels = tauVals,
#                         ordered = TRUE, labels=tauLabs)
#     propBiasMu$Nlab <- factor(propBiasMu$N, levels = Nvals,
#                          ordered = TRUE, labels=Nlabs)
#     
#     
#     # # print(
#     #   plot1=ggplot(filter(propBiasMu,tau == aTau & N==Nval),aes(method,proportionalBias))+
#     #         geom_violin()+
#     #         stat_summary(fun=mean, geom="point",color="cyan")+
#     #         ggtitle(paste("mu, ","true tau=",aTau,sep=""))+
#     #         geom_hline(yintercept = 0)+
#     #     theme(axis.text.x=element_text(angle = 45,hjust=1))
#     # # )
# 
#     plot1=ggplot(filter(propBiasMu, N %in% Nvals),aes(method,proportionalBias))+
#         geom_violin()+
#         stat_summary(fun=mean, geom="point",color="cyan")+
#         ggtitle(expression(mu))+
#         geom_hline(yintercept = 0)+
#         facet_wrap(~Nlab+tau,labeller = labeller(.cols = label_parsed, .multi_line = FALSE),ncol=2)+
#         theme(axis.text.x=element_text(angle = 45,hjust=1))
#     ######### For tau
#     propBias=final.CIout %>% group_by(method,parameter,tau,N) %>%
#       filter(parameter=="tau")%>%
#       mutate(proportionalBias=(estimate-tau)/tau*100)
#     propBias$N = as.factor(propBias$N)
#     
#     propBias$tau <- factor(propBias$tau, levels = tauVals,
#                              ordered = TRUE, labels=tauLabs)
#     propBias$Nlab <- factor(propBias$N, levels = Nvals,
#                               ordered = TRUE, labels=Nlabs)
#     
# 
#     # # print(
#     #   plot2=ggplot(filter(propBias,tau == aTau & N==Nval),aes(method,proportionalBias))+
#     #         geom_violin()+
#     #         stat_summary(fun=mean, geom="point",color="cyan")+
#     #         ggtitle(paste("tau, ","true tau=",aTau,sep=""))+
#     #         geom_hline(yintercept = 0)+
#     #     theme(axis.text.x=element_text(angle = 45,hjust=1))
#     # # )
#     plot2=ggplot(filter(propBias, N %in% Nvals),aes(method,proportionalBias))+
#       geom_violin()+
#       stat_summary(fun=mean, geom="point",color="cyan")+
#       ggtitle(expression(tau))+
#       geom_hline(yintercept = 0)+
#       facet_wrap(~Nlab+tau,labeller = labeller(.cols = label_parsed, .multi_line = FALSE),ncol=2)+
#       theme(axis.text.x=element_text(angle = 45,hjust=1))
#     print(grid.arrange(plot1, plot2, ncol=2))
#     
#   # }
#   
#   dev.off()
# 
# }
# 
# GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
#                            draw_group = function(self, data, ..., draw_quantiles = NULL) {
#                              data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
#                              grp <- data[1, "group"]
#                              newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
#                              newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
#                              newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
#                              
#                              if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
#                                stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
#                                                                          1))
#                                quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
#                                aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
#                                aesthetics$alpha <- rep(1, nrow(quantiles))
#                                both <- cbind(quantiles, aesthetics)
#                                quantile_grob <- GeomPath$draw_panel(both, ...)
#                                ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
#                              }
#                              else {
#                                ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
#                              }
#                            })
# 
# geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
#                               draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
#                               show.legend = NA, inherit.aes = TRUE) {
#   layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
#         position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
#         params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
# }
# propBiasViolinPlotsForSelectN_twoCIout=function(CIout,CIout2,Nvals,Nlabs,tauVals,tauLabs,dets,plotht,plotwd){
#   CIout=filter(CIout,method != "LaplaceLaplace")
#   CIout2=filter(CIout2,method != "LaplaceLaplace")
#   
#   CIout=CIout %>% group_by(method,parameter,tau,N) %>% 
#     mutate(Repeat=1:n())
#   
#   CIoutRepVal=max(CIout$Repeat)
#   
#   bad.CIout=filter(CIout,min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0)
#   final.CIout=filter(CIout,is.na(min.bfmi) | !(min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0))
#   
#   final.CIout$method = factor(final.CIout$method, 
#                               levels=c("DL","SJ","REML","PM",
#                                        "LaplaceLaplace","GaussGauss","SkewtGaussian","LaplaceGaussian"),
#                               labels=c("DL","SJ","REML","MP",
#                                        "Laplace+Laplace","Gauss+Gauss","Skew Student+Gauss","Laplace+Gauss"))
#   ##############second CIout
#   CIout2=CIout2 %>% group_by(method,parameter,tau,N) %>% 
#     mutate(Repeat=1:n())
#   
#   CIout2RepVal=max(CIout2$Repeat)
#   
#   bad.CIout2=filter(CIout2,min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0)
#   final.CIout2=filter(CIout2,is.na(min.bfmi) | !(min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0))
#   
#   final.CIout2$method = factor(final.CIout2$method, 
#                                levels=c("DL","SJ","REML","PM",
#                                         "LaplaceLaplace","GaussGauss","SkewtGaussian","LaplaceGaussian"),
#                                labels=c("DL","SJ","REML","MP",
#                                         "Laplace+Laplace","Gauss+Gauss","Skew Student+Gauss","Laplace+Gauss"))
#   #############################
#   
#   outfilepath=myfilepath
#   
#   pdf(paste(outfilepath,"TWOpropBiasViolinPlotsForNis",paste(Nvals,collapse = "_"),"and",final.CIout$simulation[1],
#             "SimDataWith",max(final.CIout$Repeat),"Repeats",dets,".pdf",sep=""),width = plotwd,height = plotht)
#   
#   ######### For mu
#   propBiasMu1=final.CIout %>% group_by(method,parameter,tau,N,simulation) %>%
#     filter(parameter=="mu")%>%
#     mutate(proportionalBias=(estimate-mu)/mu*100)
#   propBiasMu1$N = as.factor(propBiasMu1$N)
#   
#   propBiasMu2=final.CIout2 %>% group_by(method,parameter,tau,N,simulation) %>%
#     filter(parameter=="mu")%>%
#     mutate(proportionalBias=(estimate-mu)/mu*100)
#   propBiasMu2$N = as.factor(propBiasMu2$N)
#   
#   propBiasMu=rbind(propBiasMu1,propBiasMu2)
#   
#   propBiasMu$tau <- factor(propBiasMu$tau, levels = tauVals,
#                            ordered = TRUE, labels=tauLabs)
#   propBiasMu$Nlab <- factor(propBiasMu$N, levels = Nvals,
#                             ordered = TRUE, labels=Nlabs)
#   
#   # dodge <- position_dodge(width = .8)
#   plot1=ggplot(filter(propBiasMu, N %in% Nvals),aes(method,proportionalBias,col=simulation,fill=simulation))+
#     # geom_violin(position = dodge)+
#     geom_split_violin()+
#     # stat_summary(fun=mean, geom="point",aes(color=simulation),position = dodge)+
#     ggtitle(expression(mu))+
#     geom_hline(yintercept = 0)+
#     facet_wrap(~Nlab+tau,labeller = labeller(.cols = label_parsed, .multi_line = FALSE),ncol=2)+
#     theme(axis.text.x=element_text(angle = 45,hjust=1),legend.position="bottom")+
#     labs(x = "", y = "Proportional Bias", fill = "Random Effect Simulated From: ",color = "Random Effect Simulated From: ")+
#     # scale_fill_manual(labels = c("normal", "skewt"), values = c("blue", "red"))+
#     # scale_color_manual(labels = c("normal", "skewt"), values = c("blue", "red")) 
#     scale_color_hue(labels = c("Normal", "Skew Student"))+
#     scale_fill_hue(labels = c("Normal", "Skew Student"))
#   ######### For tau
#   propBias1=final.CIout %>% group_by(method,parameter,tau,N,simulation) %>%
#     filter(parameter=="tau")%>%
#     mutate(proportionalBias=(estimate-tau)/tau*100)
#   propBias1$N = as.factor(propBias1$N)
#   
#   propBias2=final.CIout2 %>% group_by(method,parameter,tau,N,simulation) %>%
#     filter(parameter=="tau")%>%
#     mutate(proportionalBias=(estimate-tau)/tau*100)
#   propBias2$N = as.factor(propBias2$N)
#   
#   propBias=rbind(propBias1,propBias2)
#   
#   propBias$tau <- factor(propBias$tau, levels = tauVals,
#                          ordered = TRUE, labels=tauLabs)
#   propBias$Nlab <- factor(propBias$N, levels = Nvals,
#                           ordered = TRUE, labels=Nlabs)
#   
#   
#   # # print(
#   #   plot2=ggplot(filter(propBias,tau == aTau & N==Nval),aes(method,proportionalBias))+
#   #         geom_violin()+
#   #         stat_summary(fun=mean, geom="point",color="cyan")+
#   #         ggtitle(paste("tau, ","true tau=",aTau,sep=""))+
#   #         geom_hline(yintercept = 0)+
#   #     theme(axis.text.x=element_text(angle = 45,hjust=1))
#   # # )
# 
#   my_labeller <- label_bquote(
#     rows = .(Nvals) ,
#     cols = .(tau) 
#   )
#   
#   plot2=ggplot(filter(propBias, N %in% Nvals),aes(method,proportionalBias,fill=simulation,col=simulation))+
#     geom_split_violin()+
#     # stat_summary(fun=mean, geom="point",color="cyan")+
#     ggtitle(expression(tau))+
#     geom_hline(yintercept = 0)+
#     labs(x = "", y = " ")+
#     facet_wrap(~Nlab+tau,labeller = labeller(.cols = label_parsed, .multi_line = FALSE),ncol=2)+#label_wrap_gen(multi_line = F,width = 50))+
#     theme(axis.text.x=element_text(angle = 45,hjust=1))
#   
#   
#   #extract legend
#   #https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
#   g_legend<-function(a.gplot){
#     tmp <- ggplot_gtable(ggplot_build(a.gplot))
#     leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#     legend <- tmp$grobs[[leg]]
#     return(legend)}
#   
#   mylegend<-g_legend(plot1)
#   
#   # print(grid.arrange(arrangeGrob(plot1+ theme(legend.position="none"), 
#   #                    plot2+ theme(legend.position="none"),
#   #                    mylegend, ncol=3,widths = c(7,7,1.5))))
#   grid.arrange(arrangeGrob(plot1 + theme(legend.position="none"),
#                                  plot2 + theme(legend.position="none"),
#                                  nrow=1, bottom = textGrob("Method",vjust = -1)),
#                      mylegend, nrow=2,heights=c(9, 1))
#   # }
#   
#   dev.off()
#   
# }
# 
# 
# 
# CPplotForSmallN=function(CIout,CIout2,Nvals,Nlabs,tauVals,tauLabs,dets,plotht,plotwd){
#   
#   CIout=filter(CIout,method != "LaplaceLaplace")
#   CIout2=filter(CIout2,method != "LaplaceLaplace")
#   
#   CIout=CIout %>% group_by(method,parameter,tau,N) %>% 
#     mutate(Repeat=1:n())
#   
#   CIoutRepVal=max(CIout$Repeat)
#   
#   bad.CIout=filter(CIout,min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0)
#   final.CIout=filter(CIout,is.na(min.bfmi) | !(min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0))
#   
#   final.CIout$method = factor(final.CIout$method, 
#                               levels=c("DL","SJ","REML","PM",
#                                        "LaplaceLaplace","GaussGauss","SkewtGaussian","LaplaceGaussian"),
#                               labels=c("DL","SJ","REML","MP",
#                                         "Laplace+Laplace","Gauss+Gauss","Skew Student+Gauss","Laplace+Gauss"))
# 
#   ##############second CIout
#   CIout2=CIout2 %>% group_by(method,parameter,tau,N) %>% 
#     mutate(Repeat=1:n())
#   
#   CIout2RepVal=max(CIout2$Repeat)
#   
#   bad.CIout2=filter(CIout2,min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0)
#   final.CIout2=filter(CIout2,is.na(min.bfmi) | !(min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0))
#   
#   final.CIout2$method = factor(final.CIout2$method, 
#                               levels=c("DL","SJ","REML","PM",
#                                        "LaplaceLaplace","GaussGauss","SkewtGaussian","LaplaceGaussian"),
#                               labels=c("DL","SJ","REML","MP",
#                                        "Laplace+Laplace","Gauss+Gauss","Skew Student+Gauss","Laplace+Gauss"))
#   
#   
#   #############################
#   outfilepath=myfilepath
#   
#   ######## coverage probs
#   
#   CPsum1=final.CIout %>% group_by(method,parameter,tau,N,simulation) %>% 
#     summarise(numberOfSims=n(),Coverage=sum(inInt)/numberOfSims)
# 
#   CPsum2=final.CIout2 %>% group_by(method,parameter,tau,N,simulation) %>% 
#     summarise(numberOfSims=n(),Coverage=sum(inInt)/numberOfSims)
#   
#   CPsum=rbind(CPsum1,CPsum2)
#   
#   CPsum$tau <- factor(CPsum$tau, levels = tauVals,
#                       ordered = TRUE, labels=tauLabs)
#   CPsum$Nlab <- factor(CPsum$N, levels = Nvals,
#                        ordered = TRUE, labels=Nlabs)
#   # CPsum$parameter <- factor(CPsum$parameter, levels = c("mu","tau"),
#   #                     ordered = TRUE, labels=c(expression(mu), expression(tau)))
#   # 
#   pdf(paste(outfilepath,"coverageProbsForBoth",
#             "SimDataWithNvals",paste(Nvals,collapse = "_"),"and",max(final.CIout$Repeat),"Repeats",dets,".pdf",sep=""),width = plotwd,height = plotht)
#   # 
#   # print(
#   #   ggplot(filter(CPsum,N==Nval),aes(parameter,Coverage,col=method))+
#   #     geom_point()+
#   #     geom_hline(aes(yintercept=.95))+
#   #     facet_wrap(~tau,labeller = label_parsed)+
#   #     scale_x_discrete(labels = c('mu' = expression(mu),
#   #                                 'tau'   = expression(tau)))+
#   #     xlab("Parameter")
#   #     # ggtitle(paste("tau=",aTau,sep=""))+
#   #     # facet_wrap(~tau,labeller = labeller(tau = c("1"="expression(tau=1)","5"="ugh2")))
#   # )
#   # 
#   
#   print(
#     # ggplot(filter(CPsum,N %in% Nvals),aes(method,Coverage,col=parameter))+
#     ggplot(filter(CPsum,!is.na(method) & N %in% Nvals),aes(method,Coverage,shape=parameter,col=simulation,group=simulation))+
#       geom_point(position = position_dodge(0.5))+
#       geom_hline(aes(yintercept=.95))+
#       facet_wrap(~Nlab+tau,labeller = labeller(.cols = label_parsed, .multi_line = FALSE),ncol=2)+
#       theme(axis.text.x=element_text(angle = 45,hjust=1))+
#       xlab("Method")+
#       labs(shape = "Parameter",color = "Random Effect")+
#       # scale_fill_manual(labels = c("normal", "skewt"), values = c("blue", "red"))+
#       # scale_color_manual(labels = c("normal", "skewt"), values = c("blue", "red")) 
#       scale_color_hue(labels = c("Normal", "Skew Student"))+
#       scale_shape(labels = c(expression(mu), expression(tau)))
#     
#     
#     # ggtitle(paste("tau=",aTau,sep=""))+
#     # facet_wrap(~tau,labeller = labeller(tau = c("1"="expression(tau=1)","5"="ugh2")))
#   )
#   
#   dev.off()   
# }
# 
# 
# shortProcessResults=function(CIout){
#   
#   checkBadRuns(CIout)
#   
#   CIout=CIout %>% group_by(method,parameter,tau,N) %>% 
#     mutate(Repeat=1:n(),totalIters=n(),
#            good=is.na(min.bfmi) | !(min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0))
#   View(CIout)
# 
#   keptPercent=CIout %>% group_by(method,parameter,tau,N) %>% 
#     summarise(total=n(),goodTot=sum(good),kept=sum(good)/total)
#   View(keptPercent)
#   
# 
#   bad.CIout=filter(CIout,min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0)
#   final.CIout=filter(CIout,is.na(min.bfmi) | !(min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0))
#   
#   the.N.vals=sort(unique(final.CIout$N))
#   the.tau.vals=sort(unique(final.CIout$tau))
# 
#   CPsum=final.CIout %>% group_by(method,parameter,tau,N) %>% 
#     summarise(numberOfSims=n(),coverage=sum(inInt)/numberOfSims)
#   View(CPsum)
#   
# }
# 
# ################081821 runs, changed the prior for tau to df = 4, sd =2 for N=5 and tau=1, reran all
# CIout=collectRuns(filepath = myfilepath,runDate = "251208",howDatSim = "normal") 
# 
# shortProcessResults(CIout)
# 
# 
# normalRes=processResults(CIout = CIout,"test")
# 
# 
# 
# tauValz=c(1)
# tauLabz=c(expression(paste(tau,"=1")))
# 
# Nvalz=c(5,10)
# Nlabz=c(expression(paste("N=5")), expression(paste("N=10")))
# CPplotForSmallN(CIout,CIout,Nvals = Nvalz,Nlabs = Nlabz,
#                 tauVals = tauValz,tauLabs = tauLabz,"210818NEWPRIOR",plotht=7,plotwd=6)
# 
# propBiasViolinPlotsForSelectN(CIout,Nvals = Nvalz,
#                               Nlabs = Nlabz,tauVals = tauValz,tauLabs = tauLabz,"210818NEWPRIOR",
#                               plotht=5,plotwd=8)
# 
