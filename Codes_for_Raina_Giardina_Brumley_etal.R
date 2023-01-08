library(plyr)
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(rstatix)

set.seed(123)

######################################################
# Analysing co-growth data with repeated measure ANOVA
######################################################

Cogrowth <- read.csv("cogrowth.csv")

# Reorder the dataset
Cogrowth$time <- ordered(Cogrowth$time,levels = c("T0", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12"))

#==================
# Check assumptions
#==================

# Outliers
Cogrowth %>%
  group_by(treatment, time) %>%
  identify_outliers(Syne_ml)

# Normality
Cogrowth %>%
  group_by(treatment, time) %>%
  shapiro_test(Syne_ml)

ggqqplot(Cogrowth, "Syne_ml", ggtheme = theme_bw()) +
  facet_grid(time ~ treatment, labeller = "label_both")

#======================
# Compute the main test
#======================

res.aov <- anova_test(
  data = Cogrowth, dv = Syne_ml, wid = id, within = c(treatment, time))

# The following function produces an ANOVA table that is automatically corrected for eventual deviation from the sphericity assumption. It applies the Greenhouse-Geisser sphericity correction by default but this can be changed if needed.
get_anova_table(res.aov)
write.csv(res.aov, "ANOVA.csv")

#============================
# (1) Simple Main effect test
#============================

# This will test the effect of the treatment at each time point
one.way <- Cogrowth %>%
  group_by(time) %>%
  anova_test(dv = Syne_ml, wid = id, within = treatment) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way

write.csv(one.way, "SMET.csv")

#==================================================
# (2) Pairwise comparisons between treatment groups
#==================================================
pwc <- Cogrowth %>%
  group_by(time) %>%
  pairwise_t_test(
    Syne_ml ~ treatment, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

write.csv(pwc, "pwc.csv")

############################################
# Graph growth data
############################################

#Co-growth
Cogrowth <- read.csv("cogrowth.csv")

Cogrowth$time <- ordered(Cogrowth$time,levels = c("T0", "T1", "T2", "T3", "T4", "T5","T6", "T7", "T8", "T9", "T10", "T11", "T12"))
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

dsum <- data_summary(Cogrowth, varname="Syne_log", groupnames=c("treatment", "time"))

# Convert time to a factor variable
dsum$time=as.factor(dsum$time)

# Graph using mean + SD
ggplot(dsum, aes(x=time, y=Syne_log, group=treatment, color=treatment)) + 
  geom_line()+
  geom_pointrange(aes(ymin=Syne_log-sd, ymax=Syne_log+sd))

#############################################################################

#Growth curves
Growth_curves <- read.csv("curves.csv")

Growth_curves$time <- ordered(Growth_curves$time,levels = c("T0", "T2", "T4", "T6", "T8", "T10","T12", "T14", "T16", "T18", "T20", "T22", "T24"))
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

dsum <- data_summary(Growth_curves, varname="Cells_ml", 
                     groupnames=c("treatment", "time"))

# Convert time to a factor variable
dsum$time=as.factor(dsum$time)

# Quick graph using mean + SD
ggplot(dsum, aes(x=time, y=Cells_ml, group=treatment, color=treatment)) + 
  geom_line()+
  geom_pointrange(aes(ymin=Cells_ml-sd, ymax=Cells_ml+sd))

###################################################################
# Analysing NanoSIMS data with non parametric tests
###################################################################

my_data <- read.csv("Fxnet_N.csv")

####################################################################
# Kruskal-Wallis N
kruskal.test(nitrogen ~ treatment, data = my_data)

pairwise.wilcox.test(my_data$nitrogen, my_data$treatment,
                     p.adjust.method = "BH")


# Kruskal-Wallis C
kruskal.test(carbon ~ treatment, data = my_data)

pairwise.wilcox.test(my_data$carbon, my_data$treatment,
                     p.adjust.method = "BH")

#############################################################################

#Nitrogen
N <- ggbarplot(my_data, x = "concentration", y = "nitrogen", 
          add = c("mean_se"),
          color = "treatment", palette = "jco",
          position = position_dodge(0.8))
#Set the y axis to range
N + scale_y_continuous(expand = c(0, 0))
N + coord_cartesian(ylim=c(0.37,1.2))

#########################################

my_data <- read.csv("Fxnet_C.csv")

#Carbon
C <- ggbarplot(my_data, x = "concentration", y = "carbon", 
               add = c("mean_se"),
               color = "treatment", palette = "jco",
               position = position_dodge(0.8))
#Set the y axis to range 
C + scale_y_continuous(expand = c(0, 0))

########################################

#Violin plots
dN <- ggplot(my_data, aes(x=concentration, y=nitrogen, fill=treatment)) + 
  geom_violin()+
  labs(x="Synechococcus concentration", y = "Assimilation")
dN + theme_classic()

dC <- ggplot(my_data, aes(x=concentration, y=carbon, fill=treatment)) + 
  geom_violin()+
  labs(x="Synechococcus concentration", y = "Assimilation")
dC + theme_classic()

########################################

#Jitter plots

# Nitrogen
my_data <- read.csv("Fxnet_N.csv")
# 0.2 : degree of jitter in x direction
n<-ggplot(my_data, aes(x=concentration, y=nitrogen, color=treatment)) + 
  geom_jitter(position=position_jitter(0.4))

n + stat_summary(fun.y=mean, geom="point", shape=20,
                 size=3, color="black")

# Carbon
my_data <- read.csv("Fxnet_C_figS10.csv")
# 0.2 : degree of jitter in x direction
c<-ggplot(my_data, aes(x=concentration, y=carbon, color=treatment)) + 
  geom_jitter(position=position_jitter(0.4))

c + stat_summary(fun.y=mean, geom="point", shape=20,
                 size=3, color="black")

########################################
# Graph Chemotaxis data
my_data <- read.csv("Chemotaxis_2022.csv")

#Chemotaxis
I <- ggbarplot(my_data, x = "treatment", y = "index", 
               add = c("mean_se", "jitter"),
               position = position_dodge(0.8))
#Set the y axis to range 
I + scale_y_continuous(expand = c(0, 0),limits = c(0, 5.2))




