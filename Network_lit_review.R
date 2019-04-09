### clean the memory to avoid unnecessary errors:
rm(list = ls())

### set directory to the folder of analytic data

# Get the directory of the current R script
curWD <- dirname(rstudioapi::getSourceEditorContext()$path)

# Set the working directory to the directory where this script is 
setwd(curWD)

# load the packages needed, if not exist, download from cran
if (!require(tidyverse)) {install.packages("tidyverse",repos = "http://cran.us.r-project.org"); require(tidyverse)}
if (!require(dplyr)) {install.packages("dplyr",repos = "http://cran.us.r-project.org"); require(dplyr)}
if (!require(psych)) {install.packages("psych",repos = "http://cran.us.r-project.org"); require(psych)}
library("psych")
library("dplyr")
library("tidyverse")

#load data
d <- read.csv("/Users/apple/Desktop/Dataset4.csv", header = TRUE)
#drop NAs
d <- drop_na(d)
#summary data
dim(d)
summary(d)
#reverse score
data_reverse <- d %>% 
  dplyr::mutate(A6 = recode( A6 , '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                A7 = recode( A7 , '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                A8 = recode( A8 , '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                A9 = recode( A9 , '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                C15 = recode(C15, '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                C16 = recode(C16, '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                C17 = recode(C17, '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                C18 = recode(C18, '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                N24 = recode(N24, '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                N25 = recode(N25, '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                N26 = recode(N26, '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                O35 = recode(O35, '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                O36 = recode(O36, '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                E42 = recode(E42, '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                E43 = recode(E43, '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1),
                E44 = recode(E44, '1' = 5, '2' = 4, '3' = 4, '4' = 2, '5' = 1)) %>%
  dplyr::select(A1, A2, A3, A4, A5, A6, A7, A8, A9,
                C10, C11, C12, C13, C14, C15, C16, C17, C18,
                N19, N20, N21, N22, N23, N24, N25, N26, 
                O27, O28, O29, O30, O31, O32, O33, O34, O35, O36,
                E37, E38, E39, E40, E41, E42, E43, E44)
#check reversed data
summary(data_reverse)
#install qgraph
install.packages("qgraph")
library(qgraph)

#Reminder:
#FDR: most frequently used, node placement and distance are not meaningful
#MDS: node placement and distance are meaningful--Divide into dimensions
#PCA: extract meaningful dimensions

##########################1 FDA network for 5 traits#############################
data_zerooder <- cor(data_reverse)
qgraph(data_zerooder, layout="spring",
       groups = list("Agreeableness" = 1:9, 
                     "Conscientiousness" = 10:18,
                     "Neuroticism" = 19:26,
                     "Openness" = 27:36,
                     "Extroversion" = 37:44),
       color = c("lightblue", "lightsalmon", "antiquewhite", "coral1", "darkolivegreen3"))


##########################2 networks for 2 traits #############################
#agreement and extraversion are said interacting with each other
ae <- data_reverse %>%
  dplyr::select(A1, A2, A3, A4, A5, A6, A7, A8, A9,
                E37, E38, E39, E40, E41, E42, E43, E44)
##FDA-partial
ae_zerooder <- cor(ae)
qgraph(ae_zerooder, layout="spring",
       groups = list("Agreeableness" = 1:9, 
                     "Extroversion" = 10:17),
       color = c("lightblue", "lightsalmon"), graph = "pcor")
##FDA-glasso
ae_glasso <- EBICglasso(cor(ae), n=175) 
qgraph(ae_glasso,
       groups = list("Agreeableness" = 1:9, "Extroversion" = 10:17),
       color = c("lightblue", "lightsalmon"), vsize=4)
text(-1,-1, paste("Stress=", round(ae_MDS_mspline$stress,2)))

##MDS-Partial
install.packages("smacof")
library(smacof) 
ae_dissimilarity <-sim2diss(ae_zerooder)
ae_mds <- mds(ae_dissimilarity) 
head(round(ae_mds$conf, 2))
#determine transformation type 
#ordinal
ae_MDS_ordinal <- mds(ae_dissimilarity, type="ordinal")
plot(ae_MDS_ordinal, plot.type = "Shepard", main="Ordinal")
text(1.1,0.3, paste("Stress =", round(ae_MDS_ordinal$stress,2)))
#ratio
ae_MDS_ratio <- mds(ae_dissimilarity, type="ratio")
plot(ae_MDS_ratio, plot.type = "Shepard", main="Ratio")
text(1.1,0.3, paste("Stress =", round(ae_MDS_ratio$stress,2)))
#interval
ae_MDS_interval <- mds(ae_dissimilarity, type="interval")
plot(ae_MDS_interval, plot.type = "Shepard", main="Interval")
text(1.1,0.3, paste("Stress =", round(ae_MDS_interval$stress,2)))
#mspline
ae_MDS_mspline <- mds(ae_dissimilarity, type="mspline")
plot(ae_MDS_mspline, plot.type = "Shepard", main="Spline")
text(1,0.3, paste("Stress =", round(ae_MDS_mspline$stress,2)))#choose

qgraph(ae_zerooder, layout=ae_MDS_mspline$conf,
       groups = list("Agreeableness" = 1:9, "Extroversion" = 10:17), 
       color = c("lightblue","lightsalmon"),graph = "pcor", 
       vsize=4)
text(-1,-1, paste("Stress=",round(ae_MDS_mspline$stress,2)))

#MDS-glasso
ae_glasso <- EBICglasso(cor(ae), n=175) 
qgraph(ae_glasso,layout=ae_MDS_mspline$conf,
       groups = list("Agreeableness" = 1:9, "Extroversion" = 10:17),
       color = c("lightblue", "lightsalmon"), vsize=4)
text(-1,-1, paste("Stress=", round(ae_MDS_mspline$stress,2)))

#principle components
library("psych")
ae_PCA <- principal(cor(ae), nfactors = 2) 
qgraph(ae_glasso, layout=ae_PCA$loadings,
       groups =list("Agreeableness" = 1:9, 
                    "Extroversion" = 10:17),color = c("lightblue", "lightsalmon"), title="ae", layoutOffset=c(.3,.1), vsize=1)
text(1.5,-.8, paste("% var=", round(sum(ae_PCA$values[1:2]/ length(ae_PCA$values)),2)))
title(xlab="Component 1", ylab= "Component 2")

##########################3 networks for 3 traits #############################
#agreeableness, extroversion and neuroticism (neuroticism is said to be relatively independent)
aen <- data_reverse %>%
  dplyr::select(A1, A2, A3, A4, A5, A6, A7, A8, A9,
                N19, N20, N21, N22, N23, N24, N25, N26, 
                E37, E38, E39, E40, E41, E42, E43, E44)

##FDA-partial
aen_zerooder <- cor(aen)
qgraph(aen_zerooder, layout="spring",
       groups = list("Agreeableness" = 1:9, 
                     "neurotism" = 10:17,
                     "Extroversion" = 18:25),
       color = c("lightblue", "lightsalmon","antiquewhite"), graph = "pcor")
##FDA-glasso
qgraph(aen_zerooder, layout="spring",
       groups = list("Agreeableness" = 1:9, 
                     "neurotism" = 10:17,
                     "Extroversion" = 18:25),
       color = c("lightblue", "lightsalmon","antiquewhite"), sampleSize= 175, graph = "glasso")
text(-1,1, paste("Stress=", round(ae_MDS_mspline$stress,3)))

#MDS-partial
aen_dissimilarity <-sim2diss(aen_zerooder)
aen_mds <- mds(aen_dissimilarity) 
head(round(ae_mds$conf, 3))
#determine transformation type 
#ordinal
aen_MDS_ordinal <- mds(aen_dissimilarity, type="ordinal")
plot(aen_MDS_ordinal, plot.type = "Shepard", main="Ordinal")
text(1.1,0.3, paste("Stress =", round(aen_MDS_ordinal$stress,2)))
#ratio
aen_MDS_ratio <- mds(aen_dissimilarity, type="ratio")
plot(aen_MDS_ratio, plot.type = "Shepard", main="Ratio")
text(1.1,0.3, paste("Stress =", round(aen_MDS_ratio$stress,2)))
#interval
aen_MDS_interval <- mds(aen_dissimilarity, type="interval")
plot(aen_MDS_interval, plot.type = "Shepard", main="Interval")
text(1.1,0.3, paste("Stress =", round(aen_MDS_interval$stress,2)))
#mspline
aen_MDS_mspline <- mds(aen_dissimilarity, type="mspline")
plot(aen_MDS_mspline, plot.type = "Shepard", main="Spline")
text(1,0.3, paste("Stress =", round(aen_MDS_mspline$stress,2)))#choose

qgraph(aen_zerooder, layout=aen_MDS_mspline$conf,
       groups = list("Agreeableness" = 1:9, 
                     "neurotism" = 10:17,
                     "Extroversion" = 18:25), 
       color = c("lightblue","lightsalmon", "antiquewhite"), graph = "pcor", 
       vsize=4)
text(-1,-1, paste("Stress=",round(aen_MDS_mspline$stress,3)))

#glasso
aen_glasso <- EBICglasso(cor(aen), n=175) 

qgraph(aen_glasso, layout="spring",
       groups = list("Agreeableness" = 1:9, 
                     "neurotism" = 10:17,
                     "Extroversion" = 18:25),
       color = c("lightblue", "lightsalmon","antiquewhite"))

qgraph(aen_glasso,layout=aen_MDS_mspline$conf,
       groups = list("Agreeableness" = 1:9, 
                     "neurotism" = 10:17,
                     "Extroversion" = 18:25),
       color = c("lightblue", "lightsalmon","antiquewhite"), vsize=4)
text(-1,-1, paste("Stress=", round(aen_MDS_mspline$stress,3)))
#simplify the graph to avoid overlapping
library("wordcloud")
qgraph(aen_glasso,layout=aen_MDS_mspline$conf,
       groups = list("Agreeableness" = 1:9, 
                     "neurotism" = 10:17,
                     "Extroversion" = 18:25),
       color = c("lightblue", "lightsalmon","antiquewhite"), vsize=0, rescale=FALSE, labels=FALSE)
points(aen_MDS_mspline$conf, pch=16) 
textplot(aen_MDS_mspline$conf[,1]+.03, aen_MDS_mspline$conf[,2]+.03, colnames(aen), new=F)

#principle components
library("psych")
aen_PCA <- principal(cor(aen), nfactors = 2) 
qgraph(aen_glasso, layout=aen_PCA$loadings,
groups =list("Agreeableness" = 1:9, 
             "neurotism" = 10:17,
             "Extroversion" = 18:25),color = c("lightblue", "lightsalmon","antiquewhite"), title="aen", layoutOffset=c(.3,.1), vsize=4)
text(1.5,-.8, paste("% var=", round(sum(aen_PCA$values[1:2]/ length(aen_PCA$values)),2)))
title(xlab="Component 1", ylab= "Component 2")