setwd("A:/MSc/KB_Project/Knowledge-based-Project-5")
seed <- set.seed(42)
options(warn = -1)
Sys.setenv(JAVA_HOME='')
# dependencies
# BiocManager::install("clusterProfiler", force=TRUE)
# BiocManager::install("org.Hs.eg.db", force = TRUE)
# install.packages("DT")
library('ggfortify')
library('dplyr')
library('rmcfs')
library('rJava')
library('R.ROSETTA')
library('devtools')
library('kableExtra')
library('devtools')
library('clusterProfiler')
library('org.Hs.eg.db')
library('wordcloud') 
library('RColorBrewer')
library('tm')
library('DT')
# install_github("komorowskilab/VisuNet", force=TRUE)
library(VisuNet)

# -------------------------------------
#   DATA EXPLORATION/PREPROCESSING
# -------------------------------------


## load data
df <- read.csv("Project5.csv", sep = "\t")
df.head <- head(df, 100)
kbl(df) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center")
dim(df)

tab <- as.data.frame(table(df$Host))
colnames(tab) <- c("Host", "Frequency")
kbl(tab) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center", font_size = 80)




# initial read shows that some T AAs were misread as "TRUE" and F as "FALSE". Let's fix that.
df[] <- lapply(df, function(x) ifelse(x == "TRUE", "T", x))
df[] <- lapply(df, function(x) ifelse(x == "FALSE", "F", x))
df.head <- head(df, 10)
kbl(df.head) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center")

# Remove null columns
df <- df %>%
dplyr::select(where(~ !all(. == "?")))
df.head <- head(df, 10)
kbl(df.head) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center")
dim(df)


df$P772  # Some values in here, but sparsity is extremely high. Should consider just removing these columns
df <- df[, sapply(df, function(x) mean(x == "?") < 0.8)]
df.head <- head(df, 10)
kbl(df.head) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center")
dim(df)



# -------------------------
#           MCFS
# -------------------------

result <- mcfs(Host ~ ., df, projections = 3000, projectionSize = 0.1, splits = 5, cutoffPermutations = 20, threadsNumber = 16, seed = seed)
plot(result, type="distances", legend=TRUE)
result2 <- result$RI[1:result$cutoff_value,]
result2
df$Host
df[1, c(1: 25, dim(df)[2])]
df$P48

kbl(result$RI) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center")

dim(result$data)

kbl(result$data) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center")

df <- df[2:dim(df)[2]]
df.head <- head(df, 10)
kbl(df.head) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center")

table(df$Host)

dim(result$data)

tab <- as.data.frame(table(result$data$P48))
colnames(tab) <- c("AA", "Freq")
kbl(tab) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center", font_size = 80)




# # prepare id graph plot
gid <- build.idgraph(result, orphan_nodes = T)
plot.idgraph(gid, label_dist = 1, cex=2, seed=seed)
gid <- build.idgraph(result, size=10, orphan_nodes = T)
plot.idgraph(gid, label_dist = 1, cex=2, seed=seed)

# ?build.idgraph
# ?plot.idgraph

#------------------
#   R. ROSETTA
#------------------
#Documentation
?rosetta


# ------- Run Rosetta classification
# AVIAN
rosAvGenetic <- rosetta(result$data, discrete=T, roc = TRUE, reducer="Genetic", clroc = "Avian")
rosAvJohnson <- rosetta(result$data, discrete=T, roc = TRUE, clroc = "Avian")

# HUMAN
rosHumGenetic <- rosetta(result$data, discrete=T, roc = TRUE, reducer="Genetic", clroc = "Human")
rosHumJohnson <- rosetta(result$data, discrete=T, roc = TRUE, clroc = "Human")

# roc
rosAvGenetic$ROCstats

plotRule(df, recAutconJohnson, type="heatmap", discrete=T, ind=topRuleInd)
plotMeanROC(rosAvJohnson)
plotMeanROC(rosAvGenetic)



# ------- Build Rulesets
# AVIAN
rlsAvGenetic <- rosAvGenetic$main
rlsAvJohnson <- rosAvJohnson$main
res_rlsAvGenetic <- viewRules(rlsAvGenetic[rlsAvGenetic$decision == "Avian", ])
res_rlsAvJohnson <- viewRules(rlsAvJohnson[rlsAvJohnson$decision == "Avian", ])


# HUMAN
rlsHumGenetic <- rosHumGenetic$main
rlsHumJohnson <- rosHumJohnson$main
res_rlsHumGenetic <- viewRules(rlsHumGenetic[rlsHumGenetic$decision == "Human", ])
res_rlsHumJohnson <- viewRules(rlsHumJohnson[rlsHumJohnson$decision == "Human", ])

# ANALYZE rule based model
# significant rules in the model
tabS <- table(rlsAvGenetic[rules$pValue < 0.05,]$decision)
tabS

# fraction of significant rules, in [%]
tabS[1]/length(rlsAvGenetic$decision=="Avian")*100
tabS[2]/length(rlsAvGenetic$decision=="Human")*100

featsJohnson <- getFeatures(rlsAvGenetic, filter = T, filterType = "pvalue", thr = 0.05)
featsGenetic <- getFeatures(rlsAvGenetic, filter = T, filterType = "pvalue", thr = 0.05)
featsJohnson
featsGenetic


featsJohnsonAv <- featsJohnson$features$Avian
featsJohnsonAv
featsJohnsonHum <- featsJohnson$features$Human
featsJohnsonHum

topRuleIndGA <- which(sapply(rosAvGenetic$main$features, function(x) length(strsplit(x, ',')[[1]])) == 2 & rlsAvGenetic$decision == "Avian")[1]
topRuleIndGA
topRuleIndGH <- which(sapply(rosAvGenetic$main$features, function(x) length(strsplit(x, ',')[[1]])) == 2 & rlsAvGenetic$decision == "Human")[1]
topRuleIndGH




rlsAvGenetic$pValue[topRuleIndG]
rlsAvGenetic$pValue[topRuleIndG]
rlsAvJohnson$pValue[1]
rlsAvJohnson$pValue[1]



# ------- View Rules

# JOHNSON
kbl(rlsAvJohnson) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center")
kbl(rlsHumJohnson) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center")


# GENETIC
kbl(rlsAvGenetic) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center")
kbl(rlsHumGenetic) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center")

# ------- Visualize the results in VisuNet
# JOHNSON
visAvJohnson <- VisuNet::visunet(rlsAvJohnson)  # visuNet - Avian - Johnson
visHumJohnson <- VisuNet::visunet(rlsHumJohnson)  # visunet - Human - Johnson

# GENETIC
visAvGenetic <- VisuNet::visunet(rlsAvGenetic)  # visunet - Avian - Genetic
visHumGenetic <- VisuNet::visunet(rlsHumGenetic)  # visunet - Avian - Genetic


table(unlist(strsplit(as.character(recAutconJohnson$levels), ",")))

Hums <- visAvGenetic$Human
rlsHum <- rlsHumGenetic$RO
 
Avs <- visAvGenetic$Avian


Hums


visAvGenetic$all$nodes


# View results
visAvJohnson
rlsAvJohnson
visAvGenetic
rlsAvGenetic

rlsSigGen <- rlsAvGenetic[rlsAvGenetic$pValue < 0.05 & rlsAvGenetic$decision == "Avian", ]
dim(rlsSigGen)
?visunet
?rosetta


