seed <- set.seed(42)
options(warn = -1)
Sys.setenv(JAVA_HOME='')
# dependencies
# BiocManager::install("clusterProfiler", force=TRUE)
# BiocManager::install("org.Hs.eg.db", force = TRUE)
# install.packages("DT")
library(kableExtra)
library(ggfortify)
library(dplyr)
library(rmcfs)
library(rJava)
library(R.ROSETTA)
library(devtools)
library(kableExtra)
library(devtools)
library('VisuNet')
library('clusterProfiler')
library('org.Hs.eg.db')
library('wordcloud') 
library('RColorBrewer')
library('tm')
library('DT')

df <- read.csv("Project5.csv", sep = "\t")
df.head <- head(df, 10)
kbl(df.head) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center")
dim(df)

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


df$P772  # Some values in here, but sparsity is extremely high. Should consider just removing these columns?

result <- mcfs(Host ~ ., df, projections = 1500, projectionSize = 0.1, splits = 5, splitSetSize = 500, cutoffPermutations = 6, threadsNumber = 8, seed = seed)
plot(result, type="distances", legend=TRUE)
result2 <- result$RI[1:result$cutoff_value,]
result2

df$Host
df[1, c(1: 25, dim(df)[2])]
df$P48
# 
# rownames(df) <- df[1]
# df <- df[2:dim(df)[2]]
# df.head <- head(df, 10)
kbl(df.head) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center")

table(df$Host)

dim(result$data)

# # prepare id graph plot
# gid <- build.idgraph(result, size = 20)
# plot.idgraph(gid, label_dist = 0.5)
# 
# # 
rosAv <- rosetta(result$data, discrete=T, roc = TRUE, clroc = "Avian")
rosAv

rules <- rosAv$main
rlsAv <- viewRules(rules[rules$decision == "Avian", ])
rlsAv
rosHum <- rosetta(result$data, discrete=T, roc = TRUE, clroc = "Avian")
rosHum
rules <- rosHum$main
rlsHum <- viewRules(rules[rules$decision == "Human", ])
rlsHum


kbl(result$data) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center")

visAv <- VisuNet::visunet(rules)
visAv


