#load libraries
library(data.table)
library(tidyverse)
library(data.table)
library(ggpubr)
library(dplyr)
library(data.table)

#set working directory
#you need to reset this or figure out a way to pipe a variable in from realpath or something
setwd('~/Documents/Documents - Ellie’s MacBook Pro (2)/Mutation_Rate/scripts_ratestools/')

###load data
#get list of dpgq txt files in directory
dpqg_files <- list.files(pattern = "*.txt")
for(filename in dpqg_files){
  print(filename)
}
#?????????????????????????????????because this makes any sense
df <- dpqg_files %>%
  set_names(.) %>%
  map_df(read.delim, header=FALSE, sep = " ", .id = "FileName")

#clean up empty columns and column 1 for plotting
df_clean <- subset(df, select = -c(V3,V5))
df_clean <- df_clean %>% 
  rename(chr=V1, position=V2, depth=V4, qual=V6, individual=FileName)
#these filters assume that your IDs are before either a '_' or a '.'
df_clean$individual <- sub("_\\S+","", df_clean$individual, perl=TRUE)
df_clean$individual <- gsub("\\..*","",df_clean$individual)
df_clean$depth_log<- log(as.numeric(as.character(df_clean$depth)))
df_clean$qual_log<- log(as.numeric(as.character(df_clean$qual)))
df_clean$depth<- as.numeric(as.character(df_clean$depth))
df_clean$qual<- as.numeric(as.character(df_clean$qual))


#generate plots 
depth <- ggplot(df_clean, aes(x=depth, color=individual)) + 
  geom_histogram(fill="white", alpha=0.5, position="identity") +
  xlab("") + ylab("log depth") +
  theme(legend.position="top", legend.title = element_blank())
qual <- ggplot(df_clean, aes(x=qual, color=individual)) + 
  geom_histogram(fill="white", alpha=0.5, position="identity") +
  xlab("") + ylab("log GQ") +
  theme(legend.position="top", legend.title = element_blank())

figure <- ggarrange(depth, qual,
                    ncol = 1, nrow = 2)
#need to add something that produces png
jpeg(file="log_depth_qual.jpg", width = 550, height = 700)
figure
dev.off()

#generate table of quantiles per individual

depth_quant <- quantile(df_clean$depth, na.rm=TRUE)
depth_quant$mean <- mean(df_clean$depth, na.rm=TRUE)
depth_table <- data.frame(matrix(unlist(depth_quant), nrow=1, byrow=TRUE),stringsAsFactors=FALSE)
colnames(depth_table) = c("0%","25%","50%","75%","100%","mean")

qual_quant <- quantile(df_clean$qual, na.rm=TRUE)
qual_quant$mean <- mean(df_clean$qual, na.rm=TRUE)
qual_table <- data.frame(matrix(unlist(qual_quant), nrow=1, byrow=TRUE),stringsAsFactors=FALSE)
colnames(qual_table) = c("0%","25%","50%","75%","100%","mean")


write.csv(depth_table, file="depth_ratestools.csv", row.names=FALSE)
write.csv(qual_table, file="qual_ratestools.csv", row.names=FALSE)
