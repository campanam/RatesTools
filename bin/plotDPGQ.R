#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------------
# plotDPGQ version 0.3.0
# Ellie E. Armstrong & Michael G. Campana, 2022
# Stanford University and Smithsonian Institution

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to RatesTools;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.
 
# We politely request that this work be cited as:
# Campana, M.G. & E.E. Armstrong. 2020. RatesTools: Pipeline to calculate de novo
# mutation rates from parent-offspring trios. Smithsonian Institution and Stanford
# University. <https://github.com/campanam/RatesTools>.
#----------------------------------------------------------------------------------------

#load libraries
library(data.table)
library(tidyverse) # Also loads dplyr library

stem = commandArgs(trailingOnly=TRUE) # Get filestem for output

#get list of dpgq txt files in directory
dpgq_files <- list.files(pattern = "*.txt")

df <- dpgq_files %>%
  set_names(.) %>%
  map_df(read.delim, header=FALSE, sep = " ", .id = "FileName")

#clean up empty columns and column 1 for plotting
df_clean <- subset(df, select = -c(V3,V5))
df_clean <- df_clean %>% 
  rename(chr=V1, position=V2, depth=V4, qual=V6, individual=FileName)
#Remove extraneous bits of file name
df_clean$individual <- sub(".*_offspring","", df_clean$individual, perl=TRUE)
df_clean$individual <- sub(".variants.txt","",df_clean$individual)
df_clean$depth_log <- log(as.numeric(as.character(df_clean$depth)))
df_clean$qual_log <- log(as.numeric(as.character(df_clean$qual)))
df_clean$depth <- as.numeric(as.character(df_clean$depth))
df_clean$qual <- as.numeric(as.character(df_clean$qual))


#generate plots
png(file=paste(stem,"_log_depth.png",sep = ""), width = 550, height = 700)
ggplot(df_clean, aes(x=depth_log, color=individual)) + 
  geom_histogram(fill="white", alpha=0.5, position="identity") +
  xlab("") + ylab("log depth") +
  theme(legend.position="top", legend.title = element_blank())
dev.off()
png(file=paste(stem,"_depth.png",sep = ""), width = 550, height = 700)
ggplot(df_clean, aes(x=depth, color=individual)) + 
  geom_histogram(fill="white", alpha=0.5, position="identity") +
  xlab("") + ylab("depth") +
  theme(legend.position="top", legend.title = element_blank())
dev.off()
png(file=paste(stem,"_log_qual.png",sep = ""), width = 550, height = 700)
ggplot(df_clean, aes(x=qual_log, color=individual)) + 
  geom_histogram(fill="white", alpha=0.5, position="identity") +
  xlab("") + ylab("log GQ") +
  theme(legend.position="top", legend.title = element_blank())
dev.off()
png(file=paste(stem,"_qual.png",sep = ""), width = 550, height = 700)
ggplot(df_clean, aes(x=qual, color=individual)) + 
  geom_histogram(fill="white", alpha=0.5, position="identity") +
  xlab("") + ylab("log GQ") +
  theme(legend.position="top", legend.title = element_blank())
dev.off()

samples <- c(unique(df_clean$individual))
depth_table <- data.frame(matrix(ncol = 7, nrow = 0))
qual_table <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(depth_table) = c("Sample","0%","25%","50%","75%","100%","mean")
colnames(qual_table) = c("Sample","0%","25%","50%","75%","100%","mean")

for (ind in samples) {
#generate table of quantiles per individual
depth_quant <- quantile(df_clean[df_clean$individual == ind,]$depth, na.rm=TRUE)
depth_quant$mean <- mean(df_clean[df_clean$individual == ind,]$depth, na.rm=TRUE)
tmp <- c(ind,unlist(depth_quant))
depth_table[nrow(depth_table)+1,] <- tmp

qual_quant <- quantile(df_clean[df_clean$individual == ind,]$qual, na.rm=TRUE)
qual_quant$mean <- mean(df_clean[df_clean$individual == ind,]$qual, na.rm=TRUE)
tmp <-c(ind,unlist(qual_quant))
qual_table[nrow(qual_table)+1,] <- tmp
}
#Calculate values for all individuals together
depth_quant <- quantile(df_clean$depth, na.rm=TRUE)
depth_quant$mean <- mean(df_clean$depth, na.rm=TRUE)
tmp <- c("AllCombined",unlist(depth_quant))
depth_table[nrow(depth_table)+1,] <- tmp

qual_quant <- quantile(df_clean$qual, na.rm=TRUE)
qual_quant$mean <- mean(df_clean$qual, na.rm=TRUE)
tmp <- c("AllCombined",unlist(qual_quant))
qual_table[nrow(qual_table)+1,] <- tmp

write.csv(depth_table, file=paste(stem,"_depth_ratestools.csv", sep = ""), row.names=FALSE)
write.csv(qual_table, file=paste(stem,"_qual_ratestools.csv", sep = ""), row.names=FALSE)
