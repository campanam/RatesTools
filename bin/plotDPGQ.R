#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------------
# plotDPGQ version 0.2
# Ellie E. Armstrong & Michael G. Campana, 2021
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
#library(ggpubr)

stem = commandArgs(trailingOnly=TRUE) # Get filestem for output

#get list of dpgq txt files in directory
dpqg_files <- list.files(pattern = "*.txt")

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
df_clean$depth_log <- log(as.numeric(as.character(df_clean$depth)))
df_clean$qual_log <- log(as.numeric(as.character(df_clean$qual)))
df_clean$depth <- as.numeric(as.character(df_clean$depth))
df_clean$qual <- as.numeric(as.character(df_clean$qual))


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
jpeg(file=paste(stem,"log_depth_qual.jpg",sep = ""), width = 550, height = 700)
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


write.csv(depth_table, file=paste(stem,"depth_ratestools.csv", sep = ""), row.names=FALSE)
write.csv(qual_table, file=paste(stem,"qual_ratestools.csv", sep = ""), row.names=FALSE)
