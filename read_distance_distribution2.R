#this script takes an annotated probe report from seqmonk, expected to contain read counts
#in some windowed interval (usually 50kb from TrAEL-seq)
#columns 13 onwards contain read counts

#it also takes a list of point features from which to calculate distances, expecting
#chromosome in colume 1
#position in column 16

#it will then calculate the mean, 95% CI and N for windows of size specified by distance step
#away from the nearest feature

#it writes out a single file containing the distance categories in the first column
#then sets of 3 columns (mean, ci, n) for each dataset in the original annotated probe report

# Load necessary libraries 
library(dplyr)

# Define input file paths
seq_read_file <- "Input and output/Mouse stuff/abdulkadir mesc trael 50kb windows.txt"
features_file <- "../../Origin calling with OKseqHMM/E14__HMMsegments_IZ.txt"
output_file<-"Input and output/Mouse stuff/abdulkadir mesc trael 50kb windows in 100kb distance.txt"
distance_step<-100000

#note there are more variables defined lower down so that we can process multiple datasets from a single file

# Load input files
seq_read_data_source <- read.table(seq_read_file, header = TRUE,sep = "\t")
features_data_source <- read.table(features_file, header = FALSE, sep = "\t",skip=1)

# Extract relevant columns
seq_read_data <- seq_read_data_source %>%
  select(Chromosome, Start, End)
colnames(seq_read_data)<-c("chromosome","start","end")
features_data <- features_data_source %>%
  select(chromosome = V1, feature_position = V16)

# Calculate mid-points of windows
seq_read_data <- seq_read_data %>%
  mutate(mid_point = (start + end) / 2)

# Initialize a vector to store distances
distances <- numeric(nrow(seq_read_data))

# Calculate distance to the nearest feature for each window
for (i in 1:nrow(seq_read_data)) {
  chr <- seq_read_data$chromosome[i]
  mid_point <- seq_read_data$mid_point[i]
  
  # Get feature positions for the same chromosome
  feature_positions <- features_data %>%
    filter(chromosome == chr) %>%
    pull(feature_position)
  
  # Calculate the distance to the nearest feature
  distances[i] <- min(abs(mid_point - feature_positions))
}

# Add distances to the data frame
seq_read_data <- seq_read_data %>%
  mutate(distance_to_feature = distances)

# Define distance categories
distance_categories <- seq(0, 10000000-distance_step, by = distance_step)+1

output_table_full<-data.frame(Distance_Category = distance_categories)

for(i in 13:ncol(seq_read_data_source)-12) {
  counts_column<-i+12
#  output_file <- paste(output_path,colnames(seq_read_data_source)[counts_column],".txt", sep="")
#  print(output_file)

seq_read_data$read_count<-seq_read_data_source[,counts_column]

# Initialize a list to store output
output_list <- list()

# Calculate statistics for each distance category
for (d in distance_categories) {
  windows_in_category <- seq_read_data %>%
    filter(distance_to_feature>d & distance_to_feature<d+distance_step-1)
  
  mean_read_count <- mean(windows_in_category$read_count, na.rm = TRUE)
  sd_read_count <- sd(windows_in_category$read_count, na.rm = TRUE)
  n_windows <- nrow(windows_in_category)
  ci <- qt(0.975, df = n_windows - 1) * (sd_read_count / sqrt(n_windows))
  
  output_list[[as.character(d)]] <- c(mean_read_count, sd_read_count, ci, n_windows)
}

# Convert output list to data frame
output_table <- do.call(rbind, output_list)
colnames(output_table) <- c("Mean_Read_Count", "SD_Read_Count", "CI", "Num_Windows")
output_table <- as.data.frame(output_table)
output_table$Distance_Category <- distance_categories

#output_table<-output_table %>% select("Distance_Category","Mean_Read_Count", "CI", "Num_Windows","SD_Read_Count")
output_table<-output_table %>% select("Mean_Read_Count", "CI", "Num_Windows")
colnames(output_table)<-c(colnames(seq_read_data_source)[counts_column],colnames(seq_read_data_source)[counts_column],colnames(seq_read_data_source)[counts_column])

output_table_full<-cbind(output_table_full,output_table)

# Write output table to file
}
write.table(output_table_full, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

