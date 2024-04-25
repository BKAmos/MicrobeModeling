# Load Packages of interest
library(bnstruct)
library(caret)
library(data.table)
library(rsample)
library(dbnR)
library(dplyr)

###### Below is isolate level DBN learning structure across all time points at the same time
#load the time points and the layer files that are associated with them at the isolate level as well as loading the header and data files in as dataframes for validation # nolint
timepoint <- BNDataset('all_data.txt', '/isolate_complete_header.txt')


# Load Layer File
layersTP <- scan(file = "BNStructData/isolate_complete_layers.txt")


# Load matrix for layer definition in Network creation
layerStruct2 <- matrix(c(0,0,1,0),2,2)


# Load the data as a matrix
tpData = read.table('/all_data_index.txt',sep=' ')


# Add a timewindow column with initial values set to 1
tpData <- tpData %>%
  mutate(timeWindow = 1)

# Use cumsum() to increment the new column value each time the same value is seen in the other column
tpData <- tpData %>%
  group_by(V1) %>%
  mutate(timeWindow = cumsum(timeWindow))


# Get unique integers from the column
unique_integers <- unique(tpData$V1)

# Remove the first four batches
unique_integers <- tail(unique_integers, length(unique_integers) - 4)

list_of_train_dfs <- list()
list_of_test_dfs <- list()

for (i in 1:10) {
# Step 2: Calculate Number of Integers for Training and Testing
  n_train <- ceiling(0.8 * length(unique_integers))
  n_test <- length(unique_integers) - n_train


# Step 3: Randomly Sample Integers for Training and Testing
  # set.seed(123)  # for reproducibility
  train_integers <- sample(unique_integers, size = n_train, replace = FALSE)
  test_integers <- sample(setdiff(unique_integers, train_integers), size = n_test, replace = FALSE)


# Filter rows where the values in the specified column match the integers in the list
  test_df <- tpData[tpData$V1 %in% test_integers, ]
  train_df <- tpData[tpData$V1 %in% train_integers, ]
  
  list_of_test_dfs <- append(list_of_test_dfs, list(test_df))
  list_of_train_dfs <- append(list_of_train_dfs, list(train_df))
  
}

# Read in data header file
tpDataHeader = read.table('isolate_complete_header.txt', sep=' ')

# Format data approapriately for the BNDataset constructor based on the header file
names = as.character(tpDataHeader[1,])
cardinality = as.numeric(tpDataHeader[2,])
discreteVec = as.character(tpDataHeader[3,])

acc1 <- c()
list_of_matrices <- list()

for (i in seq_along(list_of_train_dfs)) {
  train_data <- list_of_train_dfs[[i]]
  # print(dim(train_data))
  train_data <- train_data[, !(names(train_data) %in% c("V1", "timeWindow"))]
  # print(dim(train_data))
  trainData <- BNDataset(train_data, discreteVec, names, cardinality)
  dbntrain <- learn.network(trainData, algo='mmhc', layering=layersTP, layer.struct=layerStruct2, alpha = .05, max.parents = 3)
  list_of_matrices[[i]] <- dbntrain@dag
 
  for (q in 1:100) {
    for (j in seq_along(list_of_test_dfs)) {
      test_data <- list_of_test_dfs[[i]]
      test_data_tp <- test_data[test_data$timeWindow %in% j,]
      # print(test_data_tp$V1)
      # print(test_data_tp$timeWindow)
      test_data_tp <- test_data_tp[, !(names(test_data_tp) %in% c("V1", "timeWindow"))]
      test_data_tp <- BNDataset(test_data_tp, discreteVec, names, cardinality)
      dbn1test = learn.params(dbntrain, test_data_tp, use.imputed.data = F)
      # print(test_data_tp)
      enginetest1 <- InferenceEngine(dbn1test)
      enginetest1 <- belief.propagation(enginetest1)
      new.net.test1 <- updated.bn(enginetest1)
      mpv1TestInf <- get.most.probable.values(new.net.test1)
      acc1 <- c(acc1, mpv1TestInf)
  
  
    }
  }
}



# Example list
# dbntraincpts <- dbntrain@cpts
# testCPTS <- melt(dbntraincpts)
# 
# write.table(testCPTS, file = "", row.names = FALSE, col.names = TRUE, quote = TRUE, sep = ',')


output_directory <- '/graphs/'
# file_name <- paste0(output_directory, "batch234_maxP2_smallalpha.csv")
# write.csv(dbntrain@dag, file=file_name)



write.table(acc1, 'all_10foldCV_inf_std_noObs.txt', row.names = F, col.names = F)

for (i in seq_along(list_of_matrices)) {
# Generate a unique file name
  file_name <- paste0(output_directory, "matrix_", i, ".txt")

  # Write the matrix to the file
  write.table(list_of_matrices[[i]], file = file_name, sep = ",", row.names = FALSE, col.names = FALSE)
}
