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
layersTP <- scan(file = "/isolate_complete_layers.txt")


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
# unique_integers <- tail(unique_integers, length(unique_integers) - 4)

list_of_train_dfs <- list()
list_of_test_dfs <- list()

# for (i in 1:10) {
# Step 2: Calculate Number of Integers for Training and Testing
  # n_train <- ceiling(0.8 * length(unique_integers))
  # n_test <- length(unique_integers) - n_train


# Step 3: Randomly Sample Integers for Training and Testing
  # set.seed(123)  # for reproducibility
  # train_integers <- sample(unique_integers, size = n_train, replace = FALSE)
  # test_integers <- sample(setdiff(unique_integers, train_integers), size = n_test, replace = FALSE)


# Filter rows where the values in the specified column match the integers in the list
b1 <- tpData[tpData$V1 == 1,]
b2 <- tpData[tpData$V1 == 2,]
b3 <- tpData[tpData$V1 == 3,]
b4 <- tpData[tpData$V1 == 4,]

test_b1_g <- tpData[tpData$V1 >= 5 & tpData$V1 <= 14, ]
test_b2_g <- tpData[tpData$V1 >= 45 & tpData$V1 <= 54, ]
test_b3_g <- tpData[tpData$V1 >= 85 & tpData$V1 <= 94, ]
test_b4_g <- tpData[tpData$V1 >= 125 & tpData$V1 <= 134, ]

train_b234_g <- tpData[tpData$V1 >= 15 & tpData$V1 <= 44, ]
train_b134_g <- tpData[tpData$V1 >= 55 & tpData$V1 <= 84, ]
train_b124_g <- tpData[tpData$V1 >= 95 & tpData$V1 <= 124, ]
train_b123_g <- tpData[tpData$V1 >= 135 & tpData$V1 <= 164, ]

test_b1_c <- rbind(b1, test_b1_g)
test_b2_c <- rbind(b2, test_b2_g)
test_b3_c <- rbind(b3, test_b3_g)
test_b4_c <- rbind(b4, test_b4_g)

train_b234_c <- rbind(b2, b3, b4, train_b234_g)
train_b134_c <- rbind(b1, b3, b4, train_b134_g)
train_b124_c <- rbind(b1, b2, b4, train_b124_g)
train_b123_c <- rbind(b1, b2, b3, train_b123_g)

list_of_test_dfs <- append(list_of_test_dfs, list(test_b1_c))
list_of_test_dfs <- append(list_of_test_dfs, list(test_b2_c))
list_of_test_dfs <- append(list_of_test_dfs, list(test_b3_c))
list_of_test_dfs <- append(list_of_test_dfs, list(test_b4_c))

list_of_train_dfs <- append(list_of_train_dfs, list(train_b234_c))
list_of_train_dfs <- append(list_of_train_dfs, list(train_b134_c))
list_of_train_dfs <- append(list_of_train_dfs, list(train_b124_c))
list_of_train_dfs <- append(list_of_train_dfs, list(train_b123_c))

# }

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
    for (j in 1:10) {
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
# write.table(testCPTS, file = "/Users/bkamos/Documents/GitHub/InRoot/timeSeriesxLotus/RTesting/isolate_testing_95/cptsForSmallGraph/dbntraincptsTest.csv", row.names = FALSE, col.names = TRUE, quote = TRUE, sep = ',')


output_directory <- 'graphs/'
#
write.table(acc1, 'Inference/all_10foldCV_inf_std_noObs.txt', row.names = F, col.names = F)

for (i in seq_along(list_of_matrices)) {
# Generate a unique file name
  file_name <- paste0(output_directory, "matrix_", i, ".txt")

  # Write the matrix to the file
  write.table(list_of_matrices[[i]], file = file_name, sep = ",", row.names = FALSE, col.names = FALSE)
}
