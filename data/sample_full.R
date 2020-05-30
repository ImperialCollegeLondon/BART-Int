populationData <- read.csv("data/ss16pil.csv")
# Employment","Sex","Education","Disability","Health_insurance","Own_child","Race"
cols <- c("ESR", "SEX", "SCHL", "DIS", "HICOV", "OC", "RAC1P", "PINCP")
income_col <- "PINCP"
Total_person_income <- populationData[, income_col]
populationData <- populationData[,cols]
populationData[,1:7] <- sapply(populationData[, 1:7], as.factor)
colnames(populationData) <- c("Employment","Sex","Education","Disability","Health_insurance","Own_child","Race", "total_income")
populationData <- populationData[complete.cases(populationData),]
num_training <- 10000
set.seed(123)
training_indices <- sample.int(nrow(populationData), num_training)
train <- populationData[training_indices,]
candidate <- populationData[-training_indices,]
write.csv(train, "data/train2.csv")
write.csv(candidate, "data/candidate2.csv")

