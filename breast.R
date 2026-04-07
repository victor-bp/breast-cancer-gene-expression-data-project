# breast cancer type prediction from micro-array gene samples
# we use a library for reading the csv file because read_csv from the r standard library cant read csv files this large
library("data.table") 
data <- fread("Breast_GSE45827.csv")


x <- model.matrix(type ~. - 1, data[, 1:500])
# we z-score normalization to standardize our data 
x <- scale(x)

# the target vector y (data$type) is a categorical variable over the possible cancer types
#  we therefore convert the possible cancer types to a multivariate uniform variable (one-hot-encoding), so that we can predict each type separately 
y <- as.data.frame(data$type)
y$basal <- ifelse(y$`data$type` == "basal", 1, 0)
y$HER <- ifelse(y$`data$type` == "HER", 1, 0)
y$cell_line <- ifelse(y$`data$type` == "cell_line", 1, 0)
y$normal <- ifelse(y$`data$type` == "normal", 1, 0)
y$luminal_A <- ifelse(y$`data$type` == "luminal_A", 1, 0)
y$luminal_B <- ifelse(y$`data$type` == "luminal_B", 1, 0)

#  the type "normal" describes if breast cancer is present or not, so we use this as our target variable

# we split the data into 80% used for training and 20% used for testing
set.seed(1)
train <- sample(1:nrow(x), round(nrow(x) * 0.8))
test <- (-train)

# to start off we do least squares regression as a baseline, we will evaluate it using k-fold cross validation since we have a very low amount of samples 
library(glmnet)

# we do 10-fold cross validation to estimate the out-of-sample error
k <- 10
errors <- numeric(10)
for (i in 0:(k-1)) {
  current_fold.test <- train[floor(length(train)*i/k+1):floor(length(train)*(i+1)/k)]
  current_fold.train <- (-current_fold.test)
  least_squares <- lm(y$normal ~ x, subset=current_fold.train)
  
  current_fold.test.x = x[current_fold.test, ]
  current_fold.test.y = y$normal[current_fold.test]
  
  current_fold.test.pred <- predict(least_squares, newdata = as.data.frame(current_fold.test.x))
  error <- mean((current_fold.test.pred[current_fold.test] - current_fold.test.y)^2)
  errors[i + 1] <- error 
}
mean_cve <- mean(errors)
print(mean_cve)

# next we will do ridge regression 
#  to select our tuning parameter lambda for the prior l2 norm term we will use cv.glmnet to calculate 
#  the mean cross validation error for different values of lambda, we will then choose the lambda that gives us the 
#  lowest cross validation error 
set.seed(1)
wide_grid <- 10^seq(1, -2, length=300)
cv.out <- cv.glmnet(x[train, ], y$normal[train], nfolds=10, alpha=0, lambda=wide_grid)
plot(cv.out)
 
min_lambda <- cv.out$lambda.min

# next we will do lasso regression
#  to select our tuning parameter lambda for the prior l1 norm term we will use adaptive validation
#  we thus do not need a test set, and we can use our entire training set to hopefully get lower bias.



# finally we calculate the out-of-sample error for the different methods on the test set we withheld at the start.



























