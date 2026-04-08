# breast cancer type prediction from micro-array gene samples
# we use a library for reading the csv file because read_csv from the r standard library cant read csv files this large
library("data.table") 
#https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida
data <- fread("D:\\Projects\\uni\\masters_notes\\2026_spring\\Statistical Inference for High Dimensional\\project\\Breast_GSE45827.csv") #https://sbcb.inf.ufrgs.br/data/cumida/Genes/Breast/GSE45827/Breast_GSE45827.csv
x <- as.matrix(data[, 3:54677])

# we use z-score normalization to standardize our data 
x <- scale(x)

# the target vector y (data$type) is a categorical variable over the possible breast cancer subtypes
#  we therefore convert the possible cancer types to a multivariate uniform variable (one-hot-encoding), so that we can predict each type separately 
y <- as.data.frame(data$type)
y$basal <- ifelse(y$`data$type` == "basal", 1, 0)
y$HER <- ifelse(y$`data$type` == "HER", 1, 0)
y$cell_line <- ifelse(y$`data$type` == "cell_line", 1, 0)
y$normal <- ifelse(y$`data$type` == "normal", 1, 0)
y$luminal_A <- ifelse(y$`data$type` == "luminal_A", 1, 0)
y$luminal_B <- ifelse(y$`data$type` == "luminal_B", 1, 0)

# the basal breast cancer subtype has the most occurrences (and thus makes our dataset the least unbalanced) so we use this as our target variable
y$target = y$normal

# we split the data into 80% used for training and 20% used for testing
set.seed(1)
train <- sample(1:nrow(x), round(nrow(x) * 0.8))
test <- (-train)

# to start off we do least squares regression as a baseline, we will evaluate it using k-fold cross validation since we have a very low amount of samples 
library(glmnet)

# we do 5-fold cross validation to estimate the out-of-sample error
#k <- 5
#errors <- numeric(k)
#for (i in 0:(k-1)) {
#  current_fold.test <- train[floor(length(train)*i/k+1):floor(length(train)*(i+1)/k)]
#  current_fold.train <- (-current_fold.test)
#  least_squares <- lm(y$normal ~ x, subset=current_fold.train)
#  
#  current_fold.test.x = x[current_fold.test, ]
#  current_fold.test.y = y$normal[current_fold.test]
#  
#  current_fold.test.pred <- predict(least_squares, newdata = as.data.frame(current_fold.test.x))
#  error <- mean((current_fold.test.pred[current_fold.test] - current_fold.test.y)^2)
#  errors[i + 1] <- error 
#}
#mean_cve <- mean(errors)
#print(mean_cve)

# we can then fit on our entire training data
# df <- data.frame(y=y$target, x)
# lsquares.fit <- lm(y ~ ., data=df[train])
lsquares.fit <- glmnet(x[train, ], y$target[train], alpha=0, lambda=0)

# next we will do ridge regression 
#  to select our tuning parameter lambda for the prior l2 norm term we will use cv.glmnet to calculate 
#  the mean cross validation error for different values of lambda, we will then choose the lambda that gives us the 
#  lowest cross validation error 
set.seed(1)
wide_grid <- 10^seq(-5, 1, length=200)
cv.out <- cv.glmnet(x[train, ], y$target[train], nfolds=5, alpha=0, lambda=wide_grid)
plot(cv.out)
ridge.min_lambda <- cv.out$lambda.min
ridge.min_lambda

# we can then fit on all our data
ridge.fit <- glmnet(x[train, ], y$target[train], alpha=0, lambda=ridge.min_lambda)


# next we will do lasso regression
#  to select our tuning parameter lambda for the prior l1 norm term we will use adaptive validation 
#  this should give us better support recovery, but for prediction cross validation should perform better
l <- function (x) max(abs(x[-1, ])) # the maximum norm is used for lasso adaptive validation
                                    #  this is because we are trying to do support recovery and 
                                    #  because we use the maximum we can capture all "reasonably large predictors"
                                    #  see fundamentals of high dimensional statistics example 4.4.1 & chapter 7.4
max_lambda <- l(t(x) %*% y$target)/dim(x)[1] # we choose the maximum lambda to be the smallest lambda 
                                             #  for which all predictors are 0, 
                                             #  see "A Practical Scheme and Fast Algorithm to Tune the Lasso With Optimality Guarantees" page 3
wide_grid <- 10^seq(-5, log(max_lambda, 10), length=200) # set of potential tuning parameters
beta_hats <- glmnet(x[train, ], y$target[train], alpha=1, lambda=wide_grid) # compute all betas ahead of time 

r <- max(wide_grid)
factor <- 3/(4 * dim(x)[1]) # 3/(4n) is a good approximation of the factor for the lasso error bound
                            #  according to fundamentals of high dimensional statistics example 4.4.1 (bottom of page 126)
r_av_found <- FALSE
r_test <- 0

while ((r != min(wide_grid)) && (!r_av_found)) {
  #beta_hat_r <- glmnet(x[train, ], y$target[train], alpha=1, lambda=r)
  beta_hat_r <- coef(beta_hats, s=r)
  r_prime <- max(wide_grid)
  message(sprintf("r: %s\n", r))
  while ((r_prime > r) && (!r_av_found)) {
    #beta_hat_r_prime <- glmnet(x[train, ], y$target[train], alpha=1, lambda=r_prime)
    beta_hat_r_prime <- coef(beta_hats, s=r_prime)
    
    #message(sprintf("r: %s\n", r))
    #message(sprintf("r_prime: %s\n", r_prime))
    #message(sprintf("l(coef(beta_hat_r_prime) - coef(beta_hat_r)): %s\n", l(beta_hat_r_prime - beta_hat_r)))
    #message(sprintf("factor * r_prime + factor * r: %s\n", factor * r_prime + factor * r))
    if (l(beta_hat_r_prime - beta_hat_r) > factor * r_prime + factor * r) {
      r <- min(wide_grid[wide_grid > r])
      r_av_found <- TRUE #exit from both loops
      r_test <- r
    }
    r_prime <- max(wide_grid[wide_grid < r_prime])
  }
  r <- max(wide_grid[wide_grid < r])
}
r_av <- r
r_av

# now that we have our tuning parameter r_av we can fit on all our training data
lasso.av.fit <- glmnet(x[train, ], y$target[train], alpha=1, lambda=r_av)


# let's try doing cross validation too to compare, this should give us better prediction 
set.seed(1)
wide_grid <- 10^seq(-5, 1, length=200)
cv.out <- cv.glmnet(x[train, ], y$target[train], nfolds=5, alpha=1, lambda=wide_grid)
plot(cv.out)

lasso.cv.min_lambda <- cv.out$lambda.min
lasso.cv.min_lambda

# we can then fit on all our data
lasso.cv.fit <- glmnet(x[train, ], y$target[train], alpha=1, lambda=lasso.cv.min_lambda )


# --- prediction ---
# finally we calculate the out-of-sample error for the different methods on the test set we withheld at the start.
# least squares
lsquares.pred <- predict(lsquares.fit, s=0, newx=x[test, ])
mean((lsquares.pred - y$target[test])^2)

# ridge regression
ridge.pred <- predict(ridge.fit, s=ridge.min_lambda, newx=x[test, ])
mean((ridge.pred - y$target[test])^2)

# lasso adaptive validation
lasso.av.pred <- predict(lasso.av.fit, s=r_av, newx=x[test, ])
mean((lasso.av.pred - y$target[test])^2)

# lasso cross validation 
lasso.cv.pred <- predict(lasso.cv.fit, s=lasso.cv.min_lambda, newx=x[test, ])
mean((lasso.cv.pred - y$target[test])^2)

# --- support recovery ---
# let's look at our support recovery, which predictors did the different methods find and are these well known predictors?

lsquares.predictors <- which(abs(coef(lsquares.fit)) > 0.05)[-1]
ridge.predictors <- which(abs(coef(ridge.fit)) > 0.001)[-1]
lasso.cv.predictors <- which(coef(lasso.cv.fit) != 0)[-1]
lasso.av.predictors <- which(coef(lasso.av.fit) != 0)[-1]

lsquares.predictors
ridge.predictors
lasso.cv.predictors
lasso.av.predictors

# and the overlapping predictors
Reduce(intersect, list(v1 = lasso.cv.predictors, v2=lasso.av.predictors))


# let's try doing least squares using only the predictors recovered from lasso with tuning parameter chosen by adaptive validation 

lsquares.lav.fit <- glmnet(x[train, lasso.av.predictors], y$target[train], alpha=0, lambda=0)
lsquares.lav.pred <- predict(lsquares.lav.fit, s=0, newx=x[test, lasso.av.predictors])
mean((lsquares.lav.pred - y$target[test])^2)

recovered_predictors <- names(data[1, 3:54677])[lasso.av.predictors]
recovered_predictors

# to find out which genes these probe identifiers correspond to we have to crossreference the annotation table
gpl <- fread("D:\\Projects\\uni\\masters_notes\\2026_spring\\Statistical Inference for High Dimensional\\project\\GPL570_limpo.txt") #https://sbcb.inf.ufrgs.br/cumida

genes <- gpl[ID %in% recovered_predictors, c("Gene Symbol", "Gene Title")]
genes
# thus we find a correlation between the presence of the genes MOCS2, AKIRIN1 and the presence of basal-subtype breast cancer
# thus we find a correlation between the presence of the genes SLC34A3 , LOC100505609, CETP and breast cancer





