set.seed(102725)
library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')

#######################################################
# Investigating the reason for log-transformation of  #
# the protein levels in biomarker-raw.csv             #
#######################################################

# Load biomarker-raw dataset
biomarker_raw <- read.csv("data/biomarker-raw.csv")
ncol(biomarker_raw) 

# there are 1320 proteins. Let's make a histogram to visualize 
# the distribution of some of their raw values.
protein_sample <- sample(colnames(biomarker_raw), 3)

biomarker_sample1 <- biomarker_raw |>
  select(protein_sample[1]) |>
  slice(-1) |>
  pull() |>
  as.numeric()

biomarker_sample2 <- biomarker_raw |>
  select(protein_sample[2]) |>
  slice(-1) |>
  pull() |>
  as.numeric()

biomarker_sample3 <- biomarker_raw |>
  select(protein_sample[3]) |>
  slice(-1) |>
  pull() |>
  as.numeric()

# make histograms of 3 randomly selected proteins
hist(biomarker_sample1, breaks = 15,
     xlab = protein_sample[1],
     main = "Distribution of Randomly Sampled Protein 1")
hist(biomarker_sample2, breaks = 15,
     xlab = protein_sample[2],
     main = "Distribution of Randomly Sampled Protein 2")
hist(biomarker_sample3, breaks = 15,
     xlab = protein_sample[3],
     main = "Distribution of Randomly Sampled Protein 3")

# All three histograms reveal that proteins are skewed to the right
# that is, they have clusters of values on the lower end of measurement
# with a long tail of smaller frequency high value measurements
# log transformation of this will rein in outliers
# and make the distribution of sampled proteins more Gaussian
# which is critical for the analysis procedures such as the t-test.


#######################################################
# Investigating outlier trends in biomarker-raw.csv   #
# and seeing outlier frequency in ASD and TD groups   #
#######################################################

# This bit of code transforms the raw data up until the trimming
biomarker_outliers <- read_csv('data/biomarker-raw.csv', 
         skip = 2,
         col_select = -2L,
         col_names = c('group', 
                       'empty',
                       pull(var_names, abbreviation),
                       'ados'),
         na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # log transform, center and scale, and trim
  mutate(across(.cols = -c(group, ados), 
                ~ scale(log10(.x))[, 1])) %>%
  # reorder columns
  select(group, ados, everything())

# comparison function for filtering the data
temp <- function(datacol) {
  return(datacol > 3)
}

# Filter any rows that don't have at least 1 outlier
# then calculate each individual's total no. of outliers
biomarker_outliers <- biomarker_outliers |>
  filter(if_any(c(-group, -ados), temp)) |>
  rowwise() |>
  mutate(
    n_indv_outliers = sum(
      across(c(-group, -ados)) > 3
    )
  )

# Sorting from most to least outliers
biomarker_outliers_arranged <- biomarker_outliers |>
  arrange(-n_indv_outliers) |>
  select(group, ados, n_indv_outliers)

# Make histograms for the distribution of outliers in ASD and TD groups
biomarker_outliers_arranged |>
  filter(group == "ASD") |>
  ggplot(aes(x = n_indv_outliers)) +
  geom_histogram(color = "white", bins = 20) +
  labs(x = "No. of outliers", y = "Count",
       title = "Distribution of outliers for ASD group") +
  theme_minimal()

biomarker_outliers_arranged |>
  filter(group == "TD") |>
  ggplot(aes(x = n_indv_outliers)) +
  geom_histogram(color = "white", bins = 20) +
  labs(x = "No. of outliers", y = "Count",
       title = "Distribution of outliers for TD group") +
  theme_minimal()

# It seems like for both the ASD and TD groups, the distribution of 
# outliers is around the same.
# The TD group in particular has a higher frequency of large outliers > 100
# than the ASD group

# Quick investigation of the arranged outliers
head(biomarker_outliers_arranged, 15)
# It seems like most of the outliers are in the TD group.
# individuals with large numbers of outliers could have had a 
# misread protein panel as well, or some sort of lab equipment malfunction
# especially the individual with 126 outliers?



# methods
# use a fuzzy intersection instead of a hard intersection 
# to combine the sets of top predictive proteins across selection methods

## MULTIPLE TESTING
####################
# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out_fuzzy <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out_fuzzy$confusion

# compute importance scores
proteins_s2_fuzzy <- rf_out_fuzzy$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out_fuzzy$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# MAJOR CHANGE: instead of intersection, we will use union
# for the "fuzzy intersection"
proteins_sstar_fuzzy <- union(proteins_s1, proteins_s2_fuzzy)

biomarker_sstar_fuzzy <- biomarker_clean %>%
  select(group, any_of(proteins_sstar_fuzzy)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split_fuzzy <- biomarker_sstar_fuzzy %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split_fuzzy), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

testing(biomarker_split_fuzzy) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(pred.fac = as.factor(pred > 0.5),
         truth.fac = as.factor(class)) |>
  class_metrics(estimate = pred.fac,
                truth = truth.fac, pred,
                event_level = 'second')

# We obtain pretty bad metric results for the regular logistic regression
# compared to the in-class analysis.
# We will use cv.glmnet to cross-validate the model and apply elastic net
# since perhaps the larger number of predictive variables is harming
# the model's fit.

# Obtain training X & Y matrices for cv.glmnet
biomarker_train_x <- biomarker_split_fuzzy |>
  training()|>
  select(-class) |>
  as.matrix()
biomarker_train_y <- biomarker_split_fuzzy |>
  training() |>
  pull(class) |>
  as.factor()

# Obtain testing X matrix for cv.glmnet prediction
biomarker_test_x <- testing(biomarker_split_fuzzy) |>
  select(-class) |>
  as.matrix()

# Fit logistic regression model using cv.glmnet
fit_fuzzy_en <- glmnet::cv.glmnet(x = biomarker_train_x,
               y = biomarker_train_y,
               family = "binomial",
               alpha = 0.1)

predicted_fuzzy_en <- predict(fit_fuzzy_en,
        newx = biomarker_test_x,
        type = "response",
        s = "lambda.min") 

testing(biomarker_split_fuzzy) |>
  mutate(pred = as.numeric(predicted_fuzzy_en),
         pred.fac = as.factor(pred > 0.5),
         truth.fac = as.factor(class)) |>
  class_metrics(estimate = pred.fac,
                truth = truth.fac, pred,
                event_level = "second")
# these estimates are better than the initial model with all variables.
# but are still worse than the in-class analysis.

# results
# Use any method to find either:
  # a simpler panel that achieves comparable classification accuracy
  # an alternative panel that achieves improved classification accuracy
# Benchmark your results against the in-class analysis.