# 5.1
# set up table 5.4
tp_rate <- 0.99
fp_rate <- 0.05
p_disease <- 0.001
mat <- matrix(c(tp_rate * p_disease,       fp_rate * (1 - p_disease), 
                (1 - tp_rate) * p_disease, (1 - fp_rate) * (1 - p_disease)), 
              ncol = 2, byrow = T)
colnames(mat) <- c("disease", "no_disease")
rownames(mat) <- c("pos", "neg")

mat
sum(mat)
apply(mat, 1, sum)
apply(mat, 2, sum)

# probability of disease given a positive test
# p(disease | pos) = p(disease, pos) / p(pos)
# or               = p(pos | disease) * p(disease) / p(pos)
post <- mat["pos", "disease"] / rowSums(mat)["pos"]
post

# if retested, what is the probability of disease if the test is negative?

# set up second table, given that first test was positive
# p(disease | neg, previous pos) = p(pos | disease) * p(disease*) / p(pos)
mat2 <- matrix(c(tp_rate * post,       fp_rate * (1 - post), 
                 (1 - tp_rate) * post, (1 - fp_rate) * (1 - post)), 
               ncol = 2, byrow = T)
colnames(mat2) <- colnames(mat)
rownames(mat2) <- rownames(mat)

mat2
sum(mat2)
apply(mat2, 1, sum)
apply(mat2, 2, sum)

# probability of disease given one positive test and one negative test
post2 <- mat2["neg", "disease"] / apply(mat2, 1, sum)["neg"]
post2

# what about just one negative test
# very, very small
mat["neg", "disease"] / apply(mat, 1, sum)["neg"]

# probability of disease given two positive tests
mat2["pos", "disease"] / apply(mat2, 1, sum)["pos"]

#-------------------------------------------------------------------------------

# 5.2

# A. frequency table
N <- 100000
freq <- mat * N
freq

# B. proportion of people that have the disease given a positive test
99 / (99 + 4995)

# C. retest all the positives out of 10M tested
freq2 <- mat2 * apply(10000000 * mat, 1, sum)["pos"]
freq2
apply(freq2, 1, sum)

# D. proportion of disease-havers out of negative re-testers
1 / (1 + 94905)

#-------------------------------------------------------------------------------

# 5.3

# A. probability of disease given negative test
neg_post <- mat["neg", "disease"] / rowSums(mat)["neg"]
neg_post

# B. probability of disease given positive test following negative test

mat3 <- matrix(c(tp_rate * neg_post,       fp_rate * (1 - neg_post), 
                (1 - tp_rate) * neg_post, (1 - fp_rate) * (1 - neg_post)), 
              ncol = 2, byrow = T)
colnames(mat3) <- colnames(mat)
rownames(mat3) <- rownames(mat)

mat3

mat3["pos", "disease"] / rowSums(mat3)["pos"]

# same as before; order doesn't matter
post2

#-------------------------------------------------------------------------------
