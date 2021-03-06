# Doing Bayesian Data Analysis Ch 11

## Null hypothesis significance testing

"In NHST, the goal of inference is to decide whether a particular value of a 
parameter can be rejected."

"We might want to know whether a coin is fair -- whether we can reject the 
"null" hypothesis that the bias of the coin has the specific value 
$\theta = 0.50$."


## Testing procedure

1. Find the probabilities of all possible outcomes; i.e. what could have 
happened

2. Find the probability of outcomes at least as extreme as observed outcome -- 
the p-value

e.g. if you fipped a coin 4 times, your possible outcomes are HHHH, HHHT, ..., 
TTTT; $2^4 = 16$ possible outcomes.

If you got HHHH, there is 1 outcome at least as extreme (HHHH), so has p-value 
1/16  = 0.0625

If you got HHHT, there are 4 + 1 outcomes at least as extreme; p = 5/16 = 0.3125


## Dependence on data collection procedure

The set of possible outcomes depends on how the data was collected, the testing 
and stopping criteria; Kruschke calls it 'experimenter intention'.

E.g. flip the coin N times, or for O minutes, or until Pth head appeared, ...

Each of these could lead to different p-values.

**BAD: A stopping criterion based on the observed results will bias the data.**


## Confidence intervals

Confidence intervals contain the range of parameter values that would not be 
rejected. (note: this is one definition)

"The CI is simply the smallest and largest values of $\theta$ that yield 
$p\leq2.5\%$."

CI is not a distribution.

  + "Some theorists have explored normalized p curves, which do integrate to 1, and 
are called confidence distributions."


## Bayesian intervals

The Bayesian highest density interval (HDI) is the smallest set of values of 
$\theta$ that have $P(\theta|D)\geq95\%$.

"those values of $\theta$ that have at least some minimal level of posterior 
credibility"

Direct interpretation wrt probability of parameters given data.

Do not depend on experimenter intentions.


## Multiple testing

"When comparing multiple conditions,a key goal in NHST is to keep the overall 
false alarm rate down to a desired maximum such as 5%. Abiding by this 
constraint depends on the number of comparisons that are to be made."

"In a Bayesian analysis, however, there is just one posterior distribution over
the parameters that describe the conditions."


## Kruschke on intentions vs priors

"The experimenterâs stopping and testing intentions are mysterious, capricious, 
idiosyncratic"

"Prior beliefs are overt, explicitly debated, and founded on publicly accessible
previous research."


## What sampling distributions are good for

"sampling distributions tell us the probabilities of possible data if we run an 
intended experiment given a particular hypothesis"

* planning an experiment

* posterior predictive check
    + "If the posterior parameter values really are good descriptions of the 
    data, then the predicted data from the model should actually look like real 
    data"
