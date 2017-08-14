
# funnel plots for rates

# required packages / functions
library(tidyverse)
library(broom)
source("./roundNiceNumber.r")

# usage - input is a df
# i.e. your call should look like this ...
# funl_Data(df.data, "practice", "locality", "O", "n", "rt"
#          , target = NULL , smoothness = 100, fnlMinEvents = NULL, fnlMaxEvents = NULL)

# @param unit - institution or geography
# @param group - groups
# @param O - observed count
# @param n - denominator
# @param rt - rate (observed count (O) / denominator (n) )
# @param target - specify a mean/target value (optional)
# @param smoothness - higher value returns smoother funnels (optional - defaults to 100)
# @param fnlMinEvents - lower limit for x-axis extent of funnels (optional)
# @param fnlMaxEvents - upper limit for y-axis extent of funnels (optional)


funl_Data <- function(df.data, col.unit, col.group, col.O, col.n, col.rt
                      , target = NULL, smoothness = 100, fnlMinEvents = NULL, fnlMaxEvents = NULL, ... ) {
  
  # as a minimum col.unit, col.O and one of col.n or col.rt are required
  # should throw an error if this happens
  # this needs some work!
  if(missing(col.rt)) {
    stop('col.rt is missing from your input df');
  }
  
  # if col.rt is not supplied then caculate it
  # col.rt = col.O / col.n
  
  # if a value for target is not supplied then use global (weighted) mean
  if (is.null(target)) {
      target <- sum(df.data[[col.O]], na.rm = TRUE) / sum(df.data[[col.O]] / df.data[[col.rt]], na.rm = TRUE)
  } else {
    target = target
  }
  
  # implement method for selecting smoothness factor if none supplied?

  # determines x-axis extent of funnels
  # calculate values for fnlMinEvents and fnlMaxEvents (if none provided in call)
  if (is.null(fnlMinEvents)) {
    fnlMinEvents <- rounddown_Nice(min(df.data[[col.O]], na.rm = TRUE))
  }
  if (is.null(fnlMaxEvents)) {
    fnlMaxEvents <- roundup_Nice(max(df.data[[col.O]], na.rm = TRUE))
  }
  
  # produce funnels
  if (smoothness > diff(range(fnlMinEvents, fnlMaxEvents))) {
      
    # if value for smoothness supplied in call to funl_Data (default = 100)
    # is greater than maxEvents - minEvents
    # then recalibrate smoothness to maxEvents - minEvents
    events <- seq.int(fnlMinEvents, fnlMaxEvents)
    
    } else {
    
    # else initialise empty vector to hold events and set first element to minEvents
    events    <- rep(NA_real_, smoothness)
    events[1] <- fnlMinEvents
    
    # for loop to fill events vector where the difference to the preceeding elements value gradually increases
    for (i in seq_len(smoothness)[-1]) {
      
      events[i] <- max(round((fnlMaxEvents / events[i - 1])^(1 / ((smoothness + 1) - i)) * events[i - 1]), events[i - 1])
    
    }
    
  }
  
  # derive denominator n
  n <- events / target

  # nested df (events series and confidence levels)
  fnlInput <- expand.grid(fnlLimit = c("threeSigma", "twoSigma")
                   , events = events
                   , stringsAsFactors = FALSE) %>% 
    mutate(conf.level = ifelse(fnlLimit == "threeSigma", 0.99, 0.95)) %>% 
    arrange(fnlLimit) %>% 
    group_by(fnlLimit) %>%
    nest(events, conf.level) %>% 
    ungroup()
  
  # calculate confidence intervals with exact Poisson test
  # poisson.test - x must be x must be finite, non-negative, and integer  
  # apply poisson.test to each combinaton of event values and confidence levels in nested df
  # use broom::tidy to store results in dfs
  fnlOutput <- fnlInput %>%
    group_by(fnlLimit) %>%
    mutate(
      pois = map(data, ~ map2_df(.$events, .$conf.level, function(x, y)
        tidy(poisson.test(x, alternative = "two.sided", conf.level = y))))
    ) %>%
    ungroup() %>% 
    unnest(pois)
 
  # tidy fnl df
  funls <- fnlOutput %>%
    select(estimate, conf.low, conf.high, fnlLimit) %>%
    mutate(
      target = target
      , n = estimate / target
      , fnlLow = conf.low / n 
      , fnlHigh = conf.high /n
    ) %>%
    select(-conf.low, -conf.high) %>% 
    rename(
      events = estimate
    )

  # create df for units and in-control T/F
  units <- as_tibble(df.data[[col.unit]])
  names(units) <- col.unit
  
  # apply poisson.test to observed event values in input df
  units$threeSigmaLow <- sapply(df.data[[col.O]], function(x) if(is.na(x)) {NA}
                               else{poisson.test(x, alternative = "two.sided", conf.level = 0.99)$conf.int[1]}) / (df.data[[col.O]] /target)
  
  units$threeSigmaHigh <- sapply(df.data[[col.O]], function(x) if(is.na(x)) {NA}
                                else {poisson.test(x, alternative = "two.sided", conf.level = 0.99)$conf.int[2]} ) / (df.data[[col.O]] / target)
  
  units$status <- ifelse(df.data[[col.rt]] < units$threeSigmaLow, "low", ifelse(df.data[[col.rt]] > units$threeSigmaHigh, "high", "in"))
    
  # return list of fnl & units
  funls <- as_tibble(funls)
  results <- list(funls, units)
  return(results)
  
  }

