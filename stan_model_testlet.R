rm(list=ls())
library(rstan)
library(reshape2)
require(data.table)
library(dplyr)
load("/home/crweber/Dropbox/Supreme Court Data/SCDB_2018_02_justiceCentered_Citation.Rdata")
dat<-SCDB_2018_02_justiceCentered_Citation
dat$y<-car::recode(dat$vote, "1=1; 2=0; 3:4=1; else=NA")
### Structure Data --  
dat$term<-tstrsplit(dat$caseId, "-", fixed=TRUE)[[1]]
dat$case_term_id<-tstrsplit(dat$caseId, "-", fixed=TRUE)[[2]]
### This creates an ID matrix, I don't think justices have an ID, so I gave them one ####
### This needs to be used after retrieving posteriors ####
JusticeID = dat %>% arrange(term, justiceName) %>% 
  mutate(justiceID=as.numeric(as.factor(justiceName))) %>% 
  group_by(caseId) %>% 
  mutate(un=mean(y)) %>%
  filter(un<1) %>%
  ungroup() %>%
  arrange(term, justiceName) %>% 
  group_by(justiceName, justiceID) %>%
  summarise()

# Misc cats 13 and 14 

Vote =  dat %>% arrange(justiceName) %>% 
  mutate(justiceID=as.numeric(as.factor(justiceName))) %>% 
  group_by(caseId) %>% 
  mutate(un=mean(y)) %>%
  filter(un<1) %>%
  ungroup() %>%
  group_by(caseId) %>% 
  mutate(un=mean(y)) %>%
  filter(un<1) %>%
  ungroup() %>%
  arrange(justiceName) %>%
  #mutate(Type=ifelse((issueArea==14 | issueArea==13), 6, issueArea )) %>%
  # Here is the recode -- we'll need to turn it into something less silly.
  mutate(Type= recode(issueArea, `1`=1, `2`=2,
                      `3`=3, `4`=4, `5`=5,
                      `6`=6, `7`=7, `8`=8,
                      `9`=9, `10`=10, `11`=11, `12`=12,
                      `13`=13, `14`=13)) %>%
  mutate(JustType=paste0(as.character(justiceID),"-",as.character(Type))) %>%
  mutate(JustTypeID=as.numeric(as.factor(JustType))) %>%
  mutate(Case=caseId) %>%
  mutate(CaseID=as.numeric(as.factor(caseId))) %>%
  arrange(justiceID, CaseID, JustTypeID) %>%
  select(c("y", "JustTypeID", "CaseID", "justiceID"))


Vote %>% 
  arrange(CaseID, justiceID, JustTypeID)

stan.dat<- list(y=Vote$y,
                Justice=Vote$justiceID,
                Type=Vote$JustTypeID,
                Case=Vote$CaseID,
                ljust=max(Vote$justiceID),
                ltype=max(Vote$JustTypeID),
                lcase=max(Vote$CaseID),
                N=length(Vote$y)
)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

### This is just a 2PL IRT model, with an effect for testlets
two.pl<-"
data {
  int<lower=1> ljust;               // Number of justices
  int<lower=1> ltype;               // Number of type-justice combinations
  int<lower=1> lcase;               // Number of cases
  int<lower=1> N;                   // Number of observations Case x Justice
  int<lower=1, upper=N> Justice[N]; // Justice Indicator
  int<lower=1, upper=N> Type[N];    // Justice-Type Indicator
  int<lower=1, upper=N> Case[N];    // Case indicator
  int<lower=0, upper=1> y[N];       // Data
}
parameters {
  vector[lcase]          alpha;   // discrimination/loading/slope
  vector[lcase]          beta_free;    // difficulty/intercept
  vector[ljust]          theta;   // judicial ideology
  vector[ltype]          testlet; // testlet effect
  real<lower=0> sigma_beta;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_testlet[ltype];
  real mu;
}

transformed parameters{
 vector[lcase] beta;
 for(i in 1:(lcase-1)) beta[i]=beta_free[i];
 beta[lcase] = -1*sum(beta_free);         //This constrains the difficulty parameters to sum to zero
}

model {
  alpha      ~ normal(0, sigma_alpha);
  beta_free  ~ normal(0,sigma_beta);
  sigma_alpha~ cauchy(0, 5);
  sigma_beta ~ cauchy(0, 5);
for (i in 1:ltype){
  testlet[i] ~ normal(0, sigma_testlet[i]);
  sigma_testlet[i] ~ cauchy(0, 5);
}
  theta ~ normal(mu, 1);

for (i in 1:N) y[i] ~ bernoulli_logit(alpha[Case[i]]*theta[Justice[i]]-beta[Case[i]]+testlet[Type[i]]);
}
generated quantities{
   vector[N] log_lik;  // Likelihood
   vector[N] r_post;
   vector[N] p_post;
   vector[N] prior_pred;
 for (i in 1:N) {
  r_post[i]=binomial_rng(y[i], inv_logit(alpha[Case[i]]*theta[Justice[i]]-beta[Case[i]]+testlet[Type[i]]));
  log_lik[i]=bernoulli_logit_lpmf(y[i] | alpha[Case[i]]*theta[Justice[i]]-beta[Case[i]]+testlet[Type[i]]);
  }
 }
 "



### This takes some time.
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fit <- stan(model_code =two.pl, 
            data = stan.dat, iter=2000, chains=1,
            pars=c("theta", "alpha", "beta",  "testlet"),
            control=list("max_treedepth"=15, "adapt_delta"=.99))

model<-stan_model(model_code=two.pl)
fit<-vb(model,  stan.dat, output_samples=2000, tol_rel_obj=.001, iter=10000)


# Sampling - memory issues I don't understand. May want to try HPC
#fit <- stan(model_code =two.pl, 
#            data = stan.dat, iter=2000, chains=3,
#            control=list("max_treedepth"=15, "adapt_delta"=.99))






