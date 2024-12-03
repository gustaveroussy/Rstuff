## FROM : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6415926/

## R code for estimating Bliss independence interaction index and 95% confidence interval
## Read in original data of Gami+BKM120 (full data are shown in table 1):
## --
##  rep1: Replication 1; rep2: Replication 2
##  DoseA: dose of drug A, DoseB: dose of drug B, y: response

# rep1 <- read.table(header = TRUE,
#                    sep = ',',
#                    text = 
#                      'DoseA,DoseB,y
#                    2.000,2.000,0.949
#                    0.667,2.000,0.506
#                    0.025,0.000,0.012
#                    0.008,0.000,0.053')

# rep2 <- read.table(header = TRUE,
#                    sep = ',',
#                    text = 
#                      'DoseA,DoseB,y
#                    2.000,2.000,0.780
#                    0.667,2.000,0.458
#                    0.025,0.000,-0.020
#                    0.008,0.000,-0.039')

# data1 <- rbind(rep1, rep2)
# attach(data1)

## FROM THE PAPER
# data1 <- read.table(header = TRUE,
#                    sep = ',',
#                    text = 
#                      'DoseA,DoseB,y
#                    2,     2,    .949
#                    .667,  2,    .506
#                    .222,  2,    .499
#                    .074,  2,    .5
#                    .025,  2,    .483
#                    .008,  2,    .505
#                    2,     .667, .659
#                    2,     .222, .510
#                    2,     .074, .435
#                    2,     .025, .434
#                    2,     .008, .401')
# attach(data1)

## OUR DATA
data1 <- readxl::read_xlsx('/home/job/WORKSPACE/B21052_HELE_01/drug_synergy/20210628_WST-1_04_formatTed_right.xlsx')
attach(data1)
## Core function to estimate Interaction Index, modified based on Harbron’s algorithm (2010)
## --
## d1: observed dose of drug A; d2: observed dose of drug B
## m1: Hill slope estimated from stage I for single drug A
## m2: Hill slope estimated from stage I for single drug B
## logIC50.1: loge(IC50) estimated from stage I for single drug A
## logIC50.2: loge(IC50) estimated from stage I for single drug B
## U: the upper asymptote of single drugs with the default value of 1
## L: the lower asymptote of single drugs with the default value of 0
## tau: interaction index
## gp: grouping parameter for taus, default is the overall tau for the surface model
## niterations: the iteration number on each call of tau.model, default 50
## upper & lower: vector of range to estimate response of combined drug at each combination dose level
## Y: the estimated response based on current tau value & single drug pamameters from stage I.

tau.model <- function(d1, d2, m1, m2, logIC50.1, logIC50.2,
                      U=1, L=0, ..., gp=rep(1,sum(d1>0 & d2>0)),
                      niterations=50) {
  Ic50.1 <- exp(logIC50.1) 
  Ic50.2 <- exp(logIC50.2)
  tau <- exp(unlist(list(...)))
  mono <- (d1==0) | (d2==0)
  if(length(gp) != sum(!mono)) {
    stop('Incorrect Number Of Grouping Parameters')
  }
  if(!all(sort(unique(gp))==1:length(tau))) {
    stop('Numbers Of Parameters & Grouping Parameters Disagree')
  }
  upper <- rep(U, length(d1))
  lower <- G <- rep(L, length(d1))
  for(i in 1:niterations) {
    Y <- (upper + lower)/ 2.0
    Part1 <- Ic50.1 * ((Y-L)/(U-Y))^(1.0/m1)
    Part2 <- Ic50.2 * ((Y-L)/(U-Y))^(1.0/m2)
    dPart1 <- d1/Part1
    dPart2 <- d2/Part2
    dPart3 <- dPart1 * dPart2
    dPart1[!mono] <- dPart1[!mono]/ tau[gp]
    dPart2[!mono] <- dPart2[!mono]/ tau[gp]
    dPart3[!mono] <- dPart3[!mono]/ tau[gp]
    G <- dPart1 + dPart2 + dPart3 - 1
    if(m1<0) {
      lower[G==0] <- upper[G==0] <- Y[G==0]
      lower[G<0] <- Y[G<0]
      upper[G>0] <- Y[G>0]
    } else {
      lower[G==0] <- upper[G==0] <- Y[G==0]
      upper[G<0] <- Y[G<0]
      lower[G>0] <- Y[G>0]
    }
  }
  Y
}


## FIRST STAGE : estimate the Hill Slope and IC50 for drugs A and B, respectively
## --
## uppers, lowers: after rescaling, the upper & lower asymptote of both single drugs
## Ic50.1 & Ic50.2: selected initial value of IC50 for drug A & B
## drugA_est & drugB_est: store the nlsLM estimated Hill slope and IC50 for each single drug

# require(minpack.lm)
uppers <- 1
lowers <- 0
Ic50.1 <- DoseA[DoseB==0][which.min(abs(y[DoseB == 0] - .5))] #pick the initial value for IC50, which is the observed single drug dosage produces the response closest to 0.5
drugA_est <- minpack.lm::nlsLM(y ~ (uppers-lowers) / (1 + exp(m * logIC50) / DoseA^m) + lowers,
                               start = list(m = 1, logIC50 = log(Ic50.1)),
                               subset = DoseB == 0 & DoseA != 0)
Ic50.2 <- DoseB[DoseA==0][which.min(abs(y[DoseA == 0] - .5))]
drugB_est <- minpack.lm::nlsLM(y ~ (uppers-lowers) / (1 + exp(m * logIC50) / DoseB^m) + lowers,
                               start = list(m = 1, logIC50 = log(Ic50.2)),
                               subset = DoseA == 0 & DoseB != 0)


## Bootstrapping the parameters estimated from stage I to estimate the 95% CI of tau estimated from stage II
## --
## drugA.para & drug B.para: store each drug’s estimated Hill slope & IC50
## drugA.var & drugB.var: store each drug’s estimated variance & covariance of Hill slope & IC50
## N: Bootstrapping size
## mu1p & mu2p: random samples of Hill slope & IC50 drawn from bivariate-normal distribution for each drug with corresponding means & variance-covariance estimated from stage I
## pars: store the estimated taus based on each mu1p & mu2p.

drugA.para <- summary(drugA_est)$coef[,1] 
drugA.var <- vcov(drugA_est) 
drugB.para <- summary(drugB_est)$coef[,1]
drugB.var <- vcov(drugB_est)
# require(mvtnorm)
N <- 1000
mu1p <- mu2p <- NULL
mu1p <- mvtnorm::rmvnorm(N, drugA.para, drugA.var) 
mu2p <- mvtnorm::rmvnorm(N, drugB.para, drugB.var) 


## SECOND STAGE : estimate Interaction Index

### To estimate an overall interaction index

pars <- NULL
for (i in 1:N) {
  print(i)
  overallTau <- minpack.lm::nlsLM(y ~ tau.model(d1 = DoseA, d2 = DoseB,
                                                m1 = mu1p[i,1], m2=mu2p[i,1],
                                                logIC50.1 = mu1p[i,2], logIC50.2 = mu2p[i,2],
                                                logtau1 = logtau1),
                                  start = list(logtau1 = 0)) 
  pars <- rbind(pars, summary(overallTau)$coef[,1])
} 
log_tau <- mean(pars)
tau <- unlist(data.frame(est=exp(log_tau), LCL= exp(log_tau-1.96*sd(pars)), HCL= exp(log_tau+1.96*sd(pars))))
print(tau)

### To estimate separated interaction indexes at each dose of drug A combining with various doses of drug B

pars <- NULL
gp <- as.numeric(as.factor(DoseA[DoseA*DoseB > 0]))
for (i in 1:N) {
  sepTauA <- minpack.lm::nlsLM(y ~ tau.model(d1=DoseA, d2=DoseB, m1=mu1p[i,1], m2=mu2p[i,1], 
                                 logIC50.1=mu1p[i,2], logIC50.2=mu2p[i,2], 
                                 logtau1=logtau1, logtau2=logtau2, logtau3=logtau3, logtau4=logtau4, logtau5=logtau5, logtau6=logtau6, gp=gp), start=c(logtau1=0, logtau2=0, logtau3=0, logtau4=0, logtau5=0, logtau6=0)) 
  pars <- rbind(pars, summary(sepTauA)$coef[,1])
}
Ataus <- data.frame(cbind(est=exp(apply(pars,2,mean)),LCL=exp(apply(pars,2,mean)-1.96*apply(pars,2,sd)), HCL=exp(apply(pars,2,mean)+1.96*apply(pars,2,sd))))
cat('The follwing output has been exponentiated')
print(Ataus)
