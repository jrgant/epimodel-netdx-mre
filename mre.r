library(EpiModel)
library(data.table)

set.seed(1971)

## set up the network
N <- 10000
nw <- network.initialize(N, directed = FALSE)
nw <- set.vertex.attribute(nw, "group", sample(1:3, N, replace = TRUE))

## set up the target stats as a named vector
target.stats <- c(2000, 100, 200, 300, 400, 500)
names(target.stats) <- c("edges",
                           "mix.group.1.2", "mix.group.1.3",
                           "mix.group.2.2", "mix.group.2.3",
                           "mix.group.3.3")

## shuffle the target stats, retaining the names
target.stats.shuffle <- target.stats[sample(1:length(target.stats),
                                            size = length(target.stats),
                                            replace = FALSE)]

## set up the formation and dissolution formulas
formation <- ~ edges + nodemix("group")
dissolution <- dissolution_coefs(~ offset(edges), duration = 50, d.rate = 0)

## estimate netork models (set seed to account for stochastic model fitting)
set.seed(888888)
est <- netest(nw, formation, target.stats, coef.diss = dissolution)

set.seed(888888)
est.shuffle <- netest(nw, formation, target.stats.shuffle, coef.diss = dissolution)

## demonstrate that statnet estimates the same parameters for each
## target statistics, regardless of the order of the target stats
all(est$coef.form - est.shuffle$coef.form == 0)

## model diagnostics
dx <- netdx(est, nsims = 10, nsteps = 50)

## show diagnostic outputs
dx$target.stats
dx$stats.table.formation

## ERROR
## rownames and target.stat names are not aligned
##   - plot.netdx() and percent diffs. in dx$stats.table.formation both affected
cbind(rownames(dx$target.stats), dx$target.stats$names)
plot(dx)

## align the simulated output with the original target.stats and
## recalculate the percent differences
tstats_match <- merge(
  as.data.table(target.stats,
                keep.rownames = TRUE)[, .(rn, target = target.stats)],
  as.data.table(dx$stats.table.formation,
                keep.rownames = TRUE)[, .(rn, `Sim Mean`)],
  by = "rn"
)

tstats_match[, corrected_pct_diff := (`Sim Mean` - target) / target * 100]

tstats_match
