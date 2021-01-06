##############################################################################################
##############################################################################################
# Check how robust/classical SPCA behave when the mean of the centres is shifted.

rm(list = ls())
source("./utils.R")
load("./scripts/robust_sim_study_crazy_means_eps_1_2_5.RData")

Eps <- length(values[[1]])

robust.meanMRE <- matrix(0, nrow = length(values), ncol = Eps)
classical.meanMRE <- matrix(0, nrow = length(values), ncol = Eps)
MREdiffs <- matrix(0, nrow = length(values), ncol = Eps)
rownames(robust.meanMRE) <- rownames(classical.meanMRE) <- names(values)
for (model in 1:length(values)) {
    for (e in 1:Eps) {
        robust.meanMRE[model, e] <- mean(values[[model]][[e]]$robust.MRE)
        classical.meanMRE[model, e] <- mean(values[[model]][[e]]$classical.MRE)
        MREdiffs[model, e] <- robust.meanMRE[model, e] - classical.meanMRE[model, e]
    }
}

# print(robust.meanMRE)
writeLines(to.tex.numeric.table(round(robust.meanMRE, 3)))
# print(classical.meanMRE)
writeLines(to.tex.numeric.table(round(classical.meanMRE, 3)))
print(mean(robust.meanMRE[MREdiffs>0]/classical.meanMRE[MREdiffs>0]))
print(mean(classical.meanMRE[MREdiffs<0]/robust.meanMRE[MREdiffs<0]))
print(MREdiffs > 0)

##############################################################################################

robust.meanACV <- matrix(0, nrow = length(values), ncol = Eps)
classical.meanACV <- matrix(0, nrow = length(values), ncol = Eps)
ACVdiffs <- matrix(0, nrow = length(values), ncol = Eps)
rownames(robust.meanACV) <- rownames(classical.meanACV) <- names(values)
for (model in 1:length(values)) {
    for (e in 1:Eps) {
        robust.meanACV[model, e] <- 1 - mean(values[[model]][[e]]$robust.ACV)
        classical.meanACV[model, e] <- 1 - mean(values[[model]][[e]]$classical.ACV)
        ACVdiffs[model, e] <- robust.meanACV[model, e] - classical.meanACV[model, e]
    }
}

# print(robust.meanACV)
writeLines(to.tex.numeric.table(round(1000*robust.meanACV, 3)))
# print(classical.meanACV)
writeLines(to.tex.numeric.table(round(1000*classical.meanACV, 3)))
print(mean(robust.meanACV[ACVdiffs>0]/classical.meanACV[ACVdiffs>0]))
print(mean(classical.meanACV[ACVdiffs<0]/robust.meanACV[ACVdiffs<0]))
print(ACVdiffs > 0)

##############################################################################################
##############################################################################################