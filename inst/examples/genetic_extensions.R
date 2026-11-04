# Small diploid example: supplied effects, no downloads or optional packages.
library(gsim)
W <- matrix(c(0,1,2,1, 2,0,1,2, 0,2,1,0), 4, 3,
            dimnames=list(paste0("animal",1:4), paste0("marker",1:3)))
additive <- matrix(c(1,-2,.5),3,1,dimnames=list(colnames(W),"D1"))
dominance <- additive * .3
interaction <- additive * -.2
environment <- setNames(c(-1,0,2,1),rownames(W))
sim <- gsim(W=W, architecture="fixed", beta=additive,
            dominance=dominance, interaction=interaction,
            environment=environment, standardize_W=FALSE,
            scale_effects=FALSE, seed=7)
stopifnot(isTRUE(all.equal(sim$G,sim$G_additive+sim$G_dominance+sim$G_interaction)))
# Freeze this model for future BED cohorts without recentering or rescaling.
state <- gsim_trait_state(sim,center=0)
# For repeated visits, pass baseline values and raw-dosage slope values separately.
records <- gsim_records(sim$G_additive+sim$G_dominance,
                        interaction_values=sim$G_slope,
                        environment=1, seed=11)
# Phased-parent workflows accept mutation_rate, for example:
# cohort <- gsim_simulate_cohort(parents,offspring,seed=11,output="cohort",
#                               mutation_rate=1e-5)
