#!/usr/bin/env Rscript

library(openxlsx)

setwd("H:/metabolomics_for_nathan/PedsObesity_R24_Serum_Targeted_Metabolites/POMMS_Interim_Cohort/")

load("POMMS_Interim_Cohort_Data.Rdata")

data        = pomms_interim_cohort_data
metabolites = colnames(data)[20:93]
formula     = "CaseControl ~ `%s` + Age + Sex + Race + Ethnicity"

results     = data.frame(Metabolite=NA, Number=1, Beta=1, SE=1, Tval=1, Pval=1)
for(metabo in metabolites) {
	data[,metabo] = log(data[,metabo])
	data[,metabo] = ifelse(abs(data[,metabo]) == Inf, NA, data[,metabo])
	f             = sprintf(formula, metabo)
	model         = tryCatch(glm(f, data=data, family="binomial"), error=function(e) { NA })
	n             = tryCatch(nobs(model),                          error=function(e) { NA })
	s             = tryCatch(summary(model)$coef[2,],              error=function(e) { rep(NA, 4) })
	results       = rbind(results, c(metabo, n, s))
}
results          = results[-1,]
results$Number   = as.numeric(results$Number)
results$Beta     = as.numeric(results$Beta)
results$SE       = as.numeric(results$SE)
results$Tval     = as.numeric(results$Tval)
results$Pval     = as.numeric(results$Pval)
results$adj.Pval = p.adjust(results$Pval, method="fdr")
results          = results[order(results$adj.Pval),]

pomms_interim_cohort_analysis_results = results
save(      pomms_interim_cohort_analysis_results, file="POMMS_Interim_Cohort_Analysis_Results.Rdata")
write.xlsx(pomms_interim_cohort_analysis_results,      "POMMS_Interim_Cohort_Analysis_Results.xlsx")
