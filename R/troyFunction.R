# Copyright 2023 Observational Health Data Sciences and Informatics
#
# This file is part of TROY
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
library(dplyr)
library(FeatureExtraction)

troyFunction <- function(cohortDatabaseSchema=cohortDatabaseSchema,
                         cohortTable=cohortTable,
                         connectionDetails=connectionDetails,
                         outputFolder=outputFolder) {
  # ParallelLogger::addDefaultFileLogger(file.path(outputFolder, "replicateLog.txt"))
  # ParallelLogger::addDefaultErrorReportLogger(file.path(outputFolder, "errorReportR.txt"))
  # on.exit(ParallelLogger::unregisterLogger("DEFAULT_FILE_LOGGER", silent = TRUE))
  # on.exit(ParallelLogger::unregisterLogger("DEFAULT_ERRORREPORT_LOGGER", silent = TRUE), add = TRUE)

  message("Running TROY...")
  fileName <- file.path('~/git/Troy/inst/csv/', "PrespecAnalyses.csv")
  prespecAnalyses <- readr::read_csv(fileName)
  referenceRct <- read.csv(file.path('~/git/Troy/inst/settings/', "TargetComparatorList.csv"))
  colnames(referenceRct) = c("targetId", "comparatorId", "design", "description", "trialName")
  trialNames <- unique(referenceRct$trialName)

  # for loop
  for (i in 1:length(trialNames)) {
    trial <- unique(trialNames)[i]
    trialIdx <- which(referenceRct[,"trialName"] == trial)
    message(trial)

    # specificationRct <- read.csv(file.path(getwd(), "inst", "csv", paste0(trial, ".csv")))
    specificationRct <- read.csv(file.path('~/git/Troy/inst/csv/', paste0(trial, ".csv")))
    specificationRct$n <- NA
    colnames(specificationRct) <- c("characteristics", "target", "targetSd", "comparator", "comparatorSd", "tag", "isNa", "targetSize", "comparatorSize", "statistics", "analysisId", "conceptIds", "hasValue", "summary", "n")
    specificationRct <- subset(specificationRct, specificationRct$isNa != 'Y')
    specificationRct[is.na(specificationRct)]<-""

    for (z in 1:2){
      connection <- DatabaseConnector::connect(connectionDetails)
      cohort <- Andromeda::andromeda()
      sql <- "select * from @cohortDatabaseSchema.@cohortTable"
      sql <- SqlRender::render(sql, cohortDatabaseSchema=cohortDatabaseSchema, cohortTable=cohortTable)
      DatabaseConnector::querySqlToAndromeda(connection, sql, andromeda = cohort, andromedaTableName = "cohort")

      if(z == 1){
        message("Eligible criteria cohorts")
      } else {
        message("Indication cohorts")
      }

      tc <- referenceRct[trialIdx[z],c(1,2)]
      targetNumber <- tc[1,1]
      comparatorNumber <- tc[1,2]
      targetCohort <- cohort$cohort %>% filter(.data$COHORT_DEFINITION_ID == targetNumber) %>% mutate(treatment = 1)
      comparatorCohort <- cohort$cohort %>% filter(.data$COHORT_DEFINITION_ID == comparatorNumber) %>% mutate(treatment = 0)
      population <- targetCohort %>% union(comparatorCohort)
      keepFirst  <- population %>%
        group_by(subject_id) %>%
        filter(cohort_start_date == min(cohort_start_date))

      sizeOE <- data.frame(treatment = -1, n = keepFirst %>% tally() %>% summarise(n = n()))
      sizeOE <- rbind(sizeOE, data.frame(population %>% group_by(.data$treatment) %>% tally()))

      if (nrow(sizeOE) < 2) {
        message(paste0("There are no eligible patients for the following trials: ", trial))
        next
      }

      covariateSettings <- troyCreateCovariateSettings(specificationRct, prespecAnalyses=prespecAnalyses)
      covariates <- getDbCovariateData(connectionDetails = connectionDetails,
                                       cdmDatabaseSchema = cdmDatabaseSchema,
                                       cohortDatabaseSchema = cohortDatabaseSchema,
                                       cohortTable = cohortTable,
                                       cohortId = c(targetNumber, comparatorNumber),
                                       covariateSettings = covariateSettings)

      # covariates$covariates <- covariates$covariates %>% left_join(population, by=("rowId"), copy = TRUE)
      covariates$covariates <- covariates$covariates %>% left_join(covariates$covariateRef)

      statisticsPooled <- covariates$covariates %>%
        group_by(.data$covariateId) %>%
        summarise(conceptId = .data$conceptId,
                  sum = sum(as.numeric(.data$covariateValue), na.rm = TRUE),
                  mean = mean(as.numeric(.data$covariateValue), na.rm = TRUE),
                  sumSqr = sum(as.numeric(.data$covariateValue)^2, na.rm = TRUE),
                  median = median(as.numeric(.data$covariateValue), na.rm = TRUE),
                  n = n(),
                  min = min(as.numeric(.data$covariateValue), na.rm = TRUE),
                  max = max(as.numeric(.data$covariateValue), na.rm = TRUE)
        ) %>%
        mutate(sd = sqrt((.data$sumSqr - (.data$sum^2 / .data$n)) / .data$n))

      df <- data.frame()
      df2 <- data.frame()
      for (idx in 1:nrow(specificationRct)) {
        print(idx)
        print(specificationRct$characteristics[idx])
        conceptIds <- as.numeric(strsplit(as.character(specificationRct$conceptIds[idx]), ";")[[1]])
        analysisId <- as.numeric(specificationRct$analysisId[idx])
        covariateIds <- paste0(conceptIds, sprintf("%03d", analysisId))
        # summary = Y
        if (specificationRct$summary[idx] == "Y") {
          if (substr(analysisId, 1, 1) == 7){
            t <- data.frame(statisticsPooled %>% filter(conceptId %in% conceptIds))
            maxCovariateId <- t[which.max(t$n),1]
            t <- data.frame(statisticsPooled %>% filter(covariateId %in% maxCovariateId))
            if (nrow(t)==0) {
              t <- 0
              sd <- 0
              n <- 0
            } else {
              n <- t$n
              n <- as.numeric(n)
              t <- t %>% select(specificationRct$statistics[idx])
              t <- as.numeric(t)
              if(specificationRct$statistics[idx]=="median") {
                sd <- summary(
                  data.frame(
                    covariates$covariates %>%
                      filter(.data$covariateId %in% maxCovariateId))$covariateValue
                )
                sd <- as.numeric(sd)
                sd <- (paste0(sd[2], "–", sd[5]))} else {sd <- data.frame(statisticsPooled %>% filter(covariateId %in% maxCovariateId))[, c("sd")]}}
          } else {
            t <- data.frame(statisticsPooled %>% filter(covariateId %in% covariateIds))[,c(specificationRct$statistics[idx])]
            n <- data.frame(statisticsPooled %>% filter(covariateId %in% covariateIds))[,c("n")]
            if(specificationRct$statistics[idx]=="median"){
              sd <- summary(
                data.frame(
                  covariates$covariates %>%
                    filter(.data$covariateId %in% covariateIds))$covariateValue
              )
              sd <- as.numeric(sd)
              sd <- (paste0(sd[2], "–", sd[5]))} else {sd <- data.frame(statisticsPooled %>% filter(covariateId %in% covariateIds))[,"sd"]}
          }
        }

        #t <- data.frame(statisticsPooled %>% filter(conceptId %in% conceptIds))[,c(specificationRct$statistics[idx])]
        #t <- as.numeric(t)

        # summary NA, hasValue NA, eg. count
        if (!specificationRct$summary[idx]=="Y"){
          if (specificationRct$hasValue[idx]==""){
            t <- data.frame(covariates$covariates %>%
                              filter(.data$covariateId %in% covariateIds) %>%
                              summarise(n = n_distinct(.data$rowId)))
            t <- as.numeric(t)
            n <- as.numeric(t)
            sd <- (t / sizeOE[sizeOE[,"treatment"] < 0,"n"]) * 100
          } else {
            #if (nchar(covariateIds) > 3) {
            if (substr(analysisId, 1, 1) == 7) {
              t <- data.frame(covariates$covariates %>%
                                filter(.data$conceptId %in% conceptIds) %>%
                                filter(rlang::parse_expr(paste0(specificationRct[idx,"hasValue"]))) %>%
                                summarise(n = n_distinct(.data$rowId)))
              t <- as.numeric(t)
              sd <- (t / sizeOE[sizeOE[,"treatment"] < 0,"n"]) * 100
            } else {
              t <- data.frame(covariates$covariates %>%
                                filter(.data$covariateId %in% covariateIds) %>%
                                filter(rlang::parse_expr(paste0(specificationRct[idx,"hasValue"]))) %>%
                                summarise(n = n_distinct(.data$rowId)))
              t <- as.numeric(t)
              sd <- (t / sizeOE[sizeOE[,"treatment"] < 0,"n"]) * 100
            }
          }
        }

        # Only IQR
        if (specificationRct$statistics[idx]=="IQR") {
          #covariateIds <- c(1002)
          t <- summary(
            data.frame(
              covariates$covariates %>%
                filter(.data$covariateId %in% covariateIds))$covariateValue
          )
          t <- as.numeric(t)
          sd <- NA
          if(is.na(t[1] > 0)){
            t <- summary(
              data.frame(
                covariates$covariates %>%
                  filter(.data$conceptId %in% conceptIds))$covariateValue
            )
            t <- as.numeric(t)
          }
          #paste0(t[2], "–", t[5])
          t <- (paste0(t[2], "–", t[5]))
          sd <- NA
        }
        print(t)
        # df <- rbind(df, t)
        # df2 <- rbind(df2, sd)
        specificationRct$pooledTroy[idx] <- t
        specificationRct$pooledTroySd[idx] <- sd
        specificationRct$n[idx] <- n
        t <- NA
        sd <- NA
        n <- NA
      }
      # specificationRct$pooledTroy <- df[,1]
      # specificationRct$pooledTroySd <- df2[,1]

      if(z == 1){
        outputCsv <- file.path(outputFolder, paste0(trial, "TroyCriteria.csv"))
      } else {
        outputCsv <- file.path(outputFolder, paste0(trial, "TroyIndication.csv"))
      }


      output <- read.csv(file.path('~/git/Troy/inst/csv/', paste0(trial, ".csv")))
      output <- output %>% left_join(select(specificationRct, characteristics, pooledTroy, pooledTroySd, n), by = c('characteristics'='characteristics'))

      if(z == 1){
        outputCsvSimple <- file.path(outputFolder, paste0(trial, "TroyCriteriaSimple.csv"))
      } else {
        outputCsvSimple <- file.path(outputFolder, paste0(trial, "TroyIndicationSimple.csv"))
      }


      for (i in 1:nrow(output)){
        if(output$statistics[i] == "n"){
          output$targetSd[i] <- paste0(round((as.numeric(output$target[i]) / as.numeric(output$targetSize)) * 100, 2), '%')
          output$comparatorSd[i] <- paste0(round((as.numeric(output$comparator[i]) / as.numeric(output$targetSize)) * 100, 2), '%')
          output$pooledTroySd[i] <- paste0(round(as.numeric(output$pooledTroySd[i]), 2), '%')
        } else {
          next
        }
      }

      for (i in 1:nrow(output)){
        if(!is.na(output$pooledTroySd[i])){
          if(output$pooledTroySd[i] == "NA%"){
            output$pooledTroySd[i] <- NA
          } else {
            next
          }
        }
      }

      output$pooledTroySize <- NA
      output$targetTroySize <- NA
      output$comparatorTroySize <- NA
      output$pooledTroySize[1] <- sizeOE %>% filter(treatment < 0) %>% summarize(n = if_else(n() > 0, first(n), 0)) %>% pull(n)
      output$targetTroySize[1] <- sizeOE %>% filter(treatment == 1) %>% summarize(n = if_else(n() > 0, first(n), 0)) %>% pull(n)
      output$comparatorTroySize[1] <- sizeOE %>% filter(treatment == 0) %>% summarize(n = if_else(n() > 0, first(n), 0)) %>% pull(n)

      write.csv(output, outputCsv)
      write.csv(output[,c(1,2,3,4,5,16,17,8,9,18,19,20)], outputCsvSimple)
      disconnect(connection)
    }
  }
}
###########
# Troy function
troyCreateCovariateSettings <- function (specificationRct, prespecAnalyses, useDemographicsGender = FALSE, useDemographicsAge = FALSE,
                                         useDemographicsAgeGroup = FALSE, useDemographicsRace = FALSE,
                                         useDemographicsEthnicity = FALSE, useDemographicsIndexYear = FALSE,
                                         useDemographicsIndexMonth = FALSE, useDemographicsPriorObservationTime = FALSE,
                                         useDemographicsPostObservationTime = FALSE, useDemographicsTimeInCohort = FALSE,
                                         useDemographicsIndexYearMonth = FALSE, useConditionOccurrenceAnyTimePrior = FALSE,
                                         useConditionOccurrenceLongTerm = FALSE, useConditionOccurrenceMediumTerm = FALSE,
                                         useConditionOccurrenceShortTerm = FALSE, useConditionOccurrencePrimaryInpatientAnyTimePrior = FALSE,
                                         useConditionOccurrencePrimaryInpatientLongTerm = FALSE,
                                         useConditionOccurrencePrimaryInpatientMediumTerm = FALSE,
                                         useConditionOccurrencePrimaryInpatientShortTerm = FALSE,
                                         useConditionEraAnyTimePrior = FALSE, useConditionEraLongTerm = FALSE,
                                         useConditionEraMediumTerm = FALSE, useConditionEraShortTerm = FALSE,
                                         useConditionEraOverlapping = FALSE, useConditionEraStartLongTerm = FALSE,
                                         useConditionEraStartMediumTerm = FALSE, useConditionEraStartShortTerm = FALSE,
                                         useConditionGroupEraAnyTimePrior = FALSE, useConditionGroupEraLongTerm = FALSE,
                                         useConditionGroupEraMediumTerm = FALSE, useConditionGroupEraShortTerm = FALSE,
                                         useConditionGroupEraOverlapping = FALSE, useConditionGroupEraStartLongTerm = FALSE,
                                         useConditionGroupEraStartMediumTerm = FALSE, useConditionGroupEraStartShortTerm = FALSE,
                                         useDrugExposureAnyTimePrior = FALSE, useDrugExposureLongTerm = FALSE,
                                         useDrugExposureMediumTerm = FALSE, useDrugExposureShortTerm = FALSE,
                                         useDrugEraAnyTimePrior = FALSE, useDrugEraLongTerm = FALSE,
                                         useDrugEraMediumTerm = FALSE, useDrugEraShortTerm = FALSE,
                                         useDrugEraOverlapping = FALSE, useDrugEraStartLongTerm = FALSE,
                                         useDrugEraStartMediumTerm = FALSE, useDrugEraStartShortTerm = FALSE,
                                         useDrugGroupEraAnyTimePrior = FALSE, useDrugGroupEraLongTerm = FALSE,
                                         useDrugGroupEraMediumTerm = FALSE, useDrugGroupEraShortTerm = FALSE,
                                         useDrugGroupEraOverlapping = FALSE, useDrugGroupEraStartLongTerm = FALSE,
                                         useDrugGroupEraStartMediumTerm = FALSE, useDrugGroupEraStartShortTerm = FALSE,
                                         useProcedureOccurrenceAnyTimePrior = FALSE, useProcedureOccurrenceLongTerm = FALSE,
                                         useProcedureOccurrenceMediumTerm = FALSE, useProcedureOccurrenceShortTerm = FALSE,
                                         useDeviceExposureAnyTimePrior = FALSE, useDeviceExposureLongTerm = FALSE,
                                         useDeviceExposureMediumTerm = FALSE, useDeviceExposureShortTerm = FALSE,
                                         useMeasurementAnyTimePrior = FALSE, useMeasurementLongTerm = FALSE,
                                         useMeasurementMediumTerm = FALSE, useMeasurementShortTerm = FALSE,
                                         useMeasurementValueAnyTimePrior = FALSE, useMeasurementValueLongTerm = FALSE,
                                         useMeasurementValueMediumTerm = FALSE, useMeasurementValueShortTerm = FALSE,
                                         useMeasurementRangeGroupAnyTimePrior = FALSE, useMeasurementRangeGroupLongTerm = FALSE,
                                         useMeasurementRangeGroupMediumTerm = FALSE, useMeasurementRangeGroupShortTerm = FALSE,
                                         useObservationAnyTimePrior = FALSE, useObservationLongTerm = FALSE,
                                         useObservationMediumTerm = FALSE, useObservationShortTerm = FALSE,
                                         useCharlsonIndex = FALSE, useDcsi = FALSE, useChads2 = FALSE,
                                         useChads2Vasc = FALSE, useHfrs = FALSE, useDistinctConditionCountLongTerm = FALSE,
                                         useDistinctConditionCountMediumTerm = FALSE, useDistinctConditionCountShortTerm = FALSE,
                                         useDistinctIngredientCountLongTerm = FALSE, useDistinctIngredientCountMediumTerm = FALSE,
                                         useDistinctIngredientCountShortTerm = FALSE, useDistinctProcedureCountLongTerm = FALSE,
                                         useDistinctProcedureCountMediumTerm = FALSE, useDistinctProcedureCountShortTerm = FALSE,
                                         useDistinctMeasurementCountLongTerm = FALSE, useDistinctMeasurementCountMediumTerm = FALSE,
                                         useDistinctMeasurementCountShortTerm = FALSE, useDistinctObservationCountLongTerm = FALSE,
                                         useDistinctObservationCountMediumTerm = FALSE, useDistinctObservationCountShortTerm = FALSE,
                                         useVisitCountLongTerm = FALSE, useVisitCountMediumTerm = FALSE,
                                         useVisitCountShortTerm = FALSE, useVisitConceptCountLongTerm = FALSE,
                                         useVisitConceptCountMediumTerm = FALSE, useVisitConceptCountShortTerm = FALSE,
                                         longTermStartDays = -365, mediumTermStartDays = -180, shortTermStartDays = -30,
                                         endDays = 0, includedCovariateConceptIds = c(), addDescendantsToInclude = FALSE,
                                         excludedCovariateConceptIds = c(), addDescendantsToExclude = FALSE,
                                         includedCovariateIds = c())
{
  covariateSettings <- list(temporal = FALSE, temporalSequence = FALSE)
  analysisNames <- paste0("use", as.vector(unlist(prespecAnalyses[prespecAnalyses$analysisId %in% unique(specificationRct$analysisId),"analysisName"])))
  formalNames <- names(formals(createCovariateSettings))
  anyUseTrue <- FALSE
  for (name in formalNames) {
    #value <- formalNames[name] %in% analysisNames
    value <- name %in% analysisNames
    if (is.null(value)) {
      value <- vector()
    }
    if (is.numeric(get(name))) {
      value <- get(name)
    }
    if (grepl("use.*", name)) {
      if (value) {
        covariateSettings[[sub("use", "", name)]] <- value
        anyUseTrue <- TRUE
      }
    }
    else {
      covariateSettings[[name]] <- value
    }
  }
  if (!anyUseTrue) {
    stop("No covariate analysis selected. Must select at least one")
  }
  attr(covariateSettings, "fun") <- "getDbDefaultCovariateData"
  class(covariateSettings) <- "covariateSettings"
  return(covariateSettings)
}

