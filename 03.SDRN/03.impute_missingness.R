# Author: pcardoso
###############################################################################

# Imputation code

###############################################################################


# Function for different imputation methods ----
imputation_methods <- function(data, method, mice.ignore.vars = NULL, mice.m = 5, mice.maxit = 10) {
	
	# load libraries
	require(tidyverse)
	
	# check whether method is one of the available options
	if (!(tolower(method) %in% c("group", "mice", "dpmm"))) {stop("The only methods available are: group / mice / DPMM")}
	
	
	# Imputation methods
	
	if (tolower(method) == "group") {
		
		## Group patients by selected variables
		### sex, ethnicity, hba1c by 10 mmol/mol
		## Impute with mean value
		
		imputation <- data %>%
				mutate(
						hba1c_groups = cut(prehba1c, breaks = seq(50, 110, 10), include.lowest = TRUE)
				) %>%
				group_by(sex, ethnicity, hba1c_groups) %>% # reomved age (?)
				mutate(
						prebmi_impute = mean(prebmi, na.rm = TRUE),
						preegfr_impute = mean(preegfr, na.rm = TRUE),
						pretotalcholesterol_impute = mean(pretotalcholesterol, na.rm = TRUE),
						prehdl_impute = mean(prehdl, na.rm = TRUE),
						prealt_impute = mean(prealt, na.rm = TRUE)
				) %>%
				ungroup() %>%
				mutate(
						prebmi_group_impute = ifelse(is.na(prebmi), prebmi_impute, prebmi),
						preegfr_group_impute = ifelse(is.na(preegfr), preegfr_impute, preegfr),
						pretotalcholesterol_group_impute = ifelse(is.na(pretotalcholesterol), pretotalcholesterol_impute, pretotalcholesterol),
						prehdl_group_impute = ifelse(is.na(prehdl), prehdl_impute, prehdl),
						prealt_group_impute = ifelse(is.na(prealt), prealt_impute, prealt)
				) %>%
				select(-c(prebmi_impute, preegfr_impute, pretotalcholesterol_impute, prehdl_impute, prealt_impute)) %>%
				as.data.frame()
		
		# output object
		output <- list(
				type = "group",
				imputation = imputation
		)
		
	} else if (tolower(method) == "mice") {
		
		# MICE
		
		# load libraries
		require(mice)
		
		if (!is.null(mice.ignore.vars)) {
			
			prediction_matrix <- mice::mice(
					data = data,
					maxit = 0
			)$predictorMatrix
			
			# change ignore variables to not be used
			prediction_matrix[mice.ignore.vars,] <- 0
			
		}
		
		# fit mice model and impute
		mice_output <- mice::mice(
				data = data,
				m = mice.m,
				maxit = mice.maxit,
				predictorMatrix = prediction_matrix
		)
		
		# datasets with imputations
		imputed_datasets <- complete(mice_output, action = "all")
		
		# list of imputed variables
		selected_list <- lapply(imputed_datasets, \(df) select(df, all_of(names(Filter(\(x) x != "", mice_output$method)))))
		
		# get average of imputed values
		imputed_variables <- Reduce("+", selected_list) / mice_output$m
		colnames(imputed_variables) <- paste0(colnames(imputed_variables), "_mice_impute")
		
		# final dataset
		imputation <- data %>%
				cbind(imputed_variables) %>%
				as.data.frame ()
		
		# output object
		output <- list(
				type = "mice",
				imputation = imputation
		)
		
	} else if (tolower(method) == "dpmm") {
		
		# DPMM
		## Fit the model (complete data)
		
		## Make the predictions
		
		stop("DPMM not coded yet")
		
	} else {stop("This is an error, the function shouldn't be here.")}
	
	class(output) = "imputation method"
	
	return(output)
	
}








