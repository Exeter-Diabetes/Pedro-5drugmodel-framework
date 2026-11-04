# Author: pcardoso
###############################################################################

# Prediction code

###############################################################################

# Predict from model ----
predict_5drugmodel <- function(data, model, drug_var = NULL, drugs = NULL, pred_col = "pred.") {
	
	# check whether objects is provided
	if (is.null(drugs)) {stop("drugs needs to be provided")}
	if (is.null(drug_var)) {stop("drug_var needs to be provided")}
	
	# check class of object (method used if class == "imputation method")
	if ("imputation method" %in% class(data)) {
		
		# check the imputation method
		if (data$type == "group") {
			# if group used
			
			# output object
			output_dataset <- data$imputation
			
			# check whether drugs provided are included in the dataset
			if (!all(drugs %in% unique(output_dataset %>% select(all_of(drug_var)) %>% unlist() %>% unique()))) {stop("some drugs provided are not in drug_var")}
			if (is.null(drugs)) {stop("drugs needs to be provided")}
			
			# name of variables imputated
			cols_imputed <- output_dataset %>% select(contains("_impute")) %>% colnames()
			orig_cols_imputed <- sub("_group_impute$", "", cols_imputed)
			
			# remove the original vars with missingness
			imputed_dataset <- output_dataset %>%
					select(-all_of(orig_cols_imputed))
			
			# rename imputed vars
			colnames(imputed_dataset) <- sub("_group_impute$", "", colnames(imputed_dataset))
			
			# predict from model
			for (drug_name in drugs) {
				
				# change drug_var to the specific drug being predicted
				current_prediction_dataset <- imputed_dataset %>%
						rename("target_var" = drug_var) %>%
						mutate(target_var = drug_name) %>%
						rename({{drug_var}} := "target_var")
				
				# predict from model
				prediction_vector <- predict(model, current_prediction_dataset)
				
				# attach this var to original dataset
				output_dataset <- cbind(output_dataset, prediction_vector) %>%
						as.data.frame() %>%
						rename(!!paste0("pred.group.", drug_name) := "prediction_vector")
				
			}
			
			
		} else if (data$type == "mice") {
			# if mice used
			
			# output object
			output_dataset <- data$imputation
			
			# check whether drugs provided are included in the dataset
			# if (!all(drugs %in% unique(output_dataset %>% select(all_of(drug_var)) %>% unlist() %>% unique()))) {stop("some drugs provided are not in drug_var")}
			if (is.null(drugs)) {stop("drugs needs to be provided")}
			
			# name of variables imputated
			cols_imputed <- output_dataset %>% select(contains("_mice_impute")) %>% colnames()
			orig_cols_imputed <- sub("_mice_impute$", "", cols_imputed)
			
			# remove the original vars with missingness
			imputed_dataset <- output_dataset %>%
					select(-all_of(orig_cols_imputed))
			
			# rename imputed vars
			colnames(imputed_dataset) <- sub("_mice_impute$", "", colnames(imputed_dataset))
			
			# predict from model
			for (drug_name in drugs) {
				
				# change drug_var to the specific drug being predicted
				current_prediction_dataset <- imputed_dataset %>%
						rename("target_var" = drug_var) %>%
						mutate(target_var = drug_name) %>%
						rename({{drug_var}} := "target_var")
				
				# predict from model
				prediction_vector <- predict(model, current_prediction_dataset)
				
				# attach this var to original dataset
				output_dataset <- cbind(output_dataset, prediction_vector) %>%
						as.data.frame() %>%
						rename(!!paste0("pred.mice.", drug_name) := "prediction_vector")
				
			}
			
		} else if (data$type == "dpmm") {
			# if DPMM used
			stop("DPMM not codded yet")
		} else {stop("This is an error, you shouldn't be here")}
		
	} else {
		# no imputation method / prediction from model
		
		# output object
		output_dataset <- data
		
		# predict from model
		for (drug_name in drugs) {
			
			# change drug_var to the specific drug being predicted
			current_prediction_dataset <- output_dataset %>%
					rename("target_var" = drug_var) %>%
					mutate(target_var = drug_name) %>%
					rename({{drug_var}} := "target_var")
			
			# predict from model
			prediction_vector <- predict(model, current_prediction_dataset)
			
			# attach this var to original dataset
			output_dataset <- cbind(output_dataset, prediction_vector) %>%
					as.data.frame() %>%
					rename(!!paste0(pred_col, drug_name) := "prediction_vector")
			
		}
		
	}
	
	return(output_dataset)
	
}
