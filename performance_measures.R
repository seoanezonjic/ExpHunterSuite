get_TP <- function(cont_table) {
	cont_table["TRUE","TRUE"]
}

get_FP  <- function(cont_table) {
	cont_table["TRUE","FALSE"]
}

get_TN <- function(cont_table) {
	cont_table["FALSE","FALSE"]
}

get_FN <- function(cont_table) {
	cont_table["FALSE","TRUE"]
}

# JRP fr future: More than happy to change the input to the contingence table directly. Or we can use method dispatch to call different methods depending on the class of the input.
# For now it's like this to avoid people making the table the wrong way round - they have to give a test_res and ground truth results as it stands.
calc_prec <- function(test_res_logical, ground_truth_logical) {
	conf_mat <- table(test_res_logical, ground_truth_logical)
	TP <- get_TP(conf_mat)
	FP <- get_FP(conf_mat)
	TP / (TP + FP)
}

calc_rec <- function(test_res_logical, ground_truth_logical) {
	conf_mat <- table(test_res_logical, ground_truth_logical)
	TP <- get_TP(conf_mat)
	FN <- get_FN(conf_mat)
	TP / (TP + FN)
}

calc_spec <- function(test_res_logical, ground_truth_logical) {
	conf_mat <- table(test_res_logical, ground_truth_logical)
	TN <- get_TN(conf_mat)
	FP <- get_FP(conf_mat)
	TN / (TN + FP)
}

calc_acc<- function(test_res_logical, ground_truth_logical) {
	conf_mat <- table(test_res_logical, ground_truth_logical)
	TP <- get_TP(conf_mat)
	TN <- get_TN(conf_mat)
	(TP + TN) / sum(conf_mat)
}

calc_f1 <- function(test_res_logical, ground_truth_logical) {
	precision_val <- calc_prec(test_res_logical, ground_truth_logical)
	recall_val <- calc_rec(test_res_logical, ground_truth_logical)
	2 * ((precision_val * recall_val) / (precision_val + recall_val))
}
