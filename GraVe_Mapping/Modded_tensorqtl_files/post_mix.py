# ## modified by Yanyu Liang
# def trc_calc(genotypes_t, log_counts_t, raw_counts_t, covariates0_t,
#              count_threshold=100, select_covariates=True):
#     mask_t = raw_counts_t >= count_threshold
#     mask_cov = raw_counts_t != 0
#
#     if select_covariates:
#         covariates_t = select_covariates(covariates0_t[mask_cov], log_counts_t[mask_cov])
#         # b_t, b_se_t = linreg_robust(covariates0_t[mask_cov, :], log_counts_t[mask_cov])
#         # tstat_t = b_t / b_se_t
#         # m = tstat_t.abs() > 2
#         # m[0] = True
#         # covariates_t = covariates0_t[:, m]
#     else:
#         covariates_t = covariates0_t
#
#     M = torch.unsqueeze(mask_t, 1).float()
#     M_cov = torch.unsqueeze(mask_cov, 1).float()
#     Y = log_counts_t.reshape(1,-1).T
#     Y[M_cov == False] = 0
#     if covariates_t.shape[1] != 0:
#         Y[M_cov == True] = regress_out(covariates_t[mask_cov, :], log_counts_t[mask_cov])
#     res = wrapper_nominal_algo1_matrixLS(Y, genotypes_t.T / 2, M, numpy = False)
#     dof = M.sum() - 2
#     return res, int(mask_t.sum()), dof


# def trc_calc2(genotypes_t, log_counts_t, raw_counts_t, covariates0_t,
#              count_threshold=100, select_covariates=True):
#     """
#     Inputs:
#       genotypes_t:  genotype dosages (variants x samples)
#       log_counts_t: log-transformed, normalized counts (e.g., log(counts/size_factors))
#       # log_counts_t: log(counts/(2*libsize)) --> TODO: use better normalization. CPM/TMM vs size factors?
#       # log_counts_t: log(counts/(2*libsize)) --> TODO: use size factor normalized counts instead
#       raw_counts_t: raw RNA-seq counts
#       covariates0_t: covariates matrix (samples x covariates)
#                      including genotype PCs, PEER factors, etc.
#                      ***with intercept in first column***
#       count_threshold: minimum read count to include a sample
#     """
#     # only use samples that pass count threshold
#
# ## modified by Yanyu Liang
#
#     mask_t = raw_counts_t >= count_threshold
#     mask_cov = raw_counts_t != 0
#
#     if select_covariates:
#         b_t, b_se_t = linreg_robust(covariates0_t[mask_cov, :], log_counts_t[mask_cov])
#         tstat_t = b_t / b_se_t
#         m = tstat_t.abs() > 2
#         m[0] = True
#         covariates_t = covariates0_t[:, m]
#     else:
#         covariates_t = covariates0_t
#
#     M = torch.unsqueeze(mask_t, 1).float()
#     M_cov = torch.unsqueeze(mask_cov, 1).float()
#     Y = log_counts_t.reshape(1,-1).T
#     Y[M_cov == False] = 0
#     if covariates_t.shape[1] != 0:
#         Y[M_cov == True] = regress_out(covariates_t[mask_cov, :], log_counts_t[mask_cov])
#     res = wrapper_nominal_algo1_matrixLS(Y, genotypes_t.T / 2, M, numpy = False)
#     dof = M.sum() - 2
#     return res, int(mask_t.sum()), dof
#
#     # my version
#
#     # mask_t = raw_counts_t >= count_threshold
#     # if select_covariates:
#     #     covariates_t = filter_covariates(covariates0_t[mask_t], log_counts_t[mask_t])
#     # else:
#     #     covariates_t = covariates0_t[mask_t, 1:]
#     #
#     # residualizer = tensorqtl.Residualizer(covariates_t)
#     # res = cis.calculate_cis_nominal(genotypes_t[:, mask_t] / 2, log_counts_t[mask_t].reshape(1,-1), residualizer)
#     # # [tstat, beta, beta_se, maf, ma_samples, ma_count], samples, dof
#     # return res, covariates_t.shape[0], residualizer.dof
#
#     if select_covariates:
#         covariates_t = filter_covariates(covariates0_t, log_counts_t)
#     else:
#         covariates_t = covariates0_t[:, 1:]
#
#     mask_t = raw_counts_t >= count_threshold
#     residualizer = tensorqtl.Residualizer(covariates_t[mask_t])
#     res = cis.calculate_cis_nominal(genotypes_t[:, mask_t] / 2, log_counts_t[mask_t].reshape(1,-1), residualizer, return)
#     # [tstat, beta, beta_se, maf, ma_samples, ma_count], samples
#     return res, int(mask_t.sum()), residualizer.dof

