import pandas as pds
import numpy as np
import scipy
import tqdm
from random import choices
import defs

from bartpy.sklearnmodel import SklearnModel
from bartpy.features.featureselection import SelectNullDistributionThreshold, SelectSplitProportionThreshold
from bartpy.diagnostics.features import *
from bartpy.extensions.baseestimator import ResidualBART

def Simulation_Gibbs_Sampler(params):
    
    DF = params["DF"]
    identifiers = params["identifiers"]
    aaa = params["aaa"]
    methods = params["methods"]
    covariates = params["covariates"]
    nbr_in_common = params["nbr_in_common"]
    nbr_iter = params["nbr_iter"]
    approx_integral = params["approx_integral"]
    a_sigma = params["a_sigma"]
    b_sigma = params["b_sigma"]
    a_sigma2 = params["a_sigma2"]
    b_sigma2 = params["b_sigma2"]
    a = params["a"]
    b = params["b"]
    alpha_pi = params["alpha_pi"]
    beta_pi = params["beta_pi"]
    montecarlo_size_CI = params["montecarlo_size_CI"]
    coverage = params["coverage"]
    burn_in = params["burn_in"]

    common_records = DF.sample(n = nbr_in_common)

    B = pds.concat([DF.sample(n = 250), common_records]).drop(['Y'], axis = 1)
    B = B.reset_index(drop=True)
    n_B = B.shape[0]

    A = pds.concat([DF.sample(n = 430), common_records])[list(identifiers.keys())+['Y']]
    A = A.reset_index(drop=True)
    n_A = A.shape[0]

    B.loc[np.random.choice(np.arange(n_B), size=5), 'id1'] = ''
    A.loc[np.random.choice(np.arange(n_A), size=5), 'id3'] = ''
    B.loc[np.random.choice(np.arange(n_B), size=5), 'id2'] = ''

    B['propensity_score'] = defs.propensity_score(B, covariates, None, False)

    ate_common_records = aaa * common_records['X2'].mean()

    cartesian_product_AB = B.merge(A, how='cross', suffixes=("_B", "_A")) 
    cartesian_product_AB["source_index_B"] = np.repeat(B.index, n_A)
    cartesian_product_AB["source_index_A"] = np.tile(A.index, n_B)

    for linking_var in identifiers.keys():
        method = methods[identifiers[linking_var]]
        df = cartesian_product_AB.filter(regex=linking_var)
        cartesian_product_AB[linking_var+"_comparison"] = np.array([method(a, b) for a,b in zip(df.iloc[:,0], df.iloc[:,1])]).astype(int).reshape(-1,1)

    comparison_vectors = cartesian_product_AB.filter(regex="comparison")

    true_linkage_z = -np.ones(n_B)
    true_linkage_z[B.iloc[-nbr_in_common:,:].index] = A.iloc[-nbr_in_common:,:].index

    intercept = np.ones(cartesian_product_AB.shape[0])
    records_treatment =  np.array(cartesian_product_AB.treatment)
    records_propensity_score = np.array(cartesian_product_AB.propensity_score)
    interaction_term = records_treatment * records_propensity_score
    records_covariates = np.array(cartesian_product_AB[covariates])
    #X = np.concatenate([intercept.reshape(-1,1), records_covariates, records_propensity_score.reshape(-1,1), interaction_term.reshape(-1,1)], axis=1)
    X = np.concatenate([intercept.reshape(-1,1), records_propensity_score.reshape(-1,1), interaction_term.reshape(-1,1)], axis=1)

    records_treatment_miss =  1 - np.array(cartesian_product_AB.treatment)
    interaction_term_miss = records_treatment_miss * records_propensity_score
    #X_miss = np.concatenate([intercept.reshape(-1,1), records_covariates, records_propensity_score.reshape(-1,1), interaction_term_miss.reshape(-1,1)], axis=1)
    X_miss = np.concatenate([intercept.reshape(-1,1), records_propensity_score.reshape(-1,1), interaction_term_miss.reshape(-1,1)], axis=1)


    new_z = -np.ones(n_B)
    new_z[::2] = np.arange(len(new_z[::2]))

    idx_match_A = new_z[new_z>=0]
    idx_match_B = np.nonzero(new_z>=0)[0]
    links = pds.MultiIndex.from_tuples(zip(idx_match_A,idx_match_B))
    pairs = pds.MultiIndex.from_frame(cartesian_product_AB[["source_index_A", "source_index_B"]])
    if not(coverage):
        dict_empirical_count_linkage = {pair:0 for pair in pairs}
    dict_params = {"unmatch":[], "match":[], "sigma_square":[], "betas":[], "sigma2_square":[], "mu2":[], "atel":[]}

    Betas = scipy.stats.multivariate_normal.rvs(np.zeros(X.shape[1]), np.eye(X.shape[1]))
    mu2 = scipy.stats.norm.rvs(0, 1)

    dev = 0.25 * np.std(cartesian_product_AB.Y)
    low_val = int(np.min(cartesian_product_AB.Y) - dev)
    high_val = int(np.max(cartesian_product_AB.Y) + dev)
    nbr_pts = int(high_val - low_val)

    for _ in tqdm.tqdm(range(nbr_iter)):

        ### --  UPDATE THETA --> match, unmatch
        comparison_vectors_for_non_matches = comparison_vectors[~pairs.isin(links)]
        unmatch = scipy.stats.beta.rvs(comparison_vectors_for_non_matches.sum(axis=0) + a, (1-comparison_vectors_for_non_matches).sum(axis=0) + b)
        comparison_vectors_for_matches = comparison_vectors[pairs.isin(links)]
        match = scipy.stats.beta.rvs(comparison_vectors_for_matches.sum(axis=0) + a, (1-comparison_vectors_for_matches).sum(axis=0) + b)
        dict_params["unmatch"].append(unmatch)
        dict_params["match"].append(match)

        ### --  UPDATE COEF OUTCOME MODEL (MATCHES AND NON MATCHES) --> Betas --> mu2, sigma2_square --> outcome model distribution
        data_for_matches = cartesian_product_AB[pairs.isin(links)]
        outcome_for_matches =  np.array(data_for_matches.Y)
        intercept = np.ones(len(outcome_for_matches))
        treatment_for_matches =  np.array(data_for_matches.treatment)
        linked_records_propensity_score = np.array(data_for_matches.propensity_score)
        interaction_term = treatment_for_matches * linked_records_propensity_score
        linked_records_covariates = np.array(data_for_matches[covariates])
        #K = np.concatenate([intercept.reshape(-1,1), linked_records_covariates, linked_records_propensity_score.reshape(-1,1), interaction_term.reshape(-1,1)], axis=1)
        K = np.concatenate([intercept.reshape(-1,1), linked_records_propensity_score.reshape(-1,1), interaction_term.reshape(-1,1)], axis=1)
        n_AB = sum(new_z>=0)
        sigma_square = scipy.stats.invgamma.rvs(a_sigma + n_AB/2, loc=0, scale=np.linalg.norm(outcome_for_matches - K @ np.array(Betas))**2 / 2 + b_sigma) 
        Sigma_beta = np.linalg.inv( (K.T @ K)/sigma_square + np.eye(K.shape[1]) )
        mu_beta = Sigma_beta @ (K.T @ outcome_for_matches)/sigma_square
        Betas = scipy.stats.multivariate_normal.rvs(mu_beta, Sigma_beta)
        residuals = cartesian_product_AB.Y - X @ Betas.T
        estimated_variance = residuals.T @ residuals / (len(residuals) - (X.shape[1]+1))
        distr_linked_outcomes = scipy.stats.norm.pdf(residuals, 0, np.sqrt(estimated_variance))
        dict_params["sigma_square"].append(sigma_square)
        dict_params["betas"].append(Betas)

        data_for_non_matches = A[~A.index.isin(idx_match_A)]
        outcome_for_non_matches =  np.array(data_for_non_matches.Y)
        sigma2_square = scipy.stats.invgamma.rvs(a_sigma2 + (n_A - n_AB)/2, loc=0, scale=b_sigma2 + sum((outcome_for_non_matches - mu2)**2)/2) 
        sigma_mu2_square = 1 / ((n_A - n_AB)/sigma2_square + 1)
        m_mu2 = sigma_mu2_square * ((outcome_for_non_matches).sum() / sigma2_square)
        mu2 = scipy.stats.norm.rvs(m_mu2, np.sqrt(sigma_mu2_square))
        distr_non_linked_outcomes = scipy.stats.norm.pdf(cartesian_product_AB.Y, mu2, np.sqrt(sigma2_square))
        dict_params["sigma2_square"].append(sigma2_square)
        dict_params["mu2"].append(mu2)
        
        ### --  UPDATE PROBABILITIES
        w1 = ( np.multiply( comparison_vectors, np.log(match/unmatch) ) + np.multiply( 1-comparison_vectors, np.log((1-match)/(1-unmatch)) ) ).sum(axis=1)
        w2 = np.log(distr_linked_outcomes / distr_non_linked_outcomes)
        probabilities = np.array(np.exp(w1+w2))
        probabilities = probabilities.reshape(n_B, n_A)
        
        for j in range(n_B):
            n_AB_ = (np.delete(new_z, j)>=0).sum()
            proba_for_unmatch = (n_A - n_AB_) * (n_B - n_AB_ - 1 + beta_pi) / (n_AB_ + alpha_pi)
            not_possible_values = list( set(np.delete(new_z, j).astype(int)) - set([-1]) )
            proba = probabilities[j,:].copy()
            proba[not_possible_values]  = 0
            proba = np.append(proba, proba_for_unmatch)
            possible_values = np.arange(n_A)
            possible_values = np.append(possible_values, -1)
            new_z[j] = choices(possible_values, weights = proba)[0]

        ### ATEL
        y_miss_new = np.random.uniform( low_val, high_val, size = (data_for_matches.shape[0], nbr_pts))
        X_miss_new = X_miss[pairs.isin(links)]
        Betas_miss_new = scipy.stats.multivariate_normal.rvs(mu_beta, Sigma_beta, size=(nbr_pts, approx_integral))
        y_miss_new = np.repeat(y_miss_new[:, np.newaxis, :], approx_integral, axis=1)
        ITS = np.mean( y_miss_new - ( Betas_miss_new @ X_miss_new.T ).T, axis = 1)
        missing_outcomes_new = np.array([ defs.inverse_transform_sampling(ITS[indiv], len(ITS[indiv])) for indiv in range(ITS.shape[0]) ])
        Y0_new = (data_for_matches.treatment==0) * data_for_matches.Y + (data_for_matches.treatment==1) * missing_outcomes_new.flatten()
        Y1_new = (data_for_matches.treatment==1) * data_for_matches.Y + (data_for_matches.treatment==0) * missing_outcomes_new.flatten()
        atel_new = np.mean(Y1_new - Y0_new)
        dict_params["atel"].append(atel_new)

        ### _ UPDATE INDICES FOR MATCHES
        idx_match_A = new_z[new_z>=0]
        idx_match_B = np.nonzero(new_z>=0)[0]
        links = pds.MultiIndex.from_tuples(zip(idx_match_A,idx_match_B))

        if not(coverage):
            for link in links:
                dict_empirical_count_linkage[link] += 1

    if burn_in:
        starting_point = 4*int(nbr_iter/5)
        ate_post_burn_in = dict_params["atel"][starting_point:]
    else:
        starting_point = 0
        ate_post_burn_in = dict_params["atel"]

    if not(coverage):
        v = np.array(list(dict_empirical_count_linkage.values())) / nbr_iter
        k = np.array(list(dict_empirical_count_linkage.keys()))
        v = v[starting_point:]
        k = k[starting_point:]

    nbr_bins = 10 * int( np.max(ate_post_burn_in) - np.min(ate_post_burn_in) )
    ate_estimator = np.array([ defs.inverse_transform_sampling(ate_post_burn_in, nbr_bins) for _ in range(montecarlo_size_CI) ])
    lower_bound, upper_bound = np.quantile(ate_estimator, [0.025, 0.975])
    mean = np.mean(ate_estimator)

    if not(coverage):
        return dict_params, v, k, true_linkage_z, mean, (mean - lower_bound, upper_bound - mean), ate_common_records, mean - ate_common_records
    else:
        return mean, (mean - lower_bound, upper_bound - mean), ate_common_records, mean - ate_common_records