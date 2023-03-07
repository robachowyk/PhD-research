import pandas as pds
import numpy as np
import scipy
import tqdm
from random import choices
import defs

def do_Gibbs_iteration(A, n_A, n_B, cartesian_product_AB, comparison_vectors, new_z, X, X_miss, covariates, pairs, links, low_val, high_val, nbr_pts, approx_integral, dict_params, a, b, a_sigma, b_sigma, a_sigma2, b_sigma2, alpha_pi, beta_pi):

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
    K = np.concatenate([intercept.reshape(-1,1), linked_records_covariates, linked_records_propensity_score.reshape(-1,1), interaction_term.reshape(-1,1)], axis=1)
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
        val = choices(possible_values, weights = proba)[0]
        new_z[j] = val

    ### ATEL
    y_miss_new = np.random.uniform( low_val, high_val, size = (data_for_matches.shape[0], nbr_pts))
    X_miss_new = X_miss[pairs.isin(links)]
    Betas_miss_new = scipy.stats.multivariate_normal.rvs(mu_beta, Sigma_beta, size=(nbr_pts, approx_integral))
    y_miss_new = np.repeat(y_miss_new[:, np.newaxis, :], approx_integral, axis=1)
    ITS = np.mean( y_miss_new - ( Betas_miss_new @ X_miss_new.T ).T, axis = 1)
    missing_outcomes_new = np.array([ defs.inverse_transform_sampling(ITS[indiv], len(ITS[indiv])) for indiv in range(ITS.shape[0]) ])
    new_Y0 = (data_for_matches.treatment==0) * data_for_matches.Y + (data_for_matches.treatment==1) * missing_outcomes_new.flatten()
    new_Y1 = (data_for_matches.treatment==1) * data_for_matches.Y + (data_for_matches.treatment==0) * missing_outcomes_new.flatten()
    new_atel = np.mean(new_Y1 - new_Y0)
    dict_params["atel"].append(new_atel)

    ### _ UPDATE INDICES FOR MATCHES
    idx_match_A = new_z[new_z>=0]
    idx_match_B = np.nonzero(new_z>=0)[0]
    links = pds.MultiIndex.from_tuples(zip(idx_match_A,idx_match_B))

    return links, new_z, dict_params

