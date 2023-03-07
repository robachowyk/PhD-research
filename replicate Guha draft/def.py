import pandas as pds
import numpy as np
import statsmodels.api as sm
import scipy
import tqdm
from random import choices

def logit(p):
    return np.log(p/(1-p))

def minmaxscaler(X):
    return (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    
def propensity_score(DF, covariates, scaler, convert_to_logit):
    
    """ Compute propensity score estimates: the probability (logistic regression) that an observation is treated or not conditioned on some covariates.
        These estimates are built conditionaly on covariates passed using a logit after transformation by scaler (when one is specified).
        Estimated probabilities can be converted into logit (convert_to_logit parameter).

        DF:                dataframe,
        covariates:        list of strings for covariates variable in DF,
        scaler:            sklearn.preprocessing function scaler for exemple,
        convert_to_logit:  boolean for converting probabilities to logit when building the propensity score estimates based on a logistic regression
    """
    exog = np.array(DF[covariates])
    if scaler != None:
        exog = scaler(exog)
    intercept = np.ones(DF.shape[0]).reshape(-1,1)
    exog = np.concatenate((exog, intercept), axis=1)
    model = sm.Logit(DF.treatment, exog).fit(disp=0)
    predictions = model.predict(exog)
    if convert_to_logit:
        return logit(predictions)
    else: 
        return predictions

def mcmc_func( params ):
    
    """   params = {'cartesian_product_AB':
                    'z_init':  
                    'comparison_vectors':  
                    'covariates':
                    'n_A':
                    'n_B':
                    'X':}
    """

    cartesian_product_AB = params['cartesian_product_AB']
    new_z = params['z_init']
    comparison_vectors = params['comparison_vectors']
    covariates = params['covariates']
    n_A = params['n_A']
    n_B = params['n_B']
    X = params['X']

    nbr_iter = 10
    a_sigma, b_sigma, a_sigma2, b_sigma2, a, b, alpha_pi, beta_pi = 1, 1, 1, 1, 1, 1, 1, 1
    Betas = scipy.stats.multivariate_normal.rvs(np.array([0,0,0]), np.eye(3))
    mu2 = scipy.stats.norm.rvs(0, 1)

    pairs = pds.MultiIndex.from_frame(cartesian_product_AB[["source_index_A", "source_index_B"]])
    dict_empirical_count_linkage = {pair:0 for pair in pairs}

    for _ in tqdm.tqdm(range(nbr_iter)):

        ### _ UPDATE INDICES FOR MATCHES
        idx_match_A = new_z[new_z>=0]
        idx_match_B = np.nonzero(new_z>=0)[0]
        links = pds.MultiIndex.from_tuples(zip(idx_match_A,idx_match_B))
        pairs = pds.MultiIndex.from_frame(cartesian_product_AB[["source_index_A", "source_index_B"]])

        ### --  UPDATE THETA --> match, unmatch
        comparison_vectors_for_non_matches = comparison_vectors[~pairs.isin(links)]
        unmatch = scipy.stats.beta.rvs(comparison_vectors_for_non_matches.sum(axis=0) + a, (1-comparison_vectors_for_non_matches).sum(axis=0) + b)
        comparison_vectors_for_matches = comparison_vectors[pairs.isin(links)]
        match = scipy.stats.beta.rvs(comparison_vectors_for_matches.sum(axis=0) + a, (1-comparison_vectors_for_matches).sum(axis=0) + b)

        ### --  UPDATE COEF OUTCOME MODEL (MATCHES AND NON MATCHES) --> Betas --> mu2, sigma2_square --> outcome model distribution
        data_for_matches = cartesian_product_AB[pairs.isin(links)]
        outcome_for_matches =  np.array(data_for_matches.Y)
        treatment_for_matches =  np.array(data_for_matches.treatment)
        linked_records_propensity_score = propensity_score(data_for_matches, covariates, None, False)
        interaction_term = treatment_for_matches * linked_records_propensity_score
        K = np.concatenate([np.ones(len(outcome_for_matches)).reshape(-1,1), linked_records_propensity_score.reshape(-1,1), interaction_term.reshape(-1,1)], axis=1)
        n_AB = sum(new_z>=0)
        sigma_square = scipy.stats.invgauss.rvs(a_sigma + n_AB/2, np.linalg.norm(outcome_for_matches - K @ np.array(Betas))**2 / 2 + b_sigma)
        Sigma_beta = np.linalg.inv( (K.T @ K)/sigma_square + np.eye(K.shape[1]) )
        mu_beta = Sigma_beta @ (K.T @ outcome_for_matches)/sigma_square
        Betas = scipy.stats.multivariate_normal.rvs(mu_beta, Sigma_beta)
        data_for_non_matches = cartesian_product_AB[~pairs.isin(links)]
        outcome_for_non_matches =  np.array(data_for_non_matches.Y)
        sigma2_square = scipy.stats.invgauss.rvs(a_sigma2 + (n_A - n_AB)/2, b_sigma2 + sum((outcome_for_non_matches - mu2)**2)/2)
        sigma_mu2_square = 1/((n_A - n_AB)/sigma2_square + 1)
        m_mu2 = sigma_mu2_square * (outcome_for_non_matches).sum() / sigma2_square
        mu2 = scipy.stats.norm.rvs(m_mu2, np.sqrt(sigma_mu2_square))
        residuals = cartesian_product_AB.Y - X @ Betas.T
        estimated_variance = residuals.T @ residuals / (len(residuals) - (X.shape[1]+1))
        distr_linked_outcomes = scipy.stats.norm.pdf(residuals, 0, np.sqrt(estimated_variance))
        distr_non_linked_outcomes = scipy.stats.norm.pdf(cartesian_product_AB.Y, mu2, np.sqrt(sigma2_square)) # scipy.stats.norm.rvs(0, 1), 10000)

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
            idx_max_proba = np.argpartition(proba, -5)[-5:]
            possible_values = np.arange(n_A)
            possible_values = np.append(possible_values, -1)
            val = choices(possible_values[idx_max_proba], weights = proba[idx_max_proba])[0]
            new_z[j] = val

        for link in links: # check burn-in
            dict_empirical_count_linkage[link] += 1
    
    v = np.array(list(dict_empirical_count_linkage.values())) / nbr_iter
    k = np.array(list(dict_empirical_count_linkage.keys()))

    return v, k