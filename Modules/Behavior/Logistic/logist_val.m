function p_ = logist_val(dat, gamma, lambda, bs)
% function p_ = logist_val(dat, gamma, lambda, bs)
%
% Uses the logistic model to compute p (probability of
% making a particular decision) for the given data (m rows
% of trials x n columns of data categories) and parameters
% (nx1). Assumes logit(p) (that is, ln(p/(1-p)) is a linear
% function of the data.
% Gamma is the guess rate, lambda is the lapse rate described
% by Abbott's law (see Strasburger 2001, or Finney 1974)

% mxn times nx1 gives mx1
p_ = gamma + (1 - gamma - lambda)./(1+exp(-(dat*bs(:))));
