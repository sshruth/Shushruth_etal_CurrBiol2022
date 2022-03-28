function [fits_,llik_,pred_,sems_] = logist_fit(data, gamma, lambda)
%
% LOGIST_FIT fits a logistic function to data using maximum likelihood
%   maximization under binomial assumptions.  It uses logit_err for
%   error calculation. The logistic assumes a linear (on the coefficients)
%   model of log-odds-ratio (ln(p/(1-p))) = SUM(B_i*x_i) = f(X), where B_i are
%   the coefficients to fit and x_i is the data. Also incorporates a guess
%   rate (gamma) and lapse rate (lambda). Thus, 
%         p = gamma + (1 - gamma - lambda) * 1/(1 + exp(-f(x))
% 
% Usage: [fits,llik,pred,sems] = logist_fit(data,gamma,lambda)
% Input values are:
%	  data,   Must be organized as rows of observations.
%  	  		   The 1st column should be a 1 if you want to estimate a 
%				   'bias' term, columns 2 thru n-1 are the values for the 
%				   independent variables, and the last column contains 1 or 
%				   0 for the choice.
%	  gamma,  Specifies guess rate. Must be on the interval of [0,1]; if
%              not (or not given), it's a free parameter.
%				   Note: set to 0.5 for unbiased forced choice
%	  lambda, Specifies lapse rate. Must be on the interval of [0,1]; if
%              not (or not given), it's a free parameter.
%	Returns values are:
%	  fits,   [gamma; lambda; betas]
%               gamma is guess rate; if given as input, returns as is
%               lambda is guess rate; if given as input, returns as is
%               betas are linear coefficients for columns 1:end-1 of data.
%	  llik,    The log likelihood of the fit
%	  pred,    A vector representing the prob(y=1) under the fit.
%    sems,    Standard errors of the fits. [gamma_sem lasem betas_sems]

% 	 9/5/2001 jig spawned from logistfit
%   4/6/2001 standard errors added (JD)
%   3/22/2001 update for new matlab (replace fmins)
%   9/17/94  M Shadlen

global Data

if nargin < 1 | isempty(data)
  return;
end
Data = data;

% Generate a first guess. For the betas, simply convert all the 1's to .99 and
% all the 0's to .01's and do a quick regression on the logits.  This is
% roughly converting the 1's to logit of 3 and 0's to logit of -3
y      = (6 * data(:,end)) - 3;
guess  = [0; 0; data(:,1:end-1)\y]';  % this is a 1st guess
nbetas = length(guess) - 2;           % useful variable

% set lower, upper bounds for the constrained fit
lb = [0 0 -inf*ones(1,nbetas)];
ub = [1 1  inf*ones(1,nbetas)];

% if gamma/lambda given, constrain to use those values
if nargin >= 2 & ~isempty(gamma) & gamma >= 0 & gamma <= 1
  guess(1) = gamma;
  lb(1)    = gamma;
  ub(1)    = gamma;
  gamma_fit = 0; % flag for later
else
  gamma_fit = 1;
end
if nargin >= 3 & ~isempty(lambda) & lambda >= 0 & lambda <= 1
  guess(2) = lambda;
  lb(2)    = lambda;
  ub(2)    = lambda;
  lambda_fit = 0; % flag for later
else
  lambda_fit = 1;
end

% Do the fit
fits_ = fmincon('logit_err', guess, [], [], [], [], lb, ub)';

% to check the hessian...
% [X,F,E,O,L,G,HESSIAN]=fmincon('logit_err',guess,[],[],[],[],lb,ub);

% return the log likelihood
if nargout > 1
  llik_ = -logit_err(fits_);
end

% return the predicted probability (of a correct response) for each observation
if nargout > 2
  pred_ = logist_val(Data(:,1:end-1), fits_(1), fits_(2), fits_(3:end));
end

% standard errors
% The covariance matrix is the negative of the inverse of the 
% hessian of the natural logarithm of the probability of observing 
% the data set given the optimal parameters.
if nargout > 3

  % the hessian is a mess. A couple of facts to keep in mind:
  %  - "the probability of observing the data set" means the product
  %     of the predicted probabilities of observing the result on each
  %     trial, which is just logist_val for correct trials but 1 - 
  %     logist_val for incorrect trials.
  %  - We take the natural logarithm of this product, which means that we
  %     can instead compute the sum of the individual logs... which is
  %     just P = SUM(ln(p(Xs, Data, correct))) + SUM(ln(1-p(Xs, Data, incorrect)))
  %  - We compute the matrix of second derivatives of P with respect to each pair
  %     of parameters (i.e., gamma (g), lambda (l), and each beta (bs))
  %  - it is mirror-symmetric, since d^2P/dxdy = d^2P/dydx
  %  - the diagonal is, of course, d^2P/dx^2
  %  - We will build it in three sections:
  %     1) a vector of [d^2P/(dg)(dbs_i) d^2P/(dg)^2   d^2P/(dg)(dl)]
  %     2) a vector of [d^2P/(dl)(dbs_i) d^2P/(dl)(dg) d^2P/(dl)^2  ]
  %     3) a matrix of d^2P/(dbs_i)(dbs_j)

  % compute correct/incorrect differently (p vs. 1-p)
  Lc = Data(:,end) == 1;
  Li = Data(:,end) == 0;
  % some useful variables to make the computations easier
  ga = fits_(1);
  la = fits_(2);
  bs = fits_(3:end);
  zc = exp(-1*(Data(Lc,1:end-1) * bs));
  zi = exp(-1*(Data(Li,1:end-1) * bs));
  kc = (ga * zc + 1 - la).^2;
  ki = (-zi + ga*zi - la).^2;
  
  % first build matrix #3, above, which is wrt the betas
  Hbs = zeros(nbetas, nbetas);
  Ac  = -((ga*zc.^2+la-1).*(-1+ga+la).*zc)./(kc.*(1+zc).^2);
  Ai  = -((-zi.^2+ga.*zi.^2+la).*(-1+ga+la).*zi)./(ki.*(1+zi).^2);
  for i = 1:nbetas
	 for j = 1:i
		Hbs(i,j) = sum(Ac.*Data(Lc,i).*Data(Lc,j)) + ...
			 sum(Ai.*Data(Li,i).* Data(Li,j));
		Hbs(j,i) = Hbs(i,j); 
	 end
  end
	
  % build the Hessian, whose size depends on whether we are considering
  % ga and/or la
  if gamma_fit == 1 & lambda_fit == 1
	 H = zeros(nbetas+2,nbetas+2);
	 H(3:end,3:end) = Hbs;
  elseif gamma_fit == 1 | lambda_fit == 1
	 H = zeros(nbetas+1,nbetas+1);
	 H(2:end,2:end) = Hbs;
  else
	 H = Hbs;
  end
  
  % now fill in H for #1 (second derivatives wrt gamma)
  % if we've fit gamma, this is always the first row & col
  % of the matrix
  if gamma_fit == 1
	 H(1,1) = -sum(zc.^2./kc) - sum(zi.^2./ki);
	 
	 if lambda_fit == 1
		H(1,2) = sum(zc./kc) + sum(zi./ki);
		H(2,1) = H(1,2);
	 end
	 
	 bind = lambda_fit + 1;
	 for i = 1:nbetas
		H(bind+i, 1) = (la - 1)*sum(Data(Lc,i).*zc./kc) + ...
			 la*sum(Data(Li,i).*zi./ki);
		H(1, bind+i) = H(bind+i, 1);
	 end
  end
  
  lind = gamma_fit + 1;
  if lambda_fit == 1
	 H(lind,lind) = -sum(1./kc) - sum(1./ki);
	 for i = 1:nbetas
		H(lind+i,lind) = -ga*sum(zc.*Data(Lc,i)./kc) + ...
					 (1-ga)*sum(zi.* Data(Li,i)./ki);
		H(lind,lind+i) = H(lind+i,lind);
	 end
  end

  % now compute the negative of the inverse, and take the diagonal
  cov_mat=-H^(-1);
  se=sqrt(diag(cov_mat));
  
  % put it into the output array
  if gamma_fit == 1 & lambda_fit == 1
	 sems_ = se;
  elseif gamma_fit == 1
	 sems_ = [se(1); 0; se(2:end)];
  elseif lambda_fit == 1
	 sems_ = [0; se];
  else
	 sems_ = [0; 0; se];
  end
end
