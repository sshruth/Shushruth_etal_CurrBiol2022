function err_ = logit_err(Xs)
% function err_ = logist_err(Xs)
% 
%	Computes the negative of the log likelihood of obtaining Data 
%    in the global arrays given a logistic model with
%    parameters in vector Xs.
%	The data is sent in as a global variable (Data) and is in columns 
%    such that Data(:,1) is ones (which must be
%    specified; they represent a bias); Data(:,[2:end-1]) are the values 
%    (linear coefficients) for the other observed variables; Data(:,m) 
%	  is a 1 (correct) or a 0 (incorrect).
%
%  Input:
%		Xs ... array of parameters:
%			Xs(1) 	 is gamma, the guess rate
%			Xs(2) 	 is lambda, the lapse rate
%			Xs(3:end) are the coefficients corresponding to
%							[3:end-1] columns in Data
%
%	Returns: err_, which is -log likelihood

global Data

ps = logist_val(Data(:,1:end-1), Xs(1), Xs(2), Xs(3:end));

% for each data entry. 
% Now calculate the joint probability of obtaining the data set conceived as
% a list of Bernoulli trials.  This is just ps for trials = 1 and 1-ps for
% trials of 0.
p    = [ps(Data(:,end)==1); 1-ps(Data(:,end)==0)];
err_ = -sum(log(p));			% minus log likelihood
