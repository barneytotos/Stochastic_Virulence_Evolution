## Morgan

* Only test that is missing at this point is stochastic - tuning. 
	* exists for deterministic - tuning, but unclear if we want tuning tests in the main test file because we may not run the tuning model anyway
* make sure metadata (beta, gamma values, beta/gamma scaling, etc.; maybe even parameter values?) get attached to sim output (attributes, or make an S3 object)
	* Currently I have the tidying done in run_sim() and parameters returned as a separate list entry. Could consider returning as one object then tidying outside of run_sim(). 
* Code for all 6 models (stochastic or deterministic for no tradeoff, tradeoff only, and efficiency) running on sherlock
	* got an error of negative values returned in abundance for one of the ode so need to check that...
	* when returned check that the sims are returning what we need to plot what we want to talk about in the paper

## Someone

* Different levels of detail: what do we need to get ensembles (rather than summaries) back? Or more quantiles? (maybe a "quantiles" argument, i.e. a vector of quantiles of traits to return [marginal, i.e. for each trait: what about 2D summaries?])
* house of cards mutation model for no-tradeoff case?

* package structure:
    * top level - analyses for the paper
	    * R files here are the specific runs we want to do
		* maybe separate directory for analyses/, outputs/ (however you prefer)
	*  in package `hpevosim/`:  http://r-pkgs.had.co.nz/
	    * document & export functions as appropriate
		* ~~import functions as appropriate (stats::rbinom or ##' @importFrom stats rbinom)~~
		* ~~NSE confuses the checker (`mean_postrait <- sd_postrait <- NULL` or `globalVariables()`?~~ https://www.google.com/search?client=firefox-b-d&q=NSE+R+%22no+visible+binding%22
https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
