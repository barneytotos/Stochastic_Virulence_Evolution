## Morgan

* Only test that is missing at this point is stochastic - tuning. 
	* deterministic - tuning exists, but it is unclear if we want tests for the tuning model in the main test file because we may not run the tuning model anyway
* make sure metadata (beta, gamma values, beta/gamma scaling, etc.; maybe even parameter values?) get attached to sim output (attributes, or make an S3 object)
	* Currently I have tidying done in run_sim() for deterministic and parameters returned as a separate list entry. Could consider returning as one object then tidying outside of run_sim(). 

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
