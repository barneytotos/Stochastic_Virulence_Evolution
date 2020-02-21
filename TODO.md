## Morgan

* Deterministic cleaned up a bit (see commit Feb 20), but still have some lingering code to clean
	* tuning won't work currently
	* some lingering "tuning" places where it shouldn't be
	* Still need to check ability for the deterministic model to run tradeoff and efficiency 
	* make sure metadata (beta, gamma values, beta/gamma scaling, etc.; maybe even parameter values?) get attached to sim output (attributes, or make an S3 object)
		* Some progress here but still need to update bookkeeping for the deterministic model (this stuff is returned for stochastic == T)
* Brainstorm what plots we need and what code updates are needed to get us there
	* Create slurm job array bash script for sending these jobs to Stanford cluster
* Some residual bad names that keep popping up. Cleaning these as I go


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
