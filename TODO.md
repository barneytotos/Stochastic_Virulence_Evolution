## Morgan

* Deterministic cleaned up a bit (see commit Feb 20), but still have some lingering code to clean
	* Tuning model fixed to the point that it is now "working" except that it is a bit odd in that parasites are evolving directly in efficiency and gamma which jointly determine beta as in the efficiency model. In the tuning model efficiency is defined by two traits (call them phi and psi), which are pre-processed to produce efficiency, which causes the tuning model to more or less just be a slightly modified version of efficiency instead of having a parasite evolving in three traits
	* make sure metadata (beta, gamma values, beta/gamma scaling, etc.; maybe even parameter values?) get attached to sim output (attributes, or make an S3 object)
		* Could consider returning as one object then tidying outside of run_sim(). Currently I have the tidying done in run_sim() and parameters returned as a separate list entry
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
