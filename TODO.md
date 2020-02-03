## Morgan

* ~~get model "2B" (efficiency, zero correlation) running on Stanford cluster~~
* Update "vir_evo_stochas_eff_bulk.R" so that it can run in parallel on Stanford cluster
* ~~Combine get_mut and update_mut_pt~~
* ~~List of naming updates (I like the positive_trait, negative_trait idea, but should integrate it throughout)~~
  - ~~alpha and gamma depending on the model~~
  - SIR vs host_mort vs host resistance and tolerance evolution or no

## Someone

* test deterministic cases
* what do we need to get ensembles (rather than summaries) back? Or more quantiles?
* house of cards mutation model for no-tradeoff case?
* package structure:
    * top level - analyses for the paper
	    * R files here are the specific runs we want to do
		* maybe separate directory for analyses/, outputs/ (however you prefer)
	* pkg/ (hpevosim/ or whatever we want to call the package): http://r-pkgs.had.co.nz/
	    * document & export functions as appropriate
		* import functions as appropriate (stats::rbinom or ##' @importFrom stats rbinom)
		* NSE confuses the checker (`mean_postrait <- sd_postrait <- NULL` or `globalVariables()` ? https://www.google.com/search?client=firefox-b-d&q=NSE+R+%22no+visible+binding%22
https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
