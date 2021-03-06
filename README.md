# Studying the effects of the North Atlantic Subtropical High on United States Rainfall

Codes for analyzing how the statistical properties of daily rainfall over the conterminous United States depened on the position of the North Atantic Subtropical High (NASH) during the summer season.

The software can be used to download daily rainfall data from the NOAA USHCN dataset and fit Bayesian statistical Models describing rainfall occurrence, probability distribution and extremes, and their dependence on the NASH position.

The codes are in R and use Stan (https://mc-stan.org/) for fitting the Bayesian models and other commonly used R packages. A quick description of the main scripts:

```
read_ushcn_daily.R
```
To download and clean the dataset (1218 daily rainfall records over the CONUS)

```
main_nash_single_station.R
```
To fit the model to a single station dataset and make some nice plots.

```
main_nash_cluster.R
```
Main script to process the entire dataset in a linux cluster. The script ```stats.sh``` can be used to submit the jobs using Slurm.

```
main_nash_plot_results.R 
```
To plot the results for the entire dataset (produce figures in the paper)


