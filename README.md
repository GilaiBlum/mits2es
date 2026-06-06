README
================

# its2es (Extended for Multiple Disruptions)

This package implements interrupted time series (ITS) analysis for both
continuous and count outcomes, and quantifies the associated effect
size, as described in *Effect size quantification for interrupted time
series analysis: Implementation in R and analysis for Covid-19 research*. 

**Version 0.2.0 Update:** The package has been substantially extended to support **multiple disruption events**, **zero-inflated count models**, **seasonal adjustments via Fourier terms**, and **mixed-effects architectures**. 

The main functions fit an ITS regression model, and then use the fitted
values and the model-based counterfactual values to quantify the effect
size. Effect sizes are reported as Cohen’s *d* (via parametric bootstrapping) for continuous outcomes and Relative Risk (RR) for count outcomes.

An example describing how to install and use this package for a multi-disruption scenario is
described below. A more detailed tutorial, including the data analysis
described in the original paper, is also available with this package (Rmd + pdf
file).

## Installation

You can install the package from its [GitHub repository](https://github.com/GilaiBlum/mits2es). You first
need to install the [devtools](https://github.com/r-lib/devtools) package.

``` r
install.packages("devtools", repos = "[http://cran.us.r-project.org](http://cran.us.r-project.org)")
```

Then install its2es using the install_github function in the devtools package.

```
library(devtools)
install_github("GilaiBlum/mits2es")
```

## Example: Multiple Disruptions

This example demonstrates how to evaluate multiple disruptions (e.g., the onset of the COVID-19 pandemic and its subsequent cessation) using the multiple-interruption architecture.

1. Load the library and the Israel all-cause mortality data.

```
library(its2es)
data <- Israel_mortality
```

2. Define the formula and a vector of intervention start indices. In this example, we set two disruptions: the onset of COVID-19 (March 2020) and a subsequent shift (e.g., January 2021).

```
form <- as.formula("percent ~ time")
start_indices <- c(which(data$Year==2020 & data$Month==3), 
                   which(data$Year==2021 & data$Month==1))
```

3. Fit a linear regression ITS model to the mortality percent. We use multi_model_type = "multiple" to evaluate the distinct effect size of each disruption independently.

```
fit <- its_lm(data = data, 
              form = form, 
              time_name = "time",
              intervention_start_indices = start_indices, 
              multi_model_type = "multiple",
              freq = 12, 
              seasonality = "full", 
              impact_model = "full",
              counterfactual = TRUE, 
              print_summary = FALSE)
```

Console Output:

| Disruption   | Cohen's d | 2.5% CI  | 97.5% CI | P-value |
|--------------|----------:|---------:|----------:|--------:|
| Disruption_1 | 1.038391  | 0.332192 | 1.715101  | 0.0025  |
| Disruption_2 | -0.421530 | -0.985100 | 0.152300 | 0.1250  |

4. Plot the predicted values and counterfactual values. The updated plotting function will automatically handle the multiple boundaries and draw the appropriate disconnected regression lines.

```
p <- plot_its_lm(data = fit$data, 
                 intervention_start_indices = start_indices, 
                 y_lab = "All-cause mortality percent", 
                 response = "percent", 
                 date_name = "Date")
p
```


