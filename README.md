# portfolio_optimization
Portfolio Optimization in MATLAB
By HX for Reference Number  16-067

The repo contains code for portfolio optimization and portfolio performance modelling. See files for specific comments.

&nbsp;

#basic_requirements
There must be a data directory containing a csv file with the tickers of assets that is been optimized. The name of the csv file must be in _name_\_assets.csv where italics is the name of the dataset. If the csv file has multiple columns and rows, the ticker names must be in the first column.

Example ticker files for the Dow Jones, S&P 500, and the TSX are provided: `data/djia_assets.csv`, `data/sp500_assets.csv`, `data/sptsx_assets.csv`

If there is any issue download or creating an asdata file, an example is available in the data directory. Run command:
```
copyfile('data/djia_asdata_example.mat','data/djia_asdata.mat');
```

&nbsp;

# workflow
See `example_script_01.m`, which goes through all the steps and performs three portfolio optimizations.

__Create Data__: `download_data` &raquo; `ascsv_to_asdata`

__Portfolio Modeling__: `mean_cvar_optimization` or `mean_variance_optimization` &raquo; `calc_portfolio_return`

&nbsp;

# sample_results
These are results from running first model in `example_script_01.m` in chart form

![Efficient Frontier and Efficient Portfolios](/images/1-1.png)

![Portfolio Performance](/images/1-3.png)

&nbsp;

# disclaimer
This is an experimental codebase, do not use in production settings.

Past performance is not indicative of future results.

&nbsp;

#task_list
This repo is a work in progress, possible future tasks:
- Test cases
- Better comments
- Additional example scripts
- Black-Litterman model
- Refactor to classes if desired
- Upload graphics, chart functions
