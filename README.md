# CCMMshadowrateVAR-code

# Code for „Forecasting with Shadow-Rate VARs“ by Carriero, Clark, Marcellino and Mertens (2021)

Andrea Carriero (Queen Mary University of London), Todd Clark (Federal Reserve Bank of Cleveland), Massimiliano Marcellino (Bocconi University, IGIER and CEPR) and Elmar Mertens (Deutsche Bundesbank)

The usual disclaimers apply; the views and results conveyed in our research are solely those of the authors and do not necessarily reflect the views of the Federal Reserve Bank of Cleveland, the Federal Reserve System, the Eurosystem, or the Deutsche Bundesbank.

The working paper and supplementary appendices are available here: https://doi.org/10.26509/frbc-wp-202109

*WORK IN PROGRESS*

## Directories
All core scripts are in the main directory. In addition, there are the following subdirectories:
- `data` for data construction
- `matlabtoolbox` for general utilities (also available at https://github.com/elmarmertens/em-matlabbox)

The default data file is `fredMD20baa-2021-07.csv` (based on the 2021-07 vintage of FRED-MD, available at https://research.stlouisfed.org/econ/mccracken/fred-databases/).

## General notes

- All scripts set the MATLAB path to point to toolboxes in `matlabtoolbox`. In addition, most scripts collect output in a temporary directory, which is by default created as subfolder `foo` within the main directory. Edit `localtemp.m` to change the location of this temp directory. Output is collected in a LaTeX file, which is also compiled at the end of each script (provided a LaTeX installation can be found on the path). To control the compilation of output, please edit `finishwrap.m`. To avoid collecting output files, comment out the call to `initwrap` in each script (and make sure to define instead a variable called `wrap` that is set to empty).

- The code requires a recent version of Matlab (we used Matlab versions 2019a-2021a) including access to Matlab’s Statistics and Machine Learning Toolbox. The codes employ `parfor` loops that are executed in parallel when a `parpool` has been created in Matlab, which requires availability of the Matlab Parallel Computing Toolbox (otherwise the loops will be executed sequentially).

- The main directory contains a bash script `gobatch.sh` that can be used to launch a sequence of multiple Matlab scripts from the shell. The Matlab scripts are executed in sequence *and in separate Matlab sessions*. Each Matlab session opens a parallel pool. For example, the shell command `sh gobatch.sh goVAR.m goVARshadowrate.m` will launch a command line session of Matlab, start a parallel pool, and then execute `goVAR.m`; once `goVAR.m` has been executed, the Matlab sessions classes, a new one is reopened for execution of `goVARshadowrate.m`. (The shell script supports as many command line arguments as supported by bash and has been written for use on macOS and Linux.) Alternatively, Matlab scripts can, of course, equally be called interactively on the Matlab GUI’s command line.

- To execute individual models for a given sample, use scripts called `do*.m`. To launch out-of-sample runs for a given model, consider scripts called `go*.m`. After computing the out-of-sample runs each of these `go*.m` scripts stores results for further post-processing in a `*.mat` file. To produce comparison tables of various formats, please use Matlab scripts called `oosEvaluationTables*.m`, which assume that previously stored `*.mat` files are stored in a directory indicated at the beginning of each of theses scripts.  

## Model names
- XXX
