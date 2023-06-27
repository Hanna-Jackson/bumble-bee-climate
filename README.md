# Bumble Bee Climate Instruction Manual
![BC_Bombus occidentalis male_Goldenrod](https://github.com/Hanna-Jackson/bumble-bee-climate/assets/71409828/988e7d0b-cf03-4afc-a4e3-a34cb0fe88c0)
Photo Credits: Sarah Johnson, coauthor 

Welcome to the code repository for Jackson et al 2022: ["Climate change winners and losers among North American bumblebees"](https://royalsocietypublishing.org/doi/10.1098/rsbl.2021.0551)

Thank you so much for taking the time to look through my work! 
Any feedback on the code and analysis is appreciated.

If you have issues, feel free to email me at `hmj2@sfu.ca` and I can try to help! 

Good luck!
<br>
\- Hanna Jackson

## Using `Main.R` 
[`Main.R`](https://github.com/Hanna-Jackson/bumble-bee-climate/blob/main/Main.R) contains everything needed to load and prepare the data and run the analysis. 

I have set up the file structure of this repository to mirror the file structure needed to run and save everything such that the code will work and all the necessary folders will be already made for objects to be saved to. 

### Setting Parameters of the Analysis
All "options" for this analysis are specified in one object called `args`. Everything from the site resolution to filtering options for the data prep to which model we're going to use and how long it will run for. Everything done in [`Main.R`](https://github.com/Hanna-Jackson/bumble-bee-climate/blob/main/Main.R) hereafter will use these parameters that you choose. 

### Data Prep
This section of [`Main.R`](https://github.com/Hanna-Jackson/bumble-bee-climate/blob/main/Main.R) sources 10 files that go through each step of the data processing (you can find those files in the [data_prep](https://github.com/Hanna-Jackson/bumble-bee-climate/tree/main/data_prep) folder. When you source a file all it does is **load in the function** that that specific file has in it. Then all we have to do is **call that function** with the arguments (from the `args` object) we want and all the code runs, doing one step in our data processing. At the end of each file the important objects get saved to your working directory and then when a different function later on in the process needs that data it just loads it!
<br>
RData objects containing the spatial object of all the sites are given already in the [data_prep/sites](https://github.com/Hanna-Jackson/bumble-bee-climate/tree/main/data_prep/sites) folder (But you can also make them yourself if you want to in the data prep workflow). 

#### Data Download
You will need to download the **bumble bee data** from the [Data Dryad repository](https://datadryad.org/stash/dataset/doi:10.5061%2Fdryad.c59zw3r8f) as well as all the **climate data** from the [CHELSAcruts database](https://chelsa-climate.org/chelsacruts/). This data should go in the [data_prep/data](https://github.com/Hanna-Jackson/bumble-bee-climate/tree/main/data_prep/data) folder to make the code I have written work with it.

If you don't want to do that (the climate files are very large and there is many of them, so automating the download process is basically required), you can skip this section and move right on to model running and instead load the pre-processed data that I provide in the correct directory. The only issue that arises with that method is that the files for the 50x50km and 100x100km site resolutions are too large to be uploaded, so only 250x250km is avaliable. 

**Floral resource data** is already available in the [data_prep/data](https://github.com/Hanna-Jackson/bumble-bee-climate/tree/main/data_prep/data) folder. 

### Model Runs
To run the model we first source a file [`11.Prep_For_Run_JAGS.R`](https://github.com/Hanna-Jackson/bumble-bee-climate/blob/main/data_prep/11.Prep_For_Run_JAGS.R) that loads the packages and data (from the [data_prep/saved/for_analysis](https://github.com/Hanna-Jackson/bumble-bee-climate/tree/main/data_prep/saved/for_analysis) file) that we need to run our model, scales the data and packages it into list called `my.data` that JAGS will use. 

Next we **source the model file** (which model file is specified in `args` by `args$model`). The models themselves are stored in the [models](https://github.com/Hanna-Jackson/bumble-bee-climate/tree/main/models) folder. 

**IMPORTANT NOTE:**
*I have written comments in our two main models [`env_model.R`](https://github.com/Hanna-Jackson/bumble-bee-climate/blob/main/models/env_model.R) and [`era_model.R`](https://github.com/Hanna-Jackson/bumble-bee-climate/blob/main/models/era_model.R) that I hope will make understanding the model relatively easy for someone unfamiliar with JAGS or with occupancy modeling. I encourage you to look at these files, as they are the core of our analysis!*

Next we **run the model** using the call to the `jags()` function. This requires that you have [JAGS](https://mcmc-jags.sourceforge.io) (the software we use to run our model) downloaded onto your computer. For context, all models take about 15 hours to run on my 2020 13" M1 Macbook Pro. 

Running the save code after that saves the model output to the [outputs](https://github.com/Hanna-Jackson/bumble-bee-climate/tree/main/output) directory. 


## Making Figures
After you have your model output saved you can then go to the [`make_figures`](https://github.com/Hanna-Jackson/bumble-bee-climate/tree/main/make_figures) folder and make your plots! Make sure that before you run the figure code you tell the `pdf()` funtion the file path so it knows where to save the figure to. 

**Main Text figures:**
<br>
**Figure 1** is made by: [`Plot_Main_Effects.R`](https://github.com/Hanna-Jackson/bumble-bee-climate/blob/main/make_figures/Plot_Main_Effects.R)
<br>
**Figure 2** is made by: (appologies, code not yet cleaned and uploaded) 
<br>
**Table 1** is made by: [`Plot_Species_Table_Heatplot.R`](https://github.com/Hanna-Jackson/bumble-bee-climate/blob/main/make_figures/Plot_Species_Table_Heatplot.R)

**Supplementary Figures:**
<br>
Figures showing the change in environmental variables through space and time is in [`Plot_Climate_Maps.R`](https://github.com/Hanna-Jackson/bumble-bee-climate/blob/main/make_figures/Plot_Climate_Maps.R)
<br>
Figures showing the same information as in Figure 1 but at multiple spatial resolutions plotted together is in [`Plot_Main_Effects_Resolutions.R`](https://github.com/Hanna-Jackson/bumble-bee-climate/blob/main/make_figures/Plot_Main_Effects_Resolutions.R)
<br>
Figures showing the distribution of the records, the number of visits and the number of species across space and time is in [`Plot_Records_Maps.R`](https://github.com/Hanna-Jackson/bumble-bee-climate/blob/main/make_figures/Plot_Records_Maps.R)

