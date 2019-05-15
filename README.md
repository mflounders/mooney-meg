# Mooney MEG Project, Overview
Publication: Neural dynamics of visual ambiguity resolution by perceptual prior, [link](https://elifesciences.org/articles/41861)
https://elifesciences.org/articles/41861 
Author: Matt Flounders [[email](mailto:flounders.matt@gmail.com)]

## Directory Structure 
______________________

### analysis
Path: /data/disk3/Matt/analysis
- preprocessed data for all analyses
- subdirectories for behavioral, commonality, decoding, TGM, behavioral indicies, rsa
- analysis README
-------------------------------------

### eprime_task_data
Path: /data/disk3/Matt/eprime_task_data
- practice, Mooney paradigm eprime task files
- Mooney paradigm eprime task files
- Mooney paradigm eprime and eprime data files, in subject specific directories
- all images used for Mooney paradigm (must be in same directory as eprime task to run)
-------------------------------------

### manuscript
Path: /data/disk3/Matt/manuscript
- subdirectories for behavioral, commonality, decoding, meg-fmri, rsa, tgm figures and pngs used for manuscript submission
- symlink (aka soft link) to eLife submitted zip files, containing scripts and data for Fig. 4-6
    - eLife_Final
    - eLife_Finalreduced
- manuscript README, links to dropbox for FINALIZED submission files (text and figures)
-------------------------------------

### misc
Path: /data/disk3/Matt/misc
- contains miscellaneous script files or shared items
- two subdirectories:
    - /collabs subdirectory for loose collaboration scripts with other lab members
    - /FromIsilon_Nov5 subdirectory contains recovered files that were lost from a random deletion of items located in my data directory on Nov 5 2016. Primiarily a back up reference point, files were redistributed to appropriate folders after recovery
-------------------------------------

### rawdata
Path: /data/disk3/Matt/rawdata
- contains all neuroimaging raw data, behavioral data, and eyetracking data collected for the Mooney MEG project
- subdirectories:
    - /[date] corresponds to MEG data collected on that date, for mapping to subject directory
    - /3TB_2015-2016 corresponds to all 3T fMRI data collected for potential source modeling
    - /behavioral_data 
    - /eyetracking_data
    - /rejecteddata correspond to rejected MEG data based on movement or incomplete datasets
- for mapping MEG data to subjects, all preprocessing scripts rely on a "naming.mat" file, located in /data/disk3/matt/scripts/utilities. This is common (when appropriate) throughout most scripts to use a static, reliable subject list and references
- more information on subject mapping can be found in the spreadsheet stored here: [Mooney Project Subject Info link](https://docs.google.com/spreadsheets/d/1xewvbsa7LNacsljWXldgkFEQM2svEJODBg478IiVO9E/edit?usp=sharing)
-------------------------------------

### scripts
Path: /data/disk3/Matt/scripts
- subdirectories for behavioral, commonality, decoding, tgm, preprocessing, rsa analyses
- /utilities contains various dependencies, plotting helpers, and .mat files used in the scripts
- /test_sensor_stats reflects early assessment of event related field and univariate exploration
- scripts README.md
-------------------------------------
______________________
