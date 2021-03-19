# prfVista
A small respository to call pRF functions in vistasoft without the need to set up vistasoft-specific data structures and session files. The assumption is that the fMRI time series are in NIFTI files, and the stimulus file is in either a NIFTI file or MATLAB file.

This requires the vistasoft tools on the paths. The code for this repository and all of its dependencies can be added to the Matlab path using the [ToolboxToolox](https://github.com/ToolboxHub/ToolboxToolbox) configuration file: [prfVista](https://github.com/WinawerLab/ToolboxRegistry/blob/master/configurations/prfVista.json).

See [rmMain.m](https://github.com/vistalab/vistasoft/blob/master/mrBOLD/Analysis/retinotopyModel/rmMain.m) for more details.

Examples:

```
%% Solve a pRF model with one scan using a small synthetic dataset
datafiles = fullfile(prfRootPath, 'local', 'vista2test_Data.nii.gz');
stimfiles = fullfile(prfRootPath, 'local', 'vista2test_Stim.nii.gz');
stimradius = 10;
results = prfVistasoft(stimfiles, datafiles, stimradius);
figure, scatter(results.model{1}.x0, results.model{1}.y0)
```

```
%% Same as above, but for coarse fit only
datafiles = fullfile(prfRootPath, 'local', 'vista2test_Data.nii.gz');
stimfiles = fullfile(prfRootPath, 'local', 'vista2test_Stim.nii.gz');
stimradius = 10;
results = prfVistasoft(stimfiles, datafiles, stimradius, 'wsearch', 'coarse');
figure, scatter(results.model{1}.x0, results.model{1}.y0)
```

```
%% Same as above, but for css model
datafiles = fullfile(prfRootPath, 'local', 'vista2test_Data.nii.gz');
stimfiles = fullfile(prfRootPath, 'local', 'vista2test_Stim.nii.gz');
stimradius = 10;
results = prfVistasoft(stimfiles, datafiles, stimradius, 'wsearch', 'coarse', 'css');
figure, scatter(results.model{1}.x0, results.model{1}.y0)
```
