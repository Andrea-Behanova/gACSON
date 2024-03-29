# gACSON

gACSON is a freely available Matlab-based software, developed for visualization, segmentation, assessment, and morphology analysis of myelinated axons in 3D-EM volumes of brain tissue samples.

If you use gACSON in your research, please cite it as

[![DOI](https://zenodo.org/badge/214974720.svg)](https://zenodo.org/badge/latestdoi/214974720)

The segmentation algorithm used in gACSON is based on ACSON pipeline described in:

Abdollahzadeh, A., Belevich, I., Jokitalo, E., Tohka, J. & Sierra, A. Automated 3D Axonal Morphometry of WhiteMatter.Sci. Reports9, 6084 (2019).

This software uses several external packages as follows:

- Bio-Formats package: https://www.openmicroscopy.org/bio-formats/downloads/
- Block-matching and 4D filtering (BM4D) algorithm for image denoising: http://www.cs.tut.fi/~foi/GCF-BM3D/
- SLIC supervoxels https://github.com/fk128/SLICSupervoxels
- Accurate fast marching and skeletonization: https://www.mathworks.com/matlabcentral/fileexchange/24531-accurate-fast-marching?s_tid=prof_contriblnk

Please download the Bio-Formats and BM4D packages and place them into the gACSON directory. SLIC supervoxels and Accurate fast marching and skeletonization code is included in gACSON package, but they are distributed under a different open-source licence. Check the licence information in the respective directories.  

This version of the software has been implemented and tested in Matlab R2020b.

A manual is available [here](https://docs.google.com/presentation/d/1qeAwYf-yLEzeY5ch08ZSZvE08q1JwcFN/edit?usp=sharing&ouid=110327177656302699585&rtpof=true&sd=true) - posted separately due to the file size. Please note that many of the slides contain videos - if they are not visible please download the file.    

___________________________________________________________________________________________________________________

**Displaying an EM image with its overlaid segmentation:** 

<img src="fig/disp_seg.gif" width="585" height="375" />


**Machine learning-based semantic segmentation of myelin:**

<img src="fig/ML_myelin_seg.gif" width="585" height="375" />


**Evaluation panel:**

<img src="fig/evaluation.gif" width="585" height="375" />
