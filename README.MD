# AVIR.R
![AVIR 1.0](https://img.shields.io/badge/AVIR.R-1.0-brightgreen)
[![Maintained by Huan Lab](https://img.shields.io/badge/Maintainer-Huan%20Lab-blue)](github.com/HuanLab)

<br/>
AVIR.R is a program developed to recognize computational variation among metabolic features in samples, thereby identifying problems in the data -- such as faulty peak integration or feature misalignment -- caused by metabolomic processing. 

<br/>
Check out our detailed user guide as well as video guide on [Youtube](put link here).

<br/>

## Quick Start Guide

<br/>
Once you download AVIR.R and unzip it, first install the necessary packages for it. Simply run the following code in R Studio:

<br/>`install.packages('e1071')`
<br/>`install.packages('caret')`

<br/>
Change the working directory for AVIR.R to the directory in which you have downloaded and unzipped it into.

<br/>
Run the script by clicking the "Source" button at the top right of the script. Alternatively press CTRL + SHIFT + ENTER.

<br/>

## File Formatting

<br/>
When adding your own datasets to be analyzed using AVIR.R, please split it into two .csv files, one with the height, and the other with the area, of the chromatographic feature peaks. This is as shown below, alternatively in the PeakArea_Demo.csv and PeakHeight_Demo.csv files in the repository.

<img src = "https://github.com/ewongma/AVIR_Demo/blob/main/images/FileFormattingInput.png" width = "850" >
