# A short description
The R program developed in this study supports the estimation of hyperparameters for the negative power transformation of the generalized extreme value (NPT-GEV) distribution. 
It implements multiple hyperparameter estimation techniques, including resampling (Resam), adding small noise (AddNoise), $k$-fold cross-validation (kFCV), and bootstrap sampling of kFCV (Boots), ensuring stable and accurate parameter estimation. 
The program integrates both maximum likelihood estimation (MLE) and L-moments estimation (LME) methods, optimizing hyperparameters via goodness-of-fit tests using the Cramér–von Mises (CvM) statistic. 
Real-world data applications validate the program’s effectiveness in extreme value modeling, particularly for inter-amount time (IAT) analysis in rainfall events. 
The program utilizes core R packages such as 'stats' and 'ismev' for optimization, statistical modeling, and visualization.

# R program
It requires three main files for execution:
1) main.R – The primary script that runs the entire analysis, including data loading, hyperparameter estimation, parameter fitting using Maximum Likelihood Estimation (MLE) and L-moments Estimation (LME), and return level computations.
2) RequiredFunction.R – A supplementary script containing essential functions for hyperparameter estimation, goodness-of-fit tests such as the CvM statistic across different techniques, and parameter estimation using MLE and L-moments (LM) methods.
3) IAT Data – The dataset containing time intervals between consecutive rainfall events exceeding predefined thresholds.

The dataset can be modified based on the research objective. In this study, we focus on flooding situation prediction by analyzing minimum IAT data, which represents the shortest time interval between heavy rainfall occurrences. The NPT-GEV model is applied to estimate extreme value parameters, providing insights into rainfall patterns and return levels critical for flood risk assessment. The program is structured to allow flexibility for researchers to apply it to other fields of extreme value analysis by adjusting the dataset accordingly.
