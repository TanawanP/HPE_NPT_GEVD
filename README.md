The R program developed in this study supports the estimation of hyperparameters for the negative power transformation of the generalized extreme value (NPT-GEV) distribution. 
It implements multiple hyperparameter estimation techniques, including resampling (Resam), adding small noise (AddNoise), $k$-fold cross-validation (kFCV), and bootstrap sampling of kFCV (Boots), ensuring stable and accurate parameter estimation. 
The program integrates both maximum likelihood estimation (MLE) and L-moments estimation (LME) methods, optimizing hyperparameters via goodness-of-fit tests using the Cramér–von Mises statistic. 
Real-world data applications validate the program’s effectiveness in extreme value modeling, particularly for inter-amount time (IAT) analysis in rainfall events. 
The program utilizes core R packages such as 'stats' and 'ismev' for optimization, statistical modeling, and visualization.
