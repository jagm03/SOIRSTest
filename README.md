Introduction
This repository contains the code associated with the paper:

Title: Nonparametric Testing of Covariate Influence on Spatial Point Processes Using LISA Functions

Spatial point processes provide an efficient mathematical framework for modeling random spatially referenced events within a support space. In many scenarios, these events are influenced by environmental or other spatially varying factors, known as covariates. Incorporating covariates into point process models enhances their accuracy and predictive capabilities.

Our approach analyzes spatial point processes by considering functional marks, specifically Local Indicators of Spatial Association (LISA) functions, to describe the second-order structure of the point process. We develop a nonparametric method for testing the hypothesis of independence between a marked point process and a covariate. This method allows us to test whether a given covariate influences the second-order structure of the point process without making strong parametric assumptions.

By calculating classical correlation coefficients between the covariate values at observed points and the corresponding values of the functional marks, we can discern how the local spatial association varies with respect to the covariate. This not only detects the presence of dependence but also provides an interpretation of the observed relationship, offering valuable insights into the underlying spatial processes.

The effectiveness of this testing method is demonstrated through simulation experiments and an application to a dataset involving fish locations in a water reservoir and associated covariates.
