# Low_Coherence_Regressor_SYSID

  
## Project Description
**Low coherence regressor design:** In this project, low coherence regressor design in sparse system identification is implemented and evaluated in comparison with state-of-the art sparse estimation algorithms.

## Project Structure
* **Files:** 
   * **coordinate_transformation.m:** This function implements the design of new regressor H.
   * **generate_highcorr_regrand.m:** This function generates random matrix with high mutual coherence and random sparse vector [1].  
   * **generate_observation.m:** This function generates observations by using the main regressor \Phi and sparse parameter vector. 
   * **lassoadmm.m:** Solve lasso problem via ADMM [2] .
   * **omp.m:** Contain OMP algorithm [3].
   * **run_demo.m:** The main code of the project that produce the results. 
* **Launch:** 
   * **Run run_demo.m file**
   * **See the results as plots**
   
   
## References
* **[1]:** https://se.mathworks.com/matlabcentral/fileexchange/68810-randcorr
* **[2]:** https://web.stanford.edu/~boyd/papers/admm/lasso/lasso.html
* **[3]:** https://github.com/indigits/sparse-plex
