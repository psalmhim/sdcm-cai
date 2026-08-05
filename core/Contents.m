% DCMCAI - Calcium Imaging DCM
% Version 0.92 05-Apr-2026
%
% Spectral DCM for calcium imaging data, including simulation, estimation, and visualization.
%
%   add_connectivity_bump              - Add a weak oscillatory mode around f0 (Hz) by modifying connectivity.
%   build_dcm_from_sim
%   check_subject_data
%   construct_validity
%   construct_validity_allinone
%   construct_validity_calcium
%   construct_validity_calcium_WORKING - Minimal WORKING construct validity for calcium spectral DCM
%   csd_plot
%   dcm_cai_setup_extrinsic
%   dcm_cai_visualize                  - EXTRACT AND DISPLAY RESULTS
%   dcm_compare_fit
%   dcm_evaluate_fit
%   evaluate_recovery
%   extract_zebra_data                 - Load and Process Zebrafish Data
%   generate_sim_figures               - Create paper-ready figures from simulation_results.mat
%   plot_simulation_figures
%   report_peb_result                  - Peb result analysis
%   report_sensitivity_result
%   report_sensitivity_result_new
%   report_sensitivity_result_old
%   run_calcium_simulation_monte
%   run_mc_simulation
%   run_mc_simulation_4node
%   run_peb_analysis                   - Parametric Empirical Bayes (PEB) Analysis for Zebrafish DCM
%   run_prior_sensitivity              - Prior Sensitivity Analysis for DCM
%   run_sdcm_construct_validity
%   run_sim_construct_validity
%   run_sim_order
%   run_zebra_peb_order
%   run_zebra_sdcm                     - CONFIGURATION
%   run_zebra_sdcm_all
%   run_zebra_sdcm_all_
%   run_zebra_sdcm_sensitivity         - Prior Sensitivity Analysis for DCM (FIXED)
%   run_zebra_sdcm_sensitivity_new     - Prior Sensitivity Analysis for DCM (FIXED)
%   sdcm_peb_A                         - Inputs
%   sdcm_sc_fc                         - 16 and 52 is missing
%   sensitivity_gain_sweep             - .m
%   simulate_ca_dcm
%   spm_csd_calcium_mtf
%   spm_dcm_calcium_csd
%   spm_dcm_calcium_csd_
%   spm_dcm_calcium_csd_data           - Get cross-spectral density data-features using a VAR model
%   spm_dcm_calcium_param_map
%   spm_dcm_calcium_priors
%   spm_dcm_calcium_priors__
%   spm_dcm_calcium_results            - Extended visualization for Calcium Spectral DCM
%   spm_dcm_calcium_results_summary    - Extended visualization and quantitative summary for Calcium Spectra...
%   spm_dcm_calcium_show_params
%   spm_dcm_mtf_calcium
%   spm_fx_calcium
%   spm_fx_calcium__
%   spm_gx_calcium
%   spm_gx_calcium_
%   spm_gx_calcium_hill
%   spm_gx_calcium_linearized
%   study_sdcm_validate
%   viz_zebra                          - Row vector for FC Mean
%   viz_zebra_peb
%   viz_zebra_sc
%   zebra_run_sdcm
