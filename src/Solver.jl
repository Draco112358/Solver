module Solver
using MKL, Base.Threads, AMQPClient, JSON, AWS, AWSS3, DotEnv, Oxygen, HTTP, JSON3
using GZip, CodecZlib, Serialization, MAT, SparseArrays, LinearAlgebra, FLoops, DelimitedFiles
using JLD2, Printf, SpecialFunctions, Interpolations, Base.Sort, StaticArrays, FFTW
using MLUtils: unsqueeze
using Base.Sort: searchsortedfirst


include("./lib/utility.jl") # Contiene download_json_gz, get_solverInput_from_s3, ecc.
include("./lib/format_input_output_solver_functions.jl")

#SHAREDRIS INCLUDE START

include("./lib/sharedRis/Compute_Lp_Self.jl")
include("./lib/sharedRis/Song_P_improved_Ivana_strategy.jl")
include("./lib/sharedRis/Song_improved_Ivana_strategy.jl")
include("./lib/sharedRis/distfcm.jl")
include("./lib/sharedRis/find_nodes_ports_or_le.jl")
include("./lib/sharedRis/calcola_P.jl")
include("./lib/sharedRis/calcola_Lp.jl")
#SHAREDRIS INCLUDE END

# SOLVER FFT INCLUDE START
include("./lib/solverFFT/is_material_conductor.jl")
include("./lib/solverFFT/From_3D_to_1D.jl")
include("./lib/solverFFT/bin_search.jl")
include("./lib/solverFFT/create_volumes_mapping.jl")
include("./lib/solverFFT/create_volume_centers.jl")
include("./lib/solverFFT/create_Grids_externals.jl")
include("./lib/solverFFT/create_nodes_ref.jl")
include("./lib/solverFFT/create_mapping_Gamma_no_rep.jl")
include("./lib/solverFFT/create_mapping_Ax_v2.jl")
include("./lib/solverFFT/create_mapping_Ay_v2.jl")
include("./lib/solverFFT/create_mapping_Az_v2.jl")
include("./lib/solverFFT/create_A_mats_and_find_borders_with_map_Zs.jl")
include("./lib/solverFFT/compute_Zs_with_indices.jl")
include("./lib/solverFFT/create_expansion_ind_Lp_x_grids_v2.jl")
include("./lib/solverFFT/create_expansion_ind_Lp_y_grids_v2.jl")
include("./lib/solverFFT/create_expansion_ind_Lp_z_grids_v2.jl")
include("./lib/solverFFT/compute_diagonals.jl")
include("./lib/solverFFT/compute_row_P_sup.jl")
include("./lib/solverFFT/find_voxels_port_pp.jl")
include("./lib/solverFFT/find_voxels_le_pp.jl")
include("./lib/solverFFT/build_center_P_Voxels.jl")
include("./lib/solverFFT/build_centers_Lp_with_air.jl")
include("./lib/solverFFT/build_center_air_Voxels.jl")
include("./lib/solverFFT/mesher_FFT.jl")
include("./lib/solverFFT/compute_Lp_Voxels.jl")
include("./lib/solverFFT/compute_P_vox_Rcc.jl")
include("./lib/solverFFT/compute_Voxels_Rcc.jl")
include("./lib/solverFFT/compute_Lp_Voxels_QS.jl")
include("./lib/solverFFT/compute_rows_Rcc_P.jl")
include("./lib/solverFFT/compute_rows_Rcc_Lp.jl")
include("./lib/solverFFT/compute_FFT_mutual_coupling_mats.jl")
include("./lib/solverFFT/build_row_vox_Rcc.jl")
include("./lib/solverFFT/compute_Circulant_Lp_Rcc.jl")
include("./lib/solverFFT/compute_Circulant_P_sup_Rcc.jl")
include("./lib/solverFFT/build_Yle_FFT.jl")
include("./lib/solverFFT/compute_Z_self.jl")
include("./lib/solverFFT/compute_Matrix_vector_fft.jl")
include("./lib/solverFFT/gmres_custom.jl")
include("./lib/solverFFT/FFT_solver_QS_S_type.jl")
include("./lib/solverFFT/do_solving_fft.jl")
# SOLVER FFT INCLUDE END

#SOLVER RIS INCLUDE START
include("./lib/solverRis/build_Yle_S.jl")
include("./lib/solverRis/compute_matrix_vector.jl")
include("./lib/solverRis/gmres_custom.jl")
include("./lib/solverRis/iter_solver_QS_S_type.jl")
include("./lib/solverRis/do_solving_ris.jl")
#SOLVER RIS INCLUDE END

# SOLVER RSI CAMPI INCLUDE START
include("./lib/solverRisCampi/get_signal.jl")
include("./lib/solverRisCampi/build_trapezoidal_pulse.jl")
include("./lib/solverRisCampi/genera_segnale_esponenziale.jl")
include("./lib/solverRisCampi/genera_segnale_Gaussiano_modulato.jl")
include("./lib/solverRisCampi/genera_segnale_sinusoidale.jl")
include("./lib/solverRisCampi/crea_freqs.jl")
include("./lib/solverRisCampi/compute_fields_components.jl")
include("./lib/solverRisCampi/computeVs.jl")
include("./lib/solverRisCampi/fft_UAq.jl")
include("./lib/solverRisCampi/complex_matrix_to_float_array_matrix.jl")
include("./lib/solverRisCampi/genera_punti_circonferenza.jl")
include("./lib/solverRisCampi/get_punti_oss_3D.jl")
include("./lib/solverRisCampi/build_Yle.jl")
include("./lib/solverRisCampi/compute_Matrix_vector.jl")
include("./lib/solverRisCampi/gmres_custom.jl")
include("./lib/solverRisCampi/compute_Ec_Gauss.jl")
include("./lib/solverRisCampi/compute_Ar_Gauss.jl")
include("./lib/solverRisCampi/compute_lambda_numeric.jl")
include("./lib/solverRisCampi/compute_E_field_Gauss.jl")
include("./lib/solverRisCampi/compute_H_field_Gauss.jl")
include("./lib/solverRisCampi/iter_solver_E_Gaussian_Is_type.jl")
include("./lib/solverRisCampi/do_solving_electric_fields.jl")
# SOLVER RSI CAMPI INCLUDE END
include("solver_start.jl")
export julia_main
end 
