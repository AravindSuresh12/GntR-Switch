include("Include.jl")
# mean center -
function mean_center_array(results_array::Array{Float64,2})::Array{Float64,2}

    # get the size -
    (NR,NC) = size(results_array)
    scaled_array = zeros(NR,NC)

    for col_index = 1:NC

        data_col = results_array[:,col_index]
        mu_value = mean(data_col)
        std_value = std(data_col)

        for row_index = 1:NR
            scaled_array[row_index,col_index] = (data_col[row_index] - mu_value)/(std_value)
        end
    end

    return scaled_array
end

# computes the model performance -
function model_performance(parameter_guess_array,index)

    # what is the host_type?
    host_type = :cell_free

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 16.0
    time_step_size = 0.01

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"
    #path_to_data_dir = "$(pwd())/data"

    # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # Set some values for the weights, gene lengths etc -
    model_data_dictionary = deepcopy(default_data_dictionary)

   # Phase 1: parameter update =========================================================================== #
   # update the paramaters in the model data dictionary -
   # for now - lets only search over dG's -
   R = model_data_dictionary["R"]
   T_K = model_data_dictionary["T_K"]

   # compute W -
   tmp_W_array = Float64[]
   for index = 1:9
       parameter_guess = parameter_guess_array[index]
       value = exp(-1*(parameter_guess))
       push!(tmp_W_array,value)
   end

   # update the control W's -
   control_parameter_dictionary = model_data_dictionary["control_parameter_dictionary"]
   control_parameter_dictionary["W_RNAP_P70"] = tmp_W_array[1]
   control_parameter_dictionary["W_RNAP_P28"] = tmp_W_array[2]
   control_parameter_dictionary["W_RNAP_mP70"] = tmp_W_array[3]
   control_parameter_dictionary["W_S70_RNAP_P70"] = tmp_W_array[4]
   control_parameter_dictionary["W_S70_RNAP_mP70"] = tmp_W_array[5]
   control_parameter_dictionary["W_S28_RNAP_P28"] = tmp_W_array[6]
   control_parameter_dictionary["W_GntR_mP70"] = tmp_W_array[7]
   control_parameter_dictionary["W_GntR_gluconate_protein"] = tmp_W_array[8]
   control_parameter_dictionary["W_AS28_S28"] = tmp_W_array[9]
   model_data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary

   #print(control_parameter_dictionary)

   binding_parameter_dictionary = model_data_dictionary["binding_parameter_dictionary"]
    binding_parameter_dictionary["n_S70_RNAP_GntR"]=parameter_guess_array[10]
	binding_parameter_dictionary["K_S70_RNAP_GntR"]=parameter_guess_array[11]
	binding_parameter_dictionary["n_S70_RNAP_S28"]=parameter_guess_array[12]
	binding_parameter_dictionary["K_S70_RNAP_S28"]=parameter_guess_array[13]
	binding_parameter_dictionary["n_S70_RNAP_AS28"]=parameter_guess_array[14]
    binding_parameter_dictionary["K_S70_RNAP_AS28"]=parameter_guess_array[15]
	binding_parameter_dictionary["n_S70_RNAP_Venus"]=parameter_guess_array[16]
    binding_parameter_dictionary["K_S70_RNAP_Venus"]=parameter_guess_array[17]
	binding_parameter_dictionary["n_GntR_mP70_AS28"]=parameter_guess_array[18]
	binding_parameter_dictionary["K_GntR_mP70_AS28"]=parameter_guess_array[19]
	binding_parameter_dictionary["n_GntR_mP70_Venus"]=parameter_guess_array[20]
	binding_parameter_dictionary["K_GntR_mP70_Venus"]=parameter_guess_array[21]
	binding_parameter_dictionary["n_S28_RNAP_BFP"]=parameter_guess_array[22]
	binding_parameter_dictionary["K_S28_RNAP_BFP"]=parameter_guess_array[23]
	binding_parameter_dictionary["n_AS28_S28_BFP"]=parameter_guess_array[24]
	binding_parameter_dictionary["K_AS28_S28_BFP"]=parameter_guess_array[25]


   time_constant_modifier_array = [
    0.0							;	# 1	GntR
    0.0							;	# 2	S28
    0.0							;	# 3	AS28
    0.0                         ;   # 4 Venus 
    0.0                         ;   # 5 BFP
    parameter_guess_array[26]	;	# 6	mRNA_GntR
    parameter_guess_array[27]	;	# 7	mRNA_S28
    parameter_guess_array[28]	;	# 8	mRNA_AS28
    parameter_guess_array[29]	;	# 9	mRNA_Venus							;
    parameter_guess_array[30]   ;   # 10 mRNA_BFP
    parameter_guess_array[31]	;	# 11	protein_GntR
    parameter_guess_array[32]	;	# 12	protein_S28
    parameter_guess_array[33]	;	# 13	protein_AS28
    parameter_guess_array[34]	;	# 14	protein_Venus
    parameter_guess_array[35]	;	# 15	protein_BFP
]


   model_data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

   # setup degradation_modifier_array -

    degradation_modifier_array = [
        0.0							;	# 1	GntR
        0.0							;	# 2	S28
        0.0							;	# 3	AS28
        0.0                         ;   # 4 Venus 
        0.0                         ;   # 5 BFP
    parameter_guess_array[36]	;	# 6	mRNA_GntR
    parameter_guess_array[37]	;	# 7	mRNA_S28
    parameter_guess_array[38]   ;   # 8 mRNA_AS28
    parameter_guess_array[39]   ;   # 9 mRNA_Venus
    parameter_guess_array[40]   ;   # 10 mRNA_BFP
    parameter_guess_array[41]	;	# 11	protein_GntR
    parameter_guess_array[42]	;	# 12	protein_Venus
    parameter_guess_array[43]	;	# 13	protein_sigma_70
    parameter_guess_array[44]	;	# 14	protein_Venus
    parameter_guess_array[45]	;	# 13	protein_BFP
]

   model_data_dictionary["degradation_modifier_array"] = degradation_modifier_array

   # update the translation time -
   model_data_dictionary["half_life_translation_capacity"] = parameter_guess_array[46]

   # lastly, update KL -
   biophysical_constants_dictionary = model_data_dictionary["biophysical_constants_dictionary"]
   biophysical_constants_dictionary["translation_saturation_constant"] = parameter_guess_array[47]
   model_data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary

   # gluconate GntR binding parameters
   gluconate_parameter_dictionary = model_data_dictionary["gluconate_parameter_dictionary"]
   gluconate_parameter_dictionary["n_gluconate_GntR"] = parameter_guess_array[48]
   gluconate_parameter_dictionary["K_gluconate_GntR"] = parameter_guess_array[49]
   model_data_dictionary["gluconate_parameter_dictionary"] = gluconate_parameter_dictionary

   # update the transcription capacity parameters
   model_data_dictionary["transcription_capacity_delay"] = parameter_guess_array[50]
   model_data_dictionary["transcription_capacity_slope"] = parameter_guess_array[51]
   
   # update the translation capacity parameters
   model_data_dictionary["translation_capacity_delay"] = parameter_guess_array[52]
   model_data_dictionary["translation_capacity_slope"] = parameter_guess_array[53]


   # grab defaults -
   species_symbol_type_array = model_data_dictionary["species_symbol_type_array"]
   protein_coding_length_array = model_data_dictionary["protein_coding_length_array"]
   gene_coding_length_array = model_data_dictionary["gene_coding_length_array"]
   time_constant_modifier_array = model_data_dictionary["time_constant_modifier_array"]
   initial_condition_array = model_data_dictionary["initial_condition_array"]

   # # get gene IC -
   idx_gene = findall(x->x==:gene,species_symbol_type_array)
   gene_abundance_array = initial_condition_array[idx_gene]

   # Precompute the translation parameters -
   translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)
   model_data_dictionary["translation_parameter_array"] = translation_parameter_array

   # Precompute the kinetic limit of transcription -
   transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)
   model_data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array


   # Dilution degrdation matrix -
   dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary, species_symbol_type_array, degradation_modifier_array)
   model_data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix

   # ===================================================================================================== #
   #print(model_data_dictionary)
   # Phase 2:  solve model equations ===================================================================== #
   # solve the balance equations -
   (TSIM,XSIM) = SolveBalances(time_start,time_stop,time_step_size,model_data_dictionary)
    # ===================================================================================================== #

    # Phase 3: compute the model performance metrics ====================================================== #
    p_Model_AUC = integrate(TSIM,XSIM[:,index]) 
    # ===================================================================================================== #

    # return the performance_array -
    return p_Model_AUC
end

function main(path_to_ensemble_file::String,index)

    # setup the sensitivity function -
    SF(P) = model_performance(P,index)

    # setup ranges -
    sample_bounds_array = Array{Tuple,1}()
    ensemble_array = readdlm(path_to_ensemble_file)
    (number_of_parameters,number_of_trials) = size(ensemble_array)
    for parameter_index = 1:(number_of_parameters)

        # get row of parameters -
        parameter_row = ensemble_array[parameter_index,:]
        min_value = minimum(parameter_row)
        max_value = maximum(parameter_row)

        # create the tuple -
        tmp_tuple = (min_value,max_value)

        # cache -
        push!(sample_bounds_array,tmp_tuple)
    end

    @show index

    # do the global sensitivity analysis -
    sensitivity_results = GlobalSensitivity.gsa(SF,Morris(total_num_trajectory=10000,num_trajectory=1000),sample_bounds_array)

    # return -
    return sensitivity_results
end

# setup paths -
path_to_ensemble_file = "./poets_ensemble_W_test/PC_T8.dat"
 
# compute a sensitivity array for the AUC of each species -
species_index_array = [9 14 15] #only doing for Venus and BFP- just try the 3 species ie mRNA Venus Protein Venus and Protein BFP, because we only have experimental data for those three. 
number_of_species = length(species_index_array)
number_of_parameters = 53
results_array = zeros(number_of_parameters,1)
for species_index in species_index_array

    global results_array

    # conduct senstivity analysis -
    sensitivity_results = main(path_to_ensemble_file,species_index)

    # get the μ and σ^2
    mu = sensitivity_results.means_star
    var = sensitivity_results.variances #can also get skew

    #@show mu, var

    results_array = [results_array transpose(mu) transpose(var)]

end

results_array = results_array[:,2:end]

fname = "./sensitivity_results/Sensitivity_matrix_gluconate.dat"
writedlm(fname,results_array)
