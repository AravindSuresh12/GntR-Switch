##

include("Include.jl")
##
function objective_function(parameter_guess_array,time_start,time_step_size,time_stop,model_data_dictionary, exp_data_dictionary)  

    # what is the host_type?
    host_type = :cell_free

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





    model_data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary

	#print(binding_parameter_dictionary)
    # time constant modifier -

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
		0.0	;	# 1	GntR
		0.0	;	# 2	Venus
		0.0	;	# 3	sigma_70
        0.0 ;   # 4 Venus
        0.0 ;   # 5 BFP
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
     
    (TSIM,XSIM)= SolveBalances(time_start, time_stop, time_step_size, model_data_dictionary)

    # ===================================================================================================== #
	#print(TSIM, XSIM)

    # Phase 3:  compute simulation error ================================================================== #
    # # compute the error - we need to do a bunch of interpolation -

	 tsim_exp_protein = exp_data_dictionary["prot_data_array"][:,1]

    # # Venus mRNA -
    # itp_Venus_mRNA =  Interpolations.LinearInterpolation(TSIM, (1000)*XSIM[:,5]);
    # mRNA_Venus_sim = itp_Venus_mRNA[tsim_exp_protein]  # convert to muM from nM

	# tsim_exp_mRNA = exp_data_dictionary["mRNA_data_array"][:,1]

	# # GntR mRNA -
	# itp_GntR_mRNA =  Interpolations.LinearInterpolation(TSIM, (1000)*XSIM[:,4]);
	# # mRNA_GntR_sim = itp_GntR_mRNA[tsim_exp_mRNA]  # convert to muM from nM
    # mRNA_GntR_sim = itp_GntR_mRNA[tsim_exp_protein]  # convert to muM from nM

    #Doing simulations for all 5 protein formations

    # GntR

    itp_GntR_protein =  Interpolations.LinearInterpolation(TSIM, XSIM[:,11]);
    protein_GntR_sim = itp_GntR_protein[tsim_exp_protein] #also an error point

    #S28

    itp_S28_protein =  Interpolations.LinearInterpolation(TSIM, XSIM[:,12]);
    protein_S28_sim = itp_S28_protein[tsim_exp_protein]

    #AS28

    itp_AS28_protein =  Interpolations.LinearInterpolation(TSIM, XSIM[:,13]);
    protein_AS28_sim = itp_AS28_protein[tsim_exp_protein]

    #Venus -

    itp_Venus_protein =  Interpolations.LinearInterpolation(TSIM, XSIM[:,14]);
    protein_Venus_sim = itp_Venus_protein[tsim_exp_protein]

    # # BFP

    itp_BFP_protein =  Interpolations.LinearInterpolation(TSIM, XSIM[:,15]);
    protein_BFP_sim = itp_BFP_protein[tsim_exp_protein]

    # # get experimental data Venus-
    # mRNA_Venus_exp = exp_data_dictionary["mRNA_data_array"][:,2]          # mean is col 2 nM
	# mRNA_Venus_std_exp = exp_data_dictionary["mRNA_data_array"][:,3]      # stdev is col 3 nM

    #only for BFP and Venus do we have experimental data

    protein_Venus_exp = exp_data_dictionary["prot_data_array"][:,2]       # mean is col 2 muM
	protein_Venus_std_exp = exp_data_dictionary["prot_data_array"][:,3]   # stdev is col 3 muM
 
    protein_BFP_exp = exp_data_dictionary["prot_data_array"][:,4]       # mean is col 2 muM
	protein_BFP_std_exp = exp_data_dictionary["prot_data_array"][:,5]   # stdev is col 3 muM
 
	# # get experimental data GntR-
    # mRNA_GntR_exp = exp_data_dictionary["mRNA_data_array"][[1,2,4,5],4]  # ABHI: ignore 4h point
	# mRNA_GntR_std_exp = exp_data_dictionary["mRNA_data_array"][:,5]      # stdev is col 3 nM

    # tsim_exp_mRNA = exp_data_dictionary["mRNA_data_array"][[1,2,4,5],1]
	# A = DataInterpolations.CubicSpline(mRNA_GntR_exp,tsim_exp_mRNA)
	# mRNA_GntR_exp = A.(tsim_exp_protein)


	# mRNA_Venus_exp = exp_data_dictionary["mRNA_data_array"][[1,2,4,5],2]
	# tsim_exp_mRNA = exp_data_dictionary["mRNA_data_array"][[1,2,4,5],1]
	# B = DataInterpolations.CubicSpline(mRNA_Venus_exp,tsim_exp_mRNA)
	# mRNA_Venus_exp = B.(tsim_exp_protein)

    #  compute error terms -
    error_term_array = zeros(5,1) 
    # # error_term_array = zeros(4,1)

    # # mRNA Venus -
    # # mRNA_Venus_std_exp[1] = 1.0   # we have 0 ic
    # # tmp_arr = 1.0./((mRNA_Venus_std_exp).^2)
    # # W_mRNA = transpose([0.001,100,1,1000,100])
    # # error_vector_1 = W_mRNA * (mRNA_Venus_exp .- mRNA_Venus_sim)
	
    # # error_vector_1a = ((mRNA_Venus_exp.-minimum(mRNA_Venus_exp))./(maximum(mRNA_Venus_exp).-minimum(mRNA_Venus_exp))) .- ((mRNA_Venus_sim.-minimum(mRNA_Venus_sim))./(maximum(mRNA_Venus_sim).-minimum(mRNA_Venus_sim)))
    # # error_vector_1b = (maximum(mRNA_Venus_exp).-maximum(mRNA_Venus_sim))./(maximum(mRNA_Venus_exp))
    # # error_term_array[1] = 100*(transpose(error_vector_1a)*error_vector_1a + transpose(error_vector_1b)*error_vector_1b)
    # error_vector_1 = mRNA_Venus_exp .- mRNA_Venus_sim
    # error_term_array[1] = (transpose(error_vector_1)*error_vector_1)


    # protein Venus -
    protein_Venus_std_exp[1] = 1.0
    tmp_arr1 = 1.0./((protein_Venus_std_exp).^2)
    W_prot1 = diagm(tmp_arr1)

    #protein BFP- THINK ABOUT THESE TWO in detail
    protein_BFP_std_exp[1] = 1.0
    tmp_arr2 = 1.0./((protein_BFP_std_exp).^2)
    W_prot2 = diagm(tmp_arr2)


    error_vector_1= protein_Venus_exp .- protein_GntR_sim
    error_vector_2= protein_Venus_exp .- protein_S28_sim
    error_vector_3= protein_Venus_exp .- protein_AS28_sim
    error_vector_4 = protein_Venus_exp .- protein_Venus_sim
    error_vector_5 = protein_BFP_exp .- protein_BFP_sim



   
    # error_term_array[2] = transpose(error_vector_2a)*error_vector_2a + transpose(error_vector_2b)*error_vector_2b this includes mRNA also- for now im not doin it- the product takes care of the square 
    error_term_array[1] = transpose(error_vector_1)*error_vector_1 #is this even allowed lol
    error_term_array[2] = transpose(error_vector_2)*error_vector_2
    error_term_array[3] = transpose(error_vector_3)*error_vector_3
    error_term_array[4] = transpose(error_vector_4)*error_vector_4
    error_term_array[5] = transpose(error_vector_5)*error_vector_5

	
    # # #create fake GntR protein data = venus protein data
    # # # Venus protein -
	# # tsim_exp_protein = exp_data_dictionary["prot_data_array"][:,1]

    # # itp_GntR_protein =  Interpolations.LinearInterpolation(TSIM, XSIM[:,7]);
    # # protein_GntR_sim = itp_GntR_protein[tsim_exp_protein]

    # # error_vector_4= protein_Venus_exp .- protein_GntR_sim
    # # error_term_array[4] = (transpose(error_vector_4)*error_vector_4)


    # # ===================================================================================================== #
	#error_total = sum(error_term_array)

    # # return -
    return error_term_array
end

# # Evaluates the objective function values -
function local_refinement_step(path_to_data_dir, parameter_array; sigma=0.05, iteration_max=100)

    #inner functions - computes a value total for each error term 
    function _compute_error_total(objective_array,W)
        value = transpose(objective_array)*W*objective_array
        return value[1] #this should always be a float value, irrespective of the fact whether it is an array or not
    end

    # initialize -
    
    l= length(parameter_array)  #should be 41- parameter array is pvec initial
    number_of_parameters=l #should be 41, should be 41
    BIG = 1e10

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 16.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"

#     # wght array -  #what does this array do?
    W = diagm(ones(5)) #this is needed for making it one by one


#     # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

#     # Set some values for the weights, gene lengths etc -
#     # model_data_dictionary = customize_data_dictionary(default_data_dictionary,host_type)
	model_data_dictionary = default_data_dictionary

#     # load the experimental data -
    exp_data_dictionary = load_experimental_data_dictionary(path_to_data_dir)

#     # setup the functions -
    OF(P) = objective_function(P,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)

    # calculate the starting error -
    parameter_array_best = parameter_array
    error_array = BIG*ones(4) #for now, it should be 3+1, you are creating something new
    error_array[1]=_compute_error_total(OF(parameter_array_best), W)
    
    # main refinement loop -
    iteration_counter = 1
    while (iteration_counter<iteration_max)

#         # take a step up -
        parameter_up = parameter_array_best.*(1 .+ sigma*rand(number_of_parameters))  
        parameter_up = check_parameter_bounds(parameter_up)

        # take a step down -
        parameter_down = parameter_array_best.*(1 .- sigma*rand(number_of_parameters))
        parameter_down = check_parameter_bounds(parameter_down)

        # Evaluate the obj function -
        error_array[2] = _compute_error_total(OF(parameter_up),W)
        error_array[3] = _compute_error_total(OF(parameter_down),W)

#         # Calculate a correction factor -
        a = error_array[2]+error_array[3] -2.0* error_array[1] 
        parameter_corrected = parameter_array_best
        if (a>0.0)
            amda = -0.5*(error_array[3]- error_array[2])/a 
            parameter_corrected = parameter_array_best .+ amda*rand(number_of_parameters)
            parameter_corrected = check_parameter_bounds(parameter_corrected)
            error_array[4] = _compute_error_total(OF(parameter_corrected), W) #so far okay
        end

        # Which step has the min error?
        min_index = argmin(error_array)
        if (min_index == 1)
            parameter_array_best = parameter_array_best
        elseif (min_index == 2)
            parameter_array_best = parameter_up
        elseif (min_index == 3)
            parameter_array_best = parameter_down
        elseif (min_index == 4)
            parameter_array_best = parameter_corrected
        end

        # Update the local error
        error_array[1] = error_array[min_index]

        @show iteration_counter,error_array[min_index]

        # update local counter -
        iteration_counter = iteration_counter + 1
    end

    return parameter_array_best
end

function check_parameter_bounds(parameter_array)

    # setup paramter bounds - User is supposed to give Parameter_bounds array
    pvec_bounds = [

        # dG's -
        0.01  5.0    	;   # 1     W_RNAP_P70
        0.01  5.0    	;   # 2     W_RNAP_P28
        0.01  5.0    	;   # 3     W_RNAP_mP70
        0.01  5.0       ;   # 4     W_S70_RNAP_P70
        0.01  5.0       ;   # 5     W_S70_RNAP_mP70
        0.01  5.0       ;   # 6     W_S28_RNAP_P28
        0.01  5.0       ;   # 7     W_GntR_mP70
        0.01  5.0       ;   #8 new  W_GntR_gluconate_protein
        0.01  5.0       ;   # 8     W_AS28_S28_P28
       
       



    

        # binding parameters -
        1 10.0            ;   # 9     n_S70_RNAP_GntR
        0.001 100.0         ;   # 10     K_S70_RNAP_GntR
		1 10.0            ;   # 11     n_S70_RNAP_S28
        0.001 100.0         ;   # 12     K_S70_RNAP_S28
        1 10.0            ;   # 13     n_S70_RNAP_AS28
        0.001 100.0         ;   # 14     K_S70_RNAP_AS28
        1 10.0            ;   # 15     n_S70_RNAP_Venus
        0.001 100.0         ;   # 16     K_S70_RNAP_Venus        
		1 10.0            ;   # 17    n_GntR_mP70_AS28
        0.001 100.0         ;   # 18    K_GntR_mP70_AS28
        1 10.0            ;   # 19    n_GntR_mP70_Venus
        0.001 100.0         ;   # 20    K_GntR_mP70_Venus
        1 10.0            ;   # 21     n_S28_RNAP_BFP
        0.001 100.0         ;   # 22     K_S28_RNAP_BFP
        1 10.0            ;   # 23    n_AS28_S28_BFP
        0.001 100.0         ;   # 24    K_AS28_S28_BFP


        # time constants -
		0.001 100.0         ;	# 25	    mRNA_GntR
		0.001 100.0         ;	# 26	    mRNA_S28
        0.001 100.0         ;	# 27	    mRNA_AS28
		0.001 100.0         ;	# 28	    mRNA_Venus
        0.001 100.0         ;	# 29	    mRNA_BFP
		0.001 100.0         ;	# 30	    protein_GntR
		0.001 100.0         ;	# 31	    protein_S28
        0.001 100.0         ;	# 32	    protein_AS28
		0.001 100.0         ;	# 33	    protein_Venus
        0.001 100.0         ;	# 34	    protein_BFP
       
        # degradation mods -

        0.001 100.0         ;	# 35	    mRNA_GntR
		0.001 100.0         ;	# 36	    mRNA_S28
        0.001 100.0         ;	# 37	    mRNA_AS28
		0.001 100.0         ;	# 38	    mRNA_Venus
        0.001 100.0         ;	# 39	    mRNA_BFP
		0.001 100.0         ;	# 40	    protein_GntR
		0.001 100.0         ;	# 41	    protein_S28
        0.001 100.0         ;	# 42	    protein_AS28
		0.001 100.0         ;	# 43	    protein_Venus
        0.001 100.0         ;	# 44	    protein_BFP
	

         # w -
        6.0 12.0           ;   # 45    translation capacity half-life

        # KL value -
        10.0 1000.0         ;   # 46    KL in muM- dilution_degradation_matrix constant

		# n_gluconate_GntR-
		1.0 5.0					; # 47

		# K_gluconate_GntR-
		0.1 100.0				; # 48 mM

        # transcription capacity terms-
        0.0 3.0; #49 decay (hours)
        0.1 0.9; # 50 slope

        #translation capacity terms-
        0.0 3.0; # 51 decay (hours)
        0.001 0.9; #52 slope
    ];
	
    #parameter_array = pvec_bounds

    # check bounds -
  
    number_of_parameters = length(pvec_initial) 

    for parameter_index = 1:number_of_parameters

        # what is the parameter value?
        p_i = pvec_initial[parameter_index]

        # is p_i outside of the bounds?
        lb_value = pvec_bounds[parameter_index,1]
        ub_value = pvec_bounds[parameter_index,2]

        if (p_i<lb_value)
            pvec_initial[parameter_index,1] = lb_value
        end

        if (p_i>ub_value)
            pvec_initial[parameter_index,1] = ub_value
        end
    end

    # return -
    return parameter_array
end

function neighbor_function(parameter_array; sigma=0.05)

    # setup -
    number_of_parameters = length(parameter_array)

    # calculate new parameter array -
    new_parameter_array = parameter_array.*(1 .+ sigma*randn(number_of_parameters))

    # check the bounds and return -
    return check_parameter_bounds(new_parameter_array)
end

function cooling_function(temperature)

  # define my new temperature -
  alpha = 0.9
  return alpha*temperature
end

function acceptance_probability_function(rank_array,temperature)
    return (exp(-rank_array[end]/temperature))
end

function main(path_to_data_dir::String, initial_parameter_array::Array{Float64,1}; rank_cutoff::Int64=4, maximum_number_of_iterations::Int64=100) 

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 16.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free
    number_of_parameters = length(initial_parameter_array)

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"

    # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # Set some values for the weights, gene lengths etc -
   
	model_data_dictionary = deepcopy(default_data_dictionary)
    # load the experimental data -
    exp_data_dictionary = load_experimental_data_dictionary(path_to_data_dir)

    # setup the functions -
    OF(P) = objective_function(P,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)
    NF(P) = neighbor_function(P;sigma=0.01)

    # make call to POETs -
    (EC,PC,RA) = estimate_ensemble(OF,NF,acceptance_probability_function,cooling_function,initial_parameter_array;rank_cutoff=rank_cutoff,maximum_number_of_iterations=maximum_number_of_iterations)

    # return -
    return (EC,PC,RA)
end

##

 # setup initial condition vector - #I'll have to actually search for values online to set this up properly... For now, just see if the code is running- actual parameters to optimize would be numbered differently later
pvec_initial = [

	# dG's - Kj/mol
	2.5    	;   # 1     W_RNAP_P70- From Abhi Work of GntR Model- 
	1.5    	;   # 2     W_RNAP_P28 
    -2.5    ;   # 3     W_RNAP_mP70 
    -2.5    ;   # 4     W_S70_RNAP_P70 
    -2.0    ;   # 5     W_S70_RNAP_mP70
    -2.25   ;   # 6     W_S28_RNAP_P28
    -1.0    ;   # 7     W_GntR_mP70 - from Abhi Work of GntR model
    -2.5    ;   # 8     new W_GntR_gluconate_protein. #originally was -2.5
	-1.0    ;   # 8     W_AS28_S28_P28 -directly responsible for BFP

	# binding parameters -
	1            ;   # 9     n_S70_RNAP_GntR
	30         ;   # 10     K_S70_RNAP_GntR
	1            ;   # 11     n_S70_RNAP_S28
	10         ;   # 12     K_70_RNAP_S28
    1            ;   # 13     n_S70_RNAP_AS28
	30         ;   # 14     K_S70_RNAP_AS28
	1            ;   # 15     n_S70_RNAP_Venus
	10         ;   # 16     K_S70_RNAP_Venus
	1            ;   # 17     n_GntR_mP70_AS28
	1         ;   # 18    K_GntR_mP70_AS28
    1            ;  # 19     n_GntR_mP70_Venus
	1         ;   # 20    K_GntR_mP70_Venus
    2.0       ;   #21 n_S28_RNAP_BFP
    2.0      ;   #22 K_S28_RNAP_BFP
    2.0       ;   #23 n_AS28_S28_BFP
    2.0      ;   #24 K_AS28_S28_BFP

	# time constants - 
	1.0         ;	#25	        #mRNA_GntR
	1.0         ;	# 26	    mRNA_S28
    1.0         ;	# 27	    mRNA_AS28
	1.0         ;	# 28	    mRNA_Venus
    1.0         ;   # 29        mRNA_BFP

	15          ;	# 30	    protein_GntR 
	1.0         ;	# 31	    protein_S28
    1.0         ;   # 32	    protein_AS28
    1.0         ;	# 33	    protein_Venus
    2.0         ;   # 34	    protein_BFP
    

	# degradation mods -
	1.0         ;	# 35	    mRNA_GntR
	1.0         ;	# 36	    mRNA_S28
    1.0         ;	# 37	    mRNA_AS28
	1.0         ;	# 38	    mRNA_Venus
    1.0         ;   # 39        mRNA_BFP

    1.0          ;	# 40	    protein_GntR 
	1.0         ;	# 41	    protein_S28
    1.0         ;   # 42	    protein_AS28
    1.0         ;	# 43	    protein_Venus
    0.6         ;   # 44	    protein_BFP


	 # w -
	 10.0           ;   # 45    translation capacity half-life

	# KL value -
	100.0         ;   # 46    KL in muM

	# n_gluconate_GntR-
	1.0					; # 47

	# K_gluconate_GntR-
	1.0					; #48 mM

    # transcription capacity terms-
    0.0; # 49 decay (hours)
    0.6; # 50 slope

    # translation capacity terms-
    0.0; #51 decay (hours)
    0.1; #52 slope


];

# # pvec_initial = [
# #
# # 	# dG's -
# # 	4.4    	;   # 1     W_GntR_RNAP
# # 	-5.0    ;   # 2     W_GntR_sigma_70
# # 	3.3    	;   # 3     W_Venus_RNAP
# # 	-0.98    ;   # 4     W_Venus_sigma_70
# # 	-0.021    	;   # 5     W_Venus_GntR
# #
# # 	# binding parameters -
# # 	0.5            ;   # 6     n_GntR_sigma_70
# # 	3.7         ;   # 7     K_GntR_sigma_70
# # 	2.01            ;   # 8     n_Venus_sigma_70
# # 	2.66         ;   # 9     K_Venus_sigma_70
# # 	9.027            ;   # 10     n_Venus_GntR
# # 	74.99         ;   # 11    K_Venus_GntR
# #
# # 	# time constants -
# # 	4.44         ;	# 12	    mRNA_GntR
# # 	0.91         ;	# 13	    mRNA_Venus
# # 	1.41         ;	# 14	    protein_GntR
# # 	1.96         ;	# 14	    protein_Venus
# #
# # 	# degradation mods -
# # 	0.02 	    ;	# 16	    mRNA_GntR
# # 	0.86	    ;	# 17	    mRNA_Venus
# # 	0.24 	    ;	# 18	    protein_GntR
# # 	3.08 	    ;	# 19	    protein_Venus
# # 	27.7 	    ;	# 20	    protein_sigma_70
# #
# # 	 # w -
# # 	 9.22           ;   # 21    translation capacity half-life
# #
# # 	# KL value -
# # 	228.0         ;   # 22    KL in muM
# # ];

# setup -
path_to_data_dir = "$(pwd())"
pV = neighbor_function(pvec_initial; sigma=0.25)
EC = 0
PC = 0
RA = 0

##

# execute -
number_of_trials = 8 #anything above 8+ it crashes- occasionally crashes out midway, but working on that- but because of model parameters and not because of arithmetic
for trial_index = 1:number_of_trials

    global pV
    global EC
    global PC
    global RA
 

   # do a local step - #it gets faltered over here- 
    if (mod(trial_index,2) == 0)

        # find the lowest score pV -
        sum_error_array = sum(EC,dims=1)
        best_p_index = argmin(vec(sum_error_array))
        pV_best = PC[:,best_p_index]

        # local refine -
        pV = local_refinement_step(path_to_data_dir, pV_best; iteration_max=20)
    end


    # main -
    (EC,PC,RA) = main(path_to_data_dir, vec(pV); rank_cutoff=4,maximum_number_of_iterations=20)

    # dump results to disk -
    fname = "./poets_ensemble_W_test/RA_T$(trial_index).dat"
    writedlm(fname,RA)
    fname = "./poets_ensemble_W_test/EC_T$(trial_index).dat"
    writedlm(fname,EC)
    fname = "./poets_ensemble_W_test/PC_T$(trial_index).dat"
    writedlm(fname,PC)

    @show trial_index
end


