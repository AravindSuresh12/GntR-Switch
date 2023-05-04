# ----------------------------------------------------------------------------------- #
# Copyright (c) 2019 Varnerlab
# Robert Frederick School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #

# include -
#include("Include.jl")- #for now, no need

# Script to solve the balance equations -
time_start = 0.0
time_stop = 10.0
time_step_size = 0.01

# what is the host_type?
host_type = :cell_free

# path to parameters -
path_to_biophysical_constants_file = "./CellFree.json"

# Load the data dictionary (uses the default biophysical_constants file)
data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)
#print(data_dictionary)

R = data_dictionary["R"]
T_K = data_dictionary["T_K"]
# compute W -

PC = readdlm("./poets_ensemble_W_test/PC_T5.dat")
poets_params = PC[:,1]

# Update the data dictionary
control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
control_parameter_dictionary["W_RNAP_P70"] = exp(-1*(poets_params[1]))
control_parameter_dictionary["W_RNAP_P28"] = exp(-1*(poets_params[2]))
control_parameter_dictionary["W_RNAP_mP70"] = exp(-1*(poets_params[3]))
control_parameter_dictionary["W_S70_RNAP_P70"] = exp(-1*(poets_params[4]))
control_parameter_dictionary["W_S70_RNAP_mP70"] = exp(-1*(poets_params[5]))
control_parameter_dictionary["W_S28_RNAP_P28"] = exp(-1*(poets_params[6]))
control_parameter_dictionary["W_GntR_mP70"] = exp(-1*(poets_params[7]))
control_parameter_dictionary["W_AS28_S28_P28"] = exp(-1*(poets_params[8]))

data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary

    binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
    binding_parameter_dictionary["n_S70_RNAP_GntR"]=poets_params[9]
	binding_parameter_dictionary["K_S70_RNAP_GntR"]=poets_params[10]
	binding_parameter_dictionary["n_S70_RNAP_S28"]=poets_params[11]
	binding_parameter_dictionary["K_S70_RNAP_S28"]=poets_params[12]
	binding_parameter_dictionary["n_S70_RNAP_AS28"]=poets_params[13]
    binding_parameter_dictionary["K_S70_RNAP_AS28"]=poets_params[14]
	binding_parameter_dictionary["n_S70_RNAP_Venus"]=poets_params[15]
    binding_parameter_dictionary["K_S70_RNAP_Venus"]=poets_params[16]
	binding_parameter_dictionary["n_GntR_mP70_AS28"]=poets_params[17]
	binding_parameter_dictionary["K_GntR_mP70_AS28"]=poets_params[18]
	binding_parameter_dictionary["n_GntR_mP70_Venus"]=poets_params[19]
	binding_parameter_dictionary["K_GntR_mP70_Venus"]=poets_params[20]
	binding_parameter_dictionary["n_S28_RNAP_BFP"]=poets_params[21]
	binding_parameter_dictionary["K_S28_RNAP_BFP"]=poets_params[22]
	binding_parameter_dictionary["n_AS28_S28_BFP"]=poets_params[23]
	binding_parameter_dictionary["K_AS28_S28_BFP"]=poets_params[24]

data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary

time_constant_modifier_array = [
    0.0							;	# 1	GntR
    0.0							;	# 2	S28
    0.0							;	# 3	AS28
    0.0                         ;   # 4 Venus
    0.0                         ;   # 5 BFP
    poets_params[25]            ;  #  6 mRNA GntR
    poets_params[26]            ;  #  7 mRNA S28
    poets_params[27]            ;  #  8 mRNA AS28
    poets_params[28]            ;  #9 mRNA Venus
    poets_params[29]            ; #10 mRNA BFP
    poets_params[30]            ;  #11 Protein GntR
    poets_params[31]            ; #12 Protein S28
    poets_params[32]            ; #13 Protein AS28
    poets_params[33]            ; #14 Protein Venus
    poets_params[34]            ; #15 Protein BFP
]

data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

degradation_modifier_array = [
    0.0							;	# 1	GntR
    0.0							;	# 2	S28
    0.0							;	# 3	AS28
    0.0                         ;   # 4 Venus
    0.0                         ;   # 5 BFP
    poets_params[35]            ;  #  6 mRNA GntR
    poets_params[36]            ;  #  7 mRNA S28
    poets_params[37]            ;  #  8 mRNA AS28
    poets_params[38]            ;  #9 mRNA Venus
    poets_params[39]            ; #10 mRNA BFP
    poets_params[40]            ;  #11 Protein GntR
    poets_params[41]            ; #12 Protein S28
    poets_params[42]            ; #13 Protein AS28
    poets_params[43]        ;     #14 Protein Venus
    poets_params[44]            ; #15 Protein BFP
    
]

data_dictionary["degradation_modifier_array"] = degradation_modifier_array

# update the translation time -
data_dictionary["half_life_translation_capacity"] = poets_params[45]

biophysical_constants_dictionary = data_dictionary["biophysical_constants_dictionary"]
biophysical_constants_dictionary["translation_saturation_constant"] = poets_params[46]
data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary

# gluconate GntR binding parameters
gluconate_parameter_dictionary = data_dictionary["gluconate_parameter_dictionary"]
gluconate_parameter_dictionary["n_gluconate_GntR"] = poets_params[47]
gluconate_parameter_dictionary["K_gluconate_GntR"] = poets_params[48]
data_dictionary["gluconate_parameter_dictionary"] = gluconate_parameter_dictionary


species_symbol_type_array = data_dictionary["species_symbol_type_array"]
protein_coding_length_array = data_dictionary["protein_coding_length_array"]
gene_coding_length_array = data_dictionary["gene_coding_length_array"]
time_constant_modifier_array = data_dictionary["time_constant_modifier_array"]
initial_condition_array = data_dictionary["initial_condition_array"]

# # get gene IC -
idx_gene = findall(x->x==:gene,species_symbol_type_array)
gene_abundance_array = initial_condition_array[idx_gene]

# Precompute the translation parameters -
translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)
data_dictionary["translation_parameter_array"] = translation_parameter_array

# Precompute the kinetic limit of transcription -
transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)
data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array

# Dilution degrdation matrix -
dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary, species_symbol_type_array, degradation_modifier_array)
data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix

# Solve the model equations -

gluc_conc_array=[1e-5,1e-3,1e-2,1e-1,1e0,1e1]

new_array_BFP=[]
new_array_Venus=[]

for element in gluc_conc_array
    local gluconate_parameter_dictionary["gluconate_concentration"]=element   
    local (T,X) = SolveBalances(time_start,time_stop,time_step_size,data_dictionary)
    Venus_ele= X[end,14]
    BFP_ele=  X[end,15]
    push!(new_array_Venus,Venus_ele)
    push!(new_array_BFP,BFP_ele)
end

new_array_BFP #contains the simulated values of BFP_array
new_array_Venus #contains simulated values of BFP


p1= Plots.plot(log10.(gluc_conc_array),new_array_Venus, xlabel="Gluc_Conc (log(mM))", ylabel="Protein_concentration (uM)", label="Venus_simulated", legend= :topright)

# plot experimental data Venus
prot_data = CSV.read("./data/FINAL_FULL_DOSE_RESPONSE_POOLED.csv",DataFrame)
Venus = prot_data[!,"Venus_uM"]
stdev_Venus = prot_data[!,"Venus_uM_STDERR"]
BFP=prot_data[!,"BFP_uM"]
stdev_BFP = prot_data[!,"BFP_uM_STDERR"]

p2 = Plots.scatter!(log10.(gluc_conc_array),Venus,label = "Venus_Exp_Prot",legend = :topright,markercolor = "black", markersize = 3,yerror=stdev_Venus)

p3= Plots.plot!(log10.(gluc_conc_array),new_array_BFP, label="BFP_simulated", legend= :topright)

p4= Plots.scatter!(log10.(gluc_conc_array),Venus,label = "Venus_Exp_Prot",legend = :topright,markercolor = "black", markersize = 3,yerror=stdev_BFP)
