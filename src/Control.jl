# ----------------------------------------------------------------------------------- # #Note for myself- The control array values are robust 
using PyCall
np = pyimport("numpy")
# Copyright (c) 2021 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
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
#
# ----------------------------------------------------------------------------------- #
# Function: calculate_transcription_control_array
# Description: Calculate the transcriptional control array at time t
# Generated on: 2021-04-29T20:40:54.294
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => 
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Transcriptional control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_transcription_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})

	# initialize the control 
	control_array = zeros(5) #5 genes, 5 CFS

	# Alias the species - 
	GntR = x[1]
	S28 = x[2]
	AS28 = x[3]
    Venus=  x[4]
	BFP= x[5]
	mRNA_GntR = x[6]
	mRNA_S28 = x[7]
	mRNA_AS28 = x[8]
	mRNA_Venus = x[9]
	mRNA_BFP = x[10]
	protein_GntR = x[11]
	protein_S28 = x[12]
	protein_AS28 = x[13]
	protein_Venus = x[14]
	protein_BFP = x[15]	

	# Alias the binding parameters -
	binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]

	n_S70_RNAP_GntR=binding_parameter_dictionary["n_S70_RNAP_GntR"]
	K_S70_RNAP_GntR=binding_parameter_dictionary["K_S70_RNAP_GntR"]
	n_S70_RNAP_S28=binding_parameter_dictionary["n_S70_RNAP_S28"]
	K_S70_RNAP_S28=binding_parameter_dictionary["K_S70_RNAP_S28"]
	n_S70_RNAP_AS28=binding_parameter_dictionary["n_S70_RNAP_AS28"]
	K_S70_RNAP_AS28=binding_parameter_dictionary["K_S70_RNAP_AS28"]
	n_S70_RNAP_Venus=binding_parameter_dictionary["n_S70_RNAP_Venus"]
	K_S70_RNAP_Venus=binding_parameter_dictionary["K_S70_RNAP_Venus"]
	n_GntR_mP70_AS28=binding_parameter_dictionary["n_GntR_mP70_AS28"]
	K_GntR_mP70_AS28=binding_parameter_dictionary["K_GntR_mP70_AS28"]
	n_GntR_mP70_Venus=binding_parameter_dictionary["n_GntR_mP70_Venus"]
	K_GntR_mP70_Venus=binding_parameter_dictionary["K_GntR_mP70_Venus"]
	n_S28_RNAP_BFP=binding_parameter_dictionary["n_S28_RNAP_BFP"] 
	K_S28_RNAP_BFP=binding_parameter_dictionary["K_S28_RNAP_BFP"]
	n_AS28_S28_BFP=binding_parameter_dictionary["n_AS28_S28_BFP"]
	K_AS28_S28_BFP=binding_parameter_dictionary["K_AS28_S28_BFP"] 




	# Alias the control function parameters -
	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
	W_RNAP_P70 = control_parameter_dictionary["W_RNAP_P70"] 
	W_RNAP_P28=control_parameter_dictionary["W_RNAP_P28"]
	W_RNAP_mP70=control_parameter_dictionary["W_RNAP_mP70"] 
	W_S70_RNAP_P70=control_parameter_dictionary["W_S70_RNAP_P70"] 
	W_S70_RNAP_mP70=control_parameter_dictionary["W_S70_RNAP_mP70"] 
	W_S28_RNAP_P28=control_parameter_dictionary["W_S28_RNAP_P28"] 
	W_GntR_mP70=control_parameter_dictionary["W_GntR_mP70"]  #binds to mP70	
	W_AS28_S28_P28=control_parameter_dictionary["W_AS28_S28_P28"] 
	

	#Since there are two GntR receptors in which binding takes place:
	gluconate_parameter_dictionary = data_dictionary["gluconate_parameter_dictionary"]
	Gluc_conc = gluconate_parameter_dictionary["gluconate_concentration"]
	n_gluconate_GntR = gluconate_parameter_dictionary["n_gluconate_GntR"]
	K_gluconate_GntR = gluconate_parameter_dictionary["K_gluconate_GntR"]
	protein_S70= gluconate_parameter_dictionary["Protein_sigma_70"]

	A1= Gluc_conc ^ n_gluconate_GntR 
	term1A=real(Gluc_conc^n_gluconate_GntR)
	term1B=real(K_gluconate_GntR^n_gluconate_GntR) 
	B1= term1A+ term1B

	f_gluc_GntR = real(A1/B1)  #equivalent to f bound

	
	protein_GntR= (1-f_gluc_GntR)*protein_GntR #Remaining GntR unbound fraction (free from gluconate) #this protein changes but not actor set 2

	

	actor_set2= [
	protein_GntR
	]

	actor=prod(protein_GntR)



 #################################################################################################################################################
	#Control Terms

	#1- For GntR

	actors= protein_S70 #Since this is assumed to be a time invariant source, no need to take prod of it ie 3.5 e-3, because k for this is in uM 

	f_S70_RNAP_GntR= (actors ^n_S70_RNAP_GntR) / ( (actors ^ n_S70_RNAP_GntR) + (K_S70_RNAP_GntR ^ n_S70_RNAP_GntR) )  #will also be a constant number always

	#print("the f value controlling actor is : $(f_S70_RNAP_GntR) \n")
	

	num1= W_RNAP_P70 + W_S70_RNAP_P70*f_S70_RNAP_GntR

	den1= 1+ W_RNAP_P70 + W_S70_RNAP_P70*f_S70_RNAP_GntR


	control_array[1]= num1/den1  


	######################################

	#2 for S28 production- 

	f_S70_RNAP_S28= (actors ^n_S70_RNAP_S28) / ( (actors ^n_S70_RNAP_S28) + (K_S70_RNAP_S28 ^n_S70_RNAP_S28) ) 

	num2= W_RNAP_P70 + W_S70_RNAP_P70*f_S70_RNAP_S28 #because S28 gene is under the control of P70 promoter
	den2= 1+ W_RNAP_P70 + W_S70_RNAP_P70*f_S70_RNAP_S28

	control_array[2]= num2/den2

	######################################

	#3 Control function for AS28 production- since its under control of mP70-GntR

	
	#f function which talks about Gluconate binding on to GntR 

	#two stuffs bind on to mP70 containing AS28 gene- One is GntR protein, second one is S70-RNAP complex

	#the f function which talks about GntR binding with mP70 in AS28 

#assume you remove all the real/complex things

#actor changes properly 

    A=real(actor^n_GntR_mP70_AS28) 
	B=real(actor^(n_GntR_mP70_AS28)) + real((K_GntR_mP70_AS28)^n_GntR_mP70_AS28)


	f_GntR_mP70_AS28= A/B #GntR binds to the mP70 promoter which is upstream of as28 gene 

	#print("For gluconate value:   $(Gluc_conc), the f value is:   $(f_GntR_mP70_AS28) \n")

	#The f function which talks about S70-RNAP complex binding with mP70 on AS28 gene

	#print("The value of actor in this case:",actor,"\n")
	f_S70_RNAP_AS28= (actors ^n_S70_RNAP_AS28) / ( actors ^n_S70_RNAP_AS28 + K_S70_RNAP_AS28 ^n_S70_RNAP_AS28) 
	

	#print("The value of f_gluc is $(f_gluc_GntR), while the value for f_GntR_mP70_AS28 is:  $(f_GntR_mP70_AS28)   for glucose concentration $(Gluc_conc) \n")
	a1=W_RNAP_mP70 + W_S70_RNAP_mP70*f_S70_RNAP_AS28
	a2=W_GntR_mP70*(f_gluc_GntR)*(f_GntR_mP70_AS28)

	num3= (a1) + (a2)

    den3=1+(a1)+(a2)

	control_array[3]=num3/den3


	#print("The a1 array is $(a1)","\n")
	#print("The a2 array is $(100*a2)","\n")
	#print("The sum of a1 and a2 only terms come to: $(a1+100*a2)", "\n")
	#print("The control array 3 term is $(control_array[3])" ,"\n")

	######################################
	
	##4 Control function for Venus Production

	#f function which talks about Gluconate binding on to GntR 

	#the f function which talks about GntR binding with mP70 

	A4= actor ^ n_GntR_mP70_Venus
	B4= actor ^ (n_GntR_mP70_Venus) + K_GntR_mP70_Venus ^ n_GntR_mP70_Venus

    f_GntR_mP70_Venus= A4/B4 #(actor ^ n_GntR_mP70_Venus) /( actor ^ n_GntR_mP70_Venus + K_GntR_mP70_Venus ^ n_GntR_mP70_Venus ) 

	#two stuffs bind on to mP70 containing Venus gene- One is GntR protein, second one is S70-RNAP complex

	#The f function which talks about S70-RNAP complex binding with mP70 for downstream Venus gene

	f_S70_RNAP_Venus= (actors ^n_S70_RNAP_Venus) / ( actors ^ n_S70_RNAP_Venus + K_S70_RNAP_Venus ^ n_S70_RNAP_Venus) 

	

	num4= W_RNAP_mP70 + W_S70_RNAP_mP70*f_S70_RNAP_Venus + f_gluc_GntR* W_GntR_mP70*f_GntR_mP70_Venus #Here, both S70 AND GntR bind to mP70 promoter sequence 
    den4= 1 + num4



	control_array[4]= num4/den4 

	#######################################

	#5 for BFP 

	#AS28 is like an inhibitor molecule for S28- like how Gluconate is to GntR, but in reverse 

	actor_set3= [protein_AS28]'  
		
	actor3= prod(actor_set3)

	A3= real(K_AS28_S28_BFP^n_AS28_S28_BFP)
	term1=real(actor3^n_AS28_S28_BFP) #potential error
	term2= real(K_AS28_S28_BFP^n_AS28_S28_BFP)
	B3= term1+ term2

	f_AS28_S28_BFP= real(A3/B3) 

	protein_S28= (1-f_AS28_S28_BFP)*protein_AS28 
	
	actor_set4= [protein_S28]'

	actor4= prod(actor_set4)

	f_S28_RNAP_BFP =(actor4^ n_S28_RNAP_BFP) /( actor4^ n_S28_RNAP_BFP + K_S28_RNAP_BFP ^ n_S28_RNAP_BFP) #this will always be real


	num5= W_RNAP_P28 + W_S28_RNAP_P28 *f_S28_RNAP_BFP 

	f_anti_S28_RNAP_BFP=(K_S28_RNAP_BFP^ n_S28_RNAP_BFP) /( actor4^ n_S28_RNAP_BFP + K_S28_RNAP_BFP ^ n_S28_RNAP_BFP)

	
	den5= 1 + num5  + 1E-2/ ((W_AS28_S28_P28)*f_anti_S28_RNAP_BFP)  #idek if this is allowed
	#the last term in denominator refers to the inhibitory role that AS28 plays in blocking BFP production

	control_array[5]= num5/den5 

	
	#estimate transcription_capacity_delay and transcription_capacity_slope parameters
	transcription_capacity_delay = data_dictionary["transcription_capacity_delay"] 
	transcription_capacity_slope = data_dictionary["transcription_capacity_slope"] 
	# build control array with these two logistic function parameters
	f(y)= (1+exp((-transcription_capacity_delay*pi)/(transcription_capacity_slope*sqrt(3))))/(1+exp(((y-transcription_capacity_delay)*pi)/(transcription_capacity_slope*sqrt(3))))
	correction_term_manual=f(t)
	#control_array=control_array*correction_term_manual #Try for now ? originally ABHI commented it out 
	control_array=control_array

	# return -
	return control_array
end

#
# ----------------------------------------------------------------------------------- #
# Function: calculate_translation_control_array
# Description: Calculate the translation control array at time t
# Generated on: 2021-04-29T20:40:54.447
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1) (Basically its a unit vector with 9 entries of 1.0 each)
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Translation control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_translation_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})

	# initialize the control -
	control_array = ones(5)

	# #manually add txtl resources decay term for now (logistic decay)
	# f(y)=(1+exp((-6*pi)/(0.6*sqrt(3))))/(1+exp(((y-6)*pi)/(0.6*sqrt(3))))
	# correction_term_manual=f(t)
	# control_array=control_array*correction_term_manual

	# estimate transcription_capacity_delay and transcription_capacity_slope parameters #these are what he added
	translation_capacity_delay = data_dictionary["translation_capacity_delay"] 
	translation_capacity_slope = data_dictionary["translation_capacity_slope"] 
	# build control array with these two logistic function parameters
	f(y)= (1+exp((-translation_capacity_delay*pi)/(translation_capacity_slope*sqrt(3))))/(1+exp(((y-translation_capacity_delay)*pi)/(translation_capacity_slope*sqrt(3))))
	correction_term_manual=f(t)
	#control_array=control_array*correction_term_manual #THIS WAS ORIGINALLY COMMENTED OUT
	
	
	correction_term = (x[16]/100.0) 
    control_array = control_array*correction_term
	# # # return -
	return control_array
end
