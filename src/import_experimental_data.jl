using CSV
using DataFrames
using PyPlot
using Printf

function it_starts_with(str, pattern)
  return length(str) >= length(pattern) && str[1:length(pattern)]==pattern
end

function import_CVtoDataFrame_path(f_name)
    df = DataFrame(U = Float32[], I = Float32[], t = Float32[])
    line_is_valid=false
    open(f_name) do file
        for ln in eachline(file)
            if line_is_valid
                push!(df,[parse(Float32,el) for el in split(ln)])
            else
                if ln=="End Comments"
                    line_is_valid=true
                end
            end
        end
    end
    df
end

#  "22-10-2020-point_2_PROKOP.txt"

function import_EIStoDataFrame_path(f_name)
    #df = DataFrame(f = Float32[], Re = Float32[], Im = Float32[])
    if length(f_name) >= 10 && f_name[end-9:end]=="PROKOP.txt"
      # position of frequency column
      freq_col_idx = parse(Int32,f_name[end-10])
      
      df = CSV.File(f_name, decimal=',') |> DataFrame            
      return DataFrame(
        f = reverse(deepcopy(df[!, freq_col_idx])), 
        Z = reverse(deepcopy(df[!, freq_col_idx+1] .- im*df[!, freq_col_idx+2]))
      )
    end
    df = DataFrame(f = Float32[], Z = Complex[])
    line_is_valid=false
    ###
    imaginary_compensation = 1.0
    ###
    open(f_name) do file
        for ln in eachline(file)
            if line_is_valid
                splitted_row = [parse(Float32,el) for el in split(ln)]
                #push!(df, (splitted_row[1], splitted_row[5], splitted_row[6] ))
                push!(df, (splitted_row[1], splitted_row[5] +  splitted_row[6]*im*imaginary_compensation ))
            else
                if ln=="End Comments"
                    line_is_valid=true
                end
            end
        end
    end
    df.f = reverse(df.f)
    df.Z = reverse(df.Z)    
    return df
end

function get_U_I_from_IVpoint(f_name)
    line_is_valid=false
    last_line = []
    open(f_name) do file
        for ln in eachline(file)
            if line_is_valid
                last_line = [parse(Float32,el) for el in split(ln)]
            else
                if ln=="End Comments"
                    line_is_valid=true
                end
            end
        end
    end
    return (last_line[1], last_line[2])
end

function import_IVtoDataFrame_folder(;TC, pO2, bias_array, folder)
  df = DataFrame(U = Float32[], I = Float32[], t = Float32[])
  for bias in bias_array
    if bias == 0
      push!(df, (0,0,0)) 
    else
      push!(df,(get_U_I_from_IVpoint(folder*"/pol$(Int(bias*1000)).cor")..., 0) ) 
    end
  end
  return df
end

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

function import_CVtoDataFrame(;TC,pO2, data_set="MONO_110")
  pO2=Int64(pO2)
  pO2_Dan = @sprintf("%02g",pO2)
  pO2_2_digits = @sprintf("%01g",pO2)
  #
  TC=Int64(TC)
  #
  if data_set=="OLD_MONO_100"
    fNAME=string("snehurka/experimental_data_PSS/YSZ_09-2019_oxygen100/100 750to850 0to100%O2/",TC,"C/100 ",TC,"C ",pO2_Dan,"% do 1V/CV.cor")
  elseif data_set == "POLY"
    fNAME=string("snehurka/experimental_data_PSS/jako asi 6/$(TC) $(pO2_Dan) 6/cv1.cor")
  elseif data_set == "MONO_110"
    fNAME=string("snehurka/experimental_data_PSS/YSZ 110/110 $(TC) $(pO2_Dan)/cv1.cor")
  #
  #
  elseif it_starts_with(data_set, "cyc_HebbWagner")
    # cyc means ->  >> bias << has the meaning of >> number with respect to time order << of experiments
    plane = data_set[16:18]
    ordering_number = data_set[25]
    
    
#     bias_range = vcat(
#                         collect(0 : 0.05 : 1),
#                         collect(0.95 : -0.05 : -0.95),
#                         collect(-1 : 0.05 : 0)
#     )
#     
    bias_range = collect(1 : 1 : 81)
    df = DataFrame(U = Float32[], I = Float32[], t = Float32[])
    for bias in bias_range
      total_order_str = ""
      if plane == "100"
        total_order_str = "09"
      elseif plane == "110"
        total_order_str = "10"
      elseif plane == "111"
        if ordering_number in ['5', '6', '7']
          total_order_str = "12"
        else              
          total_order_str = "11"
        end
      end
      
      bias_return_value = 1.0
      cell_specification = ""
      sequence_str = ""
      if plane == "110"
        sequence_str = "druhe_mereni/"
      end    
      if plane == "111"
        if ordering_number == '1'
          sequence_str = "mereni 1/"
        elseif ordering_number == '2'
          sequence_str = "mereni 2/"
        elseif ordering_number == '3'
          sequence_str = "mereni 3 otocene kabely/"
        elseif ordering_number == '4'        
          sequence_str = "mereni 4 pm2V/"
          bias_return_value = 2.0
        elseif ordering_number == '5'
          cell_specification = " LESTENA"
        elseif ordering_number == '6'
          cell_specification = " LESTENA"
          sequence_str = "mereni 2/"
        elseif ordering_number == '7'
          cell_specification = " LESTENA" 
          sequence_str = "mereni 3 do +- 2V/"
          bias_return_value = 2.0
        end
      end
      
      bias_digit_precision = 3
      bias_epsilon = 10e-5
      
      bias_idx = Int(bias)
      bias_step = 0.05
      dir = 1
      bias_value = 0.0
      bias_str = ""
      
      return_counter = 0
      for i in 1:(bias_idx - 1)    
        bias_value = round(bias_value + dir*bias_step, digits=bias_digit_precision)
        if bias_value > bias_return_value && dir > 0
          return_counter += 1
          bias_value = round(bias_value + -2*bias_step, digits=bias_digit_precision)
          dir = -1
        end
        if bias_value < -bias_return_value && dir < 0
          return_counter += 1
          bias_value = round(bias_value + 2*bias_step, digits=bias_digit_precision)
          dir = 1
        end
        if return_counter > 1 && bias_value > 0        
          bias_value = Inf
          break
        end
      end
      bias_value = round(bias_value, digits=bias_digit_precision)
      
      
      bias_str = "$(bias_value)"
      backward_string = ""
      ocp_token = ""
      if bias_idx == 1      
        ocp_token = "vyrovnavaci OCP_Rp0$(pO2_2_digits)"
      elseif (bias_value > 0 + bias_epsilon) && (dir < 0)      
        backward_string = "backward_"
      elseif bias_value == 0.0 && dir < 0      
        bias_str = "-0.0"
        ocp_token = "OCP_negativni-sken_Rp0$(pO2_2_digits)"
      elseif bias_value < 0 - bias_epsilon && dir > 0
        backward_string = "backward_"
      elseif bias_idx > 1 && bias_value == 0.0 && dir > 0
        ocp_token = "OCP_konec_Rp0$(pO2_2_digits)"
      end       
      
      if ocp_token == ""
        file_speficitaion = "$(backward_string)pot-static_$(bias_str)V v ref_Rp0$(pO2_2_digits).cor"
      else
        file_speficitaion = "$(ocp_token).cor"
      end
      
      fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/$(total_order_str) YSZ $(plane) - Au-Au$(cell_specification)/$(sequence_str)$(TC) C/$(file_speficitaion)")
      
      try  
        push!(df,(get_U_I_from_IVpoint(fNAME)..., 0) ) 
      catch e
      
      end
    end
    return df
  
  
  
  
  
  
  
  
  
  elseif data_set=="POLY_I-V"
    # this needs to be added to separate simulation !!! ... IV_simulation
    return import_IVtoDataFrame_folder(TC=TC, pO2=pO2, bias_array=vcat(collect(0 : 0.1 : 1), collect(0.9 : -0.1 : -0.9), collect(-1 : 0.1 : 0)), 
          folder="snehurka/experimental_data_PSS/jako asi 6/$(TC) $(pO2_Dan) 6/")
  elseif data_set=="MONO_110_I-V"
    return import_IVtoDataFrame_folder(TC=TC, pO2=pO2, bias_array=vcat(collect(0 : 0.1 : 1), collect(0.9 : -0.1 : -0.9), collect(-1 : 0.1 : 0)), 
          folder="snehurka/experimental_data_PSS/YSZ 110/110 $(TC) $(pO2_Dan)/")
  else
    fNAME=string("snehurka/experimental_data_PSS/individual_files/$(data_set)")
  end
  return import_CVtoDataFrame_path(fNAME)
end





function default_import_EIStoDataFrame(;TC, pO2, bias, data_set="MONO_110", extra_tokens=Dict())
  ET = extra_tokens
  #
  pO2=Int64(pO2)
  pO2_2_digits = @sprintf("%01g",pO2)
  pO2_Dan = @sprintf("%02g",pO2)
  pO2_Michal = @sprintf("%03g",pO2)
  #
  TC=Int64(TC)
  #
  bias=Float64(bias)
  bias_mV = Int32(bias*1000)
  #
  if bias == 0.0
    ocp_token = "_ocp"
  else
    ocp_token = ""
  end
  #
  if length(data_set) > 3 && (data_set[end-3 : end-1] == "OCV")
    bias_mV = "0_"*data_set[end]
  end
  #
  if length(data_set) >= 4 && data_set[1:4]=="POLY"
    if abs(bias_mV) < 10 && data_set!="DAN_kapacity"
      bias_mV = "0_1"
    end
    fNAME=string("snehurka/experimental_data_PSS/jako asi 6/$(TC) $(pO2_Dan) 6/eis_$(bias_mV).z") 
  elseif data_set=="Zahner"
    fNAME=string("snehurka/experimental_data_PSS/individual_files/TEST DRT - Zahner - dummy cell.z")
  elseif data_set=="Solartron"
    fNAME=string("snehurka/experimental_data_PSS/individual_files/TEST DRT - Solartron - dummy cell.z")
  elseif data_set=="HebbWagner"
    # TC \in (600 : 20 : 720) ... bias = 0.3 ... pO2 = nizke, temer nulove
    fNAME=string("snehurka/experimental_data_PSS/HebbWagner/$(TC) C/$(TC)_EIS $(bias)V v ref 50mV amplituda.z")
  
  elseif data_set=="NEW_MONO_111"    
    if bias_mV == 0
      bias_mV = "mid_1"
    end
    fNAME=string("snehurka/experimental_data_PSS/DAN_nove_MONO/111/$(TC) $(pO2_Dan)/eis_$(bias_mV).z")
  
  
  elseif data_set=="Mince_l-h"
    fNAME=string("snehurka/experimental_data_PSS/Mince/1. cyklus - low-high/eis_$(TC).z")
  elseif data_set=="Mince_h-l"
    fNAME=string("snehurka/experimental_data_PSS/Mince/2. cyklus - high-low/eis_$(TC).z")
  
  
  ##########################################################
  # ! ! ! ! ! pO2 used as a repetition number ! ! ! ! ! ! 
  elseif data_set=="HebbWagner_110_more_biases"
    # TC \in (600 : 20 : 720) ...
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/02 YSZ 110 Ag-Ag/700 a 800 C vice biasu/$(TC) C/EIS $(bias)V v ref 50mV amplituda_Rp0$(pO2).z")
  
  elseif data_set=="HebbWagner_111_more_biases"
    # TC \in (600 : 100 : 700) ...
    if bias < 0.3 
      backwards_tag = "backwards_"
    else  
      backwards_tag = ""  
    end
    if bias == 0
      bias="0"
    end
    
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/03 YSZ 111 Ag-Ag/$(TC) C/multibiasy/$(backwards_tag)EIS $(bias)V v ref 50mV amplituda_Rp0$(pO2).z")
  
  elseif data_set=="HebbWagner_100"
    # TC \in (600 : 20 : 720) ...
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/04 YSZ 100 Ag-Ag/$(TC) C/EIS $(bias)V v ref 50mV amplituda_Rp0$(pO2).z")
    
  elseif data_set=="HebbWagner_110"
    # TC \in (600 : 20 : 720) ...
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/02 YSZ 110 Ag-Ag/vsechny teploty a jeden bias/$(TC) C/EIS $(bias)V v ref 50mV amplituda_Rp0$(pO2).z")
    
  elseif data_set=="HebbWagner_111"
    # TC \in (600 : 20 : 720) ...
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/03 YSZ 111 Ag-Ag/$(TC) C/EIS $(bias)V v ref 50mV amplituda_Rp0$(pO2).z")
      
      
      
  elseif data_set=="HW_TEST"
    # TC \in (600 : 20 : 720) ...
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/03 YSZ 111 Ag-Ag/$(TC) C/EIS $(bias)V v ref 50mV amplituda_Rp0$(ET["rep"]).z")
      
      
  elseif data_set=="HebbWagner_111_2nd"
    # TC \in (600 : 20 : 720) ...
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/06 YSZ 111 Ag-Ag/!data/$(TC) C/EIS $(bias)V v ref 50mV amplituda_Rp0$(pO2).z")
        
  elseif data_set=="HebbWagner_110_2nd"
    # TC \in (600 : 20 : 720) ...
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/07 YSZ 110 Ag-Ag/!data/$(TC) C/EIS $(bias)V v ref 50mV amplituda_Rp0$(pO2).z")
       
  
  elseif data_set=="HebbWagner_100_GOLD"
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/09 YSZ 100 - Au-Au/$(TC) C/EIS_$(bias)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")
  
  elseif data_set=="HebbWagner_100_GOLD_backward"
    backward_string = "backward_"
    if ocp_token != ""
      ocp_token = "_ocp_konec"
      backward_string = ""
    end
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/09 YSZ 100 - Au-Au/$(TC) C/$(backward_string)EIS_$(bias)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")
    
    
  ##  
    
  elseif data_set=="HebbWagner_110_GOLD"
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/10 YSZ 110 - Au-Au/druhe mereni/$(TC) C/EIS_$(bias)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")
  
  elseif data_set=="HebbWagner_110_GOLD_backward"
    backward_string = "backward_"
    if ocp_token != ""
      ocp_token = "_ocp_konec"
      backward_string = ""
    end
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/10 YSZ 110 - Au-Au/druhe mereni/$(TC) C/$(backward_string)EIS_$(bias)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")  

  ##  
    
  elseif data_set=="HebbWagner_111_GOLD"
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/11 YSZ 111 - Au-Au/mereni 1/$(TC) C/EIS_$(bias)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")
  
  elseif data_set=="HebbWagner_111_GOLD_backward"
    backward_string = "backward_"
    if ocp_token != ""
      ocp_token = "_ocp_konec"
      backward_string = ""
    end
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/11 YSZ 111 - Au-Au/mereni 1/$(TC) C/$(backward_string)EIS_$(bias)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")  
    
  ##  
    
  elseif data_set=="HebbWagner_111_GOLD_2"
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/11 YSZ 111 - Au-Au/mereni 2/$(TC) C/EIS_$(bias)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")
  
  elseif data_set=="HebbWagner_111_GOLD_2_backward"
    backward_string = "backward_"
    if ocp_token != ""
      ocp_token = "_ocp_konec"
      backward_string = ""
    end
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/11 YSZ 111 - Au-Au/mereni 2/$(TC) C/$(backward_string)EIS_$(bias)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")    
    
  ##  
    
  elseif data_set=="HebbWagner_111_GOLD_3"
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/11 YSZ 111 - Au-Au/mereni 3 otocene kabely/$(TC) C/EIS_$(bias)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")
  
  elseif data_set=="HebbWagner_111_GOLD_3_backward"
    backward_string = "backward_"
    if ocp_token != ""
      ocp_token = "_ocp_konec"
      backward_string = ""
    end
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/11 YSZ 111 - Au-Au/mereni 3 otocene kabely/$(TC) C/$(backward_string)EIS_$(bias)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")  
    
  ##  
    
  elseif data_set=="HebbWagner_111_GOLD_4"
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/11 YSZ 111 - Au-Au/mereni 4 pm2V/$(TC) C/EIS_$(bias)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")
  
  elseif data_set=="HebbWagner_111_GOLD_4_backward"
    backward_string = "backward_"
    if ocp_token != ""
      ocp_token = "_ocp_konec"
      backward_string = ""
    end
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/11 YSZ 111 - Au-Au/mereni 4 pm2V/$(TC) C/$(backward_string)EIS_$(bias)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")    

  elseif data_set=="cyc_HebbWagner_111_GOLD_8_4_cycles"
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/12 YSZ 111 - Au-Au LESTENA/mereni 4 cyklovani +-0.5 V/EIS jen jako poradi/$(TC) C/111-Au-Au_$(TC)C_20O2_$(Int(bias)).z")
    
  ##
  
#   elseif data_set=="HebbWagner_111_GOLD_5"
#     fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/11 YSZ 111 - Au-Au/mereni 4 pm2V/$(TC) C/EIS_$(bias)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")
#   
#   elseif data_set=="HebbWagner_111_GOLD_5_backward"
#     backward_string = "backward_"
#     if ocp_token != ""
#       ocp_token = "_ocp_konec"
#       backward_string = ""
#     end
#     fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/11 YSZ 111 - Au-Au/mereni 4 pm2V/$(TC) C/$(backward_string)EIS_$(bias)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")
    
  ##########################################################  
  elseif it_starts_with(data_set, "cyc_HebbWagner")
    # cyc means ->  >> bias << has the meaning of >> number with respect to time order << of experiments
    plane = data_set[16:18]    
    
    total_order_str = ""
    if plane == "100"
      total_order_str = "09"
    elseif plane == "110"
      total_order_str = "10"
    elseif plane == "111"
      ordering_number = data_set[25]
      if ordering_number in ['5', '6', '7']
        total_order_str = "12"
      else              
        total_order_str = "11"
      end
    end
    
    bias_return_value = 1.0
    cell_specification = ""
    sequence_str = ""
    if plane == "110"
      sequence_str = "druhe mereni/"
    end    
    if plane == "111"
      if ordering_number == '1'
        sequence_str = "mereni 1/"
      elseif ordering_number == '2'
        sequence_str = "mereni 2/"
      elseif ordering_number == '3'
        sequence_str = "mereni 3 otocene kabely/"
      elseif ordering_number == '4'        
        sequence_str = "mereni 4 pm2V/"
        bias_return_value = 2.0
      elseif ordering_number == '5'
        cell_specification = " LESTENA"
      elseif ordering_number == '6'
        cell_specification = " LESTENA"
        sequence_str = "mereni 2/"
      elseif ordering_number == '7'
        cell_specification = " LESTENA" 
        sequence_str = "mereni 3 do +- 2V/"
        bias_return_value = 2.0
      end
    end
    
    bias_digit_precision = 3
    bias_epsilon = 10e-5
    
    bias_idx = Int(bias)
    bias_step = 0.05
    dir = 1
    bias_value = 0.0
    bias_str = ""
    
    return_counter = 0
    for i in 1:(bias_idx - 1)    
      bias_value = round(bias_value + dir*bias_step, digits=bias_digit_precision)
      if bias_value > bias_return_value && dir > 0
        return_counter += 1
        bias_value = round(bias_value + -2*bias_step, digits=bias_digit_precision)
        dir = -1
      end
      if bias_value < -bias_return_value && dir < 0
        return_counter += 1
        bias_value = round(bias_value + 2*bias_step, digits=bias_digit_precision)
        dir = 1
      end
      if return_counter > 1 && bias_value > 0        
        bias_value = Inf
        break
      end
    end
    bias_value = round(bias_value, digits=bias_digit_precision)
    
    
    bias_str = "$(bias_value)"
    backward_string = ""
    ocp_token = ""
    if bias_idx == 1      
      ocp_token = "_ocp"
    elseif (bias_value > 0 + bias_epsilon) && (dir < 0)      
      backward_string = "backward_"
    elseif bias_value == 0.0 && dir < 0      
      bias_str = "-0.0"
      ocp_token = "_ocp"
    elseif bias_value < 0 - bias_epsilon && dir > 0
      backward_string = "backward_"
    elseif bias_idx > 1 && bias_value == 0.0 && dir > 0
      ocp_token = "_ocp_konec"
    end                
    
    fNAME=string("snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/$(total_order_str) YSZ $(plane) - Au-Au$(cell_specification)/$(sequence_str)$(TC) C/$(backward_string)EIS_$(bias_str)DC_10ac$(ocp_token)_Rp0$(pO2_2_digits).z")    
    
  
  ##########################################################
  
  
  elseif length(data_set) >= 8 && data_set[1:8]=="MONO_110"
    if abs(bias_mV) < 10 && data_set!="DAN_kapacity"
      bias_mV = "0_1"
    end  
    fNAME=string("snehurka/experimental_data_PSS/YSZ 110/110 $(TC) $(pO2_Dan)/eis_$(bias_mV).z")
  elseif data_set=="OLD_MONO_100"
    fNAME=string("snehurka/experimental_data_PSS/YSZ_09-2019_oxygen100/100 750to850 0to100%O2/$(TC)C/100 $(TC)C $(pO2_Dan)% do 1V/is $(bias)DC 50AC.z")
  elseif data_set=="MICHAL"
    fNAME=string("snehurka/experimental_data_PSS/K4/$(TC) C/$(pO2_Michal) O2/EIS_$(bias)DC_50AC$(ocp_token).z")
  
  elseif data_set=="K06"
    fNAME=string("snehurka/experimental_data_PSS/K06/$(TC) C/$(pO2_Michal) O2/EIS_$(bias)DC_50AC$(ocp_token).z")
  
  elseif data_set=="K07"
    fNAME=string("snehurka/experimental_data_PSS/K07/$(TC) C/$(pO2_Michal) O2/EIS_$(bias)DC_50AC$(ocp_token).z")  
  
  elseif data_set=="DAN_kapacity"
    fNAME=string("snehurka/experimental_data_PSS/PolySymAgAg/$(TC) C/eis_$(bias_mV).z")  
  
  
  
  else
    fNAME=string("snehurka/experimental_data_PSS/individual_files/$(data_set)")
  end
  
  return import_EIStoDataFrame_path(fNAME)
end

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

function CV_plot(CV_df, my_label="")
#    plot(df.U, df.Ib ,"blue"    ,label="bulk")
    #plot(phi_range[cv_range].-phi0, Ibb_range[cv_range] ,label="bulk_grad")
#    plot(df.U, df.Is, "green"   ,label="surf")
#    plot(df.U, df.Ir, "red"     ,label="reac")    
    PyPlot.title("CV curve")
    PyPlot.xlabel(L"\eta \ (V)")
    PyPlot.ylabel(L"I \ (A)")
    
    PyPlot.plot(CV_df.U, CV_df.I, label=my_label)

    if !(my_label == "")
        legend(loc="best")
    end
end

function Nyquist_plot(EIS_df, my_label="")
    title("Nyquist plot")
    xlabel("Re\$(Z)\$")
    ylabel("-Im\$(Z)\$")
   
    plot(real(EIS_df.Z), -imag(EIS_df.Z), "x-", label = my_label)
    
    if !(my_label == "")
        legend(loc="best")
    end
    grid(true)
    ax = gca()
    ax.set_aspect(1.0)
end

function main()    
    CVraw = import_CVtoDataFrame(fNAME)
end
