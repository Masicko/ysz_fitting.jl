using Printf
using PyPlot
#using Plots
using DataFrames
#using LeastSquaresOptim
using Optim
using CSV

using NNLS
#using LinearAlgebra
using LsqFit
using ProgressMeter

#include("../src/DRT.jl")


# structure ############
# input:
# # # EEC_structure
# # # alpha_range
# # # list of initial guesses (maybe from a file???)
# # # data_set_string + experimental conditions
# 
# EEC_data struct
#   function, 
#   list of parameters (strings), 
#   string representing EEC_structure
#   actual_fitted_parmeters
# build EEC (EEC_structure::string)
#   from string "L-RC-RCPE-RCPE" and return its representation (as a function) and list of parameters as strings
# import 
#   via implemented function
# perform fit!(alpha_range, list_of_initial_guesses)
#   use of LM .. .probably
# export_results
# 
#
#

L_units = 1.0e-6 # inductance is in micro-Henry


mutable struct EEC_data_struct
  structure::String
  #Z::Function # Z(omega, prms)::Complex
  prms_names::Array{String}
  prms_values::Array{Float32}
  #
  standard_deviations::Array{Float32}
  lower_limits_for_fitting::Array{Float32}
  upper_limits_for_fitting::Array{Float32}
  #
  EEC_data_struct() = new()
end

function EEC_data_struct(structure::String)
  this = EEC_data_struct()
  this.structure = structure
  return this
end


# function get_omega_and_prms_tuple_string(prms_names)
#   output = "(omega, "
#   for item in prms_names
#     output = "$(output)$item, "
#   end
#   output = output[1:end-1]
#   output = "$(output))"
#   return output
# end
  
function get_EEC(EEC_structure)

  
  EEC_actual = EEC_data_struct(EEC_structure)
  EEC_structure_splitted = split(EEC_actual.structure, "-")
  prms_names = Array{String}(undef, 0)
  
  for (i, token) in enumerate(EEC_structure_splitted)
    if token=="R"
      append!(prms_names, ["R$(i)"])
    elseif token=="L"
      append!(prms_names, ["L$(i)"])
    elseif token=="C"
      append!(prms_names, ["C$(i)"])  
    elseif token=="RCPE(tilt)"
      append!(prms_names, ["R$(i)", "C$(i)", "alpha$(i)"])
    elseif token=="RCPE"
      append!(prms_names, ["R$(i)", "C$(i)", "alpha$(i)"])
    elseif token=="RC"
      append!(prms_names, ["R$(i)", "C$(i)"])
    elseif token=="RL"
      append!(prms_names, ["R$(i)", "L$(i)"])
    else
      println("ERROR: build_EEC: Incorrect input of EEC_structure $(EEC_structure)")
      throw(Exception)
    end
  end
  
  append!(prms_names, ["Err"])
  EEC_actual.prms_names = prms_names
  return  EEC_actual
end

function g(EEC::EEC_data_struct, name::String)
  #@show name
  return EEC.prms_values[findall(x -> x==name, EEC.prms_names)[1]]
end

function evaluate_EEC(EEC::EEC_data_struct, omega)

  total_Z = 0
  for (i, token) in enumerate(split(EEC.structure, "-"))
    if token=="R"
      total_Z += g(EEC, "R$(i)")
    elseif token=="L"
      total_Z += im*omega*g(EEC, "L$(i)")*L_units      
    elseif token=="C"
      total_Z += 1/(im*omega*g(EEC, "C$(i)"))
    elseif token=="RCPE(tilt)"
      total_Z += g(EEC, "R$(i)")    /( 1 + g(EEC, "R$(i)") * g(EEC, "C$(i)") * (im*omega)^(g(EEC, "alpha$(i)")) )
    elseif token=="RCPE"
      total_Z += g(EEC, "R$(i)")    /( 1 + (g(EEC, "R$(i)") * g(EEC, "C$(i)") * im*omega)^(g(EEC, "alpha$(i)")) )
    elseif token=="RC"
      total_Z += g(EEC, "R$(i)")    /( 1 + im * omega* g(EEC, "R$(i)") * g(EEC, "C$(i)") )
    elseif token=="RL"
      total_Z += im*omega* g(EEC, "L$(i)") * g(EEC, "R$(i)")       /( g(EEC, "R$(i)") + im * omega* g(EEC, "L$(i)") )
    else
      println("ERROR: evaluate_EEC: Incorrect input of EEC_structure $(EEC_structure)")
      throw(Exception)
    end
  end
  return total_Z
end

function get_EIS_from_EEC(EEC_actual::EEC_data_struct; f_range=[])
  EIS_out = DataFrame( f = [], Z = [])
  #@show EEC_actual.prms_names
  #@show EEC_actual.prms_values
  for f in f_range
    push!(
      EIS_out, 
      (f, evaluate_EEC(EEC_actual, 2*pi*f))
    )
  end
  return EIS_out
end



# function HF_LF_check(prms_values)
#   if prms_values[3]*prms_values[4] > prms_values[6]*prms_values[7]
#     return true
#   else
#     return false
#   end
# end

function HF_LF_correction!(EEC)
  
  RCPE_count = get_count_of_RCPEs(EEC)
  
  
  characteristic_frequency_df = DataFrame(char_f = [], idx = [])
  for i in 1:RCPE_count
    push!(characteristic_frequency_df, (EEC.prms_values[3 + (i-1)*3]*EEC.prms_values[4 + (i-1)*3], i))
  end
  
  sort!(characteristic_frequency_df, :char_f)
  
  copy_of_prms_values = deepcopy(EEC.prms_values)
  for i in 1:RCPE_count
    for j in 3:5
      EEC.prms_values[j + (i-1)*3] = copy_of_prms_values[j + (characteristic_frequency_df.idx[i] - 1)*3]
    end
  end
end


function get_RC_prms_from_Z(omega, Z)
  a = real(Z)
  b = imag(Z)
  return (
    (a*a + b*b)/a,
    - b/(omega*(a*a + b*b))
    )
end

function get_RC_from_EIS(df)
  output = []
  for i in collect(1 : 1 : size(df, 1))
       push!(output, ysz_fitting.get_RC_prms_from_Z(df.f[i]*2*pi, df.Z[i]))
  end
  return output
end

function plot_RC(RC_list)
  R_list = []
  C_list = []
  x_range = []
  for (i, RC_item) in enumerate(RC_list)
      push!(R_list, RC_item[1])
      push!(C_list, RC_item[2])
      push!(x_range, log(i))
  end
  figure(325)
  plot(x_range, R_list)
  figure(326)
  plot(x_range, C_list)
  
  return
end

#
#     mask can be generated via string input, e.g. every "L" should be constant >>
#
function EEC_find_fit!(EEC_actual::EEC_data_struct, EIS_exp::DataFrame; mask=Nothing, alpha_low=0.0, alpha_upp=1.0, with_errors=false, error_type)
  
  projection_plot_maximum = 0.0002
  projection_plot_minimum = 10
  
  function EEC_plot_error_projection(prms_values, prms_names, error)
    fig_num = 333
    prms_length = length(prms_names)
    plot_edge = ceil(sqrt(prms_length))

    if projection_plot_maximum < error && error < 20
      projection_plot_maximum = error
    end
    
    if error < projection_plot_minimum
      projection_plot_minimum = error
    end
    
    figure(fig_num)
    for i in 1:prms_length
      subplot(plot_edge, plot_edge, i)
      xlabel(prms_names[i])
      ylabel("err")
      ylim(projection_plot_minimum, projection_plot_maximum)
      plot(prms_values[i], error, "o")
    end
    
    return
  end
  
  function prepare_prms(mask, x0, x)
      prms = zeros(0)
      xi = 1
      for i in collect(1 : 1 :length(mask))
          if convert(Bool,mask[i])
              append!(prms, x[xi])
              xi += 1
          else
              append!(prms, x0[i])
          end
      end
      return prms
  end
  
  function prepare_masked_stuff(mask, x0, lower_bounds, upper_bounds)
    x0M = zeros(0)
    lowM = zeros(0)
    uppM = zeros(0)
    for i in collect(1 : 1 : length(mask))
        if convert(Bool,mask[i])
            append!(x0M, x0[i])
            append!(lowM, lower_bounds[i])
            append!(uppM, upper_bounds[i])
        end
    end
    return x0M, lowM, uppM
  end
  
  function get_EIS_value_from_gf(EIS_df, gf_list)
    output = []
    for gf in gf_list
      if mod(gf,2) == 0
        append!(output, real(EIS_df.Z[div(gf,2)]))
      else
        append!(output, imag(EIS_df.Z[div(gf,2)]))
      end
    end
    return output
  end

  x_previous = Nothing
  dramatic_tolerance = 1.0e-4
  function check_dramatic_change(x)
    for (i, item) in enumerate(x)
      if abs(x[i] - x_previous[i]) > dramatic_tolerance
        return true
      end
    end
    return false
  end
  
  function to_optimize(x)    
    
    EEC_actual.prms_values = prepare_prms(mask, x0, x)
    EIS_EEC = get_EIS_from_EEC(EEC_actual, f_range=EIS_exp.f)        
    
    err = fitnessFunction(SIM, EIS_EEC, EIS_exp, error_type=error_type)
    
    
    
#     ################### 
#     #pause(0.1)
#     
#     EIS_EEC_plot = get_EIS_from_EEC(EEC_actual, f_range=EIS_get_checknodes_geometrical((1, 10000, 1.4)...))
#     if check_dramatic_change(x) || true
#       #typical_plot_sim(SIM, EIS_EEC_plot)
#     end
#     
#     #PyPlot.show()
#     
#       
#     EEC_plot_error_projection(
#       take_only_masked(mask, EEC_actual.prms_values), 
#       take_only_masked(mask, EEC_actual.prms_names), 
#       err)
#     
#      println("~~~~~ LM e = $(err)\nx = $(x)")
#     ##################
#     
#     
#     if !(check_x_in(x, lowM, uppM))
# #       ##############
#        @show x
#        @show lowM
#        @show uppM
#        println("    OUT OF THE BOUNDS   \n")
# #       ##############
# #       return 1000
#     end
    
    
    return err
  end
  
  function model(gf, x)    
    EEC_actual.prms_values = prepare_prms(mask, x0, x)
    EIS_EEC = get_EIS_from_EEC(EEC_actual, f_range=EIS_exp.f)
    SIM = EIS_simulation(800, 100, 0, use_DRT=false, plot_option="Bode Nyq", plot_legend=false)[1]
    #typical_plot_sim(SIM, EIS_EEC)
    #println("e = $(fitnessFunction(SIM, EIS_EEC, EIS_exp, error_type=error_type))\nx = $(x)")
    return get_EIS_value_from_gf(EIS_EEC, gf)
  end  
  
  SIM = EIS_simulation(800, 100, 0, use_DRT=false, plot_option="Bode Nyq", plot_legend=false)[1]
  
  prms_length = length(EEC_actual.prms_values)
  if mask == Nothing
    mask = Array{Int16}(undef, prms_length)
    mask .= 1
  end
  
  
  x0 = copy(EEC_actual.prms_values)
  
  lower_bounds = Array{Float32}(undef, prms_length)
  upper_bounds = Array{Float32}(undef, prms_length)
  lower_bounds_threshold = 1e-11
  for (i, name) in enumerate(EEC_actual.prms_names)
      if occursin("alpha", name)
        lower_bounds[i] = max(alpha_low, lower_bounds_threshold)
        upper_bounds[i] = alpha_upp
      else
        lower_bounds[i] = lower_bounds_threshold
        upper_bounds[i] = Inf
      end
  end
  for (i, name) in enumerate(EEC_actual.prms_names)
    lower_bounds[i] = max(lower_bounds[i], EEC_actual.lower_limits_for_fitting[i])
    upper_bounds[i] = min(upper_bounds[i], EEC_actual.upper_limits_for_fitting[i])
  end

  
  #######################################
  #######################################
  #######################################
  #######################################
  x0M, lowM, uppM = prepare_masked_stuff(mask, x0, lower_bounds, upper_bounds)
  #@show lowM
  #@show uppM  
 
  x_previous = x0M
 
  begin # LeastSquaresOptim
#   fit_LSO = optimize(to_optimize, x0M, lower=lowM, upper=uppM, 
#     #Δ=1000000, 
#     x_tol=1.0e-14,
#     f_tol=1.0e-14, 
#     g_tol=1.0e-14, 
#     #autodiff=:central,
#     LevenbergMarquardt(),
#     #Dogleg(),
#     )
  end
 
  if with_errors
    begin # LsqFit
      gf_range = collect( 2 : 1 : 1 + 2*length(EIS_exp.f))
      y_data = get_EIS_value_from_gf(EIS_exp, gf_range)   
      fit = curve_fit(model, gf_range, y_data, x0M, lower=lowM, upper=uppM)
      EEC_actual.prms_values = prepare_prms(mask, x0, fit.param)
    end  
  else
    begin # Optim
      fit_O = optimize(to_optimize, x0M, #lower=lowM, upper=uppM, 
        #Δ=1000000, 
        #x_tol=1.0e-18,
        #f_tol=1.0e-18,
        #g_tol=1.0e-18,
        #autodiff=:central,
        Optim.Options(iterations = 2000, f_tol=1.0e-18, g_tol=1.0e-21),
        #Optim.NelderMead(
        #  initial_simplex = Optim.AffineSimplexer()        
        #)
        #LevenbergMarquardt(),                
        #Dogleg(),
      )
      EEC_actual.prms_values = prepare_prms(mask, x0, fit_O.minimizer)  
    end
  end

  
  
  #@show fit_O
  
  
  error0 = deepcopy(x0)
  error0 .= -Inf;
  try
    EEC_actual.standard_deviations = prepare_prms(mask, error0, stderror(fit))
  catch
    #println("SO EIN PECH!")
    EEC_actual.standard_deviations = deepcopy(EEC_actual.prms_values)
  end
  return
end






















function get_left_right_width_of_EIS(EIS_df, N_for_sum=2)
  real_sum_left = 0
  real_sum_right = 0  
  for i in 1:N_for_sum
    real_sum_left += real(EIS_df.Z[end-i+1])
    real_sum_right += real(EIS_df.Z[i])
  end
  left = real_sum_left/N_for_sum
  right = real_sum_right/N_for_sum
  width = right - left
  return left, right, width
end





function get_count_of_RCPEs(EEC::EEC_data_struct)
  EEC_structure_splitted = split(EEC.structure, "-")
  return sum(map(x -> (x == "RCPE" ? 1 : 0), EEC_structure_splitted))
end


function get_init_values_RC_element(EIS_exp, EEC::EEC_data_struct)
  
end

function get_init_values(EIS_exp, EEC::EEC_data_struct)
  
  function get_capacitance_spread_item(i, total_number)
    return (i - total_number/2 - 0.5)
  end
  
  #smaller_circle_resistance_ratio = 0.25  # with respect to full width
  capacitance_spread_factor = 7/8
  
  # R | L | RCPE | RCPE structure REQUIERED!
  output = Array{Float64}(undef, length(EEC.prms_names))
  
  # resistors
  left, right, width = get_left_right_width_of_EIS(EIS_exp)
  
  minimum = Inf
  minimum_idx = -1
  for (i, item) in enumerate(imag(EIS_exp.Z))
    if item < minimum
      minimum = item
      minimum_idx = i
    end
  end
  highest_freq = EIS_exp.f[minimum_idx]
  
  # R_ohm
  output[1] = left
  
  # L
  #output[2] = 1.64
  if imag(EIS_exp.Z[end]) > 0
    #@show  imag(EIS_exp.Z[end])
    output[2] = imag(EIS_exp.Z[end])/(2*pi*EIS_exp.f[end])*(1/L_units)
  else
    output[2] = 0.0
  end
  
  if EEC.structure == "R-L-RC"
    EIS_aux = deepcopy(EIS_exp)    
    
    Z_aux = real.(EIS_aux.Z) .- left
    EIS_aux.Z = Z_aux  .+ im .* imag.(EIS_aux.Z)      
    RC_prms = get_RC_from_EIS(EIS_aux)     
    
    fitted_frequency_idx = 1
    output[3] = RC_prms[fitted_frequency_idx][1]
    output[4] = RC_prms[fitted_frequency_idx][2]
    #@show RC_prms[fitted_frequency_idx]
    return output
  
  else ## THIS SUPPOSES there are RCPE components    
    RCPE_count = get_count_of_RCPEs(EEC)
    
    for i in 1:RCPE_count
      # RCPE-R
      output[3 + (i-1)*3] = width*(1.0/RCPE_count)
      
      # RCPE-C
      central_capacitance = 1/(2*pi*highest_freq*output[3])
      if RCPE_count == 1
        spread_item = 0
      else
        spread_item = capacitance_spread_factor*(
              get_capacitance_spread_item(i, RCPE_count) /
              get_capacitance_spread_item(RCPE_count, RCPE_count)
            )
      end
      output[4 + (i-1)*3] = central_capacitance*(1 + spread_item)
      
      # RCPE-alpha
      output[5 + (i-1)*3] = 0.9        


  # #   # R3, R4
  # #   output[3] = width*smaller_circle_resistance_ratio
  # #   output[6] = width*(1 - smaller_circle_resistance_ratio)
  # #   
  # #   # alphas
  # #   output[5] = 1.0
  # #   output[8] = 0.8
  # #   
  # #   # C3, C4
  # # 
  # #   
  # #   output[4] = central_capacitance*(capacitance_spread_factor)
  # #   output[7] = central_capacitance*(capacitance_spread_factor)
  # #   
  #   @show central_capacitance
  #   @show output[4]
  #   @show output[7]
    end
    
    # C
    if length(output)==2 + 3*RCPE_count + 1
      output[2 + 3*RCPE_count + 1] = max(0, -1/(2*pi*EIS_exp.f[1]*imag(EIS_exp.Z[1])))
    end
  end
  
  
  #@show imag(EIS_exp.Z[end]) 
  #@show left, right, width
  #@show output
  return output
  #return [2.0, 6.2, 0.5 , 0.001, 1.0,    0.6, 0.01, 0.8]
end

function set_fitting_limits_to_EEC_from_EIS_exp!(EEC::EEC_data_struct, EIS_exp)

  left, right, width = get_left_right_width_of_EIS(EIS_exp)
    
  prms_length = length(EEC.prms_names)
  lower_limits = Array{Float32}(undef, prms_length)
  upper_limits = Array{Float32}(undef, prms_length)
  
  # R-L-RC structure
  if EEC.structure == "R-L-RC"
    lower_limits .= Inf
    upper_limits .= Inf
  else
    
    # R | L | RCPE | RCPE structure REQUIERED!
    
    RCPE_count = get_count_of_RCPEs(EEC)
    
    # R_ohm
    #if RCPE_count > 0
    if false
      lower_limits[1] = left/2
      upper_limits[1] = right
    else
      lower_limits[1] = -Inf
      upper_limits[1] = Inf
    end
    
    # L2    
    if imag(EIS_exp.Z[end]) < 0.0
      lower_limits[2] = -0.1
      upper_limits[2] = 0.1
    else
      lower_limits[2] = max(-0.1, (imag(EIS_exp.Z[end])/(2*pi*EIS_exp.f[end])*(1/L_units))/10 )
      upper_limits[2] = max(0.1, (imag(EIS_exp.Z[end])/(2*pi*EIS_exp.f[end])*(1/L_units))*10 )  
    end
    
    
    for i in 1:RCPE_count
      #R
      lower_limits[3 + (i-1)*3] = -Inf
      upper_limits[3 + (i-1)*3] = Inf # width*2   
      
      #C
      lower_limits[4 + (i-1)*3] = 0.0
      upper_limits[4 + (i-1)*3] = 10    
      
      #alpha
      lower_limits[5 + (i-1)*3] = -Inf
      upper_limits[5 + (i-1)*3] = Inf    
    end
    
    # C[end]
    if prms_length == 2 + 3*RCPE_count + 1      
      lower_limits[end] = 0.0
      upper_limits[end] = Inf
    end
  end
  
  EEC.lower_limits_for_fitting = lower_limits
  EEC.upper_limits_for_fitting = upper_limits
  return
end


function view_EEC(;
      EEC_structure="RL-R-L-RCPE-RCPE", 
      prms_values=[1, 0.001,      1.7, 0,        1. , 0.001, 1.0,    1.0, 0.01, 0.8], 
      f_range=(0.01, 10000, 1.2),
      print_bool=false, pyplot=1, plot_legend=true, use_DRT=true, DRT_draw_semicircles=false,
      )
  
  EEC = get_EEC(EEC_structure)
  
  function recursive_view_EEC_call(output_prms, plot_names, plot_values, active_idx)
      if active_idx > size(prms_names,1)

        ###################################################################
        # for each combination of output_prms from prms_lists perform this #######
        prms_names_in=[]
        prms_values_in=[]
        
        append!(prms_names_in, prms_names)
        append!(prms_values_in, output_prms)
      
        
        if check_equal_size(prms_names, prms_values)  
          EEC.prms_values = output_prms
          EIS_EEC = get_EIS_from_EEC(EEC, f_range=EIS_get_checknodes_geometrical(f_range...))   
        end    
        if size(plot_names,1) < 1
          plot_prms_string = ""
        else
          plot_prms_string = " $(string(plot_names)) = $(plot_values)"
        end
        if print_bool
          #@show EEC.prms_names
          #@show EEC.prms_values
          
          println(get_string_EEC_prms(EEC, plot_errors=false))
        end
        if pyplot > 0
            ysz_fitting.typical_plot_exp(EIS_simulation(800, 100, 0.0, f_range=f_range, plot_legend=plot_legend, use_DRT=use_DRT, DRT_draw_semicircles=DRT_draw_semicircles,)[1], EIS_EEC, "!EEC"*plot_prms_string) 
        end 
        
        
        return
        ###################################################################
        ###################################################################
      
      end
      if size(prms_values[active_idx],1)>1
        for i in prms_values[active_idx]
          recursive_view_EEC_call(
            push!(deepcopy(output_prms),i),
            push!(deepcopy(plot_names),prms_names[active_idx]),
            append!(deepcopy(plot_values),i),
            active_idx + 1)
        end
      else
        recursive_view_EEC_call(
          push!(deepcopy(output_prms),prms_values[active_idx][1]),
          plot_names,
          plot_values,
          active_idx + 1)
      end
    end
  
  prms_names = EEC.prms_names
  recursive_view_EEC_call([], Array{String}(undef,(0)), Array{Float64}(undef,(0)), 1)
  
  return EEC
end

function get_string_EEC_prms(EEC::EEC_data_struct; plot_errors=false)
  output = ""
  for (i, item) in enumerate(EEC.prms_values)
    output *= "$(@sprintf("%8s", EEC.prms_names[i]))  =  $(@sprintf("%12.8f", item))    "
    plot_errors ? (output *="(rel_err: $(@sprintf("%12.8f", EEC.standard_deviations[i]/item)))\n") : output *= "\n"
  end
  return output[1:end-1]
end









mutable struct EEC_data_holder_struct
  TC
  pO2
  bias
  data_set
  
  prms_names
  # indexing of data: [TC, pO2, bias, data_set, prms_values]
  data::Array{Any}
  
  EEC_data_holder_struct() = new()
end

function EEC_data_holder_struct(TC, pO2, bias, data_set, prms_names)
  this = EEC_data_holder_struct()
  this.TC = TC
  this.pO2 = pO2
  this.bias = bias
  this.data_set = make_array_from_string(data_set)
  this.prms_names = prms_names
  
  this.data = Array{Any}(undef, (length(this.TC), length(this.pO2), length(this.bias), length(this.data_set), length(this.prms_names)))
  return this
end



function save_EEC_prms_item_to_file(TC, pO2, bias, data_set, prms_names, prms_values, saving_destination; append=true, save_only_R1=false)

  full_prms_names = ("TC", "pO2", "bias", "data_set", (prms_names)...)
  full_prms_values = (TC, pO2, bias, data_set, prms_values...)

  df_out = DataFrame()
    
  if save_only_R1
    if append==true      
      df_out[!, Symbol("R1")] = [prms_values[findall(x -> x=="R1", prms_names)[1]]]
    end
  else
    for (i, name) in enumerate(full_prms_names)
      df_out[!, Symbol(name)] = length(full_prms_values) == length(full_prms_names) ? [full_prms_values[i]] : []
    end
  end
  CSV.write(saving_destination, df_out, delim="\t", append=append)
end
        


function run_EEC_fitting(;TC=800, pO2=80, bias=0.0, data_set="MONO_110",
                        EIS_input=Nothing,
                        f_interval=Nothing, succes_fit_threshold = 0.002, error_type="normalized",
                        fixed_prms_names=[], fixed_prms_values=[],
                        init_values=Nothing, alpha_low=0.2, alpha_upp=1, 
                        EEC_structure="R-L-RCPE-RCPE",
                        plot_bool=false, plot_legend=true, plot_best_initial_guess=false, plot_fit=true, plot_exp=true,
                        show_all_initial_guesses=false,
                        with_errors=false, which_initial_guess="both",
                        use_DRT=false, DRT_draw_semicircles=false,
                        trim_inductance=false,
                        save_file_bool=true, save_to_folder="../data/EEC/", file_name="default_output.txt", save_R1_file=false,
                        print_headers_bool=false,
                        EIS_preprocessing_control = EIS_preprocessing_control()
#                           ,EIS_preprocessing_control = ysz_fitting.EIS_preprocessing_control(
#                                   f_interval=Nothing, 
#                                   add_inductance=0,
#                                   trim_inductance=false, 
#                                   outlayers_threshold=5.5,                                    
#                                   use_DRT=false, DRT_control=ysz_fitting.DRT_control_struct()
#                            )
#                         )
                        
        )
  ####
  ####  TODO:
  ####  [?] initial values handling (from file, probably)
  ####  [x] R_ohm estimation (+ R2, R3 estimation?)
  ####  [x] swap RCPEs in the case Hf and Lf interchanged
  ####  [x] check, if fitting converged
  ####  [x] output format
  ####  [x] view_EEC ... recursive call
  ####  [x] define default SIM for ploting Nyquists! and turn of the legend!
  ####  [x] f_interval ... to crop the Nyquist
  ####  [ ] sjednotit plot options, abych mohl zapinat DRT a menit ploty v klidecku!
  ####  [x] upravit hranice parametrickeho prostoru, aby nehledal kraviny!
  ####  [ ] fixed_prms_names ... bug for *3 a *4  ... elements can be interchanged !!!
  ####  [ ] maybe try several automatic initial guesses and pick up the best one?
  ####  [x] progress bar
  ####  [ ] plot_EEC_data_general -> vymyslet zaznamenavani stabilnich velicin
  ####  [ ] !!! 750, 100, 0.0, "MONO_110" dost blbe trefuje nizke frekvence
  ####  [ ] !!! kolem TC 800, MONO_110 se dejou divne veci v R1
  
  ####  [ ] EEC_structure by melo obsahovat strukturu obvodu
  
  ########## Questions
  ####  [x] Jak se mam zachovat, kdyz fit neni dobry? Mam proste preskocit, nic nezapisovat do souboru a jet dal? Nebo zapsat a warning?
  ####  [x] Format zapisu do souboru -> nemel by to spis byt *.csv soubor s cisly oddelenymi carkami ci tabulatory?
  ####  [x] da se orezat Nyquist pro nizke frekvence, at nejde do zapornych cisel. Chceme?
  ####  [x] teoreticky by se plotovani obrazku mohlo vzdy vypnout, kdyz by clovek chtel ukladat do souboru. Ale myslim, ze to neni nutne
  ####  [x] nazor na to, ze pocitam prumer systemu jako 1.2 cm, pricemz to je jne prumer elektrody, ale ellyt ma 2.5 cm v prumeru
  ####  [x] da se udelat nejake moudre orezani, ktere nedovoli vice hodnotam jit doprava oproti pruseciku s realnou osou
  ####  [ ] nechcete udelat treba animace nebo dalsi obrazky? (vyhledove)
  ####  [ ] chcete vypisovat L2 v mikro-Henry nebo v Henry? stejne jako zadavani a kdekoliv
  ####  [ ] jsou limity pro alphu [0.2, 1.0] ok? nebo by bylo lepsi 0.1? (protoze pro 0.01 to delalo divne veci) 750, 80, -0.5, "MONO_110"
  
  function succesful_fit(error, EIS_EEC, EIS_exp, error_type)
    #if fitnessFunction(EIS_simulation(), EIS_EEC, EIS_exp, error_type=error_type) > succes_fit_threshold
    if error > succes_fit_threshold
      return false
    else
      return true
    end
  end
  

    
  function get_fitting_mask_and_apply_fixed_prms_values!(EEC, init_values_list, fixed_prms_names, fixed_prms_values)
    if length(fixed_prms_names) != length(fixed_prms_values)
      println("ERROR: length(fixed_prms_names) != length(fixed_prms_values) ==>> $(length(fixed_prms_names)) != $(length(fixed_prms_values))")
      return throw(Exception)
    else
      output = []
      skip_bool = false
      for (i, original_name) in enumerate(EEC.prms_names)
        skip_bool = false
        for (j, name) in enumerate(fixed_prms_names)
          if name == original_name
            append!(output, 0)
            for init_values in init_values_list
              init_values[i] = fixed_prms_values[j]
            end
            skip_bool = true
            break
          end
        end
        skip_bool || append!(output, 1)
      end
      
      if length(fixed_prms_names) + sum(output) != length(EEC.prms_names)
        println("ERROR: some of fixed_prms_names $(fixed_prms_names) not found!")
        return throw(Exception)
      end
      #@show output
      return output
    end
  end
  
  EEC_actual = get_EEC(EEC_structure)
      
  data_set = make_array_from_string(data_set)  
  EEC_data_holder = EEC_data_holder_struct(TC, pO2, bias, data_set, EEC_actual.prms_names)
  
  warning_buffer = ""
  
  if save_file_bool 
    mkpath(save_to_folder)
    saving_destination = save_to_folder*file_name
    save_EEC_prms_item_to_file([], [], [], [], EEC_actual.prms_names, [], saving_destination, append=false)
    if save_R1_file
      R1_file_name = split(file_name, '.')[1]*"_R1."*split(file_name, '.')[2]
      R1_saving_destination = save_to_folder*R1_file_name
      save_EEC_prms_item_to_file([], [], [], [], EEC_actual.prms_names, [], R1_saving_destination, append=false, save_only_R1=true)
    end
  end
  
  cycle_number = 0
  cycle_number_for_progressbar = 0
  previous_bias = Inf
  
  progress_meter = Progress( 
                        length(TC)*length(pO2)*length(bias)*length(data_set), 
                        dt=0.5                        
                      )
  
  for   (TC_idx, TC_item) in enumerate(TC), 
        (pO2_idx, pO2_item) in enumerate(pO2), 
        (bias_idx, bias_item) in enumerate(bias), 
        (data_set_idx, data_set_item) in enumerate(data_set)
    cycle_number += 1
    
    
    if save_file_bool
      ProgressMeter.update!(progress_meter, cycle_number_for_progressbar)
      ProgressMeter.update!(progress_meter, cycle_number_for_progressbar)    
    end
    cycle_number_for_progressbar +=1
  
    SIM_list = ysz_fitting.EIS_simulation(TC_item, pO2_item, bias_item, 
                  use_DRT=use_DRT, DRT_draw_semicircles=DRT_draw_semicircles, 
                  plot_option="Bode Nyq DRT RC", plot_legend=plot_legend, data_set=data_set_item)
    SIM = SIM_list[1]    
    
    
    
    
    
    try
      if EIS_input != Nothing
        ########### TODO !!! very ugly solution !!!
        EIS_exp = EIS_input
      else
        EIS_exp = ysz_fitting.import_data_to_DataFrame(SIM)
      end
    catch 
      warning_buffer *= " ->->-> ERROR: file for (TC, pO2, bias, data_set_item) = ($(TC_item), $(pO2_item), $bias_item, $data_set_item) ... NOT FOUND! \n"
      missing_array = Array{Any}(undef, length(EEC_actual.prms_names))
      missing_array .= Missing
      
      EEC_data_holder.data[TC_idx, pO2_idx, bias_idx, data_set_idx, :] = deepcopy(missing_array)
      if save_file_bool
        save_EEC_prms_item_to_file(TC_item, pO2_item, bias_item, data_set_item, EEC_actual.prms_names, missing_array, save_to_folder*file_name)
        if save_R1_file
          save_EEC_prms_item_to_file(TC_item, pO2_item, bias_item, data_set_item, EEC_actual.prms_names, missing_array, save_to_folder*R1_file_name, save_only_R1=true)
        end
      end
      ProgressMeter.update!(progress_meter, cycle_number_for_progressbar)
      cycle_number += -1
      continue
    end

    
    
    
    
    EIS_exp = EIS_preprocessing(EIS_exp, EIS_preprocessing_control)

## HERE
    set_fitting_limits_to_EEC_from_EIS_exp!(EEC_actual, EIS_exp)
    
## HERE
    init_values_list = []
    if init_values == Nothing || cycle_number > 1
      if  which_initial_guess == "alg"
        push!(init_values_list, get_init_values(EIS_exp, EEC_actual))      
      elseif which_initial_guess == "both"
        push!(init_values_list, get_init_values(EIS_exp, EEC_actual))
        if cycle_number > 1
          push!(init_values_list, deepcopy(EEC_actual.prms_values))
        end
      else
        if which_initial_guess == "previous" && length(bias) > 1
          
          bias_step = bias[2] - bias[1]
          
          if abs((bias_item - bias_step) - previous_bias) < 0.000001
            push!(init_values_list, deepcopy(EEC_actual.prms_values))
          else
            push!(init_values_list, get_init_values(EIS_exp, EEC_actual))
          end
        else
          push!(init_values_list, get_init_values(EIS_exp, EEC_actual))
        end
      end
    else
      push!(init_values_list, init_values)
    end
    
    
    
    mask = get_fitting_mask_and_apply_fixed_prms_values!(EEC_actual, init_values_list, fixed_prms_names, fixed_prms_values)
    
    plot_bool && plot_exp && ysz_fitting.typical_plot_exp(SIM, EIS_exp)
            
    best_error = Inf
    best_prms_values = []
    best_init_values_idx = -1
    
    
    for (init_values_idx, init_values) in enumerate(init_values_list)
        EEC_actual.prms_values = init_values        
        
        EIS_EEC_pre = get_EIS_from_EEC(EEC_actual, f_range=EIS_exp.f)
        if show_all_initial_guesses
          plot_bool && ysz_fitting.typical_plot_sim(SIM, EIS_EEC_pre, "!EEC initial guess $(init_values_idx)")
        end
        
        EEC_find_fit!(EEC_actual, EIS_exp, mask=mask, alpha_low=alpha_low, alpha_upp=alpha_upp, with_errors=with_errors, error_type=error_type)
        
        HF_LF_correction!(EEC_actual)        
        
        EIS_EEC = get_EIS_from_EEC(EEC_actual, f_range=EIS_exp.f)
        actual_error = fitnessFunction(EIS_simulation(), EIS_EEC, EIS_exp, error_type=error_type)
        if actual_error < best_error
          best_error = actual_error
          best_prms_values = deepcopy(EEC_actual.prms_values)
          best_init_values_idx = init_values_idx
        end
        
        if show_all_initial_guesses && !(save_file_bool) 
          println(" --------- $(init_values_idx) --------- ")
                  
          output_string = "========== TC, pO2, bias: ($TC_item, $(pO2_item), $bias_item) --- data_set_item = $data_set_item) ==========\n"
          output_string *= "$(get_string_EEC_prms(EEC_actual))"
        
          println(output_string)
          println("init_values $(init_values_idx) = $(EEC_actual.prms_values)")
          println("FITTED_values $(init_values_idx) = $(EEC_actual.prms_values)   <<<<<<<<<<<<<<<<<<<<<< ")
          println(">>> error $(init_values_idx) >>> $(actual_error)\n")
        end
        
        EIS_EEC = get_EIS_from_EEC(EEC_actual, f_range=EIS_exp.f)
        if plot_fit && show_all_initial_guesses
          plot_bool && ysz_fitting.typical_plot_sim(SIM, EIS_EEC, "!EEC fit $(init_values_idx)")
        end
    end
    
    if save_file_bool
      ProgressMeter.update!(progress_meter, cycle_number_for_progressbar)
    end    
    
    EEC_actual.prms_values .= best_prms_values
    if !save_file_bool   
      if show_all_initial_guesses
        println(" ------------ BEST ---------- ")
      end
      output_string = "========== TC, pO2, bias: ($TC_item, $(pO2_item), $bias_item) --- data_set_item = $data_set_item) ==========\n"
      output_string *= "$(get_string_EEC_prms(EEC_actual))"
      
      println(output_string)
      println("best init_values = $(init_values_list[best_init_values_idx])")
      println("best FITTED_values = $(best_prms_values)")
      println(">>> best error >>> $(best_error)\n")
    end
    if print_headers_bool
      println("========== TC, pO2, bias: ($TC_item, $(pO2_item), $bias_item) --- data_set_item = $data_set_item) ========== ")
    end
    
    
    if plot_best_initial_guess
      EEC_actual.prms_values = init_values_list[best_init_values_idx]
      EIS_EEC_pre = get_EIS_from_EEC(EEC_actual, f_range=EIS_exp.f)
      plot_bool && ysz_fitting.typical_plot_sim(SIM, EIS_EEC_pre, "!EEC best initial guess")
    end
    
    EEC_actual.prms_values = deepcopy(best_prms_values)
    EIS_EEC = get_EIS_from_EEC(EEC_actual, f_range=EIS_exp.f)
    actual_error = best_error
    EEC_actual.prms_values[end] = best_error

    
    if plot_fit
      plot_bool && ysz_fitting.typical_plot_sim(SIM, EIS_EEC, "!EEC best fit")
    end
    
    if save_file_bool
      save_EEC_prms_item_to_file(TC_item, pO2_item, bias_item, data_set_item, EEC_actual.prms_names, EEC_actual.prms_values, save_to_folder*file_name)
      if save_R1_file        
        save_EEC_prms_item_to_file(TC_item, pO2_item, bias_item, data_set_item, EEC_actual.prms_names, EEC_actual.prms_values, save_to_folder*R1_file_name, save_only_R1=true)    
      end
    end
    
    EEC_data_holder.data[TC_idx, pO2_idx, bias_idx, data_set_idx, :] = deepcopy(EEC_actual.prms_values)
    
    
    if !(succesful_fit(actual_error, EIS_EEC, EIS_exp, error_type))
      warning_buffer *="WARNING: (TC, pO2, bias, data_set_item) = ($(TC_item), $(pO2_item), $bias_item, $data_set_item) maybe NOT CONVERGED !!! // error $(actual_error) > defined threshold $(succes_fit_threshold) //\n"
      # TODO ... optionally save picture of the result for humanous check
      # or this reports can go to another error_log_file.txt
    end
          
    previous_bias = bias_item
    
  end
  
  println()
  println("          =================================")
  println("          =================================")
  println("          =========== MESSAGES ============")
  if length(warning_buffer) > 0 
    println(warning_buffer)
  else
    println(" no warnings :) ")
  end
  
  return EEC_data_holder
#   T_holder = TC_holder .+ 273.15
#   
#   if length(bias_holder) > 1 || length(T_holder) > 1
#     if bias_holder[1] == bias_holder[2]
#       figure(42)
#       title("sigma VS bias")
#       plot(1000 ./ (T_holder), log.(sigma_holder .* T_holder), label="bias "*string(bias_holder[1]))
#       legend(loc="best")
#       xlabel("1000/T")
#       ylabel("log T sigma")
#     else
#       figure(44)
#       title("sigma VS \"TC\"")
#       plot(bias_holder, sigma_holder, label="TC "*string(TC_holder[1]))
#       legend(loc="best")
#       xlabel("bias")
#       ylabel("sigma")
#     end
#   end
#   
  
end


function save_EEC_data_holder(EEC_data_holder; folder="../data/EEC/", file_name="default.txt")
  mkpath(folder)
  saving_destination = folder*file_name
  save_EEC_prms_item_to_file([], [], [], [], EEC_data_holder.prms_names, [], saving_destination; append=false)
  
  for (TC_idx, TC_item) in enumerate(EEC_data_holder.TC), 
      (pO2_idx, pO2_item) in enumerate(EEC_data_holder.pO2), 
      (bias_idx, bias_item) in enumerate(EEC_data_holder.bias), 
      (data_set_idx, data_set_item) in enumerate(EEC_data_holder.data_set)
    save_EEC_prms_item_to_file(TC_item, pO2_item, bias_item, data_set_item, 
                EEC_data_holder.prms_names,
                EEC_data_holder.data[TC_idx, pO2_idx, bias_idx, data_set_idx, : ], saving_destination)
  end
end

function load_EEC_data_holder(;folder="../data/EEC/", file_name="default.txt")
  function get_unique_list_and_repetitivity_and_period(list)
    # how many items in the row are the same
    repetitivity = -1
    for (i, item) in enumerate(list)
      if item != list[1]
        repetitivity = i-1
        break
      end
    end
    if repetitivity == -1
      return list[1:1], 0, 0
    end
        
    # how long does it take to arrive at the same item
    period = -1
    for (i, item) in enumerate(list[1 : repetitivity : length(list)])
      
      if item == list[1] && i != 1
        period = repetitivity*(i-1)
        break
      end
    end    
    if period == -1
      period = length(list)
    end
    list[1 : repetitivity : period], repetitivity, period
  end
  
  saving_destination = folder*file_name
  #data_df = CSV.read(saving_destination)
  data_df = CSV.File(saving_destination) |> DataFrame
  
  TC_list, TC_repet, TC_per = get_unique_list_and_repetitivity_and_period(data_df.TC)
  pO2_list, pO2_repet, pO2_per = get_unique_list_and_repetitivity_and_period(data_df.pO2)
  bias_list, bias_repet, bias_per = get_unique_list_and_repetitivity_and_period(data_df.bias)
  data_set_list, data_set_repet, data_set_per = get_unique_list_and_repetitivity_and_period(data_df.data_set)
    
  @show names(data_df)[5 : end]  
  EEC_data_holder = EEC_data_holder_struct(TC_list, pO2_list, bias_list, data_set_list, names(data_df)[5 : end])
  
  overall_counter = 0
  for (TC_idx, TC_item) in enumerate(EEC_data_holder.TC), 
      (pO2_idx, pO2_item) in enumerate(EEC_data_holder.pO2), 
      (bias_idx, bias_item) in enumerate(EEC_data_holder.bias), 
      (data_set_idx, data_set_item) in enumerate(EEC_data_holder.data_set)
    
    overall_counter += 1
    for (prm_name_idx, prm_name) in enumerate(EEC_data_holder.prms_names)    
      token = data_df[overall_counter, Symbol(prm_name)]            
      if ((token == "Missing") || (token == "missing"))
        token = Missing        
      else
        token = parse(Float32, token)
      end
      EEC_data_holder.data[TC_idx, pO2_idx, bias_idx, data_set_idx, prm_name_idx] = token
    end
  end
  return EEC_data_holder
end


function plot_EEC_data_general(EEC_data_holder; 
                                y_name="C3",
                                x_name="bias",
                                #
                                TC_interval = [-Inf, Inf],
                                pO2_interval = [-Inf, Inf],
                                bias_interval = [-Inf, Inf],
                                data_set = Nothing,                              
                                #
                                fig_num=102, 
                                plot_legend=true, plot_all_prms=true,
                                reversed_x = false,
                                x_mapping = Nothing,
                                marker="-x"
                                )
  
  function make_range_list(prm_name, interval)
    output_list = []
    if prm_name == "data_set"
      for (i, item) in enumerate(getfield(EEC_data_holder, Symbol(prm_name)))
        if item in interval
          append!(output_list, i)
        end
      end
    else
      for (i, item) in enumerate(getfield(EEC_data_holder, Symbol(prm_name)))
        if interval[1] <= item && item <= interval[2]
          append!(output_list, i)
        end
      end
    end
    return output_list
  end
  
  function add_legend_contribution(idx, plot_idx)
    if length(plot_idx) == 1 && length(range_list[idx]) > 1
      return "$(identifier_list[idx])=$(getfield(EEC_data_holder, Symbol(identifier_list[idx]))[plot_idx]) "
    else
      return ""
    end
  end
  
  x_name_idx = -1
  identifier_list = ["TC", "pO2", "bias", "data_set"]
  range_list = Array{Any}(undef, 4)
  range_idx_list = Array{Any}(undef, 4)
  for i in 1:4
    if x_name == identifier_list[i]
      x_name_idx = i
    end
  end
    
  if y_name == "C3&C4"    
    plot_template = "1.0/(1/C3 + 1/ C4)"
  elseif (length(y_name) >= 5) && (y_name[1:5] == "(CPE)")
    if y_name[6:end] == "C3"          
      plot_template = "((C3*R3)^(1.0/alpha3))/R3"
    end
    if y_name[6:end] == "C4"          
      plot_template = "((C4*R4)^(1.0/alpha4))/R4"
    end
    if y_name[6:end] == "C3&C4"          
      plot_template = "1.0  /  (1/(((C3*R3)^(1.0/alpha3))/R3)   +    1/(((C4*R4)^(1.0/alpha4))/R4))"
    end
  elseif (length(y_name) >= 8) && (y_name[1:8] == "(CPEsin)")
    if y_name[9:end] == "C3"          
      plot_template = "(((C3*R3)^(1.0/alpha3))/R3 )*sin(alpha3*pi/2)"
    end
  else
    plot_template = y_name
  end
  
  
  range_list[1] = make_range_list("TC", TC_interval)
  range_list[2] = make_range_list("pO2", pO2_interval)
  range_list[3] = make_range_list("bias", bias_interval)
  if data_set != Nothing 
    data_set = make_array_from_string(data_set)
  else
    data_set = make_array_from_string(EEC_data_holder.data_set)
  end
  range_list[4] = make_range_list("data_set", data_set)

  figure(fig_num)
  xlabel(x_name)
  ylabel(y_name)
  
  function get_plot_list_from_template(plot_template, EEC_data_holder, TC_idx, pO2_idx, bias_idx, data_set_idx)    
    output_length = length(TC_idx)*length(pO2_idx)*length(bias_idx)*length(data_set_idx)
    substitued_template_list = Array{Any}(undef, output_length)
    substitued_template_list .= deepcopy(plot_template)
    
    for (prm_idx, prm) in enumerate(EEC_data_holder.prms_names)
      active_prm_vs_x_name_list = EEC_data_holder.data[TC_idx, pO2_idx, bias_idx, data_set_idx, prm_idx]  
      active_prm_vs_x_name_list = replace(active_prm_vs_x_name_list, Missing => missing)      
      for (value_idx, value) in enumerate(active_prm_vs_x_name_list)
        substitued_template_list[value_idx] = replace(substitued_template_list[value_idx], prm => value)
      end
    end    
    return eval.(Meta.parse.(substitued_template_list))
  end
  
  function x_mapping_default(x)
    return x
  end
  
  y_to_plot = []
  for TC_idx in (x_name_idx != 1 ? range_list[1] : [range_list[1]]),
      pO2_idx in (x_name_idx != 2 ? range_list[2] : [range_list[2]]),
      bias_idx in (x_name_idx != 3 ? range_list[3] : [range_list[3]]),
      data_set_idx in (x_name_idx != 4 ? range_list[4] : [range_list[4]])
    
    x_to_plot = getfield(
        EEC_data_holder, 
        Symbol(identifier_list[x_name_idx])
      )[range_list[x_name_idx]]
               
    if x_mapping!=Nothing
      x_to_plot = map(x_mapping, x_to_plot)
    end
               
    y_to_plot = get_plot_list_from_template(plot_template, EEC_data_holder, TC_idx, pO2_idx, bias_idx, data_set_idx)
        
    valid_idxs = map( x -> typeof(x) != Missing, y_to_plot)    
    plot(    
      reversed_x ? reverse(x_to_plot[valid_idxs]) : x_to_plot[valid_idxs] , 
      y_to_plot[valid_idxs] ,
      label= "$(add_legend_contribution(1, TC_idx))$(add_legend_contribution(2, pO2_idx))$(add_legend_contribution(3, bias_idx))$(add_legend_contribution(4, data_set_idx))"*"$(reversed_x ? "rev" : "")",
      marker
    )
  end
  
  if plot_all_prms
    full_prms_string = ""
    for (i, name) in enumerate(identifier_list)
      full_prms_string *=" $(name)=$(string(range_list[i]))"
    end
    THE_list = EEC_data_holder.bias[range_list[3]]
    if length(range_list[3])>1
      bias_step = round(THE_list[2] - THE_list[1], digits=5)
    else
      bias_step = 666
    end
    bias_aux_string = (length(THE_list) < 8 ? THE_list : "collect($(THE_list[1]) : $(bias_step) : $(THE_list[end]))")    
    title("TC=$(EEC_data_holder.TC[range_list[1]])    pO2=$(EEC_data_holder.pO2[range_list[2]])\nbias=$(bias_aux_string)    data_set=$(EEC_data_holder.data_set[range_list[4]]) ", fontsize=10)
    
#     title("TC=$(EEC_data_holder.TC[range_list[1]])    pO2=$(EEC_data_holder.pO2[range_list[2]])\nbias=$(bias_aux_string)    data_set=$(EEC_data_holder.data_set[range_list[4]]) ", fontsize=10)
  else
    title(y_name*"  vs  "*x_name)    
  end
  
  if plot_legend
    legend(loc="best")
  end

  PyPlot.show()
  return y_to_plot
end

function get_joint_EEC_data_via_bias(EEC_data_1, EEC_data_2)
  if (EEC_data_1.TC == EEC_data_2.TC &&
     EEC_data_1.pO2 == EEC_data_2.pO2 &&
     EEC_data_1.data_set == EEC_data_2.data_set &&
     EEC_data_1.prms_names == EEC_data_2.prms_names)
    
    EEC_data = Array{EEC_data_holder_struct}(undef, 2)
    EEC_data = [EEC_data_1, EEC_data_2]
    
    if EEC_data_1.bias[end] < EEC_data_2.bias[end]
      i_smaller = 1
      i_greater = 2
    else
      i_smaller = 2
      i_greater = 1
    end
    
    
    if EEC_data[i_smaller].bias[1] > EEC_data[i_smaller].bias[end]
      apply_reverse = true
    else
      apply_reverse = false
    end
    
    length_of_smaller = length(EEC_data[i_smaller].bias)
    
    # TODO ... the following condition should be MUCH !!! better
    if    EEC_data[i_smaller].bias[end] == EEC_data[i_greater].bias[1] ||
          EEC_data[i_smaller].bias[1] == EEC_data[i_greater].bias[1]
      start_of_greater = 2      
    else
      start_of_greater = 1
    end
    
    bias_OUT = [ (apply_reverse ?
                  EEC_data[i_smaller].bias[end : -1 : 1] :  
                  EEC_data[i_smaller].bias)... ,
                EEC_data[i_greater].bias[start_of_greater:end]...]
      
    EEC_data_OUT = EEC_data_holder_struct(EEC_data[1].TC, EEC_data[1].pO2, bias_OUT, EEC_data[1].data_set, EEC_data[1].prms_names)
      
    EEC_data_OUT.data[:,:, 1:length_of_smaller, :, :] = ( apply_reverse ? 
                                                          deepcopy(EEC_data[i_smaller].data[:, :, end : -1 : 1, :, :]) :
                                                          deepcopy(EEC_data[i_smaller].data) )
    EEC_data_OUT.data[:,:, length_of_smaller+1:end, :, :] = deepcopy(EEC_data[i_greater].data[:, :, start_of_greater:end, :, :])
    
    return EEC_data_OUT
  else
    prinltn("ERROR: meta_data are not the same!")
    return throw(Exception)
  end
end
  
  
function display_fit_vs_exp(EEC_data_holder;TC, pO2, bias, data_set, 
                            use_DRT=false, DRT_draw_semicircles=false, plot_legend=false,                             
                            error_type="normalized",
                            EEC_structure="R-L-RCPE-RCPE",
                            EIS_preprocessing_control = EIS_preprocessing_control()
                  )

  for (TC_idx, TC_item) in enumerate(TC), 
      (pO2_idx, pO2_item) in enumerate(pO2), 
      (bias_idx, bias_item) in enumerate(bias), 
      (data_set_idx, data_set_item) in enumerate(make_array_from_string(data_set))
    
    SIM_list = ysz_fitting.EIS_simulation(TC_item, pO2_item, bias_item, 
                    use_DRT=use_DRT, DRT_draw_semicircles=DRT_draw_semicircles, plot_option="Bode Nyq DRT RC", plot_legend=plot_legend, data_set=data_set_item)
    SIM = SIM_list[1]
    EIS_exp = ysz_fitting.import_data_to_DataFrame(SIM)    
    
    EIS_exp = EIS_preprocessing(EIS_exp, EIS_preprocessing_control)
    typical_plot_exp(SIM, EIS_exp)
    
    
    # EEC part
    function find_idx_in_list(item, list)      
      findall(x -> x==item, list)[1]
    end
    
    prms_values = deepcopy(EEC_data_holder.data[
          find_idx_in_list(TC_item, EEC_data_holder.TC), 
          find_idx_in_list(pO2_item, EEC_data_holder.pO2), 
          find_idx_in_list(bias_item, EEC_data_holder.bias), 
          find_idx_in_list(data_set_item, EEC_data_holder.data_set),
          : ])
    
    EEC_actual = get_EEC(EEC_structure)
    EEC_actual.prms_values = prms_values
    EIS_EEC = get_EIS_from_EEC(EEC_actual, f_range = EIS_exp.f)
    
    typical_plot_sim(SIM, EIS_EEC)
    
    println("--- ($(TC_item), $(pO2_item), $(bias_item), $(data_set)) --- ")
    @show prms_values
    println("fitnessFunction = $(fitnessFunction(EIS_simulation(), EIS_EEC, EIS_exp, error_type=error_type))")
    
  end  
end
  

function model_R1_from_TC(x, TC)
  return x[1]*exp((x[2]* e0)/(kB*TCtoT(TC)))
end

function model_TC_from_R1(x, R1)
  return TtoTC((x[2]/R)*(1/log(R1/x[1])))
end

function findfit_from_for_R(TC_data, R1_data)
  function to_optimize(x)
    error = 0
    for (i, TC) in enumerate(TC_data)
      error += (R1_data[i] - model_R1_from_TC(x, TC))^2
    end
    return sqrt(error)/length(TC_data)
  end
  
  
  
  fit_O = optimize(to_optimize, [1., 1.])
  
  R1_fitted_values = []
  TC_plot_range = 700 : 10 : 850
  for TC in TC_plot_range
    append!(R1_fitted_values, model_R1_from_TC(fit_O.minimizer, TC))
  end
  
  plot(TC_data, R1_data)
  plot(TC_plot_range, R1_fitted_values)
  return fit_O
end



#############
#############
#############

function get_computed_data_MONO_110()
  cdm = load_EEC_data_holder(folder="../data/EEC/", file_name="computed_data_minus_pO2.txt")
  cdp = load_EEC_data_holder(folder="../data/EEC/", file_name="computed_data_plus_pO2.txt")
  return get_joint_EEC_data_via_bias(cdm, cdp)
end

function test_R3_stealing()
  println("TODO!!!")
end

function get_bias_from_chrono(bias_idx, bias_step, bias_return_value; num_cycles=1, take_each=1)
    bias_digit_precision = 3
    bias_epsilon = 10e-5
    
    bias_step = 0.05
    dir = 1
    bias_value = 0.0
        
    if take_each==2    
      bias_idx =   bias_idx/2 + 0.5*mod(bias_idx, 2) - div(bias_idx-1, 82)
    else
      bias_idx = bias_idx    
    end
    
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
      if return_counter > (num_cycles*2 - 1) && bias_value > 0        
        bias_value = Inf
        break
      end
    end
    bias_value = round(bias_value, digits=bias_digit_precision)
end
