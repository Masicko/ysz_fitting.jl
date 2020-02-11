module ysz_fitting
#######################
# TODO 
# [!] general global search using projection to each variable
# [x] compute EIS exacly on checknodes and therefore remove plenty of "EIS_apply_checknodes"
# [ ] put appropriate and finished stuff into "CV_fitting_supporting_stuff"
# [x] implement LM algorithm
# [ ] sjednotit, co znamena pO2, jeslti jsou to procenta a jak se prenasi do simulace!  a taky T .. Celsia a Kelvina !!!
# [ ] srovnat vodivost elektrolytu s experimentem CV i EIS naraz
# [x] vymyslet novy relevantni vektor parametru
# [ ] aplikovat masku pro fitting pro scan_2D_recursive
# [o] snehurka, velke objemy dat, maska a metafile
# ->[x] save_dir preposilany parametrem a odlisit scripted_dir
# ->[ ] sneh - zadavani i jinych parametru nez prms
# ->[ ] sneh - pouzivat jen jeden soubor pro spouteni sbatch ... v hlavicce #!(..)/julia
# ->[ ] mozna pouzit precompile
# [ ] postarat se o modularitu ysz_model_COSI, at se nemusi vytvaret znovu "experiments_COSI" -> JUERGEN
#######################


using Printf
using PyPlot
using DataFrames
using LeastSquaresOptim

include("../examples/ysz_experiments.jl")
include("../src/CV_fitting_supporting_stuff.jl")
include("../src/EIS_fitting_supporting_stuff.jl")
include("../src/import_experimental_data.jl")


###########
###########


function CV_get_shared_checknodes()
    return get_checknodes(0.06,0.95,-0.95,-0.06,0.12)
end

function EIS_get_shared_omega_range()
    # omegas = (w0, w1, w_fac)
    return omega_range = (1, 7500, 1.2)
end

function EIS_get_shared_checknodes()
    return EIS_get_checknodes_geometrical(EIS_get_shared_omega_range()[1], EIS_get_shared_omega_range()[2], EIS_get_shared_omega_range()[3])
end

function get_shared_prms() 
    #     old_prms =  [21.71975544711280, 20.606423236896422, 0.0905748, -0.708014, 0.6074566741435283, 0.1]
    # first dump fit = (A0, R0, DGA, DGR, beta, A) = (18.8, 19.2,     -0.14,      -0.5,   0.6074566741435283,   0.956)
    return (A0, R0, DGA, DGR, beta, A) = (18.8, 19.2,     -0.14,      -0.5,   0.6074566741435283,   0.956)
end

function get_shared_add_prms()
    # cap_fit + resistance fit = (DD, nu, nus, ms_par) = (4.35e-13,    0.85,      0.21,   0.05)
    return (DD, nu, nus, ms_par) = (4.35e-13,    0.85,      0.21,   0.05)
end


function filename_format_prms(; save_dir="./nouze/", prefix="", prms=Nothing, prms_names=("A0", "R0", "DGA", "DGR", "beta", "A"), scripted_tuple)

  function consistency_check()
    if (size(prms_names,1) != size(scripted_tuple,1) || 
        size(prms_names,1) != size(prms,1))
      return false
    end
    return true
  end
  
  if !consistency_check()
    println("ERROR: file_format_prms(): shape mismatch (!consistency_check())")
    return throw(Exception)
  end

  scripted_dir = ""
  out_name = prefix
  for i in 1:size(scripted_tuple, 1)
    if scripted_tuple[i] == 1
      scripted_dir = string(scripted_dir,"_",prms_names[i],@sprintf("%.2f",prms[i]))
    end
    out_name = string(out_name, "_", prms_names[i], @sprintf("%.2f",prms[i]))
  end
  
  out_path = string(save_dir,scripted_dir,"/")
  out_name = string(out_name, ".csv")

  (out_path, out_name)
end








##########
##########



function CV_apply_checknodes(CV_in, checknodes)
    DataFrame( U = checknodes[:,1], I = CV_get_I_values(CV_in, checknodes))
end

function EIS_apply_checknodes(EIS_in,checknodes)
    DataFrame( f = checknodes, Z = EIS_get_Z_values(EIS_in, checknodes))
end

function CV_load_file_prms(;save_dir, prms, prms_names=("A0", "R0", "DGA", "DGR", "beta", "A"), scripted_tuple=(0, 0, 0, 0, 0, 0))
    (out_path, out_name) = filename_format_prms( save_dir=save_dir, prefix="CV", prms=prms, prms_names=prms_names, scripted_tuple=scripted_tuple)

    CSV.read(
        string(out_path,out_name),
    )
    return
end

function CV_save_file_prms(df_out, save_dir, prms, prms_names, scripted_tuple)
    (out_path, out_name) = filename_format_prms( save_dir=save_dir, prefix="CV", prms=prms, prms_names=prms_names, scripted_tuple=scripted_tuple)
    run(`mkdir -p $out_path`)

    CSV.write(
        string(out_path,out_name),
        df_out
    )
    return
end

function EIS_load_file_prms(;save_dir, prms, prms_names=("A0", "R0", "DGA", "DGR", "beta", "A"), scripted_tuple=(0, 0, 0, 0, 0, 0))
    (out_path, out_name) = filename_format_prms( save_dir=save_dir, prefix="EIS", prms=prms, prms_names=prms_names, scripted_tuple=scripted_tuple)
    
    df = CSV.read(string(out_path, out_name))
    DataFrame(f = df.f, Z = (df.Re .+ df.Im.*im))
end

function EIS_save_file_prms(df_out, save_dir, prms, prms_names, scripted_tuple)
    (out_path, out_name) = filename_format_prms( save_dir=save_dir, prefix="EIS", prms=prms, prms_names=prms_names, scripted_tuple=scripted_tuple)
    run(`mkdir -p $out_path`)

    CSV.write(
        string(out_path,out_name),
        DataFrame(f = df_out.f, Re = real(df_out.Z), Im = imag(df_out.Z))
    )
    return
end

function CV_simple_run(;pyplot=false)
    old_prms = [21.71975544711280, 20.606423236896422, 0.0905748, -0.708014, 0.6074566741435283, 0.1]
    
    (A0, R0, DGA, DGR, beta, A) = get_shared_prms()
    (A0, R0, DGA, DGR, beta, A) =[23.0,   23.,    0.4,    0.4,    0.4,    0.8]
    
    (DD, nu, nus, ms_par) = get_shared_add_prms()
    
   CV_sim =  ysz_experiments.run_new(
        out_df_bool=true, voltammetry=true, sample=8, pyplot_finall=true,
        prms_in= (A0, R0, DGA, DGR, beta, A),
        add_prms_in=(DD, nu, nus, ms_par)
    )
    @show CV_sim
    checknodes =  CV_get_shared_checknodes()
    if pyplot
        figure(1)
        CV_plot(CV_apply_checknodes(CV_sim,checknodes), "s_r")
        CV_plot(CV_apply_checknodes(import_CVtoDataFrame(T=800, pO2=100),CV_get_shared_checknodes()), "exp")
    end    
    return
end

function EIS_simple_run(;pyplot=false)    
    
#     pO2_in_sim = 1.0
#     EIS_bias=-0.5
    
    (A0, R0, DGA, DGR, beta, A) = get_shared_prms()
    (A0, R0, DGA, DGR, beta, A) =[29.620786013494122,   25.34421543847122,    0.4,    0.4,    0.43268062905804937,    0.21181634115352904]
    
    (DD, nu, nus, ms_par) = get_shared_add_prms()
    
    EIS_df = ysz_experiments.run_new(
        pyplot=pyplot, EIS_IS=true, out_df_bool=true, EIS_bias=EIS_bias, omega_range=EIS_get_shared_omega_range(),
        dx_exp=-8,
        # TODO !!! T a pO2
        prms_in=[A0, R0, DGA, DGR, beta, A],
        add_prms_in=(DD, nu, nus, ms_par)
    )
    #@show EIS_df
    checknodes =  EIS_get_shared_checknodes()
    if pyplot
        figure(2)
        Nyquist_plot(EIS_apply_checknodes(EIS_df,checknodes), "s_r")
        Nyquist_plot(EIS_apply_checknodes(import_EIStoDataFrame(T=800, pO2=100, bias=0.0),EIS_get_shared_checknodes()), "exp")
    end
end

# # # function get_dfs_to_interpolate(
# # #                             EIS_bool=true, pyplot=false,
# # #                            T=800, pO2=100, EIS_bias=0.0,
# # #                             #rprm1=range(15, stop=20, length=3),
# # #                             #rprm2=range(15, stop=20, length=3),
# # #                             #rprm1=range(18.4375, stop=18.59375, length=3),
# # #                             #rprm2=range(20.59375, stop=21, length=3),
# # #                             rprm1=range(19.3, stop=21.0, length=2),
# # #                             #rprm2=range(15.0, stop=23, length=7),
# # #                             rprm2=range(-3.0, stop=-1.6, length=2),
# # #                             nx=100, ny=100,
# # #                             wp=[1, 1, 0, 0, 0, 0], #TODO !!!
# # #                             depth_remaining=1,
# # #                             approx_min=Inf,
# # #                             approx_min_position=[0, 0],
# # #                             recursive=false,
# # #                             recursive_string="",
# # #                             )
# # # 
# # #     println("nx = ",nx, " ..... ny = ",ny)
# # #     #PyPlot.close()
# # #     #PyPlot.close()
# # #     
# # #     ######################################################
# # #     #rprm1 = collect(-5.0 : 1.0 : 5.0)
# # #     #rprm1 = collect(-3. : 0.5 : 3.0)
# # #     #rprm1=range(15, stop=20, length=3)
# # #     
# # #     #rprm2 = collect(-3. : 0.5 : 3.0)
# # #     #rprm2=range(15, stop=20, length=3)
# # #     ######################################################
# # #     
# # #     wpn = ["A0","R0","DGA","DGR","beta","A"]
# # #     #A0 = 21.71975544711280
# # #     #R0 = 19.53
# # #     A0 = 19.50
# # #     R0 = 19.85
# # #     DGA = 0.0905748
# # #     DGR = -0.708014
# # #     beta = 0.6074566741435283
# # #     A = -0.356
# # #     
# # #     (A0, R0, DGA, DGR, beta, A) = (21.71975544711280, 20.606423236896422, 0.0905748, -0.708014, 0.6074566741435283, 0.1)
# # # 
# # #     println(recursive_string,"computing 2D scan... ")
# # #     lrprm1 = size(rprm1,1)
# # #     lrprm2 = size(rprm2,1)
# # #     
# # #     global_min = Inf
# # #     global_min_position = [0,0]
# # #     
# # #     CV_matrix = Array{DataFrame}(undef,lrprm1, lrprm2)
# # #     ok_matrix = Array{Int32}(undef,lrprm1, lrprm2)
# # #     if EIS_bool
# # #         println(recursive_string,"Computing EISs")
# # #     else
# # #         println(recursive_string,"Computing CVs")
# # #     end
# # #     @time(
# # #         for i1 in 1:lrprm1
# # #             for i2 in 1:lrprm2
# # #                 try
# # #                     print(recursive_string,rprm1[i1], " / ", rprm2[i2])
# # #                     
# # #                     ##########################################
# # #                     R0 = rprm1[i1]
# # #                     A = rprm2[i2]
# # #                     ##########################################
# # #                     
# # #                    if EIS_bool 
# # #                         CV_matrix[i1,i2] = ysz_experiments.run_new(
# # #                             pyplot=false, EIS_IS=true, out_df_bool=true, 
# # #                             tref=0,
# # #                             prms_in=[A0, R0, DGA, DGR, beta, A],  
# # #                             nu_in=0.9, pO2_in=1.0, DD_in=9.5658146540360312e-10
# # #                         )
# # #                     else
# # #                         
# # #                         CV_matrix[i1,i2] = ysz_experiments.run_new(
# # #                             out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
# # #                             prms_in=[A0, R0, DGA, DGR, beta, A]
# # #                         )
# # #                     end
# # #                     ok_matrix[i1, i2] = 1
# # #                     println("  .. ok :)")
# # #                 catch e
# # #                     if e isa InterruptException
# # #                         rethrow(e)
# # #                     else
# # #                         ok_matrix[i1, i2] = 0
# # #                         println("  :(   ")
# # #                         continue
# # #                     end
# # #                 end
# # #             end
# # #         end
# # #     )
# # #     
# # #     #PyPlot.figure(figsize=(6,4))
# # #     #xlabel("parameter")
# # #     #zlabel("error")
# # #     #title("Fitness Function")
# # # 
# # #     if pyplot && !(recursive)
# # #         figure(figsize=(3,3))
# # #     end
# # #     #print_ok_matrix(ok_matrix)
# # #     
# # #     println(recursive_string,"Preparing for interpolation ... ")
# # #         for i1 in 1:lrprm1-1
# # #             for i2 in 1:lrprm2-1
# # #                 if all_four_nodes_ok(ok_matrix,i1,i2)
# # #                     print(recursive_string," intrp i1 / i2 : ",i1, " / ", i2, " for prms: ", rprm1[i1], " / ", rprm2[i2])
# # #                     
# # #                     if EIS_bool
# # #                             EIS_list = [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]]
# # #                             EIS_with_checknodes_list = []
# # #                             checknodes =  EIS_get_shared_checknodes()
# # #                             for i in 1:size(EIS_list,1)
# # #                                 push!(EIS_with_checknodes_list, DataFrame(f = checknodes[:], Z = EIS_get_Z_values(EIS_list[i], checknodes)))
# # #                             end
# # #                             return EIS_with_checknodes_list
# # #                     else
# # #                         # TODO !!!
# # #                         FF=CV_get_FF_interp_2D(
# # #                             [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
# # #                             rprm1[i1:i1+1],
# # #                             rprm2[i2:i2+1];
# # #                             nx=nx, ny=ny,
# # #                             T=T, pO2=pO2
# # #                         )
# # #                     end
# # #                 end 
# # #             end
# # #         end
# # # end

function EIS_view_exp(T=800, pO2=100, EIS_bias=0.0,)
    EIS_raw = import_EIStoDataFrame(;T=T,pO2=pO2,bias=EIS_bias)
    checknodes =  EIS_get_shared_checknodes()
    EIS_exp = DataFrame(f = checknodes[:], Z = EIS_get_Z_values(EIS_raw, checknodes))
    
    Nyquist_plot(EIS_exp,"exp")
end

function EIS_view_interpolation_at(Q_list, x, y)
    EIS_intrp = DataFrame(
        f = Q_list[1].f,
        Z = EIS_biliComb(
            Q_list[1],
            Q_list[2],
            Q_list[3],
            Q_list[4],
            x,y).Z
    )
    Nyquist_plot(EIS_intrp, string("x/y = ",x,"/",y))
end

function plot_all(df_list, prms_list)
    PyPlot.figure(figsize=(6,4))
    #PyPlot.figure(figsize=(10,8))
    for i in 1:size(df_list,1)
        #println(i)
        CV_plot(
            df_list[i],
            string("prms = ",prms_list[i])
        )
    end
    
    CV_orig = import_CVtoDataFrame(;T=800,pO2=100)
    #checknodes = get_checknodes(0.01,0.95,-0.95,-0.01,0.05)
    checknodes = CV_get_shared_checknodes()
    CV_exp = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_orig, checknodes))
    
    CV_plot(CV_exp, "exp")
end

function plot_CV_matrix(CV_matrix, ok_matrix, rprm1, rprm2, T, pO2)
    PyPlot.figure(figsize=(6,4))
    
    
    for i1 in 1:size(rprm1,1)
        for i2 in 1:size(rprm2,1)
            if Bool(ok_matrix[i1,i2])
                #println(" ploting i1 / i2 : ",i1, " / ", i2)
                CV_plot(
                    CV_matrix[i1,i2],
                    string("prms = ", rprm1[i1], " / ", rprm2[i2])
                )
            end
        end
    end

    CV_orig = import_CVtoDataFrame(;T=T,pO2=pO2)
    #checknodes = get_checknodes(0.01,0.95,-0.95,-0.01,0.05)
    checknodes = CV_get_shared_checknodes()
    CV_exp = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_orig, checknodes))
    
    CV_plot(CV_exp, "exp")
end



function plot_EIS_matrix(EIS_matrix, ok_matrix, rprm1, rprm2, T, pO2, EIS_bias=0.0)
    PyPlot.figure(figsize=(6,4))
    
    checknodes =  EIS_get_shared_checknodes()  
    
    for i1 in 1:size(rprm1,1)
        for i2 in 1:size(rprm2,1)
            if Bool(ok_matrix[i1,i2])
                #println(" ploting i1 / i2 : ",i1, " / ", i2)
                Nyquist_plot(
                    EIS_apply_checknodes(EIS_matrix[i1, i2], checknodes),
                    string("prms = ", rprm1[i1], " / ", rprm2[i2])
                )
            end
        end
    end

    EIS_raw = import_EIStoDataFrame(;T=T,pO2=pO2,bias=EIS_bias)
    EIS_exp = DataFrame(f = checknodes[:], Z = EIS_get_Z_values(EIS_raw, checknodes))
    
    Nyquist_plot(EIS_exp, "exp")
end


#####################
## 2D optimization ##
#####################

mutable struct FF_2D
    x_matrix
    y_matrix
    err_matrix
end



function surfplot_FF_2D(FF::FF_2D; my_title="Objective function")
    surf(FF.x_matrix, FF.y_matrix, FF.err_matrix, cstride=1, #cmap=ColorMap("gray"), 
        alpha=0.8);
    xlabel("prm1")
    ylabel("prm2")
    title(my_title)
    zlabel("error")
end


function CV_get_FF_interp_2D(CV_list, prm1_list, prm2_list; nx=20, ny=20, T=800, pO2=100)
    # ordered for prms A, B as [(A1,B1), (A1,B2), (A2,B1), (A2, B2)]
    # prms_list ordered as [A1 B1 A2 B2]
    CV_raw = import_CVtoDataFrame(;T=T,pO2=pO2)
    checknodes = get_checknodes(0.05,0.95,-0.95,-0.05,0.01)
    CV_exp = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_raw, checknodes))
    
    CV_with_checknodes_list = []
    
    for i in 1:size(CV_list,1)
        push!(CV_with_checknodes_list, DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_list[i], checknodes)))
    end
    
    # bilinear interpolation
    rx = range(0.0, stop=1.0, length=nx)
    ry = range(0.0, stop=1.0, length=ny)
    
    x_matrix = zeros(nx,ny)
    y_matrix = zeros(nx,ny)
    err_matrix = zeros(nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            x_matrix[ix, iy] = rx[ix]*(prm1_list[2]-prm1_list[1]).+prm1_list[1]
            y_matrix[ix, iy] = ry[iy]*(prm2_list[2]-prm2_list[1]).+prm2_list[1]
            err_matrix[ix, iy] = fitnessFunction(
                    CV_exp,
                    DataFrame(
                        U = checknodes[:,1], 
                        I = biliComb_CVs(
                            CV_with_checknodes_list[1],
                            CV_with_checknodes_list[2],
                            CV_with_checknodes_list[3],
                            CV_with_checknodes_list[4],
                            rx[ix],ry[iy]).I
                    )
                )
        end
    end
    
    return FF_2D(x_matrix, y_matrix, err_matrix)
end



function EIS_get_FF_interp_2D(EIS_list, prm1_list, prm2_list; nx=20, ny=20, T=800, pO2=100, EIS_bias=0.0)
    # ordered for prms A, B as [(A1,B1), (A1,B2), (A2,B1), (A2, B2)]
    # prms_list ordered as [A1 B1 A2 B2]
    EIS_raw = import_EIStoDataFrame(;T=T,pO2=pO2,bias=EIS_bias)
    checknodes =  EIS_get_shared_checknodes()
    EIS_exp = DataFrame(f = checknodes[:], Z = EIS_get_Z_values(EIS_raw, checknodes))
    
    EIS_with_checknodes_list = []
    
    for i in 1:size(EIS_list,1)
        push!(EIS_with_checknodes_list, DataFrame(f = checknodes[:], Z = EIS_get_Z_values(EIS_list[i], checknodes)))
    end
    
    # bilinear interpolation
    rx = range(0.0, stop=1.0, length=nx)
    ry = range(0.0, stop=1.0, length=ny)
    
    x_matrix = zeros(nx,ny)
    y_matrix = zeros(nx,ny)
    err_matrix = zeros(nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            x_matrix[ix, iy] = rx[ix]*(prm1_list[2]-prm1_list[1]).+prm1_list[1]
            y_matrix[ix, iy] = ry[iy]*(prm2_list[2]-prm2_list[1]).+prm2_list[1]
            err_matrix[ix, iy] = EIS_fitnessFunction(
                    EIS_exp,
                    DataFrame(
                        f = checknodes[:], 
                        Z = EIS_biliComb(
                            EIS_with_checknodes_list[1],
                            EIS_with_checknodes_list[2],
                            EIS_with_checknodes_list[3],
                            EIS_with_checknodes_list[4],
                            rx[ix],ry[iy]).Z
                    )
                )
        end
    end
    
    return FF_2D(x_matrix, y_matrix, err_matrix)
end


function all_four_nodes_ok(ok_matrix, x, y)
    for i in 0:1
        for j in 0:1
            if ok_matrix[x+i, y + j] == 0
                return false
            end
        end
    end
    return true
end

function print_ok_matrix(ok_matrix)
println("... printing ok_matrix ...")
    for i2 in 1:size(ok_matrix,2)
        for i1 in 1:size(ok_matrix,1)
        
            print("|", ok_matrix[i1, i2])
        end
        print("|\n")
    end
end

function find_mins(FF, dx; debug_print_bool=false)
    hard_min = Inf
    count_of_mins = 0
    values = []
    positions = []
    curv = []
    
    if size(FF.err_matrix,1) < 3 || size(FF.err_matrix,2) < 3
        println("ERROR: is_min_there: some dimension of FF.err_matrix is smaller than 3")
        return throw(Exception)
    end
    
    for i1 in 2:size(FF.err_matrix,1)-1
        for i2 in 2:size(FF.err_matrix,2)-1
            pot_min = FF.err_matrix[i1, i2]
            if pot_min < hard_min
                hard_min = pot_min
            end
            
            if debug_print_bool
                println(i1, " ", i2, " > ",pot_min)
                
                println(pot_min < FF.err_matrix[i1-1, i2], " ",
                    pot_min < FF.err_matrix[i1+1, i2], " ",
                    pot_min < FF.err_matrix[i1, i2 - 1], " ",
                    pot_min < FF.err_matrix[i1, i2 + 1])
            end
                
            if (pot_min < FF.err_matrix[i1 - 1, i2 - 1] &&
                pot_min < FF.err_matrix[i1 - 1, i2] &&
                pot_min < FF.err_matrix[i1 - 1, i2 + 1] &&
                pot_min < FF.err_matrix[i1, i2 - 1] &&
                pot_min < FF.err_matrix[i1, i2 + 1] &&
                pot_min < FF.err_matrix[i1 + 1, i2 - 1] &&
                pot_min < FF.err_matrix[i1 + 1, i2] &&
                pot_min < FF.err_matrix[i1 + 1, i2 + 1]
                )
                
                c1 = (FF.err_matrix[i1 - 1, i2] - 2*pot_min + FF.err_matrix[i1 + 1, i2])/(dx[1]*dx[1])
                c2 = (FF.err_matrix[i1, i2 - 1] - 2*pot_min + FF.err_matrix[i1, i2 + 1])/(dx[2]*dx[2])
                c12= (FF.err_matrix[i1 - 1, i2 - 1] + FF.err_matrix[i1 + 1, i2 + 1] 
                    - FF.err_matrix[i1 + 1, i2 - 1] - FF.err_matrix[i1 - 1, i2 + 1])/(4*dx[1]*dx[2])
                count_of_mins += 1
                append!(values,pot_min)
                push!(positions, [FF.x_matrix[i1, i2], FF.y_matrix[i1, i2]])
                push!(curv,[c1, c2,c12])
            end
        end
    end
    
    return count_of_mins, values, positions, curv, hard_min
end

function pickup_min(min_count, min_values, min_positions)
    value = Inf;
    position = [0, 0]
    for i in 1:min_count
        if min_values[i] < value
            value = min_values[i]
            position = min_positions[i]
        end
    end
    if value == Inf
        println("ERROR: pickup_min: value = -1")
        throw(Exception)
    else
        return value, position
    end
end

function scan_2D_recursive(;pyplot=false, 
                            EIS_bool=true,
                            T=800, pO2=100, EIS_bias=0.0,
                            
                            wp=[1, 1, 0, 0, 0, 0], #TODO !!!
                            rprm1=range(0.1, stop=0.9, length=10),
                           #rprm2=range(15.0, stop=23, length=7),
                            rprm2=range(0.8, stop=1.0, length=4),
                            
                            nx=30, ny=30,

                            depth_remaining=1,
                            approx_min=Inf,
                            approx_min_position=[0, 0],
                            recursive=false,
                            recursive_string="",
                            )

    println("nx = ",nx, " ..... ny = ",ny)
    #PyPlot.close()
    #PyPlot.close()
    
    ######################################################
    #rprm1 = collect(-5.0 : 1.0 : 5.0)
    #rprm1 = collect(-3. : 0.5 : 3.0)
    #rprm1=range(15, stop=20, length=3)
    
    #rprm2 = collect(-3. : 0.5 : 3.0)
    #rprm2=range(15, stop=20, length=3)
    ######################################################
    
    wpn = ["A0","R0","DGA","DGR","beta","A"]

    (A0, R0, DGA, DGR, beta, A) = get_shared_prms()
   (DD, nu, nus, ms_par)  =  get_shared_add_prms()
    
   
    println(recursive_string,"computing 2D scan... ")
    lrprm1 = size(rprm1,1)
    lrprm2 = size(rprm2,1)
    
    global_min = Inf
    global_min_position = [0,0]
    
    CV_matrix = Array{DataFrame}(undef,lrprm1, lrprm2)
    ok_matrix = Array{Int32}(undef,lrprm1, lrprm2)
    if EIS_bool
        println(recursive_string,"Computing EISs")
    else
        println(recursive_string,"Computing CVs")
    end
    @time(
        for i1 in 1:lrprm1
            for i2 in 1:lrprm2
                try
                    print(recursive_string,rprm1[i1], " / ", rprm2[i2])
                    
                    ##########################################
                    beta = rprm1[i1]
                    A = rprm2[i2]
                    ##########################################
                    
                   if EIS_bool
                        CV_matrix[i1,i2] = ysz_experiments.run_new(
                            pyplot=false, EIS_IS=true, out_df_bool=true, EIS_bias=EIS_bias, omega_range=EIS_get_shared_omega_range(),
                            dx_exp=-8,
                            # TODO !!! T a pO2
                            prms_in=[A0, R0, DGA, DGR, beta, A],
                            add_prms_in=(DD, nu, nus, ms_par)
                        )
                    else
                        CV_matrix[i1,i2] = ysz_experiments.run_new(
                            out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
                            prms_in=[A0, R0, DGA, DGR, beta, A],
                            add_prms_in=(DD, nu, nus, ms_par)
                        )
                    end
                    ok_matrix[i1, i2] = 1
                    println("  .. ok :)")
                catch e
                    if e isa InterruptException
                        rethrow(e)
                    else
                        ok_matrix[i1, i2] = 0
                        println("  :(   ")
                        continue
                    end
                end
            end
        end
    )
    
    #PyPlot.figure(figsize=(6,4))
    #xlabel("parameter")
    #zlabel("error")
    #title("Fitness Function")

    if pyplot && !(recursive)
        figure(figsize=(3,3))
    end
    #print_ok_matrix(ok_matrix)
    
    println(recursive_string,"Computing interpolation ... ")
    @time(
        for i1 in 1:lrprm1-1
            for i2 in 1:lrprm2-1
                if all_four_nodes_ok(ok_matrix,i1,i2)
                    print(recursive_string," intrp i1 / i2 : ",i1, " / ", i2, " for prms: ", rprm1[i1], " / ", rprm2[i2])
                    
                    if EIS_bool
                         FF=EIS_get_FF_interp_2D(
                            [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
                            rprm1[i1:i1+1],
                            rprm2[i2:i2+1];
                            nx=nx, ny=ny,
                            T=T, pO2=pO2, EIS_bias=EIS_bias
                        )                       
                    else
                        FF=CV_get_FF_interp_2D(
                            [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
                            rprm1[i1:i1+1],
                            rprm2[i2:i2+1];
                            nx=nx, ny=ny,
                            T=T, pO2=pO2
                        )
                    end
                    
                    min_count, min_values, min_positions, min_curv, hard_min = find_mins(
                        FF, 
                        [(rprm1[i1+1] - rprm1[i1])/nx,
                         (rprm2[i2+1] - rprm2[i2])/ny]
                    )
                    print("... h_min = ",hard_min)
                    
                    
                    if min_count > 0
                        println(" ... min_c ", min_count, "\n    val ", min_values, "\n    pos ", min_positions, "\n    curv ", min_curv)
                        #println(" ... min ", min_count)#
                        
                        approx_min, approx_min_position = pickup_min(min_count, min_values, min_positions)
                        
                        if depth_remaining > 0
                            println(recursive_string,">> Going deeper in intrp i1 / i2 : ",i1, " / ", i2)
                            new_min, new_min_position = scan_2D_recursive(;pyplot=true, 
                                T=T, pO2=pO2,
                                rprm1=range(rprm1[i1], stop=rprm1[i1+1], length=3),
                                rprm2=range(rprm2[i2], stop=rprm2[i2+1], length=3),
                                approx_min=approx_min,
                                approx_min_position=approx_min_position,
                                recursive=true,
                                depth_remaining=depth_remaining - 1,
                                recursive_string=string(recursive_string, "<", depth_remaining - 1, ">"),
                                nx = Int32(round(nx/1.2)), ny = Int32(round(ny/1.2))
                            )
                            if new_min < global_min
                                global_min = new_min
                                global_min_position =new_min_position
                            end
                        else
                            if approx_min < global_min
                                global_min = approx_min
                                global_min_position = approx_min_position
                            end
                            if pyplot
                                surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
                            end
                        end
                    else
                        if pyplot
                            surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
                        end
                        println()
                    end
                end 
            end
        end
    )
    
    
    
    if pyplot && !(recursive)
        if EIS_bool 
            plot_EIS_matrix(CV_matrix, ok_matrix, rprm1, rprm2, T, pO2, EIS_bias)
        else
            plot_CV_matrix(CV_matrix, ok_matrix, rprm1, rprm2, T, pO2)
        end
    end
    return global_min, global_min_position
end


# # # # function scan_2D_recursive(;pyplot=false, 
# # # #                             T=800, pO2=100,
# # # #                             #rprm1=range(15, stop=20, length=3),
# # # #                             #rprm2=range(15, stop=20, length=3),
# # # #                             #rprm1=range(18.4375, stop=18.59375, length=3),
# # # #                             #rprm2=range(20.59375, stop=21, length=3),
# # # #                             rprm1=range(10.0, stop=28., length=4),
# # # #                             #rprm2=range(15.0, stop=23, length=7),
# # # #                             rprm2=range(-2.5, stop=1.0, length=4),
# # # #                             nx=100, ny=100,
# # # #                             wp=[1, 1, 0, 0, 0, 0], #TODO !!!
# # # #                             depth_remaining=1,
# # # #                             approx_min=Inf,
# # # #                             approx_min_position=[0, 0],
# # # #                             recursive=false,
# # # #                             recursive_string="",
# # # #                             )
# # # # 
# # # #     println("nx = ",nx, " ..... ny = ",ny)
# # # #     #PyPlot.close()
# # # #     #PyPlot.close()
# # # #     
# # # #     ######################################################
# # # #     #rprm1 = collect(-5.0 : 1.0 : 5.0)
# # # #     #rprm1 = collect(-3. : 0.5 : 3.0)
# # # #     #rprm1=range(15, stop=20, length=3)
# # # #     
# # # #     #rprm2 = collect(-3. : 0.5 : 3.0)
# # # #     #rprm2=range(15, stop=20, length=3)
# # # #     ######################################################
# # # #     
# # # #     wpn = ["A0","R0","DGA","DGR","beta","A"]
# # # #     #A0 = 21.71975544711280
# # # #     #R0 = 19.53
# # # #     A0 = 19.50
# # # #     R0 = 19.85
# # # #     DGA = 0.0905748
# # # #     DGR = -0.708014
# # # #     beta = 0.6074566741435283
# # # #     A = -0.356
# # # # 
# # # #     println(recursive_string,"computing 2D scan... ")
# # # #     lrprm1 = size(rprm1,1)
# # # #     lrprm2 = size(rprm2,1)
# # # #     
# # # #     global_min = Inf
# # # #     global_min_position = [0,0]
# # # #     
# # # #     CV_matrix = Array{DataFrame}(undef,lrprm1, lrprm2)
# # # #     ok_matrix = Array{Int32}(undef,lrprm1, lrprm2)
# # # #     println(recursive_string,"Computing CVs")
# # # #     @time(
# # # #         for i1 in 1:lrprm1
# # # #             for i2 in 1:lrprm2
# # # #                 try
# # # #                     print(recursive_string,rprm1[i1], " / ", rprm2[i2])
# # # #                     
# # # #                     ##########################################
# # # #                     R0 = rprm1[i1]
# # # #                     A = rprm2[i2]
# # # #                     ##########################################
# # # #                     
# # # #                    if true
# # # #                         CV_matrix[i1,i2] = ysz_experiments.run_new(
# # # #                             out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
# # # #                             prms_in=[A0, R0, DGA, DGR, beta, A]
# # # #                         )
# # # #                     else
# # # #                         CV_matrix[i1,i2] = ysz_experiments.run_new(
# # # #                             pyplot=false, EIS_IS=true, out_df_bool=false, 
# # # #                             tref=0,
# # # #                             prms_in=[A0, R0, DGA, DGR, beta, A],  
# # # #                             nu_in=0.9, pO2_in=1.0, DD_in=9.5658146540360312e-10
# # # #                         )
# # # #                     end
# # # #                     ok_matrix[i1, i2] = 1
# # # #                     println("  .. ok :)")
# # # #                 catch e
# # # #                     if e isa InterruptException
# # # #                         rethrow(e)
# # # #                     else
# # # #                         ok_matrix[i1, i2] = 0
# # # #                         println("  :(   ")
# # # #                         continue
# # # #                     end
# # # #                 end
# # # #             end
# # # #         end
# # # #     )
# # # #     
# # # #     #PyPlot.figure(figsize=(6,4))
# # # #     #xlabel("parameter")
# # # #     #zlabel("error")
# # # #     #title("Fitness Function")
# # # # 
# # # #     if pyplot && !(recursive)
# # # #         figure(figsize=(3,3))
# # # #     end
# # # #     #print_ok_matrix(ok_matrix)
# # # #     
# # # #     println(recursive_string,"Computing interpolation ... ")
# # # #     @time(
# # # #         for i1 in 1:lrprm1-1
# # # #             for i2 in 1:lrprm2-1
# # # #                 if all_four_nodes_ok(ok_matrix,i1,i2)
# # # #                     print(recursive_string," intrp i1 / i2 : ",i1, " / ", i2, " for prms: ", rprm1[i1], " / ", rprm2[i2])
# # # #                     
# # # #                     FF=get_FF_interp_2D(
# # # #                         [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
# # # #                         rprm1[i1:i1+1],
# # # #                         rprm2[i2:i2+1];
# # # #                         nx=nx, ny=ny,
# # # #                         T=T, pO2=pO2
# # # #                     )
# # # #                     
# # # #                     min_count, min_values, min_positions, min_curv, hard_min = find_mins(
# # # #                         FF, 
# # # #                         [(rprm1[i1+1] - rprm1[i1])/nx,
# # # #                          (rprm1[i2+1] - rprm1[i2])/ny]
# # # #                     )
# # # #                     print("... h_min = ",hard_min)
# # # #                     
# # # #                     
# # # #                     if min_count > 0
# # # #                         println(" ... min_c ", min_count, "\n    val ", min_values, "\n    pos ", min_positions, "\n    curv ", min_curv)
# # # #                         #println(" ... min ", min_count)#
# # # #                         
# # # #                         approx_min, approx_min_position = pickup_min(min_count, min_values, min_positions)
# # # #                         
# # # #                         if depth_remaining > 0
# # # #                             println(recursive_string,">> Going deeper in intrp i1 / i2 : ",i1, " / ", i2)
# # # #                             new_min, new_min_position = scan_2D_recursive(;pyplot=true, 
# # # #                                 T=T, pO2=pO2,
# # # #                                 rprm1=range(rprm1[i1], stop=rprm1[i1+1], length=3),
# # # #                                 rprm2=range(rprm2[i2], stop=rprm2[i2+1], length=3),
# # # #                                 approx_min=approx_min,
# # # #                                 approx_min_position=approx_min_position,
# # # #                                 recursive=true,
# # # #                                 depth_remaining=depth_remaining - 1,
# # # #                                 recursive_string=string(recursive_string, "<", depth_remaining - 1, ">"),
# # # #                                 nx = Int32(round(nx/1.2)), ny = Int32(round(ny/1.2))
# # # #                             )
# # # #                             if new_min < global_min
# # # #                                 global_min = new_min
# # # #                                 global_min_position =new_min_position
# # # #                             end
# # # #                         else
# # # #                             if approx_min < global_min
# # # #                                 global_min = approx_min
# # # #                                 global_min_position = approx_min_position
# # # #                             end
# # # #                             if pyplot
# # # #                                 surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
# # # #                             end
# # # #                         end
# # # #                     else
# # # #                         if pyplot
# # # #                             surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
# # # #                         end
# # # #                         println()
# # # #                     end
# # # #                 end 
# # # #             end
# # # #         end
# # # #     )
# # # #     
# # # #     
# # # #     
# # # #     if pyplot && !(recursive)
# # # #         plot_CV_matrix(CV_matrix, ok_matrix, rprm1, rprm2, T, pO2)
# # # #     end
# # # #     return global_min, global_min_position
# # # # end



# function scan_2D(;pyplot=false, 
#                 T=800, pO2=100,
# #                rprm1=range(15, stop=20, length=3),
# #                rprm2=range(15, stop=20, length=3),
#                 #rprm1=range(18.4375, stop=18.59375, length=3),
#                 #rprm2=range(20.59375, stop=21, length=3),
#                 rprm1=range(18.4375, stop=18.59375, length=3),
#                 rprm2=range(20.59375, stop=21, length=3),
#                 wp=[1, 1, 0, 0, 0, 0], #TODO !!!
#                 depth_remaining=2,
#                 recursive=false,
#                 nx=100, ny=100
#                 )
# 
#     
#     #PyPlot.close()
#     #PyPlot.close()
#     
#     ######################################################
#     #rprm1 = collect(-5.0 : 1.0 : 5.0)
#     #rprm1 = collect(-3. : 0.5 : 3.0)
#     #rprm1=range(15, stop=20, length=3)
#     
#     #rprm2 = collect(-3. : 0.5 : 3.0)
#     #rprm2=range(15, stop=20, length=3)
#     ######################################################
#     
#     wpn = ["A0","R0","DGA","DGR","beta","A"]
#     A0 = 21.71975544711280
#     R0 = 19.53
#     DGA = 0.0905748
#     DGR = -0.708014
#     beta = 0.6074566741435283
#     A = -0.356
# 
#     println("computing 2D scan... ")
#     lrprm1 = size(rprm1,1)
#     lrprm2 = size(rprm2,1)
#     
#     CV_matrix = Array{DataFrame}(undef,lrprm1, lrprm2)
#     ok_matrix = Array{Int32}(undef,lrprm1, lrprm2)
#     println("Computing CVs")
#     @time(
#         for i1 in 1:lrprm1
#             for i2 in 1:lrprm2
#                 try
#                     print(rprm1[i1], " / ", rprm2[i2])
#                     
#                     ##########################################
#                     A0 = rprm1[i1]
#                     R0 = rprm2[i2]
#                     ##########################################
#                     
#                     CV_matrix[i1,i2] = ysz_experiments.run_new(
#                         out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
#                         prms_in=[A0, R0, DGA, DGR, beta, A]
#                         )
#                     ok_matrix[i1, i2] = 1
#                     println("  .. ok :)")
#                 catch e
#                     if e isa InterruptException
#                         rethrow(e)
#                     else
#                         ok_matrix[i1, i2] = 0
#                         println("  :(   ")
#                         continue
#                     end
#                 end
#             end
#         end
#     )
#     
#     #PyPlot.figure(figsize=(6,4))
#     #xlabel("parameter")
#     #zlabel("error")
#     #title("Fitness Function")
# 
#     if pyplot && !(recursive)
#         figure(figsize=(3,3))
#     end
#     print_ok_matrix(ok_matrix)
#     
#     println("Computing interpolation ... ")
#     @time(
#         for i1 in 1:lrprm1-1
#             for i2 in 1:lrprm2-1
#                 if all_four_nodes_ok(ok_matrix,i1,i2)
#                     print(" intrp i1 / i2 : ",i1, " / ", i2)
#                     
#                     FF=get_FF_interp_2D(
#                         [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
#                         rprm1[i1:i1+1],
#                         rprm2[i2:i2+1];
#                         nx=nx, ny=ny,
#                         T=T, pO2=pO2
#                     )
#                     
#                     min_count, min_values, min_positions, min_curv = find_mins(
#                         FF.err_matrix, 
#                         [(rprm1[i1+1] - rprm1[i1])/nx,
#                          (rprm1[i2+1] - rprm1[i2])/ny]
#                     )
#                     
# 
#                     
#                     if pyplot
#                         println("jsem tu")
#                         surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
#                     end
#                     
#                     if min_count > 0
#                         print(" ... min ", min_count, "\n    val ", min_values, "\n    pos ", min_positions, "\n    curv ", min_curv)
#                         
#                         if depth_remaining > 0
#                             println("\n >>>>>>>> Going deeper in intrp i1 / i2 : ",i1, " / ", i2)
#                             scan_2D(;pyplot=true, 
#                                 T=T, pO2=pO2,
#                                 rprm1=range(rprm1[i1], stop=rprm1[i1+1], length=3),
#                                 rprm2=range(rprm2[i2], stop=rprm2[i2+1], length=3),
#                                 depth_remaining=depth_remaining - 1,
#                                 recursive=true
#                                 )
#                         end
#                     end
#                     
#                     print("\n")
#                 end 
#             end
#         end
#     )
#     
#     
#     
#     if pyplot && !(recursive)
#         plot_CV_matrix(CV_matrix, ok_matrix, rprm1, rprm2, T, pO2)
#     end
#     return
# end



function del_pyplot()
    for i in 1:10
        PyPlot.close()
    end
end










#####################
## 1D optimization ##
#####################
function get_FF_interp_1D(CV_list, prms_list; n=20)
    CV_raw = import_CVtoDataFrame(;T=800,pO2=100)
    checknodes = get_checknodes(0.05,0.95,-0.95,-0.05,0.01)
    CV_exp = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_raw, checknodes))

    CV1 = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_list[1],checknodes))
    CV2 = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_list[end],checknodes))
    node_list = range(0.0, stop=1.0, length=n)
    err_list_interp = [fitnessFunction(
            CV_exp,
            DataFrame(
                U = checknodes[:,1], 
                I = (linComb_CVs(CV1, CV2, i)).I
            )
        )
        for i in node_list
    ]
    prms_list_interp = node_list.*(prms_list[end]-prms_list[1]).+prms_list[1]
    FF_out = DataFrame(prms = prms_list_interp, err = err_list_interp)
end

function comparison()
    CV_list = []
    prms_list = []
    PyPlot.close()
    
    #for A in collect(-0.36: 0.1 : -0.32)
    #for A in [-0.345, -0.33]
    for A in [-0.356, -0.354]
    #for A in [-0.1, -0.6]
    #for A in collect(-0.1 : -0.01 : -0.6)
    #for A in [-0.3893]
    #for R0 in [21.4, 21.8]
    #for R0 in collect(21.4: 0.02 : 21.8)
    for R0 in [20.53]
        println(A)
        ww=DataFrame()
        push!(CV_list,ysz_experiments.run_new(
            out_df_bool=true, voltammetry=true, sample=10,
            prms_in=[21.71975544711280, R0, 0.0905748, -0.708014, 0.6074566741435283, A]
            )
        )
        push!(prms_list,A)
    end
    end
    println("jdeme dal")
    plot_all(CV_list, prms_list)
    plot_error(CV_list, prms_list)
    return
end


function scan_1D()
    PyPlot.close()
    PyPlot.close()
    #rA = [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.2]    
    #rA = collect(-0.4 : 0.02 : -0.3)
    
    #rprms = [19, 20, 21, 22]
    rprms = collect(-0.5 : 0.1 : 0.0)
    
    
    A0 = 21.71975544711280
    R0 = 20.606423236896422
    DGA = 0.0905748
    DGR = -0.708014
    beta = 0.6074566741435283
    A = -0.356
    
    println("computing ... ")
    CV_list = []
    ok_prms_list = []
    for i in 1:size(rprms,1)
        try
        print(rprms[i])
        
        ########
        A = rprms[i]
        ########
        
        push!(CV_list,ysz_experiments.run_new(
            out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
            prms_in=[A0, R0, DGA, DGR, beta, A]
            )
        )
        push!(ok_prms_list,rprms[i])
        println("  .. ok")
        catch e
            if e isa InterruptException
                rethrow(e)
            else
                println()
                continue
            end
        end
    end

    PyPlot.figure(figsize=(6,4))
    xlabel("parameter")
    ylabel("error")
    title("Fitness Function")
    FF = DataFrame()
    for i in 1:(size(ok_prms_list,1)-1)
        FF = get_FF_interp_1D([CV_list[i], CV_list[i+1]],[ok_prms_list[i], ok_prms_list[i+1]])
        plot(FF.prms, FF.err, linewidth=3.0, label=string("FF ",i))
        plot(FF.prms[1],FF.err[1], "kx", markersize=12.0, ) 
        #plot(FF.prms, FF.err)
    end
    plot(FF.prms[end],FF.err[end], "kx", markersize=12.0)
    legend(loc="best")
    grid()
    
    plot_all(CV_list, ok_prms_list)
    
    return
end



####################################
####################################
####################################
######## Levenberg-Maquard ##########
####################################

function LM_optimize(;EIS_opt_bool=false, CV_opt_bool=false, pyplot=false)
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

   function to_optimize(x) 
        #err = run_new(print_bool=false, fitting=true, voltametry=true, pyplot=false, voltrate=0.005, sample=8, bound=0.41, 
        #    prms_in=x)
        prms = prepare_prms(mask, x0, x)
        print(" >> mask = ",mask)
        print(" || prms = ",prms)
        
       
        (DD, nu, nus, ms_par) = get_shared_add_prms()
        
        err = 0.0
        err_string = ""
        if CV_opt_bool
            CV_sim = ysz_experiments.run_new(
                            out_df_bool=true, voltammetry=true, sample=8, #pyplot=true,
                            prms_in=prms,
                            add_prms_in=(DD, nu, nus, ms_par)
                        )
            if pyplot
                figure(1)
                CV_plot(CV_sim)
            end
            
            err += fitnessFunction(
                CV_apply_checknodes(CV_sim, CV_get_shared_checknodes()), 
                CV_exp
            )
            err_string = string(err_string," CV ")
        end
        if EIS_opt_bool
            EIS_sim = ysz_experiments.run_new(
                            pyplot=false, EIS_IS=true, out_df_bool=true, EIS_bias=EIS_bias, omega_range=EIS_get_shared_omega_range(),
                            dx_exp=-8,
                            # TODO !!! T a pO2
                            prms_in=prms,
                            add_prms_in=(DD, nu, nus, ms_par)
                        )

            
            if pyplot
                figure(2)
                Nyquist_plot(EIS_sim)
            end
            
            err += EIS_fitnessFunction(EIS_sim, EIS_exp)
            err_string = string(err_string," EIS ")
        end
        println(" || ",err_string,": err =", err)
        return [err]
    end
    

    
    T=800
    pO2 = 100
    EIS_bias=0.0
    
    if pyplot
        PyPlot.close()
    end
    
    if CV_opt_bool
        CV_exp = CV_apply_checknodes(
            import_CVtoDataFrame(T=T, pO2=pO2),
            CV_get_shared_checknodes()
        )
        if pyplot
            figure(1)
            CV_plot(CV_exp, "exp")
        end
    end
    
    if EIS_opt_bool
        EIS_exp = EIS_apply_checknodes(
            import_EIStoDataFrame(T=T,pO2=pO2,bias=EIS_bias),
            EIS_get_shared_checknodes()
        )

        if pyplot
            figure(2)
            Nyquist_plot(EIS_exp, "exp")
        end
    end
            
    #x0 = zeros(2)
    #optimize(rosenbrock, x0, LevenbergMarquardt())
    
    lower_bounds = [-20, -20, -1.2, -1.2, 0.1, -1.0]
    upper_bounds = [22, 22, 0.5, 0.5, 0.9, 1.0]
    
    #x0 = (18.8, 19.2,     -0.14,      -0.5,   0.6074566741435283,   0.956)
    x0 = (20.8, 19.2,     -0.14,      0.0,   0.6074566741435283,   0.156)
    x0 = [20.33593076473848, 19.03903812609013, -0.5950818911594278, 0.02485398906461962, 0.6789640127849936, 0.15628515354230042]
    x0 = [20.84794940132141, 19.479956891781775, -0.8232180847274697, -0.17334101552015208, 0.2432389397701255, -0.0625559216024488]
    x0 = [19.620786013494122,19.34421543847122,-0.6539832280099412,-0.37515172426829546,0.2326806290580494,-0.21181634115352901]
    x0 = [19.620786013494122,19.34421543847122,-0.6539832280099412,-0.37515172426829546,0.23268062905804937,-0.21181634115352904]
    
    mask = [1, 1, 1, 1, 1, 1] # determining, which parametr should be fitted
    
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
    
   
    
    #to_optimize(x0M)
    println(optimize(to_optimize, x0M, lower=lowM, upper=uppM, Δ=1000, f_tol=1.0e-14, g_tol=1.0e-14, LevenbergMarquardt()))
    
    ####optimize(to_optimize, x0M, Dogleg())
    return
end


###############################
###############################
#### SNEHURKA COMPUTATION ####
###############################



function par_study(;A0_list=17, 
                    R0_list=16, 
                    DGA_list=0, 
                    DGR_list=0, 
                    beta_list=0.4, 
                    A_list=0.2, 
                    save_dir="../snehurka/data/par_study_default/", scripted_tuple=(1, 0, 1, 0, 0, 0), prms_names=("A0", "R0", "DGA", "DGR", "beta", "A"), 
                    CV_bool="false", EIS_bool="false", mode="test")
    
  function consistency_check()
      
      
      # TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      # (ale az v pripade, ze to bude obecne napsane)
      
      return true
  end

  CV_err_counter::Int32 = 0
  CV_good_counter::Int32 = 0
  EIS_err_counter::Int32 = 0
  EIS_good_counter::Int32 = 0
  all_counter::Int32 = 0
  
  CV_bool == "true" ? (CV_bool=true) : (CV_bool=false)
  EIS_bool == "true" ? (EIS_bool=true) : (EIS_bool=false)

  ##### TODO !!! ###
  EIS_bias = 0.0
  ################

  


  if mode!="go"
        println(" >> number of sets of parameters: ",
        length(A0_list)*length(R0_list)*length(DGA_list)*length(DGR_list)*length(beta_list)*length(A_list))
        @show A0_list
        @show R0_list
        @show DGA_list 
        @show DGR_list 
        @show beta_list
        @show A_list
        @show CV_bool, EIS_bool
        return
  else
    println(" >> number of sets of parameters: ",
        length(A0_list)*length(R0_list)*length(DGA_list)*length(DGR_list)*length(beta_list)*length(A_list))
    @show A0_list
    @show R0_list
    @show DGA_list
    @show DGR_list 
    @show beta_list
    @show A_list
    @show CV_bool, EIS_bool
        
    for A0_i in A0_list
        for R0_i in R0_list
            for DGA_i in DGA_list
                for DGR_i in DGR_list
                    for beta_i in beta_list
                        for A_i in A_list
                        
                            all_counter = all_counter + 1
                            println(string(" <><><><><><><><><><><><> all_counter <><><> calculating ",all_counter," of ", 
                                length(A0_list)*length(R0_list)*length(DGA_list)*length(DGR_list)*length(beta_list)*length(A_list)))
                            print("parameters: (A0, R0, DGA, DGR, beta, A) = (",A0_i,", ",R0_i,", ",DGA_i,", ",DGR_i,", ",beta_i,", ",A_i,")")
                        
                            prms = (A0_i, R0_i, DGA_i, DGR_i, beta_i, A_i)
                            (DD, nu, nus, ms_par) = get_shared_add_prms()
                            
                            if CV_bool
                                try
                                    CV_sim = ysz_experiments.run_new(
                                                    out_df_bool=true, voltammetry=true, sample=10,
                                                    prms_in=prms,
                                                    add_prms_in=(DD, nu, nus, ms_par)
                                                )
                                    CV_save_file_prms(CV_sim, save_dir, prms, prms_names, scripted_tuple)           
                                    CV_good_counter += 1
                                    print("   CV ok :)         ")
                                 catch e
                                    if e isa InterruptException
                                        rethrow(e)
                                    else
                                        println(e)
                                        CV_err_counter += 1
                                        print("<<<<< CV FAIL! >>>>>")
                                    end
                                end
                            end
                            if EIS_bool
                                try
                                    EIS_sim = ysz_experiments.run_new(
                                                pyplot=false, EIS_IS=true, out_df_bool=true, EIS_bias=EIS_bias, omega_range=EIS_get_shared_omega_range(),
                                                dx_exp=-8,
                                                # TODO !!! T a pO2
                                                prms_in=prms,
                                                add_prms_in=(DD, nu, nus, ms_par)
                                            )
                                    EIS_save_file_prms(EIS_sim, save_dir, prms, prms_names, scripted_tuple)           
                                    EIS_good_counter += 1
                                    print("   EIS ok :)        ")
                                 catch e
                                    if e isa InterruptException
                                        rethrow(e)
                                    else
                                        println(e)
                                        EIS_err_counter += 1
                                        print("<<<<< EIS FAIL! >>>>")
                                    end
                                end
                            end
                        
                            println()
                        end
                    end
                end
            end
        end
    end
    println(string("<<<<< CV_err_count/ CV_good_count >>> ", CV_err_counter,"  /  ", CV_good_counter), 
                    string(" |||||| EIS_err_count/ EIS_good_count >>> ", EIS_err_counter,"  /  ", EIS_good_counter, " >>>>>>>"))
  end
end

function par_study_script_wrap(
                    prms_lists_string=("(17, [16, 17, 18], 0.0, 0.0, 0.4, [0.0, 0.2])"),
                    save_dir="./kadinec/", 
                    scripted_tuple_string=("(1, 0, 0, 0, 0, 0)"),
                    prms_names_string="(\"A0\", \"R0\", \"DGA\", \"DGR\", \"beta\", \"A\")", 
                    CV_bool="false", 
                    EIS_bool="false", 
                    mode="test")
  
  prms_lists = eval(Meta.parse(prms_lists_string))
  scripted_tuple = eval(Meta.parse(scripted_tuple_string))
  prms_names = eval(Meta.parse(prms_names_string))
  
  # TODO !!! more general for arbitraty number of parameters
  par_study( 
                    A0_list=prms_lists[1], 
                    R0_list=prms_lists[2], 
                    DGA_list=prms_lists[3], 
                    DGR_list=prms_lists[4], 
                    beta_list=prms_lists[5], 
                    A_list=prms_lists[6], 
                    save_dir=save_dir, 
                    scripted_tuple=scripted_tuple, 
                    prms_names=prms_names, 
                    CV_bool=CV_bool, 
                    EIS_bool=EIS_bool, 
                    mode=mode)
end

function meta_run_par_study()  
  function recursive_bash_call(output_prms_lists, active_idx)
    if active_idx > size(scripted_tuple,1)
      scripted_tuple_string = string(scripted_tuple)
      output_prms_lists_string = string(output_prms_lists)
      prms_names_string = string(prms_names)
      if bash_command == ""
        return run(`
                    $run_file_name 
                    $output_prms_lists_string 
                    $save_dir 
                    $scripted_tuple_string 
                    $prms_names_string 
                    $CV_bool 
                    $EIS_bool 
                    $mode
        `)      
      else
        return run(`
                    $bash_command  
                    $run_file_name 
                    $output_prms_lists_string 
                    $save_dir 
                    $scripted_tuple_string 
                    $prms_names_string 
                    $CV_bool 
                    $EIS_bool 
                    $mode
        `)
      end
    end
    if scripted_tuple[active_idx] == 1
      for i in prms_lists[active_idx]
        recursive_bash_call(push!(deepcopy(output_prms_lists),i), active_idx + 1)
      end
    else
      recursive_bash_call(push!(deepcopy(output_prms_lists),prms_lists[active_idx]), active_idx + 1)
    end
  end
  
  function consistency_check()
    if (size(prms_names,1) != size(scripted_tuple,1) || 
        size(prms_names,1) != size(prms_lists,1))
      return false
    end
    
    for i in 1:size(prms_lists,1)
      if size(prms_lists[i],1) < 1
        return false
      end
    end
    
    return true
  end
  
  # prms definition
  prms_names = ("A0", "R0", "DGA", "DGR", "beta", "A")
  scripted_tuple = (1, 0, 1, 0, 0, 0)
  prms_lists = (
    collect(13 : 5.0 : 14),  # A0
    collect(13 : 4.0 : 13),  # R0
    collect(-0.4 : 0.5 : 0.3), # DGA
    collect(-0.7 : 1.5 : 0.0), # DGR
    collect(0.4),  # beta
    collect(-0.2 : 1.0 : 0.8) # A
  )
  
  # preparing bash output
  bash_command = "echo"
  
  #run_file_name = "kru_run_ysz_fitting_par_study_wrap.chi"
  run_file_name = "../snehurka/run_ysz_fitting_par_study-prms-.jl"
  
  save_dir = "../snehurka/data/prvni_metavrh/"
  CV_bool = "true"
  EIS_bool = "true"
  mode = "go"
  
  if !consistency_check()
    println("ERROR: meta_run_par_study(): shape mismatch (!consistency_check())")
    return throw(Exception)
  end  
  
  # counter of output files ... TODO !!!
  recursive_bash_call([], 1)
    
  # saving metafile   TODO!!!!
  println()
  metafile_string = "METAFILE for par_study\n"
  prms_strings = [string(item) for item in prms_lists]
  for (i,str) in enumerate(prms_strings)
    metafile_string = string(metafile_string,prms_names[i],"_list=",str,"\n")
  end
  metafile_string = string(metafile_string,"prms_names=", prms_names,"\n")
  metafile_string = string(metafile_string,"scripted_tuple=", scripted_tuple,"\n")
  metafile_string = string(metafile_string,"save_dir=", save_dir,"\n")
  metafile_string = string(metafile_string,"CV_bool=", CV_bool,"\n")  
  metafile_string = string(metafile_string,"EIS_bool=", EIS_bool,"\n")
  metafile_string = string(metafile_string,"mode=", mode,"\n")
  metafile_string = string(metafile_string,"bash_command=", bash_command,"\n")
  metafile_string = string(metafile_string,"run_file_name=", run_file_name,"\n")
  metafile_string = string(metafile_string,"prms_lists=", prms_lists,"\n")
  
  println(metafile_string)
  write("$(save_dir)/_metafile_par_study.txt", metafile_string)
  
  println("ok :) ")
  
end
        


end

