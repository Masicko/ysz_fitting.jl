module YSZ_LSM_model_ARO

# TODO
# [x] surface species & reactions
#       [x] orientation and rates are OK
#       [x] OCV
#       [x] bulk-interface weights
# [ ] YSZ current impedance - old 
# [ ] compare the freq images of the 3 currents 


using Plots
using VoronoiFVM
using Printf
using ExtendableGrids
using GridVisualize
# using ForwardDiff
using Plots
using PyPlot
using DataFrames
using Base: @kwdef

using DataFramesMeta
# 
const bulk_domains = (Ω_LSM,Ω_YSZ) = (1, 2)
const surface_domains = (Γ_LSM,Γ_YSZ,Γ) = (1, 2, 3)

const bulk_species = (ie, iy, iphi) = (1, 2, 3)
const surface_species = (ies, iys, ios) = (4, 5, 6)
const species = (ie, iy, iphi, ies, iys, ios)
const species_names = ("ie", "iy", "iphi", "ies", "iys", "ios")

# basic physics constants
const kB = 1.3806503e-23
const e0 = 1.602176634e-19
const eps0 = 8.8541878128e-12 
# constant parameters
const m_par = 2.0
const ze=-1.0
const za=-2.0
const zC_LSM=1.0

mutable struct fooNode
    region
end 


mutable struct reaction_struct
    r::Float64   # surface adsorption coefficient [ m^-2 s^-1 ]
    r_A::Float64
    r_B::Float64
    r_C::Float64
    DG::Float64 # difference of gibbs free energy of adsorption  [ J ]
    DG_A::Float64
    DG_B::Float64
    DG_C::Float64
    beta::Float64 # symmetry of the adsorption    
    S::Float64 # stechiometry compensatoin of the adsorption
    exp::Bool# bool deciding if EXP should be used instead of LoMA
    #
    
  end
  

    
@kwdef mutable struct materialParameters <: VoronoiFVM.AbstractData

    # geometry parameters
    # grid
    L_LSM::Float64=1.0e-3
    L_YSZ::Float64=1.0e-3
    h_min::Float64=1.01e-12
    h_max::Float64=1.9e-1
    Γ_node_idx::Int32=-1
    #
    S_ellyt::Float64=1.0 # [m^2]
    
    # Material parameters
    # inductance
    L::Float64 = 0.0 #2.45e-6 # [H]
    
    ## Bulk
    epsY::Float64=(1.0+27.0)
    epsE::Float64=(1.0+27.0)
  
    DD::Float64=1.0e-11
    De::Float64=1.0e-6            
    T::Float64 = 1073.0
    
    ### 
    CO::Float64=1.0
    COmm::Float64=1.0
    
    
    
    #####################
    
    ### YSZ
    x_frac::Float64 = 0.13
    nu::Float64 = 0.285
    
   
    ## Interface  
    #COs::Float64=1.0e18
    # oxide adsorption from YSZ
    #                                          r,r_A,r_B,r_C    DG,DG_A,DG_B, DG_C, beta,S,  exp
    O::reaction_struct = reaction_struct( 1.0e21,  1,  1,  1,  0.0,   1,   1,    1,  0.5,1,  false)
    E::reaction_struct = reaction_struct( 1.0e25,  1,  1,  1,  0.0,   1,   1,    1,  0.5,1,  false)
    R::reaction_struct = reaction_struct( 1.0e21,  1,  1,  1,  0.0,   1,   1,    1,  0.5,1,  false)
    A::reaction_struct = reaction_struct( 1.0e21,  1,  1,  1,  0.0,   1,   1,    1,  0.5,1,  false)

    # boundary conditions
    bias::Float64=0.0
    # additional model parameters
    # - [ ] y, (1-y), y*(1-y)
    # - [ ] rate types
    separate_vacancy=true
    pO2::Float64=1.0
    ie_bulk_eqn_scaling::Float64 = 1.0e-2
      
    ### implied quantities ####
    nC_YSZ::Float64 = 2.98023223876953e28 # 1/m^3
    zC_YSZ::Float64  = 4*(1-x_frac)/(1+x_frac) + 3*2*x_frac/(1+x_frac) - 2*m_par*nu
    y_YSZ::Float64  = -zC_YSZ/(za*m_par*(1-nu))    
    numax = (2+x_frac)/m_par/(1+x_frac)
    
    ### LSM
    e_LSM::Float64 = -zC_LSM/ze 
    nC_LSM::Float64= 2.98023223876953e28 # = nC_YSZ
    
    # constants
    e0::Float64 = 1.602176634e-19
end 

function update_parameters!(this)
      @show this.y_YSZ
      @show this.zC_YSZ
      @show this.nu
    
    this.zC_YSZ = 4*(1-this.x_frac)/(1+this.x_frac) + 3*2*this.x_frac/(1+this.x_frac) - 2*m_par*this.nu
    this.y_YSZ = -this.zC_YSZ/(za*m_par*(1-this.nu))
    this.e_LSM = -zC_LSM/ze
    this.numax = (2+this.x_frac)/m_par/(1+this.x_frac)
    
    grid = makegrid(this)
    X = grid.components[XCoordinates]
    @show length(X)
    this.Γ_node_idx = Int((length(X) + 1)/2)
end

function conv_prms(prms_names, prms_values)
  prms_names_out = deepcopy(prms_names)
  prms_values_out = deepcopy(prms_values)
  
  for (idx_name, name) in enumerate(prms_names)
    if name == "DD"
      prms_names_out[idx_name] = "DD"
    end
  end
  
  push!(prms_names_out, "De"); push!(prms_values_out, "1.0e-11"); 
  push!(prms_names_out, "DD"); push!(prms_values_out, "1.0e-6");
  
  return prms_names_out, prms_values_out
end

# function set_parameters!(parameters, d::Dict)
#     for key in keys(d)
#         #key in fieldnames(typeof(parameters)) ? setfield!(parameters, Symbol(key), d[key]) : println("Field $key does not belong to Parameters.")
#         key in fieldnames(typeof(parameters)) ? setfield!(parameters, Symbol(key), d[key]) : nothing
#       @show key, fieldnames(typeof(parameters))
#       printfields(parameters)
#       @show  key in fieldnames(typeof(parameters))
#     end
#     update_parameters!(parameters)
# end

# function set_parameters!(parameters, d::Dict)
#     for key in keys(d)
#         #key in fieldnames(typeof(parameters)) ? setfield!(parameters, Symbol(key), d[key]) : println("Field $key does not belong to Parameters.")
#         
#         if occursin('.', key)
#           (key_struct, key_attribute) = split(key, '.')
#         end
#         
#           
#         key in fieldnames(typeof(parameters)) ? setfield!(parameters, Symbol(key), d[key]) : nothing            
#     end
#     update_parameters!(parameters)
# end


function set_parameters!(this, d::Dict)
  prms_names = String.(collect(keys(d)))
  prms_values = collect(values(d))
  
  found = false
  for (i,name_in) in enumerate(prms_names)
    found = false
    # supposing there is at most one '.'
    if occursin('.', name_in)
        (name_in, attribute_name_in) = split(name_in, '.')
    end
    
    # this is here for backward compatibility .... COmm = ms_par*(1 - nus)
    if name_in =="ms_par"
      setfield!(this, Symbol("COmm"), convert(Float64, prms_values[i]*(1 - 0.85)))
      continue
    elseif name_in =="OC"
      setfield!(this, Symbol("CO"), convert(Float64, prms_values[i]))
      continue
    end
    
    for name in fieldnames(typeof(this))
      if name==Symbol(name_in)          
        if name_in in ["A", "R", "O", "E"]
          actual_struct = getfield(this, name)
          for attribute_name in fieldnames(typeof(actual_struct))
            if attribute_name==Symbol(attribute_name_in)              
              if attribute_name_in in ["r", "S"]
                setfield!(actual_struct, attribute_name, Float64(10.0^prms_values[i]))
              elseif attribute_name_in in ["DG"]
                setfield!(actual_struct, attribute_name, Float64(prms_values[i]*this.e0))   #  [DGA] = eV
              else
                setfield!(actual_struct, attribute_name, Float64(prms_values[i]))
              end
            end
          end
        else
          setfield!(this, name, convert(typeof(getfield(this, name)), prms_values[i]))
        end
        found = true
        break
      end
    end
    if !(found)
      println("ERROR: set_parameters: parameter name \"$(name_in)\" not found!")
      throw(Exception)
    end
  end
  update_parameters!(this)
end


function equilibrium_voltage(sys)
    # adjust
    data = sys.physics.data
    E = (1/2)*(
          kB*data.T/e0*log(
            (
              (data.pO2)^(1/2)
              *
              (data.e_LSM)^2 
              *
              ((1-data.y_YSZ)/(data.y_YSZ))
            ) # pokus s prevracenou hodnotou neni dobry           
          )
          -
          (1/e0)*(
            data.A.DG + data.R.DG + 2*data.E.DG + data.O.DG/2
          )
        )
end

function set_bcs!(sys)
    data=sys.physics.data
    boundary_dirichlet!(sys,iphi,Γ_YSZ, 0.0)
    boundary_dirichlet!(sys,iphi,Γ_LSM, data.bias + equilibrium_voltage(sys))
    # densities
    @show data.y_YSZ
    boundary_dirichlet!(sys,iy,Γ_YSZ,data.y_YSZ)
    boundary_dirichlet!(sys,ie,Γ_LSM,data.e_LSM)    
end

function update_problem!(sys,d::Dict)
    set_parameters!(sys.physics.data, d)
    set_bcs!(sys)
    # YSZ_LSM_model_ARO.update_parameters!(sys.physics)
    # @show phibc = equilibrium_voltage(sys) + sys.physics.data.bias
    # boundary_dirichlet!(sys,iphi,Γ_LSM, phibc )
end

function simplegrid(;n=50)
	h=1.0/n
	X=collect(0:h:1); 
	grid=simplexgrid(X)
	ExtendableGrids.bfacemask!(grid,[0.5],[0.5],Γ)
	ExtendableGrids.cellmask!(grid,[0.0],[0.5],Ω_LSM)
	ExtendableGrids.cellmask!(grid,[0.5],[1.0],Ω_YSZ)
    
end

function makegrid(;hmin=1e-9, hmax=1e-1, L_YSZ=1.0, L_LSM=1.0, tol=1e-10)
    tol=min(tol, hmin/3.0, 1e-10)
	YSZ = geomspace(0.0,L_YSZ,hmin,L_YSZ*hmax, tol=tol)
	LSM = geomspace(-L_LSM,0.0,L_LSM*hmax,hmin, tol=tol)
	X=glue(LSM,YSZ, tol=tol)
	grid = simplexgrid(X)
    #tol = min(hmin, 1e-10)
	ExtendableGrids.cellmask!(grid,[0.0],[L_YSZ],Ω_YSZ,  tol=tol)
    ExtendableGrids.cellmask!(grid,[-L_LSM],[0.0],Ω_LSM, tol=tol)
    ExtendableGrids.bfacemask!(grid,[0.0],[0.0],Γ,     tol=tol)
    if num_bfaces(grid) > 3
        throw(ErrorException)
    end
    #for x in grid.components[BFaceNodes]
    #    @show grid.components[Coordinates][x]
    #end
	grid
end

makegrid(data::materialParameters) = makegrid(hmin=data.h_min, hmax=data.h_max, L_YSZ=data.L_YSZ, L_LSM=data.L_LSM)

function flux!(f,_u,edge, data)
    # sedanflux
	if edge isa fooNode
		u = _u
	else
    	u = unknowns(edge,_u)
	end
    # 
    if edge.region == Ω_LSM
	    f[iphi]= eps0*data.epsE*(u[iphi,1] - u[iphi,2])
	    gie = -e0*ze/data.T/kB*(u[iphi,2]-u[iphi,1])
        f[ie] = data.ie_bulk_eqn_scaling*data.De*sedan_part(gie, u[ie,1], u[ie,2])
    elseif edge.region == Ω_YSZ
	    f[iphi]= eps0*data.epsY*(u[iphi,1] - u[iphi,2])
	    mu1=log(1-u[iy,1])
        mu2=log(1-u[iy,2])
        giy = +(mu2-mu1)-e0*za/data.T/kB*(u[iphi,2]-u[iphi,1])
        # f[iy] = data.Di*m_par*(1-data.nu)*sedan_part(giy, u[iy,1], u[iy,2])
        f[iy] = data.DD*sedan_part(giy, u[iy,1], u[iy,2])
    end
end

function sedan_part(g, iy1, iy2)
    bp,bm=fbernoulli_pm(g)
    return -(
                  iy2*bp - iy1*bm 
             )
end

function reaction!(f,u,node, data)   
    f[ie]=0.0
    f[iy]=0.0
    if node.region == Ω_LSM
        f[iphi]= -e0*data.nC_LSM*(zC_LSM + ze*u[ie])
    elseif node.region == Ω_YSZ
        f[iphi]= -e0*data.nC_YSZ*(data.zC_YSZ + (1-data.nu)*m_par*za*u[iy])
    end
end

function storage!(f,u,node,data)
    f[iphi]= 0.0
    f[ie]= 0.0
    f[iy]= 0.0
    if node.region == Ω_LSM
        f[ie]= data.ie_bulk_eqn_scaling*u[ie]
    elseif node.region == Ω_YSZ
        f[iy]= u[iy]
    end
end

function EXP_reaction_template(this, RR::reaction_struct; PI_activites)
    # PI_activities = a_products/a_reactants
    # e^(beta*affinity) - e^((beta-1)*affinity)
    # affinity = - DG/kB/T + log(PI_activites)
    return  (            
              #(RR.r/(RR.S*RR.S))
              (RR.r/RR.S)
              #(RR.r)
              *(
                  exp(-RR.beta*RR.S*
                    (
                      RR.DG
                    )
                    /(kB*this.T)
                  )                  
                  *
                  (
                    PI_activites
                  )^(-RR.beta*RR.S)
                  -
                  exp((1 - RR.beta)*RR.S*
                  (
                    RR.DG 
                  )
                  /(kB*this.T))
                  *(
                    PI_activites
                  )^((1 - RR.beta)*RR.S)
              )    
            )
end

# surface reactions
function oxide_desorption(this, u; debug_bool=false)
    if this.A.r > 0
        #  <><><><><><><>  this is the correct direction ! <><><><><><><><>
        # O-2(s) + V(y) => O-2(y) + V(s)
        if Bool(this.A.exp)
          the_fac = 1
        else  
          # LoMA
        
          if this.separate_vacancy
            the_fac = (
                (u[iy]*(1-u[iy]))
                *
                (u[iys]*(1-u[iys]))               
              )^(this.A.S/2.0)
          else
            the_fac = (
                (u[iy]*(1-u[iy]))
                *
                (u[iyAs]*(1-u[iys]-u[iyos]))               
              )^(this.A.S/2.0)
          end
        end
        rate = the_fac*EXP_reaction_template(
          this, 
          this.A, 
          PI_activites= (this.separate_vacancy ?
                          (
                            (u[iy]/(1-u[iy]))
                            /
                            (u[iys]/(1-u[iys]))
                          )
                        :
                          (
                            (u[iy]/(1-u[iy]))
                            /
                            (u[iys]/(1-u[iys]-u[iyos]))
                          )
                        )
        )
    else
      the_fac = 0
      rate=0
    end
    if debug_bool                  
      print("  A > ")
      @show the_fac, rate, this.A
    end
    return rate
end

function electroreaction(this, u; debug_bool=false)
    if this.R.r > 0
        # O(s) + 2e-(s) => O-2(s)
        if Bool(this.R.exp)
          the_fac = 1
        else
          # LoMA
          if this.separate_vacancy
            the_fac = (
                (u[ios]*(1-u[ios]))
                *
                (u[iys]*(1-u[iys]))              
                *
                ((u[ies])^2.0)
              )^(this.R.S/2.0)
          else
            the_fac = (
                (u[iys])
                *
                (u[ios])
                *
                (u[ies])^2.0
              )^(this.R.S/2.0)
          end
        end
        rate = the_fac*EXP_reaction_template(
          this, 
          this.R, 
          PI_activites= (this.separate_vacancy ?
                          (
                            u[iys]/(1-u[iys])
                            /
                            (u[ios]/(1-u[ios]))
                            /
                            (u[ies]^2)
                          )
                        :
                          (
                            u[iys]
                            /
                            u[ios]
                            /
                            u[ies]
                          )
                        )
                    )
    else
      the_fac = 0
      rate = 0
    end
    if debug_bool
      print("  R > ")
      @show the_fac, rate
    end
    return rate
end

function oxygen_adsorption(this, u; debug_bool=false)
    if this.O.r > 0 && !(this.pO2 == 0)
        # O2(g) => 2O(s)
        if Bool(this.O.exp)
          the_fac = 1
        else  
          # LoMA
          if this.separate_vacancy
            the_fac = (
                (this.pO2)
                *
                (u[ios]*(1-u[ios]))^2                                            
              )^(this.O.S/2.0)
          else
            the_fac = (
                (u[ios]*(1-u[ios]-u[iys]))^2                                            
                *
                (this.pO2)
              )^(this.O.S/2.0)
          end
        end
        rate = the_fac*EXP_reaction_template(
          this, 
          this.O, 
          PI_activites= (this.separate_vacancy ?
                          (
                            (u[ios]/(1-u[ios]))^2
                            /
                            (this.pO2)
                          )
                        :
                          (
                            (u[ios]/(1-u[ios]-u[iys]))^2
                            /
                            (this.pO2)
                          )
                        )
        )
    else
      the_fac = 0
      rate=0
    end
    if debug_bool
      print("  O > ")
      @show the_fac, rate
    end
    return rate
end

function electron_adsorption(this, u; debug_bool=false)
    # e-(LSM) => e-(s) 
    if this.E.r > 0
        if Bool(this.E.exp)
          the_fac = 1
        else
          # LoMA
          the_fac = (
               (u[ie])
               *
               (u[ies])
            )^(this.E.S/2.0)
        end
        rate = the_fac*EXP_reaction_template(
          this, 
          this.E, 
          PI_activites= (
                            u[ies]/u[ie]
                        )
        )
    else
      the_fac = 0
      rate = 0
    end
    if debug_bool
      print("  E > ")
      @show the_fac, rate
    end
    return rate# LoMA (with beta = 0.5)
end


function bstorage!(f,u, node, data)
    # adjust
    if node.region==3
        f[iys] = data.COmm*u[iys]
        f[ios] = data.CO*u[ios]
        f[ies] = u[ies]        
    end
end


function breaction!(f,u,node, data)
    for ispec in species
        f[ispec] = 0.0
    end
    if node.region==Γ   
        r_E = electron_adsorption(data, u)
        r_A = oxide_desorption(data, u)
        r_R = electroreaction(data,u)
        r_O = oxygen_adsorption(data, u)
        # bulk
        f[ie]  = data.ie_bulk_eqn_scaling*data.nC_LSM^( -1 )*(  r_E  )
        f[iy]  = (m_par*(1-data.nu)*data.nC_YSZ)^( -1 )*( -r_A    )
        # surface
        f[ies] = data.nC_LSM^(-2/3)*(- r_E + 2*r_R)
        f[iys] = (m_par*(1-data.nu)*data.nC_YSZ)^(-2/3)*(  r_A -  r_R )
        f[ios] = (m_par*(1-data.nu)*data.nC_YSZ)^(-2/3)*(- 2*r_O + r_R )
        ## surface Poisson
        f[iphi] =  -e0*(
            (m_par*(1-data.nu)*data.nC_YSZ)^(2/3)*za*u[iys] 
            +
            data.nC_YSZ^(2/3)*data.zC_YSZ
            + 
            data.nC_LSM^(2/3)*(
                ze*u[ies] + zC_LSM
            )
        )
        #f[iphi] = 0.0
    end
end



































function get_oxide_electric_flux(sys, solution)
  dummy = fooNode(0)
  X = sys.grid.components[XCoordinates]
  
  dummy.region = Ω_YSZ
  
  uk = solution[:, end - 1]
  ul = solution[:, end]
  ukl = cat(uk, ul, dims=2)
  h = X[end] - X[end - 1]
  flux_f = zeros(sys.physics.num_species)
  
  
  flux!(flux_f, ukl, dummy, sys.physics.data)
  return sys.physics.data.S_ellyt*(flux_f[iy]/h)*e0*(-2)*sys.physics.data.nC_YSZ*m_par*(1 - sys.physics.data.nu)
end

function get_electron_electric_flux(sys, solution)
  dummy = fooNode(0)
  X = sys.grid.components[XCoordinates]
  
  dummy.region = Ω_LSM
  uk = solution[:, 1]
  ul = solution[:, 2]
  ukl = cat(uk, ul, dims=2)
  h = X[2] - X[1]
  flux_f = zeros(sys.physics.num_species)
  
  
  flux!(flux_f, ukl, dummy, sys.physics.data)
  return sys.physics.data.S_ellyt*sys.physics.data.ie_bulk_eqn_scaling^(-1)*(flux_f[ie]/h)*e0*(-1)*sys.physics.data.nC_LSM
end
  
function printfields(this)
    for name in fieldnames(typeof(this))
        @printf("%8s = ",name)
        println(getfield(this,name))
    end
end

function test_system(sys)
    data = sys.physics.data
        @show equilibrium_voltage(sys)
        @show species_names
        flux_f    =zeros(sys.physics.num_species)
        reaction_f=zeros(sys.physics.num_species)
        storage_f =zeros(sys.physics.num_species)
        breaction_f=zeros(sys.physics.num_species)
        bstorage_f =zeros(sys.physics.num_species)
        uk = [0.4, 0.5, 0.6, 0.5, 0.6, 0.1]
        ul = [0.5, 0.6, 0.4, 0.4, 0.5, 0.2]
        ukl = cat(uk,ul,dims=2)
        dummy = fooNode(0)
        # bulk_domains = (Ω_LSM,Ω_YSZ) = (1, 2)
        # surface_domains = (Γ_LSM,Γ_YSZ,Γ) = (1, 2, 3)
        @show 
        println("---------------------------------------")
        println(" Bulk regions")
        println("---------------------------------------")
        for region in bulk_domains #[Ω_LSM,Ω_YSZ]
            println("-----------------------")
            @show region
            dummy.region=region
            for test in [uk, ul]
                reaction!(reaction_f, test, dummy, sys.physics.data)
                @show reaction_f 
                storage!(storage_f, test, dummy, sys.physics.data)
                @show storage_f
            end
            flux!(flux_f, ukl, dummy, sys.physics.data)
            @show flux_f
        end
        println("---------------------------------------")
        println(" Surface regions")
        println("---------------------------------------")
        for bregion in surface_domains
            @show bregion
            dummy.region=bregion
            for test in [uk, ul]
                # R_rates =  R_rate(test, 1, sys.physics.data), R_rate(test, 2, sys.physics.data), R_rate(test, 3, sys.physics.data)
                # @show R_rates
                breaction!(breaction_f, test, dummy, sys.physics.data)
                @show breaction_f 
                bstorage!(bstorage_f, test, dummy, sys.physics.data)
                @show bstorage_f
                println("---")
            end
            println("-----------------------")
        end
end

function sys(;params_dict=:None)
    data=materialParameters()            
    if params_dict != :None
        set_parameters!(data, params_dict)
    end  
    
    grid=makegrid(hmin=data.h_min, hmax=data.h_max, L_YSZ=data.L_YSZ, L_LSM=data.L_LSM)
    physics=VoronoiFVM.Physics(data=data,num_species=6, flux=flux!, storage=storage!, reaction=reaction!, bstorage=bstorage!, breaction=breaction!)
    sys = VoronoiFVM.System(grid,physics) # ?? Sparse
    # enable species
    enable_species!(sys,iphi,[Ω_LSM,Ω_YSZ])
    enable_species!(sys,iy,[Ω_YSZ])
    enable_species!(sys,ie,[Ω_LSM])
    enable_boundary_species!(sys,iys,[Γ])
    enable_boundary_species!(sys,ies,[Γ])
    enable_boundary_species!(sys,ios,[Γ])
    # set boundary conditions consistent with parameters
    set_bcs!(sys)
    return sys
end

function equilibrium_solution(sys; testing=false)
    inival=unknowns(sys)
    inival[iphi,:] .= 0.0
    inival[iy,:] .= sys.physics.data.y_YSZ
    inival[ie,:] .= sys.physics.data.e_LSM
    inival[iys,:] .= sys.physics.data.y_YSZ
    inival[ies,:] .= sys.physics.data.e_LSM
    inival[ios,:] .= 0.5
    #
    control=VoronoiFVM.NewtonControl()
    control.tol_absolute = 1e-13
     control.max_iterations = 1000
     control.damp_initial = 1e-8
     control.damp_growth = 1.1
    testing ? control.verbose=true : control.verbose=false

    update_problem!(sys, Dict(:bias => 0.0))
    solve!(inival,inival,sys, control=control)
    if testing 
        @show LSM_current(sys,inival)
        @show YSZ_current(sys,inival)
        plotsolution(sys,inival, zoom=5.0e-9)
    end
    return inival
end

function stationary_update!(inival=nothing,sys=nothing)
    control=VoronoiFVM.NewtonControl()
    control.tol_absolute = 1e-11
    # control.max_iterations = 1000
    # control.damp_initial = 1e-8
    # control.damp_growth = 1.1
    solve!(inival,inival,sys, control=control)
end

function phi_stationary_sweep(sys, equilibrium_solution; bias_range=collect(0.01:0.01:1.0))
    eq_voltage = equilibrium_voltage(sys)
    df = DataFrame(bias=Float64[], solution=typeof(equilibrium_solution)[])
    push!(df, 
          Dict( :bias    => 0.0, 
                :solution => equilibrium_solution
            )
        )
    for dir in [1,-1]
        result = deepcopy(equilibrium_solution)
        for pbias in bias_range
            update_problem!(sys, Dict(:bias => dir*pbias))
            stationary_update!(result, sys)
            push!(df, 
                  Dict( :bias    => dir*pbias, 
                        :solution => deepcopy(result)
                    ) 
                )
        end
    end
    # 
    sort!(df, :bias)
    return df
end


function test_IV(sys=Nothing; bound=1.0, step=0.01)
    if sys==Nothing
      sys = YSZ_LSM_model_ARO.sys()
    end
    eq = equilibrium_solution(sys)
    df = phi_stationary_sweep(sys, eq) 
    p = Plots.plot()
    for cF in [LSM_current, YSZ_current, YSZ_current_neg]        
        current(x) = cF(sys, x)[1]
        curcol = Symbol(cF)
        df[!, curcol] .= current.(df.solution)
        Plots.plot!(p, df.bias, df[!,curcol], seriestype=:scatter, label=string(Symbol(cF)))
    end
    gui(p)
    return df
end

function impedance_sweep(sys,steadystate;f_range=geometric(0.9, 1.0e+5, 1.1), print_bool=false, currentF=LSM_current,excited_bc = Γ_LSM,
    excited_spec = iphi)

    # excited_bcval=sys.physics.data.bias + equilibrium_voltage(sys)
    excited_bcval=sys.boundary_values[excited_spec,excited_bc]

    function I_stdy(meas, u)
        U=reshape(u,sys)
        meas[1] = currentF(sys, U)[1]
    end
    function I_tran(meas, u)
        U=reshape(u,sys)
        meas[1] = currentF(sys, U)[2]
    end

    # Create impedance system
    isys=VoronoiFVM.ImpedanceSystem(sys,steadystate,excited_spec, excited_bc)
    dstdy=measurement_derivative(sys,I_stdy,steadystate)
    dtran=measurement_derivative(sys,I_tran,steadystate)

    df = DataFrame(f=Float64[], Z=ComplexF64[])
    #w = 2*pi*f_range[1]
    w_range = 2*pi*f_range
    for w in w_range
        print_bool && @show w
        zfreq=freqdomain_impedance(isys,w,steadystate,excited_spec,excited_bc,excited_bcval, dstdy, dtran)
        inductance = im*sys.physics.data.L*w
        push!(df, Dict(:f => w/2.0/pi, :Z => inductance + 1.0/zfreq))
        print_bool && @show zfreq
    end    
    print_bool && @show df
    return df
end

function impedance_sweep_test(sys=Nothing)
    if sys==Nothing
      sys = YSZ_LSM_model_ARO.sys()
    end
    eqsol = YSZ_LSM_model_ARO.equilibrium_solution(sys)#, testing=true)
    
    p = Plots.plot(ratio=:equal)
    for (cF, exbc) in zip([LSM_current, YSZ_current],[Γ_LSM,Γ_YSZ,Γ_LSM])
        df = impedance_sweep(sys, eqsol, currentF=cF,excited_bc = exbc)
        Plots.plot!(p, real.(df.Z), -imag.(df.Z), seriestype=:scatter, label=string(Symbol(cF)))
    end
    gui(p)
end





























### Electric currents
 ################# LSM ############
function LSM_current(sys, U)
    factory = VoronoiFVM.TestFunctionFactory(sys)
    tLSM = testfunction(factory, Γ, Γ_LSM)
    prefactor = sys.physics.data.ie_bulk_eqn_scaling^(-1)*e0*ze*sys.physics.data.nC_LSM
    currents = VoronoiFVM.integrate_stdy(sys, tLSM, U)
    return sys.physics.data.S_ellyt.*(prefactor*currents[ie], currents[iphi])
end

function LSM_transient(sys, Unew, Uold, tstep)
    factory = VoronoiFVM.TestFunctionFactory(sys)
    tLSM = testfunction(factory, Γ, Γ_LSM)
    prefactor = sys.physics.data.ie_bulk_eqn_scaling^(-1)*e0*ze*sys.physics.data.nC_LSM
    trans = VoronoiFVM.integrate(sys, tLSM, Unew, Uold, tstep)
    stdyNew = VoronoiFVM.integrate_stdy(sys, tLSM, Unew)
    stdyOld = VoronoiFVM.integrate_stdy(sys, tLSM, Uold)
    return sys.physics.data.S_ellyt.*(prefactor*trans[ie] + (stdyNew[iphi] - stdyOld[iphi])/tstep)
end

 ############ YSZ #############
function YSZ_current(sys, U)
    data = sys.physics.data
    factory = VoronoiFVM.TestFunctionFactory(sys)
    tYSZ = testfunction(factory, Γ, Γ_YSZ)
    prefactor = e0*za*(1.0 - data.nu)*m_par*data.nC_YSZ
    stdy = VoronoiFVM.integrate_stdy(sys, tYSZ, U)
    #return (-1).*(prefactor*stdy[iy], stdy[iphi])
    return sys.physics.data.S_ellyt.*(prefactor*stdy[iy], stdy[iphi])
end

function YSZ_current_neg(sys, U)
    current = YSZ_current(sys, U)
    return (-1.0).*current
end

function YSZ_trans_neg(sys, Unew, Uold, tstep)
    factory = VoronoiFVM.TestFunctionFactory(sys)
    data = sys.physics.data
    factory = VoronoiFVM.TestFunctionFactory(sys)
    tYSZ = testfunction(factory, Γ, Γ_YSZ)
    prefactor = e0*za*(1.0 - data.nu)*m_par*data.nC_YSZ
    trans = VoronoiFVM.integrate(sys, tYSZ, Unew, Uold, tstep)
    stdyNew = VoronoiFVM.integrate_stdy(sys, tYSZ, Unew)
    stdyOld = VoronoiFVM.integrate_stdy(sys, tYSZ, Uold)
    return (-1 * sys.physics.data.S_ellyt) .* (prefactor*trans[iy] + (stdyNew[iphi] - stdyOld[iphi])/tstep)
end
 
 
 ########## LEGACY ############
# return (steady, trans)
function legacy_current(sys, U)
  return ( legacy_current_steady(sys, U), legacy_current_transient(sys, U))
end
 
 # TODO ... prejmenovat funkce na spravne nazvy
function legacy_current_steady(sys, U)
  parameters = sys.physics.data
  #u_S = tpb_view(sys, U)
  u_S = U[:,parameters.Γ_node_idx]
  return parameters.S_ellyt*( -2*e0*electroreaction(parameters, u_S) )
end

# .... the following function is not needed, because Juergen implemented "integrate" in a intelligent way separating different domains
# function YSZ_reaction!(f,u,node, data)   
#     f[ie]=0.0
#     f[iy]=0.0
#     if node.region == Ω_LSM
#         f[iphi]= 0
#     elseif node.region == Ω_YSZ
#         f[iphi]= -e0*data.nC_YSZ*(data.zC_YSZ + (1-data.nu)*m_par*za*u[iy])
#     end
# end

function legacy_current_transient(sys, U)
  Qb = - integrate(sys, reaction!, U)  
  dphi_end = U[iphi, end] - U[iphi, end-1]
  X = sys.grid.components[XCoordinates]
  dx_end = X[end] - X[end-1]
  parameters = sys.physics.data
  data = parameters
  dphiB=parameters.epsY*(dphi_end/dx_end)
  
  #@show dphi_end
  #@show dx_end
  

    
  u_S = U[:,data.Γ_node_idx]
  
  #breaction_output = zeros(length(species_names))
  #breaction!(breaction_output, u_S, fooNode(Γ), parameters)
  #Qs = - breaction_output[iphi]
  
  Qs = e0*(
            (m_par*(1-data.nu)*data.nC_YSZ)^(2/3)*za*u_S[iys] 
            +
            data.nC_YSZ^(2/3)*data.zC_YSZ
            + 
            data.nC_LSM^(2/3)*(
                ze*u_S[ies] + zC_LSM
            )
        )
    
   @show Qb
   @show Qb[iphi, 2]
   @show Qs
#   @show dphiB
#   @show parameters.S_ellyt
   
  return parameters.S_ellyt*( 
                      - Qs 
                      - Qb[iphi, 2]  # YSZ naboj
                      #- Qb[iphi, 1  ] 
                      - dphiB
                    )
end
 
 
function tpb_view(sys, solution)
    subgrid_tpb = subgrid(sys.grid,[Γ], boundary=true)
    tpb_species = zeros(5)
    for i = 1:5
        tpb_species[i] = view(solution[i,:],subgrid_tpb)[1]
    end
    return tpb_species
end

function find2pos(X, tol)
    Y = abs.(X .- tol)
    val1,argminpos = findmin(Y)
    Y = abs.(X .+ tol)
    val2,argminneg = findmin(Y)
    return argminneg, argminpos
end

function phi_view(sys, solution; pos=5e-9)
    X = sys.grid.components[XCoordinates]
    ineg, ipos = find2pos(X, pos)
    subgrid_tpb = subgrid(sys.grid,[Γ], boundary=true)
    # solution[iphi,ineg]
    return solution[iphi,ineg], view(solution[iphi,:],subgrid_tpb)[1], solution[iphi,ipos]
end

function plotsolution(sys,solution;zoom=1.0)        
    @show solution[iphi,1]
    grid = deepcopy(sys.grid)
    zoom == 1.0 ? zoo=1 : zoo=3
    visualizer = GridVisualize.GridVisualizer(layout=(1,zoo),resolution=(600,300),Plotter=PyPlot,fignum=1);
    
    ysz = Ω_YSZ
    lsm = Ω_LSM

    
    subgrid_YSZ = subgrid(grid,[ysz])
    subgrid_LSM = subgrid(grid,[lsm])
    subgrid_zoom= subgrid(grid,[lsm,ysz])
    
    
	eLSM=view(solution[ie,:],subgrid_LSM)
	yYSZ=view(solution[iy,:],subgrid_YSZ)
    phi_zoom=view(solution[iphi,:],subgrid_zoom)
	scalarplot!(visualizer[1,1],subgrid_zoom, phi_zoom,clear=true,color=:red, label="phi")
	scalarplot!(visualizer[1,1],subgrid_LSM, eLSM,clear=false,show=true,color=:green, label="e")
	scalarplot!(visualizer[1,1],subgrid_YSZ, yYSZ,clear=false,show=true,color=:blue, label="y")
	
	
	
    if zoo != 1
        @show zoom
        lsm=10
        ysz=20
        ExtendableGrids.cellmask!(grid,[-1.0*zoom],[0.0],lsm, tol=1e-21)
        ExtendableGrids.cellmask!(grid,[0.0],[1.0*zoom],ysz, tol=1e-21)        
        subgrid_YSZ = subgrid(grid,[ysz])
        subgrid_LSM = subgrid(grid,[lsm])
        subgrid_zoom= subgrid(grid,[lsm,ysz])
        eLSM=view(solution[ie,:],subgrid_LSM)
	    yYSZ=view(solution[iy,:],subgrid_YSZ)
        phi_zoom=view(solution[iphi,:],subgrid_zoom)
        neg, surf, pos = phi_view(sys, solution)
        scalarplot!(visualizer[1,2],subgrid_zoom, phi_zoom,clear=true,color=:red, label="phi")
	    scalarplot!(visualizer[1,3],subgrid_LSM, eLSM,clear=true,color=:green, label="e")
	    scalarplot!(visualizer[1,3],subgrid_YSZ, yYSZ,clear=false,show=true,color=:blue, label="y")
    end
	reveal(visualizer)
    return visualizer
end

function geometric(start::Float64, stop::Float64, quotient::Float64)
    geom = []
    while start < stop
        push!(geom, start)
        start *=quotient
    end
    return geom
end

end # of the module










