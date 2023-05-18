module ysz_model_GAS_LoMA_shared

#############################################
# WHILE CREATING NEW MODEL                  #
# DO NOT FORGET TO CHECK                    #
#                                           #
# [x] boundary conditions                   #
# [x] initial conditions                    #
# [x] equilibrium quantities                #
# [x] meas_tran & meas_stdy                 #
# [x] output species, names                 #
# [x] function set_parameters               #
# [x] changed module name :)                #
#                                           #
#############################################

using Printf
using VoronoiFVM
using PyPlot
using DataFrames
using CSV
using LeastSquaresOptim

const bulk_species = (iphi, iy) = (1, 2)
const surface_species = (iyAs, iyOs) = (3, 4)
const surface_names = ("yAs", "yOs")

const index_driving_species = iphi

mutable struct reaction_struct
  r::Float64   # surface adsorption coefficient [ m^-2 s^-1 ]
  DG::Float64 # difference of gibbs free energy of adsorption  [ J ]
  beta::Float64 # symmetry of the adsorption    
  S::Float64 # stechiometry compensatoin of the adsorption
  exp::Float64 # bool deciding if EXP should be used instead of LoMA
end

mutable struct YSZParameters 

    # switches
    separate_vacancy::Bool
    weird_DD::Bool
    
    # reactions
    A::reaction_struct    # oxide adsorption from YSZ
    R::reaction_struct    # electron-transfer reaction
    O::reaction_struct    # oxygen adsorption from gas
    
    # inductance
    L::Float64  
    
    # fixed
    DD::Float64   # diffusion coefficient [m^2/s]
    conductivity::Float64 # choosen_R_ohm
    pO2::Float64 # O2 partial pressure [bar]\z
    T::Float64      # Temperature [K]
    nu::Float64    # ratio of immobile ions, \nu [1]
    nus::Float64    # ratio of immobile ions on surface, \nu_s [1]
    numax::Float64  # max value of \nu 
    nusmax::Float64  # max value of  \nu_s
    x_frac::Float64 # Y2O3 mol mixing, x [%] 
    chi::Float64    # dielectric parameter [1]
    m_par::Float64
    ms_par::Float64
  

    # known
    vL::Float64     # volume of one FCC cell, v_L [m^3]
    areaL::Float64 # area of one FCC cell, a_L [m^2]


    e0::Float64
    eps0::Float64
    kB::Float64  
    N_A::Float64 
    zA::Float64  
    mO::Float64  
    mZr::Float64 
    mY::Float64
    zL::Float64   # average charge number [1]
    yB::Float64   # electroneutral value [1]
    
    phi_eq::Float64 # equilibrium voltage [V]
    y0_eq::Float64
    yAs_eq::Float64
    yOs_eq::Float64
    
    ML::Float64   # averaged molar mass [kg]
    
    # oxygen adsorption sites coverage w.r.t. one surface YSZ cell
    OC::Float64
    
    # overvoltage influence - dimensionless number used as "supplied energy" = -e_fac*(e0)*\eta
    e_fac::Float64
    
    #
    YSZParameters()= YSZParameters( new())
end

function YSZParameters(this)
    
    #swithes
    this.separate_vacancy = true
    this.weird_DD = true
    
    # oxide adsorption from YSZ
    this.A = reaction_struct(1,1,1,1,1)
    this.A.r= 10.0^21
    this.A.DG= 0.0905748 * this.e0 # this.e0 = eV
    this.A.beta = 0.5
    this.A.S= 10^0.0
    this.A.exp= 1
    
    # electron-transfer reaction
    this.R = reaction_struct(1,1,1,1,1)
    this.R.r= 10.0^22
    this.R.DG= -0.708014 * this.e0
    this.R.beta= 0.5
    this.R.S= 10^0.0
    this.R.exp = 1
    
    # oxygen adsorption from gas
    this.O = reaction_struct(1,1,1,1,1)
    this.O.r= 10.0^20
    this.O.DG= 0.0905748 * this.e0 # this.e0 = eV
    this.O.beta = 0.5
    this.O.S= 10^0.0
    this.O.exp= 1
    
    this.L=2.3560245927364395e-6
    
    #this.DD=1.5658146540360312e-11  # [m / s^2]fitted to conductivity 0.063 S/cm ... TODO reference
    #this.DD=8.5658146540360312e-10  # random value  <<<< GOOOD hand-guess
    #this.DD=9.5658146540360312e-10  # some value  <<<< nearly the BEST hand-guess
    this.DD=4.35e-13  # testing value
    this.conductivity=-1.0  # Siemens
    this.pO2=1.0                   # O2 atmosphere 
    this.T=1073                     
    
    this.nu=0.85                     # assumption
    this.nus=0.85                    # assumption
    this.x_frac=0.13                # 13% YSZ
    this.chi=27.0                  # from relative permitivity e_r = 6 = (1 + \chi) ... TODO reference
    this.m_par =  2                 
    this.ms_par = 0.05        
    this.numax = (2+this.x_frac)/this.m_par/(1+this.x_frac)
    this.nusmax = (2+this.x_frac)/this.ms_par/(1+this.x_frac)

    # known
    this.e0   = 1.602176565e-19  #  [C]
    this.eps0 = 8.85418781762e-12 #  [As/(Vm)]
    this.kB   = 1.3806488e-23  #  [J/K]
    this.N_A  = 6.02214129e23  #  [#/mol]
    this.mO  = 16/1000/this.N_A  #[kg/#]
    this.mZr = 91.22/1000/this.N_A #  [kg/#]
    this.mY  = 88.91/1000/this.N_A #  [kg/#]

    this.vL=3.35e-29
    #this.areaL=(this.vL)^0.6666
    this.zA  = -2;
    #this.zL  = 4*(1-this.x_frac)/(1+this.x_frac) + 3*2*this.x_frac/(1+this.x_frac) - 2*this.m_par*this.nu
    #this.yB  = -this.zL/(this.zA*this.m_par*(1-this.nu))
    #this.ML  = (1-this.x_frac)/(1+this.x_frac)*this.mZr + 2*this.x_frac/(1+this.x_frac)*this.mY + this.m_par*this.nu*this.mO
    # this.zL=1.8182
    # this.yB=0.9
    # this.ML=1.77e-25   
    
    # oxygen adsorption sites coverage w.r.t. one surface YSZ cell
    this.OC = this.ms_par*(1.0-this.nus)
    
    #
    this.e_fac = 0.0
    
    return this
end

# conversions for the equilibrium case
# function y0_to_phi(this::YSZParameters, y0)
#     yB = this.yB
#     return - (this.kB * this.T / (this.zA * this.e0)) * log(y0/(1-y0) * (1-yB)/yB )
# end
# 
# function y0_activity_to_phi(this::YSZParameters, a_y0)
#   yB = this.yB
#   return - (this.kB * this.T / (this.zA * this.e0)) * log(a_y0 * (1-yB)/yB )
# end
# 
# function phi_to_y0(this::YSZParameters, phi)
#     yB = this.yB
#     X  = yB/(1-yB)*exp.(this.zA*this.e0/this.kB/this.T* (-phi))
#     return X./(1.0.+X)
# end

function equilibrium_boundary_conditions(this::YSZParameters)
    phi_eq =  (
                (this.kB*this.T*log(sqrt(this.pO2)*((1-this.yB)/this.yB)) - this.O.DG/2.0 - this.R.DG - this.A.DG)
                /
                ((+2*this.e_fac + 2)*this.e0)
              )
    a_yOs = this.pO2^(1/2.0)*exp(-this.O.DG/(2 * this.kB * this.T))
    a_yAs = a_yOs*exp( -(this.R.DG + 2*this.e_fac*this.e0*phi_eq)/(this.kB * this.T))
    a_y0 = a_yAs*exp( -this.A.DG/(this.kB * this.T))
    
#     yOs = a_yOs/(1 + a_yOs)
#     yAs = a_yAs/(1 + a_yAs)
#     y0 = a_y0/(1 + a_y0)
#     
#     @show yOs
#     @show yAs
#     @show y0    
#     
#     @show a_yOs
#     @show a_yAs
#     @show 1/(1.0 + a_yOs + a_yAs)
    
    if this.separate_vacancy
      return  phi_eq,  
              a_y0/(1 + a_y0),    
              a_yAs/(1 + a_yAs),  
              a_yOs/(1 + a_yOs)
    else
      y_V = 1/(1.0 + a_yOs + a_yAs)      
      return  phi_eq,  
              a_y0/(1 + a_y0),   
              a_yAs*(y_V),   
              a_yOs*(y_V)
    end
end

# boundary conditions
function set_typical_boundary_conditions!(sys, parameters::YSZParameters)
    sys.boundary_values[iphi,1]=parameters.phi_eq
    sys.boundary_values[iphi,2]=0.0e-3
    #
    sys.boundary_factors[iphi,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[iphi,2]=VoronoiFVM.Dirichlet
    #
    sys.boundary_values[iy,2]=parameters.yB
    sys.boundary_factors[iy,2]=VoronoiFVM.Dirichlet
end

# initial conditions
function get_typical_initial_conditions(sys, parameters::YSZParameters)   
    inival=unknowns(sys)
    inival.=0.0
    grid = sys.grid
    #
    treshold_for_linear_function = 0.6*10^(-9)
    for inode=1:size(inival,2)
        x = coordinates(grid)[inode]
        if x < treshold_for_linear_function
          inival[iphi, inode] = (((treshold_for_linear_function - x)*parameters.phi_eq)/
                                  treshold_for_linear_function)                                  
          inival[iy, inode] = (((treshold_for_linear_function - x)*parameters.y0_eq + x*parameters.yB)/
                                  treshold_for_linear_function)
        else  
          inival[iphi,inode]=0.0
          inival[iy,inode]= parameters.yB
        end
    end
    inival[iyAs,1] = parameters.yAs_eq
    inival[iyOs,1] = parameters.yOs_eq
    return inival
end

#function temperature_dependence

function YSZParameters_update!(this::YSZParameters)
    this.areaL=(this.vL)^0.6666
    this.numax = (2+this.x_frac)/this.m_par/(1+this.x_frac)
    this.nusmax = (2+this.x_frac)/this.ms_par/(1+this.x_frac)   
    
    this.zL  = 4*(1-this.x_frac)/(1+this.x_frac) + 3*2*this.x_frac/(1+this.x_frac) - 2*this.m_par*this.nu
    this.yB  = -this.zL/(this.zA*this.m_par*(1-this.nu))
    this.ML  = (1-this.x_frac)/(1+this.x_frac)*this.mZr + 2*this.x_frac/(1+this.x_frac)*this.mY + this.m_par*this.nu*this.mO
    
    
    this.phi_eq, this.y0_eq, this.yAs_eq, this.yOs_eq = equilibrium_boundary_conditions(this)

    # the following block computes the right DD according to the right conductivity
    if this.conductivity > 0.0
      testing_DD = 1.06e-11    
      this.DD = testing_DD
      conductivity_test = get_conductivity(this, perform_update=false)
      this.DD = testing_DD * (this.conductivity / conductivity_test)
    end
    
    return this
end


function set_parameters!(this::YSZParameters, prms_values, prms_names)
  found = false
  for (i,name_in) in enumerate(prms_names)
    found = false
    # supposing there is at most one '.'
    if occursin('.', name_in)
      (name_in, attribute_name_in) = split(name_in, '.')
    end
    for name in fieldnames(typeof(this))
      if name==Symbol(name_in)          
        if name_in in ["A", "R", "O"]
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
end





# time derivatives
function storage!(f,u, node, this::YSZParameters)
    f[iphi]=0
    f[iy]=this.mO*this.m_par*(1.0-this.nu)*u[iy]/this.vL
end

function bstorage!(f,u,node, this::YSZParameters)
    if  node.region==1
      if this.separate_vacancy
        f[iyAs]=this.mO*this.ms_par*(1.0-this.nus)*u[iyAs]/this.areaL
        f[iyOs]=this.mO*this.OC*u[iyOs]/this.areaL
      else
        f[iyAs]=this.mO*this.ms_par*(1.0-this.nus)*u[iyAs]/this.areaL
        f[iyOs]=this.mO*this.ms_par*(1.0-this.nus)*u[iyOs]/this.areaL
      end
    end
end

function get_conductivity(parameters; DD=parameters.DD, nu=parameters.nu, perform_update=true)
  if perform_update
    DD_orig = parameters.DD
    nu_orig = parameters.nu
    parameters.DD = DD
    parameters.nu = nu
    YSZParameters_update!(parameters)
  end
  
  f = Array{Float64}(undef, 2)
  uk = Array{Float64}(undef, 2)
  ul = Array{Float64}(undef, 2)
  
  # mental setting
  ### U = 1 V       voltage
  ### l = 1 m       YSZ_length
  ### S = 1 m^2     YSZ_cross_surface_area
  ### => d phi / dx = 1
  ### 
  ### I = abs( 2 * e0 * j^kg / mO)
  ###   = abs( 2 * e0 * f[iy] / mO)
  ###
  ### sigma = (1/R)*(l/S)    [S/m^2]
  ###       = (I/U)*(l/S)      
  ###       = I
  ###       = abs( 2 * e0 * f[iy] / mO) 
  ###        
  ###  R = (1/simga)*(l/S)      
  

  # unit gradient of electrical potential
  uk[iphi] = 0
  ul[iphi] = 1
  
  # no concentration gradient 
  # and electroneutral filling ratio
  uk[iy] = parameters.yB
  ul[iy] = parameters.yB
  
  flux_core!(f, uk, ul, parameters)
  
  if perform_update
    parameters.DD = DD_orig
    parameters.nu = nu_orig
    YSZParameters_update!(parameters)
  end
  return abs( 2 * parameters.e0 * f[iy] / parameters.mO)
end

function flux_core!(f, uk, ul, this::YSZParameters)
  f[iphi]=this.eps0*(1+this.chi)*(uk[iphi]-ul[iphi])
  
  bp,bm=fbernoulli_pm(
      (log(1-ul[iy]) - log(1-uk[iy]))
      -
      this.zA*this.e0/this.T/this.kB
      *(ul[iphi] - uk[iphi])
  )
  f[iy]= (
      this.DD
      *
      (this.weird_DD ? (1.0 + this.mO/this.ML*this.m_par*(1.0-this.nu)*0.5*(uk[iy]+ul[iy]))^2 : 1)
      *
      this.mO*this.m_par*(1.0-this.nu)/this.vL
      *
      (bm*uk[iy]-bp*ul[iy])
  )
end


# bulk flux
function flux!(f,u, edge, this::YSZParameters)
    uk=viewK(edge,u)
    ul=viewL(edge,u)
    
    flux_core!(f, uk, ul, this)
end

# sources
function reaction!(f,u, node, this::YSZParameters)
    f[iphi]=-(this.e0/this.vL)*(this.zA*this.m_par*(1-this.nu)*u[iy] + this.zL) # source term for the Poisson equation, beware of the sign
    f[iy]=0
end

function EXP_reaction_template(this::YSZParameters, RR::reaction_struct; PI_activites, overvoltage=0)
    # PI_activities = a_products/a_reactants
    
    return  (
              (RR.r/RR.S)
              *(
                  exp(-RR.beta*RR.S*
                  (
                    RR.DG + overvoltage
                  )
                  /(this.kB*this.T))
                  *(
                    PI_activites
                  )^(-RR.beta*RR.S)
                  - 
                  exp((1 - RR.beta)*RR.S*
                  (
                    RR.DG + overvoltage
                  )
                  /(this.kB*this.T))
                  *(
                    PI_activites
                  )^((1 - RR.beta)*RR.S)
              )    
            )
end

# surface reactions
function exponential_oxide_adsorption(this::YSZParameters, u; debug_bool=false)
    if this.A.r > 0
        #  <><><><><><><>  this is a correct direction ! <><><><><><><><>
        # O-2(s) + V(y) => O-2(y) + V(s)
        if Bool(this.A.exp)
          the_fac = 1
        else  
          # LoMA
        
          if this.separate_vacancy
            the_fac = (
                (u[iy]*(1-u[iy]))
                *
                (u[iyAs]*(1-u[iyAs]))               
              )^(this.A.S/2.0)
          else
            the_fac = (
                (u[iy]*(1-u[iy]))
                *
                (u[iyAs]*(1-u[iyAs]-u[iyOs]))               
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
                            (u[iyAs]/(1-u[iyAs]))
                          )
                        :
                          (
                            (u[iy]/(1-u[iy]))
                            /
                            (u[iyAs]/(1-u[iyAs]-u[iyOs]))
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

function electroreaction(this::YSZParameters, u; debug_bool=false)
    if this.R.r > 0
        # O(s) + 2e-(s) => O-2(s)
        if Bool(this.R.exp)
          the_fac = 1
        else
          # LoMA
          if this.separate_vacancy
            the_fac = (
                (u[iyOs]*(1-u[iyOs]))
                *
                (u[iyAs]*(1-u[iyAs]))               
              )^(this.R.S/2.0)
          else
            the_fac = (
                (u[iyAs])
                *
                (u[iyOs])
              )^(this.R.S/2.0)
          end
        end
        rate = the_fac*EXP_reaction_template(
          this, 
          this.R, 
          PI_activites= (this.separate_vacancy ?
                          (
                            u[iyAs]/(1-u[iyAs])
                            /
                            (u[iyOs]/(1-u[iyOs]))
                          )
                        :
                          (
                            u[iyAs]
                            /
                            u[iyOs]
                          )
                        ),
          overvoltage=this.e_fac*2*this.e0*u[iphi]
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

function exponential_gas_adsorption(this::YSZParameters, u; debug_bool=false)
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
                (u[iyOs]*(1-u[iyOs]))^2                                            
              )^(this.O.S/2.0)
          else
            the_fac = (
                (u[iyOs]*(1-u[iyOs]-u[iyAs]))^2                                            
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
                            (u[iyOs]/(1-u[iyOs]))^2
                            /
                            (this.pO2)
                          )
                        :
                          (
                            (u[iyOs]/(1-u[iyOs]-u[iyAs]))^2
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


# surface reaction + adsorption
function breaction!(f,u,node,this::YSZParameters)
    if  node.region==1
        electroR=electroreaction(this,u)
        oxide_ads = exponential_oxide_adsorption(this, u)
        gas_ads = exponential_gas_adsorption(this, u)
        f[iy]= - this.mO*oxide_ads
        # if bulk chem. pot. > surf. ch.p. then positive flux from bulk to surf
        # sign is negative bcs of the equation implementation
        f[iyAs]= this.mO*oxide_ads - this.mO*electroR 
        f[iyOs]= this.mO*electroR - this.mO*2*gas_ads
        f[iphi]=0
    else
        f[iy]=0
        f[iphi]=0
    end
end

function direct_capacitance(this::YSZParameters, PHI)
    # Clemens' analytic solution
    #printfields(this)
    
    PHI
    #PHI=collect(-1:0.01:1) # PHI = phi_B-phi_S, so domain as phi_S goes with minus
    my_eps = 0.001
    for i in collect(1:length(PHI))
        if abs(PHI[i]) < my_eps
            PHI[i]=my_eps
        end
    end
    #
    #yB = -this.zL/this.zA/this.m_par/(1-this.nu);
    yB = this.yB
    X= yB/(1-yB)*exp.(.- this.zA*this.e0/this.kB/this.T*PHI)
    y  = X./(1.0.+X)
    #
    nF = this.e0/this.vL*(this.zL .+ this.zA*this.m_par*(1-this.nu)*y)

    
    F  = - sign.(PHI).*sqrt.(
          2*this.e0/this.vL/this.eps0/(1.0+this.chi).*(
            .- this.zL.*PHI .+ this.kB*this.T/this.e0*this.m_par*(1-this.nu)*log.(
              (1-yB).*(X .+ 1.0)
             )
           )
         );
    #
    Y  = yB/(1-yB)*exp.(- this.DGA/this.kB/this.T .- this.zA*this.e0/this.kB/this.T*PHI);
    #
    CS = this.zA^2*this.e0^2/this.kB/this.T*this.ms_par/this.areaL*(1-this.nus)*Y./((1.0.+Y).^2);
    CBL  = nF./F;
    return CBL, CS, y
end



#
# Transient part of measurement functional
#

function set_meas_and_get_tran_I_contributions(meas, U, sys, parameters, AreaEllyt, X)
    Qb= - integrate(sys,reaction!,U) # \int n^F            
    dphi_end = U[iphi, end] - U[iphi, end-1]
    dx_end = X[end] - X[end-1]
    dphiB=parameters.eps0*(1+parameters.chi)*(dphi_end/dx_end)
    Qs= (parameters.e0/parameters.areaL)*parameters.zA*U[iyAs,1]*parameters.ms_par*(1-parameters.nus) # (e0*zA*nA_s)
    meas[1]= AreaEllyt*( -Qs[1] -Qb[iphi]  -dphiB)
    return ( -AreaEllyt*Qs[1], -AreaEllyt*Qb[iphi], -AreaEllyt*dphiB)
end

#
# Steady part of measurement functional
#
function set_meas_and_get_stdy_I_contributions(meas, U, sys, parameters, AreaEllyt, X)
    meas[1] = AreaEllyt*(-2*parameters.e0*electroreaction(parameters, U[:, 1]))
    return AreaEllyt*(-2*parameters.e0*electroreaction(parameters, U[:, 1]))
end


end
