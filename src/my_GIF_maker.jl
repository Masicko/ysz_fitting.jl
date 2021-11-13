module my_GIF_maker

using Plots
using Base: @kwdef
using LaTeXStrings

############# GIF PART ###############

# @userplot NyquistPlot
# @recipe function f(np::NyquistPlot)
#   x, y, i = np.args
# end

@kwdef mutable struct gif_maker_control
  x_ranges::AbstractArray=[]
  y_ranges::AbstractArray=[]
  z_ranges::AbstractArray=[]
  legend_ranges::AbstractArray=[]
  save_dir::String=""
  #
  #
  xlim::AbstractArray=[Inf, -Inf]
  ylim::AbstractArray=[Inf, -Inf]
  #
  save_name::String="gif_default"  
  fps::Float64=1.0
  title::String=""
  plot_legend::Bool=true
end


function make_GIF(gif_maker_control=gif_maker_control() )  
  
#   ENV_backup = "gksqt"
#   try
#     ENV_backup = ENV["GKSwstype"]
#   catch
#     ENV_backup = "gksqt"
#   end
#   ENV["GKSwstype"]="nul"
 
  gmc = gif_maker_control
  
  sequence_length = length(gif_maker_control.x_ranges)
  
  if (gif_maker_control.xlim[1] == Inf) && (gif_maker_control.ylim[1] == Inf)
    auto_xlim = [Inf, -Inf]
    auto_ylim = [Inf, -Inf]
    for seq_idx in collect(1:sequence_length)    
      auto_xlim[1] = min(gif_maker_control.x_ranges[seq_idx]... , auto_xlim[1])
      auto_xlim[2] = max(gif_maker_control.x_ranges[seq_idx]... , auto_xlim[2])
      #
      auto_ylim[1] = min(gif_maker_control.y_ranges[seq_idx]... , auto_xlim[1])
      auto_ylim[2] = max(gif_maker_control.y_ranges[seq_idx]... , auto_xlim[2])      
    end    
  else
    auto_xlim = gif_maker_control.xlim
    auto_ylim = gif_maker_control.ylim
  end
  
  z_range_period = sequence_length
  for i in 2:sequence_length
    if gif_maker_control.z_ranges[i] == gif_maker_control.z_ranges[1]
      z_range_period = i - 1
      break
    end
  end
  
  if z_range_period > -1
    graphs_multiplicity = Int(sequence_length/z_range_period)    
  end
  
#   @show gmc.z_ranges
#   @show graphs_multiplicity
#   @show z_range_period
#   @show sequence_length
  
  p = plot()
  anim = Animation()
  for seq_idx in collect(1:z_range_period)
      act_plot = Plots.plot(
                  title = "$(gmc.title)\nbias = $(gif_maker_control.z_ranges[seq_idx])",
                  aspect_ratio = 1.0, size=(800, 800),
                  xlim = auto_xlim, ylim = auto_ylim,
                  xlabel=L"\textrm{Re} \ Z \ [\Omega]",
                  ylabel=L"\textrm{-Im} \ Z \ [\Omega]",
                  legend=gmc.plot_legend
      )
      for rep_idx in 1:graphs_multiplicity       
        plot!(act_plot, 
                    gif_maker_control.x_ranges[seq_idx + (rep_idx-1)*z_range_period],
                    gif_maker_control.y_ranges[seq_idx + (rep_idx-1)*z_range_period],                    
                    label=gif_maker_control.legend_ranges[rep_idx + (rep_idx-1)*z_range_period]
        )
      end
      frame(anim)
  end
  mkpath(gmc.save_dir)
  gif(anim, "$(gmc.save_dir)$(gmc.save_name)", fps = gif_maker_control.fps)
  
  
#   ENV["GKSwstype"]=ENV_backup
  return
end




@kwdef mutable struct graph_3D_control
  x_ranges::AbstractArray=[]
  y_ranges::AbstractArray=[]
  z_ranges::AbstractArray=[]
  legend_ranges::AbstractArray=[]
  #
  #
  color_theme::String="nothing"
  gray_range::Tuple=(0.0, 0.7)
  #
  save_name::String=""
  save_dir::String=""
  title::String="3D plot čehosi"
  plot_legend::Bool=true
end

function plot_2D_sequence_in_3D(x_ranges, y_ranges; graph_3D_control=graph_3D_control() )
  gtc = graph_3D_control
  sequence_length = length(x_ranges)
  if sequence_length < 1
    return
  end
  if sequence_length != length(y_ranges)
    println("ERROR: length(x_ranges) != length(y_ranges): ", length(x_ranges) , " != ",  length(y_ranges))
    return throw(Exception)
  end
  
  color_controlled = false
  if graph_3D_control.color_theme == "gray"
    color_controlled = true
    color_sequence = collect(range(graph_3D_control.gray_range..., length = sequence_length))    
  end
    
  if length(graph_3D_control.z_ranges) < 1
    graph_3D_control.z_ranges = collect(range(1, sequence_length, step=1))
  end
  
  plot = Plots.plot(
          #camera=(30, -60)                    
          ylabel=L"\textrm{Re} \ Z \ [\Omega]",
          zlabel=L"\textrm{-Im} \ Z \ [\Omega]",
          xlabel=L"\textrm{bias}",
          title=gtc.title,          
          #aspect_ratio = 1.0, 1.0,
          )
  for z_idx in collect(1 : sequence_length)
    x_data = []
    y_data = []
    z_data = []
    x_z_idx_length = length(x_ranges[z_idx])
    for i in 1 : x_z_idx_length 
      if x_z_idx_length != length(y_ranges[z_idx])
        println("ERROR: length(x_ranges[z_idx]) != length(y_ranges[z_idx]): ", length(x_ranges[z_idx]) , " != ",  length(y_ranges[z_idx]), " ... z_idx = ",z_idx)
        return throw(Exception)     
      end
      y_data = x_ranges[z_idx]
      z_data = y_ranges[z_idx]
      x_data = Array{Float64}(undef, x_z_idx_length)
      x_data .= graph_3D_control.z_ranges[z_idx]            
    end
    if color_controlled
      plot!(plot, x_data, y_data, z_data, legend=false, color=RGB(color_sequence[z_idx], color_sequence[z_idx], color_sequence[z_idx]))
    else
      plot!(plot, x_data, y_data, z_data, legend=false)
    end    
  end  
  if gtc.save_name!=""
    mkpath(gtc.save_dir)
    savefig(plot, "$(gtc.save_dir)$(gtc.save_name)")
  else
    @show "JSEM  ............... TU "
    gui(plot)
  end
end

function test()
  x_ranges = []
  push!(x_ranges, collect(0 : 1 : 10))
  push!(x_ranges, collect(0 : 1 : 10))
  push!(x_ranges, collect(0 : 1 : 7))
  
  y_ranges = []
  y_ranges = []
  push!(y_ranges, collect(0 : 1 : 10))
  push!(y_ranges, collect(20 : -2 : 0)) 
  push!(y_ranges, collect(10 : -1 : 3)) 
  
  plot_2D_sequence_in_3D(x_ranges, y_ranges, graph_3D_control=graph_3D_control())
  return
end




end #of module
