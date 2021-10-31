module my_GIF_maker

using Plots
using Base: @kwdef

@kwdef mutable struct graph_3D_control
  x_ranges::AbstractArray=[]
  y_ranges::AbstractArray=[]
  z_ranges::AbstractArray=[]
  legend_ranges::AbstractArray=[]
  color_theme::String="nothing"
  gray_range::Tuple=(0.0, 0.7)
end

function plot_2D_sequence_in_3D(x_ranges, y_ranges; graph_3D_control=graph_3D_control() )
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
          xlabel="bias",
          ylabel="Re\$(\\Omega)\$",
          zlabel="-Im\$(\\Omega)\$",
          title="3D plot čehosi"
          
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
      plot!(plot, x_data, y_data, z_data, color=RGB(color_sequence[z_idx], color_sequence[z_idx], color_sequence[z_idx]))
    else
      plot!(plot, x_data, y_data, z_data)
    end
  end
  gui(plot)
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
