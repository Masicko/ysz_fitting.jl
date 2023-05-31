using CSV
using DataFrames
using PyPlot
using Printf

function export_DataFrametoCV_path(df, f_name)
  open(f_name, "w") do file
    write(file, "End Comments\n")
  end
  if size(df,2) == 3
    CSV.write(f_name, df, delim="\t", append=true)
  elseif size(df,2) == 2
    CSV.write(f_name, DataFrame(U = df.U, I = df.I, t = zeros(size(df,1))), delim="\t", append=true)
  else
    println("ERROR: inappropriate dimensions of recieved DataFrame !!!")
    return throw(Exception)
  end
  return
end

function export_DataFrametoCSV(EIS_df, f_name)
  aux_df = DataFrame(
    f = EIS_df.f, 
    ReZ = real.(EIS_df.Z),
    ImZ = imag.(EIS_df.Z)
    )
  CSV.write(f_name, aux_df)
  return 
end

function export_DataFrametoEIS_path(df, f_name)
  open(f_name, "w") do file
    write(file, "End Comments\n")
  end
  
  
  if size(df,2) == 2
    empty_col = zeros(size(df,1))
    CSV.write(f_name, DataFrame(
        f = reverse(df.f),
        A2 = empty_col,
        A3 = empty_col,
        A4 = empty_col,
        Re = reverse(real.(df.Z)),
        Im = reverse(imag.(df.Z))
        ),
        delim="\t", append=true        
    )
  else
    println("ERROR: inappropriate dimensions of recieved DataFrame !!!")
    return throw(Exception)
  end
  return
end