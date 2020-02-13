using CSV
using DataFrames

map = CSV.read("dotmap_params.txt")

input_dir = "../titan"
output_dir = "titan_new"

for (root, dirs, files) in walkdir(input_dir)
    for file in files
        if file[end-2:end] == ".py"
            file_str = read(joinpath(root,file), String)
            for row in eachrow(map)
                if !ismissing(row.old)
                    file_str = replace(file_str, "params.$(row.old)" => "params.$(row.dot)")
                end
            end
            write(joinpath(output_dir, file), file_str)
        end
    end
end
