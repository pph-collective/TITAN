using CSV
using DataFrames

cur_map = CSV.read("dotmap_params.txt")
new_params = "../dotmap_params.txt"
new_map = "dotmap_updates.txt"

input_dir = "../titan"
output_dir = "titan_new"

file_str = read(new_params, String)
for row in eachrow(cur_map)
    if !ismissing(row.old)
        global file_str = replace(file_str, "$(row.dot)," => "$(row.dot),$(row.old)")
    end
end
write(new_map, file_str)
