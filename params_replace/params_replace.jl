using CSV
using DataFrames

map = CSV.read("dotmap_updates.txt", quotechar=';')
post_map = CSV.read("post_process_map.txt")

input_dir = "../settings"
output_dir = "titan_new"

touch(joinpath(output_dir, "__init__.py"))

for (r, d, f) in walkdir(input_dir)
    for dir in d
        if isdir(joinpath(output_dir, dir))
            rm(joinpath(output_dir, dir), recursive=true)
        end
        mkdir(joinpath(output_dir, dir))
        touch(joinpath(output_dir,dir,"__init__.py"))
        for (root, dirs, files) in walkdir(joinpath(r,dir))
            for file in files
                if file[end-2:end] == ".py"
                    file_str = read(joinpath(root,file), String)

                    # prepend with imports
                    file_str = "from dotmap import DotMap\n\nparams = DotMap()\n" * file_str


                    for row in eachrow(map)
                        if !ismissing(row.old)
                            file_str = replace(file_str, "$(row.old)" => "params.$(row.dot)")
                        end
                    end

                    # and swith IDU to PWID
                    file_str = replace(file_str, "IDU" => "PWID")

                    # change some common vocab switches
                    file_str = replace(file_str, "p_value" => "prob")

                    write(joinpath(output_dir, dir, file), file_str)

                    file_base = file[1:end-3]
                    new_file = file_base * ".yml"

                    run(`poetry run python to_yaml.py $output_dir.$dir.$file_base`)


                end
            end
        end
    end
end
