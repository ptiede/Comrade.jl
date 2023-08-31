module ComradeMakieExt

using Comrade
if isdefined(Base, :get_extension)
    using Makie
else
    using ..Makie
end

end
