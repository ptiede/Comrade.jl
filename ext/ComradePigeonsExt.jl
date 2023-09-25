module ComradePigeonsExt

using Comrade

if isdefined(Base, :get_extension)
    using Pigeons
else
    using ..Pigeons
end



end
