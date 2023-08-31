module ComradeVIDAExt

using Comrade
if isdefined(Base, :get_extension)
    using VIDA
else
    using ..VIDA
end

function Comrade.template_align(div::VIDA.AbstractDivergence, m, lower, upper)

end

end
