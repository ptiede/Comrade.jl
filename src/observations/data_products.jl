abstract type DataProducts end

Base.@kwdef struct ComplexVis{T} <: DataProducts
    component::T = Stokes(:I)
end

struct Coherency end

Base.@kwdef struct Amplitudes{T} <: DataProducts
    component::T = Stokes(:I)
end

Base.@kwdef struct CPhase <: DataProducts
    diagonal::Bool = true
    minimal::Bool  = true
end


Base.@kwdef struct LogCamp <: DataProducts
    diagonal::Bool = true
    minimal::Bool  = true
end
