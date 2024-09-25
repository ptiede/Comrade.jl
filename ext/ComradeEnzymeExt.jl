module ComradeEnzymeExt

using Enzyme

function __init__()
    # We need this to ensure than Enzyme can AD through the Comrade code base
    Enzyme.API.runtimeActivity!(true)
end

end