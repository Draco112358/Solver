function isMaterialConductor(materialName::String, materials)::Bool
    material = materials[findfirst(m -> m["name"] == materialName, materials)]
    return material["conductivity"] > 0.0
end