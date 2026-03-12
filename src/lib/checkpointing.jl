module Checkpointing

using Serialization
using SHA
using JSON
using DotEnv

export verify_or_clear_checkpoints, save_checkpoint, load_checkpoint, clear_all_checkpoints


"""
Restituisce il path della cartella di checkpoint base.
"""
function get_base_checkpoint_dir()
    # Ensure .env is loaded to guarantee CHECKPOINT_DIR is picked up
    env_path = joinpath(@__DIR__, "../../.env")
    if isfile(env_path)
        DotEnv.load!(env_path)
    end
    
    path = get(ENV, "CHECKPOINT_DIR", joinpath(pwd(), "solver_checkpoints"))
    println("Using base checkpoint directory: $path")
    return path
end

"""
Restituisce il path della cartella dedicata al singolo progetto.
"""
function get_project_dir(simulation_id::String)
    return joinpath(get_base_checkpoint_dir(), simulation_id)
end

"""
Restituisce il filepath esatto per il checkpoint di quel particolare step.
"""
function get_checkpoint_filepath(simulation_id::String, step_name::String)
    project_dir = get_project_dir(simulation_id)
    if !isdir(project_dir)
        mkpath(project_dir)
    end
    return joinpath(project_dir, "$(step_name).jls")
end

"""
Verifica la validità dell'hash rispetto al current_checksum forzato.
Se non valido, cancella la cartella del progetto per resettare l'ambiente.
"""
function verify_or_clear_checkpoints(simulation_id::String, current_checksum::String)
    project_dir = get_project_dir(simulation_id)
    if !isdir(project_dir)
        mkpath(project_dir)
    end
    
    hash_file = joinpath(project_dir, "input_hash.txt")
    if isfile(hash_file)
        saved_checksum = read(hash_file, String)
        if saved_checksum != current_checksum
            println("Nuovi parametri di input rilevati per il progetto $(simulation_id) (checksum mismatch). Invalido e rimuovo i vecchi checkpoint.")
            rm(project_dir, recursive=true, force=true)
            mkpath(project_dir)
        else
            println("Parametri identici al run precedente per il progetto $(simulation_id). Checkpoint validati per l'uso.")
        end
    end
    
    # Scrittura del checksum confermato
    open(hash_file, "w") do f
        write(f, current_checksum)
    end
end

"""
Salva lo stato di una simulazione in un file isolato.
Implementa una scrittura sicura tramite file .tmp
"""
function save_checkpoint(simulation_id::String, step_name::String, state_dict::Dict)
    filepath = get_checkpoint_filepath(simulation_id, step_name)
    tmp_filepath = filepath * ".tmp"
    println("Salvando stato dello step '$(step_name)' per il progetto $(simulation_id) in: $filepath")
    
    try
        serialize(tmp_filepath, state_dict)
        mv(tmp_filepath, filepath, force=true)
    catch e
        println("Errore durante il salvataggio del checkpoint '$filepath': $e")
        if isfile(tmp_filepath)
            rm(tmp_filepath, force=true)
        end
    end
end

"""
Carica un checkpoint se esiste per lo specifico step.
Altrimenti restituisce `nothing`.
"""
function load_checkpoint(simulation_id::String, step_name::String)
    filepath = get_checkpoint_filepath(simulation_id, step_name)
    if isfile(filepath)
        println("Caricando checkpoint esistente dello step '$(step_name)' per il progetto $(simulation_id) da: $filepath")
        try
            state_dict = deserialize(filepath)
            return state_dict
        catch e
            println("Attenzione: Impossibile deserializzare il file di checkpoint '$filepath' (potrebbe essere corrotto). Verrà ricalcolato. Errore: $e")
            return nothing
        end
    end
    return nothing
end


"""
Rimuove la intera directory dei checkpoint se la simulazione è conclusa con successo.
"""
function clear_all_checkpoints(simulation_id::String)
    base_dir = get_base_checkpoint_dir()
    project_dir = joinpath(base_dir, simulation_id)
    if isdir(project_dir)
        println("Simulazione $simulation_id completata. Rimozione directory checkpoint $project_dir")
        rm(project_dir, recursive=true, force=true)
    end
end

end # module Checkpointing
