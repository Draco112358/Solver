# Questo script deve contenere le chiamate principali al tuo codice
# per "esercitare" i percorsi di codice critici.
# Ad esempio, se `solver_start.jl` definisce una funzione `start_solver()`, chiamala qui.

# Assicurati che il percorso sia corretto rispetto alla root del progetto

include("main.jl")
# Se main.jl avvia un server, potresti volerlo avviare e poi spegnerlo
# per assicurarti che le funzioni di avvio siano precompilate.
# Esempio:
# if @isdefined(start_solver)
#     solver_instance = start_solver()
#     # ... fai qualche operazione con il solver ...
#     # Se c'è una funzione di stop:
#     # stop_solver(solver_instance)
# end

# Se main.jl è un semplice script che esegue tutto, basta l'include.