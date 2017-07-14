
Pkg.add("JuMP");
Pkg.add("Cbc");
Pkg.add("Gurobi");
# Pkg.add("CPLEX"); # O Cplex pode ser instalado com esse comando.

# Pacotes para visualização de resultados:
Pkg.add("Plotly");
Pkg.add("PlotRecipes");

using JuMP
using Cbc

# Criando o modelo com o resolvedor Cbc:
model = Model(solver=CbcSolver());

# 0 <= x <= 15
@variable(model, x, lowerbound=0, upperbound=15)

# y é variável inteira não negativa:
@variable(model, y, lowerbound=0, Int)

# Problema de maximização da função x + 2y:
@objective(model, Max, x + 2y)

# Restrições:
@constraint(model, 2x + 3y <= 25)
@constraint(model, 3x + 2y <= 37)

# Impressão do modelo:
println(model)

status = solve(model)

if status == :Optimal
    println("Solução ótima encontrada!")
    println("x = $(getvalue(x)), y = $(getvalue(y))")
    println("Valor da função objetivo é $(getobjectivevalue(model))")
else
    println("Solução ótima não foi encontrada.")
end

# Posições de cada ponto a ser visitado:
citiesDict = Dict{Int, Any}()
citiesDict[1] = (523, 418)
citiesDict[2] = (527, 566)
citiesDict[3] = (435, 603)
citiesDict[4] = (386, 660)
citiesDict[5] = (346, 692)
citiesDict[6] = (431, 730)
citiesDict[7] = (419, 818)
citiesDict[8] = (347, 520)
citiesDict[9] = (332, 330)
citiesDict[10] = (165, 374)
citiesDict[11] = (196, 198)
citiesDict[12] = (187, 108)
citiesDict[13] = (210, 63)

nCities = length(citiesDict)
# Matriz para amazenar as distâncias entre os pontos:
c = zeros(nCities, nCities)

# Pkg.add("PyPlot")  # Instalação do pacote para plotagem.
using PlotRecipes  # Pacote com a função para plotar grafos

# Backend utilizado para plotagem:
pyplot()  # Melhor visualização, mas possui mais dependências.
# plotly()  # Mais simples.

# Definição da posição de cada vértice para desenho:
posX, posY = [], []

for i in sort(collect(keys(citiesDict)))
    posI = citiesDict[i]
    for j in sort(collect(keys(citiesDict)))
        posJ = citiesDict[j]
        c[i, j] = ((posI[1] - posJ[1])^2 + (posI[2] - posJ[2])^2)^0.5
    end 
    append!(posX, posI[1])
    append!(posY, posI[2])
end

# Configuração do desenho
graphplot(1:nCities, 1:nCities, names=1:nCities,
    x=posX, y=posY, fontsize=10,
    m=:white, l=:black)

using JuMP
using Gurobi
# using CPLEX

# Definição do modelo e parâmetros para a solução: 
model = Model(solver=GurobiSolver(TimeLimit=20, Threads=1, 
        Heuristics=0.0, OutputFlag=0))

# Para o CPLEX pode-se usar:
# model = Model(solver=CplexSolver(CPX_PARAM_TILIM=20, 
#         CPX_PARAM_THREADS=1, CPX_PARAM_HEURFREQ=0,
#         CPX_PARAM_SCRIND=0))

# Definição da variável matricial xij.
# Note que a variável nCities foi definida anteriormente.
@variable(model, x[i=1:nCities, j=1:nCities; i != j], Bin)

# Função objetivo:
@objective(model, Min, sum(c[i, j] * x[i, j]
        for i in 1:nCities, j in 1:nCities if i != j))

# Restrições para garantir a saída de todo vértice:
for i in 1:nCities
    @constraint(model, sum(x[i, j] for j in 1:nCities if i != j) == 1)
end 

# Restrições para garantir a chegada em cada vértice:
for j in 1:nCities
    @constraint(model, sum(x[i, j] for i in 1:nCities if i != j) == 1)
end

status = solve(model)

if status == :Optimal
    edgeOrigin = []
    edgeDest = []

    for i in keys(citiesDict)
        for j in keys(citiesDict)
            if i != j && getvalue(x[i, j]) > 0.99
                append!(edgeOrigin, i)
                append!(edgeDest, j)
            end
        end
    end
else
    println("Solução ótima não encontrada!")
end

# Visualizando a solução: 
display(
        graphplot(
            edgeOrigin, edgeDest, names=1:nCities, 
            x=posX, y=posY, fontsize=10,
            m=:white, l=:black
        )
    )

using Graphs

function lazyConstraintsCallback(cb)

    # Criando um grafo:
    g = simple_graph(nCities, is_directed=false)

    for i in 1:nCities, j in 1:nCities
        if i != j
            if getvalue(x[i, j]) > 0.01
                add_edge!(g, i, j)  # Adicionando as arestas.
            end
        end
    end

    # Encontrando os componentes conexos do grafo:
    cc = connected_components(g)

    if length(cc) > 1
        # Caso só haja uma componente conexa
        # não há subciclo e não se adiciona 
        # nenhuma restrição.
        
        minTour = sort(cc, by=length)[1]
        subtourLhs = AffExpr()
        
        # Encontrando as arestas que fazem parte do subciclo
        for i in minTour
            for j in minTour
                if i != j && getvalue(x[i, j]) > 0.01
                    subtourLhs += x[i, j]
                end 
            end
        end
        
        # Adicionando a restrição:
        @lazyconstraint(cb, subtourLhs <= length(minTour) - 1)
    end
end # End function

# [... Código com a implementação do modelo ...]

addlazycallback(model, lazyConstraintsCallback)

solve(model)

edgeOrigin = []
edgeDest = []

for i in keys(citiesDict)
    for j in keys(citiesDict)
        if i != j && getvalue(x[i, j]) > 0.01
            append!(edgeOrigin, i)
            append!(edgeDest, j)
        end
    end
end

graphplot(edgeOrigin, edgeDest, names=1:nCities,
    x=posX, y=posY, fontsize=10, m=:white, l=:black)

using Graphs

function userCutsCallback(cb)

    # Criando um grafo:
    g = simple_graph(nCities, is_directed=false)

    for i in 1:nCities, j in 1:nCities
        if i != j
            if getvalue(x[i, j]) > 0.01
                add_edge!(g, i, j)  # Adicionando as arestas.
            end
        end
    end

    # Encontrando os componentes conexos do grafo:
    cc = connected_components(g)

    minTour = sort(cc, by=length)[1]
    subtourLhs = AffExpr()

    # Encontrando as arestas que fazem parte do subciclo
    countEdges = 0
    for i in minTour
        for j in minTour
            if i != j && getvalue(x[i, j]) > 0.01
                subtourLhs += x[i, j]
                countEdges += 1
            end 
        end
    end

    # Adicionando a restrição:
    if length(cc) > 1
        @usercut(cb, subtourLhs <= length(minTour) - 1)
    elseif countEdges > length(vertices(g))
        @usercut(cb, subtourLhs == length(vertices(g)))
    end
    
end # End function

addcutcallback(model, userCutsCallback)
addlazycallback(model, lazyConstraintsCallback)

solve(model)

function userHeuristicsCallback(cb)
    for i in 1:nCities, j in 1:nCities
        if i != j
            if  getvalue(x[i, j]) > 0.1
                # Atribuindo um valor à variável:
                setsolutionvalue(cb, x[i, j], 1)
            end            
        end
    end
    
    # Submetendo solução:
    addsolution(cb)                    
    
end

addcutcallback(model, userCutsCallback)
addlazycallback(model, lazyConstraintsCallback)
addheuristiccallback(model, userHeuristicsCallback)

solve(model)
