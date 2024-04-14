"""
INPUT: Número de equações ou variáveis (n), a matriz A e o vetor coluna b, A' = [A, B].

OUTPUT: Solução x_1, x_2, ..., x_n =: x do sistema Ax = b
"""
function solve_system(A::Matrix{Float64}, b::Vector{Float64})

    n = length(b)

    dimension_failure(A, n) && return

    A' = [copy[A] copy[b]]

    for i = 1:n-1

        p = i

        for k = i:n
            if A'[k,i] != 0
                p = k
                break
            end

            if k == n
                msg_erro = true
            end
        end

        if msg_erro
            println("O sistema não admite uma única solução")
            return
        end

        a_pi = A'[p,i]

        if p != i
            v = copy(A'[p, :])
            A'[p, :] = A'[i, :]
            A'[i, :] = v
        end

        for j = i + 1:n
            m = A'[j, i] / A'[i, i]

            A'[j, :] = A'[j, :] - m*A'[i, :]
        end
    end

    if A'[n, n] == 0
        println("O sistema não admite uma única solução")
        return        
    end

    x = Vector{Float64}(undef, n)
    x[n] = A'[n, n + 1] / A'[n, n]

end

function dimension_failure(A::Matrix{Float64}, number_of_variables::Int64)

    if number_of_variables != size(A)[1]
        println("O tamanho da matriz A e do vetor b são incompatíveis")
        return true
    end

    if size(A)[1] != size(A)[2]
        println("A matriz A precisa ser quadrada!")
        return true
    end
    
    return false
end