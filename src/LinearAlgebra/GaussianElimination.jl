function solve_system(A::Matrix{Float64}, b::Vector{Float64})

    n_digits = 8

    n = length(b)

    dimension_failure(A, n) && return

    A_ext = [copy(A) copy(b)]

    for i = 1:n-1

        p = i

        msg_erro = false

        for k = i:n
            if round(A_ext[k, i], digits = n_digits) != 0
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

        a_pi = A_ext[p, i]

        if p != i
            v = copy(A_ext[p, :])
            A_ext[p, :] = A_ext[i, :]
            A_ext[i, :] = v
        end

        for j = i+1:n
            m = A_ext[j, i] / A_ext[i, i]

            A_ext[j, :] = A_ext[j, :] - m * A_ext[i, :]
        end
    end

    if round(A_ext[n, n], digits = n_digits) == 0
        println("O sistema não admite uma única solução")
        return
    end

    x = Vector{Float64}(undef, n)
    x[n] = A_ext[n, n+1] / A_ext[n, n]

    for i = n-1:-1:1
        sum = 0
        for j = i+1:n
            sum += A_ext[i, j] * x[j]
        end
        x[i] = (A_ext[i, n+1] - sum) / A_ext[i, i]
    end

    return x
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
