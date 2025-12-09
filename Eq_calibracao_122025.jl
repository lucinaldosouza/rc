using LinearAlgebra, Statistics, StatsBase
using LaTeXStrings, Plots, QuantEcon
using BenchmarkTools
using SparseArrays

# Teste de pré alocação da matriz Q
function build_Q(a_size, z_size, s_i_vals, z_chain)
    n = a_size * z_size
    Q = zeros(n, a_size, n)

    for next_s_i in 1:size(Q, 3)
        for a_i in 1:size(Q, 2)
            for s_i in 1:size(Q, 1)
                z_i = s_i_vals[s_i, 2]
                next_z_i = s_i_vals[next_s_i, 2]
                next_a_i = s_i_vals[next_s_i, 1]
                if next_a_i == a_i
                    Q[s_i, a_i, next_s_i] = z_chain.p[z_i, next_z_i]
                end
            end
        end
    end

    
    return Q
end

# Construção das matrizes de transição e utilidade
function Household(; r = 0.01,
                   w = 1.0,
                   μ = 1.0,
                   sigma = 2.0,
                   beta = 0.96,
                   ζ = 3.66,    # Parâmetro de phase-out ζ
                   m = 0.05,    # Valor base para transferência
                   y_med = 1.0,
                   τ_k = 0, 
                   τ_c = 0,
                   Tr = 0.0,    # Lump-sum transfer; p/ kd choque uma transferencia B(z)
                   taxay = 0.0,
                   trans = 0.0,
                   cmin = 0,
                   Nivel_IR = 0.825,
                   Q_pre = Q_FIX,
                   #z_chain = MarkovChain([0.9 0.1; 0.1 0.9], [0.1; 1.0]),
                   #= z_chain = MarkovChain([0.9 0.025 0.025 0.025 0.025;
                                        0.025 0.9 0.025 0.025 0.025;
                                        0.025 0.025 0.9 0.025 0.025;
                                        0.025 0.025 0.025 0.9 0.025;
                                        0.025 0.025 0.025 0.025 0.9],
                    [1, 2, 3, 4, 5]), =#
                    z_chain = MarkovChain([0.417058294170583     0.1928807119288071    0.11768823117688232   0.08209179082091791   0.055894410558944105  0.0365963403659634    0.0376962303769623    0.0265973402659734    0.0207979202079792    0.012698730126987301;
                    0.1881                0.2533                0.1382                0.0985                0.0821                0.0642                0.0618                0.0514                0.0381                0.0243;
                    0.08009199080091992   0.1151884811518848    0.25807419258074193   0.13678632136786323   0.12038796120387961   0.08789121087891212   0.06919308069193081   0.06699330066993302   0.03859614038596141   0.026797320267973206;
                    0.0524                0.0887                0.1456                0.2794                0.1366                0.0929                0.0746                0.0630                0.0431                0.0237;
                    0.035896410358964105  0.07699230076992301   0.09199080091990801   0.1272872712728727    0.21787821217878214   0.1512848715128487    0.10718928107189282   0.08929107089291072   0.06669333066693331   0.0354964503549645;
                    0.0305                0.0546                0.0813                0.0851                0.1432                0.2406                0.1487                0.1075                0.0737                0.0348;
                    0.026700000000000005  0.05500000000000001   0.05650000000000001   0.06990000000000002   0.09230000000000001   0.15360000000000001   0.23490000000000003   0.15060000000000004   0.10700000000000001   0.053500000000000006;
                    0.022097790220977905  0.04609539046095391   0.0524947505249475    0.05649435056494351   0.07939206079392061   0.09309069093090691   0.1498850114988501    0.25407459254074594   0.15638436156384364   0.08999100089991001;
                    0.0176017601760176    0.0328032803280328    0.040604060406040594  0.0477047704770477    0.059005900590058995  0.0742074207420742    0.0912091209120912    0.1549154915491549    0.3058305830583058    0.1761176117611761;
                    0.007301460292058411  0.020404080816163232  0.02090418083616723   0.027505501100220042  0.03250650130026005   0.031906381276255245  0.05251050210042008   0.08681736347269454   0.1807361472294459    0.5394078815763153],
                    #[0.17; 0.42; 0.70; 0.76; 0.83; 1.0; 1.22; 1.55; 2.20; 5.53]),
                    [0.111607; 0.275734; 0.459557; 0.498948; 0.544904; 0.65651; 0.800943; 1.01759; 1.44432; 3.6305]),
                   a_min = 1e-10,
                   a_max = 18.0,
                   a_size = 200,
                   a_vals = range(a_min, a_max, length = a_size),
                   # -Inf is the utility of dying (0 consumption)
                   u = sigma == 1 ? x -> log(x) :
                       x -> (x^(1 - sigma) - 1) / (1 - sigma))
    #

    # Create grids
    Q = Q_pre
    z_size = length(z_chain.state_values)
    # n = a_size * z_size
    n       = size(Q, 1)
    s_vals = gridmake(a_vals, z_chain.state_values)
    s_i_vals = gridmake(1:a_size, 1:z_size)

    inv_P = (z_chain.p)^1000
    inv_P = inv_P[1,:]

    # Compute the average productivity
    # y_med = sum(z_chain.state_values .* inv_P)
    y_med = 1

    # Função taxa/transferência
    T(y_e, y_i, λ, τ = 0.087) = max(y_e - λ * y_e^(1 - τ), 0) - ω(y_i; ζ=ζ, y_med=y_med)
    # T(y_e, y_i, λ, θ = 0.087, y_med = y_med) = max(exp(log(λ)*(y_e/y_med)^(-2θ)),0) - ω(y_i; ζ=ζ, y_med=y_med)

    
    # Fill in the Q and R
    # Q = spzeros(Float64, n, a_size, n) não recebe matrizes com dim > 2
    #Q = zeros(n, a_size, n)
    #= for next_s_i in 1:size(Q, 3)
        for a_i in 1:size(Q, 2)
            for s_i in 1:size(Q, 1)
                z_i = s_i_vals[s_i, 2]
                next_z_i = s_i_vals[next_s_i, 2]
                next_a_i = s_i_vals[next_s_i, 1]
                if next_a_i == a_i
                    Q[s_i, a_i, next_s_i] = z_chain.p[z_i, next_z_i]
                end
            end
        end
    end =#


    R = fill(-Inf, n, a_size)

    for new_a_i in 1:size(R, 2)
        a_new = a_vals[new_a_i]
        for s_i in 1:size(R, 1)
            a = s_vals[s_i, 1]
            z = s_vals[s_i, 2]

            # Calcula a renda do trabalho e do capital
            y_e = w * z  # Renda do trabalho (bruto)
            y_k = r * a  # Renda do capital (bruto)
            y_i = y_e + y_k  # Renda total do agente

            # Cálculo da transferência
            Tr = T(y_e, y_i, Nivel_IR)
            

            # c = w * z + (1 + r) * a - a_new
            c = (w * z + (1 + (1 - τ_k)*r)*a - Tr - a_new) / (1 + τ_c)
            # c = w * z + (1 + (1 - τ_k)*r)*a - a_new
            # tw = w * z + (1 + r) * a
            tw = w * z - Tr + (1 + (1 - τ_k) * r) * a

            if c < cmin
                tw_rel = ((cmin - c) / cmin)^(4)
            else
                tw_rel = 1
            end

            minimum = min(cmin, tw * tw_rel)

            #=
            if c > 0
                R[s_i, new_a_i] = u(c)
            end
            =#
             if c > minimum
                R[s_i, new_a_i] = u(c)
            elseif c < minimum && R[s_i, new_a_i] == -Inf && new_a_i == 1
                R[s_i, new_a_i] = u(1e-10) 
            end
        end
    end
    return (; r, w, sigma, beta, z_chain, a_min, a_max, a_size, a_vals, z_size,
            n, s_vals, s_i_vals, u, R, Q, y_med, τ_k, τ_c)
end

# Cálculo do equilíbrio (peso 0.3)
function find_equilibrium(; r0, beta, A, N, alpha, delta, a_max, cmin, Nivel_IR=0.825, tolEq = 1e-4, maxIterEq = 50, peso = 0.1)
    
    # Iniciando as variáveis
    K_supply = 0.0
    Nbar = N
    L = N
    w = 1.0
    dist = 1.0
    iter = 0
    τ_k = 0.15 # Initial guess for capital income tax rate
    #τ_c = 0.29  # Initial guess for consumption tax realmente
    τ_c = 0.0

    G_new, C, Y, G_model, T_model = 0, 0, 0, 0, 0
    c_vec, y_vec, ratio, c_obs = [], [], [], []
    y_bruto = Float64[]
    converged = false
    model_moments = [0,0,0]
    dist_prev = 0
    r_prev = 0
    resultados = Float64[]
    stationary_probs = Float64[]
    am = Float64[]
    results = Float64[]
    V0 = nothing
    total_consumption_tax = 0
    total_labor_income_tax = 0

    @show beta, cmin, Nivel_IR

    excDem = 0.0
    
    while dist > tolEq && iter < maxIterEq
        iter += 1

        # Calculando Demanda de Capital e salário inicial
        Kd = (alpha / (r0 + delta))^(1 / (1 - alpha)) * L # do matlab
        w = (1 - alpha) * (Kd / L)^alpha # do matlab

        
        # Create an instance of Household given the parameters
        am = Household(; beta, a_max, w, r = r0, cmin = cmin, τ_k = τ_k, τ_c = τ_c, Nivel_IR = Nivel_IR)
        
        # Solve the household's problem
        aiyagari_ddp = DiscreteDP(am.R, am.Q, am.beta)
        # results = QuantEcon.solve(aiyagari_ddp, VFI)

        results = isnothing(V0) ?
            QuantEcon.solve(aiyagari_ddp, VFI) :
            QuantEcon.solve(aiyagari_ddp, V0, VFI)

        V0 = results.v
        
        # Compute the stationary distribution
        stationary_probs = stationary_distributions(results.mc)[:, 1][1]

        # Compute aggregate capital K and labor supply L
        K_supply = dot(am.s_vals[:, 1], stationary_probs)
        # L = dot(am.s_vals[:,2], stationary_probs) # novo
        
        excDem = (K_supply - Kd)/((K_supply + Kd) / 2) # calculo do excesso de demanda para confirmação
        
        ### Government ###
        # Compute total consumption
        # total_consumption, consumption_vec = compute_total_consumption(results, am, stationary_probs, cmin, Nivel_IR)
        # total_consumption, c_vec, c_obs, y_vec, ratio = compute_total_consumption(results, am, stationary_probs, cmin, Nivel_IR)
        total_consumption, c_vec, y_vec, y_bruto, ratio, resultados = compute_total_consumption(results, am, stationary_probs, cmin, Nivel_IR)
        
        #total_consumption, consumption, y_vec, ratio
        C = total_consumption

        # Compute total consumption tax revenue
        total_consumption_tax = τ_c * total_consumption

        # Compute total capital income tax revenue
        total_capital_income_tax = τ_k * r0 * K_supply

        # Compute total labor income tax revenue
        total_labor_income_tax = compute_total_labor_income_tax(stationary_probs, am.s_vals, w, Nivel_IR)

        # Total tax revenue
        total_tax_revenue = total_consumption_tax + total_capital_income_tax + total_labor_income_tax

        # Compute total transfers
        total_transfers = compute_total_transfers(stationary_probs, am.s_vals, w, am.r, am.y_med)
        # total_transfers = compute_total_transfers(stationary_probs, am.s_vals, w, am.r, 1)

        G_model = total_transfers # despesa bruta do governo
        T_model = total_tax_revenue # carga tributária do governo
        
        G_new = total_tax_revenue - total_transfers
        #G_new + total_transfers = total_tax_revenue (mantendo o mesmo G_new)
        # @show G_model, G_new

        # Calculando variáveis agregadas
        Y = A * K_supply^alpha * L^(1 - alpha)
        I = delta * K_supply

        # Given that Y = C + I + G
        G_eq = Y - C - I
        
        frac_c = sum( (c_vec .< C) .* stationary_probs )

        # Compute model moments
        model_moments = [
        K_supply / Y,  # Capital-output ratio
        frac_c,  # Consumption ratio # mudar p/ target cmin ou outro
        # G_model / Y
        G_new / Y   # Debt ratio
        ]

        # Valores anteriores salvos
        dist_prev = dist
        r_prev    = r0
        
        # Novo r
        r_new = A * alpha * (L / K_supply)^(1 - alpha) - delta

        println("r_dif = $(r_new - r0) e G_dif = $(G_eq - G_new)")

        dist = maximum(abs.([r_new - r0, G_eq - G_new]))

        # Atualizando r0 como uma comb linear de r0 e r_new
        r0 = peso * r0 + (1 - peso) * r_new
        
        
        #dist = maximum(abs.([r_new - r0]))

        # Flexibilidade do peso 
        if iter > 1 && dist > dist_prev            # piorou
            peso  = min(0.99, peso + 0.05)         # aumenta amortecimento
            r0    = (r_prev + r0)/2                # retrocede meio passo
        else
            r0    = peso*r0 + (1-peso)*r_new       # passo normal (mudar r0 por r_prev?)
        end
        
        

        println("Iteration $iter: r = $r0, excDem = $excDem, dist = $dist")
        # println("Iteration $iter: r = $r0, D = $D, dist = $dist")
        
    end

    converged = dist ≤ tolEq 

    return r0, K_supply, L, w, excDem, G_new, C, Y, c_vec, y_vec, y_bruto, ratio, converged,
     iter, model_moments, resultados, stationary_probs, results, am, T_model, 
     G_model, total_consumption_tax, total_labor_income_tax, dist
    #return r0, K_supply, L, w, G_new, C, excDem, Y
end

# Transferências e T(y) ( antes m = 0.12 e ζ = 3.66)----
function ω(y_i; m = 0.05, ζ = 3.66, y_med)
    #@show m, ζ
    return (m * y_med * 2 * exp(-ζ * y_i / y_med)) / (1 + exp(-ζ * y_i / y_med))
end

function Taxa_y(y_e, λ, τ = 0.087)
    return max(y_e - λ * y_e^(1 - τ), 0)
end

function Taxa_y2(y_e, λ, y_med = 1, θ = 0.087)
    return max(exp(log(λ)*(y_e/y_med)^(-2θ)),0)
end
#-------------------------------------------------------

# Government Functions
function compute_total_consumption(results, am, stationary_probs, cmin, Nivel_IR)
    policy_indices = results.sigma  # Optimal policy indices
    #consumption = zeros(length(am.s_vals[:, 1]))
    n_states       = length(am.s_vals[:, 1])
    consumption = zeros(n_states)
    consumption_obs = zeros(n_states)
    y_vec  = zeros(n_states)
    y_bruto = zeros(n_states)
    ratio  = zeros(n_states)
    y_total = zeros(n_states)
    net_transfer = zeros(n_states)
    

    for s_i in 1:length(am.s_vals[:, 1])
        a = am.s_vals[s_i, 1]
        z = am.s_vals[s_i, 2]
        a_new = am.a_vals[policy_indices[s_i]]

        y_e = am.w * z
        y_k = am.r * a
        y_i = y_e + y_k
        y_total[s_i] = y_i
        

        T_value = Taxa_y(y_e, Nivel_IR)
        Tr_value = ω(y_i, y_med = am.y_med)
        # Tr_value = ω(y_i, y_med = 1)

        net_transfer[s_i] = Tr_value - T_value

        # Gross_income ≡ tw
        # tw = am.w * z + (1 + (1 - am.τ_k) * am.r) * a - T_value + Tr_value
        y_i = am.w * z + (1 + (1 - am.τ_k) * am.r) * a - T_value + Tr_value
        y_b = am.w * z + (1 + am.r) * a  # sem T e Tr
        c = (y_i - a_new) / (1 + am.τ_c)
        
        # gross_income = am.w * z + (1 + (1 - 0.15) * am.r) * a - T_value + Tr_value
        # c = (gross_income - a_new) / (1 + 0.29)

        #minimum = min(cmin, gross_income)

        consumption[s_i] = c
        y_vec[s_i] = y_i
        y_bruto[s_i] = y_b
        ratio[s_i] = c / y_i
    end

    total_consumption = dot(consumption, stationary_probs)
    # return total_consumption, consumption

    # Calculo do quintil
    # Ordenando por renda total
    ordem = sortperm(y_total)
    y_sorted = y_total[ordem]
    net_transfer_sorted = net_transfer[ordem]
    probs_sorted = stationary_probs[ordem]

    # Acumulando as probabilidades para os quintis
    cum_probs = cumsum(probs_sorted)
    quintis = [0.2, 0.4, 0.6, 0.8, 1.0]
    resultados = Float64[]

    start_idx = 1
    for q in quintis
        idx = findall(cum_probs .<= q)
        idx = setdiff(idx, 1:start_idx-1)
        pos_idx = idx[net_transfer_sorted[idx] .> 0.0]
        peso_quintil = sum(probs_sorted[idx])
        peso_pos = sum(probs_sorted[pos_idx])
        push!(resultados, peso_pos / peso_quintil)
        start_idx = isempty(idx) ? start_idx : maximum(idx) + 1
    end
    
    # return total_consumption, consumption, consumption_obs, y_vec, ratio
    return total_consumption, consumption, y_vec, y_bruto, ratio, resultados
end

function compute_total_labor_income_tax(stationary_probs, s_vals, w, Nivel_IR)
    total_tax = 0.0
    for s_i in 1:length(stationary_probs)
        z = s_vals[s_i, 2]
        y_e = w * z
        labor_tax = Taxa_y(y_e, Nivel_IR)
        total_tax += labor_tax * stationary_probs[s_i]
    end
    return total_tax
end

function compute_total_transfers(stationary_probs, s_vals, w, r, y_med)
    total_tr = 0.0
    for s_i in 1:length(stationary_probs)
        a = s_vals[s_i, 1]
        z = s_vals[s_i, 2]
        y_e = w * z
        y_k = r * a
        y_i = y_e + y_k
        transfer = ω(y_i; y_med=y_med)
        total_tr += transfer * stationary_probs[s_i]
    end
    return total_tr
end
#-------------------------------------------------------

function gini_weighted(y::AbstractVector, w::AbstractVector)
    idx = sortperm(y)                  # ordem crescente
    y, w  = y[idx], w[idx]             # redefinindo y e w

    w .= w ./ sum(w)                   # normalizando
    μ  = sum(y .* w)                   # média ponderada (dist. ss)
    Lcum  = cumsum(y .* w) ./ μ        # curva de Lorenz
    return 1.0 - sum(w[2:end] .* (Lcum[1:end-1] .+ Lcum[2:end]))
end

function top_share(y::AbstractVector, w::AbstractVector; p = 0.10, pr = 0.50)
    idx = sortperm(y)                  # ordenação crescente
    y, w  = y[idx], w[idx]

    w .= w ./ sum(w)                   # garante soma=1
    cumw  = cumsum(w)                  # CDF pelos pesos
    total = sum(y .* w)

    # Acumulado que não ultrapassa p (bottom 50)
    bot_idx    = findall(cumw .<= pr)
    inc_bottom = sum(y[bot_idx] .* w[bot_idx])

    # observações cujo acumulado supera 1-p (top decil, quintil etc.)
    top_idx = findall(cumw .> 1 - p)
    income_top = sum(y[top_idx] .* w[top_idx])
    # return income_top / total
    return income_top / total, inc_bottom / total
end


# Calculate supply of capital for a given r
function prices_to_capital_stock(r; beta, A, N, alpha, delta, a_max, cmin, Nivel_IR = 0.825)
    # Create an instance of Household given the parameters
    τ_k = 0.15
    τ_c = 0.29
    # Calculate the equilibrium wages
    w = A * (1 - alpha) * (A * alpha / (r + delta))^(alpha / (1 - alpha))
    
    # am = Household(; beta, a_max, w, r, cmin = cmin)
    am = Household(; beta, a_max, w, r, cmin = cmin, τ_k = τ_k, τ_c = τ_c, Nivel_IR = Nivel_IR)

    aiyagari_ddp = DiscreteDP(am.R, am.Q, am.beta)

    # Compute the optimal policy
    results = solve(aiyagari_ddp, VFI)

    # Compute the stationary distribution
    stationary_probs = stationary_distributions(results.mc)[:, 1][1]

    # Return K
    K = dot(am.s_vals[:, 1], stationary_probs)

    # Return capital
    return K
end

# Inverse Demand for capital
function r_inverse_demand(K; A, N, alpha, delta)
    return A * alpha * (N / K)^(1 - alpha) - delta
end

# Firms' parameters
A = 1
#N = 1
#alpha = 0.33
alpha = 0.45
# beta = 0.98
# delta = 0.05
delta = 0.06
# a_max = 40.0
a_max = 40 # testando grid menor após normalizar o z_grid por z_max
# cmin = 2.5

P = [0.417058294170583     0.1928807119288071    0.11768823117688232   0.08209179082091791   0.055894410558944105  0.0365963403659634    0.0376962303769623    0.0265973402659734    0.0207979202079792    0.012698730126987301;
    0.1881                0.2533                0.1382                0.0985                0.0821                0.0642                0.0618                0.0514                0.0381                0.0243;
    0.08009199080091992   0.1151884811518848    0.25807419258074193   0.13678632136786323   0.12038796120387961   0.08789121087891212   0.06919308069193081   0.06699330066993302   0.03859614038596141   0.026797320267973206;
    0.0524                0.0887                0.1456                0.2794                0.1366                0.0929                0.0746                0.0630                0.0431                0.0237;
    0.035896410358964105  0.07699230076992301   0.09199080091990801   0.1272872712728727    0.21787821217878214   0.1512848715128487    0.10718928107189282   0.08929107089291072   0.06669333066693331   0.0354964503549645;
    0.0305                0.0546                0.0813                0.0851                0.1432                0.2406                0.1487                0.1075                0.0737                0.0348;
    0.026700000000000005  0.05500000000000001   0.05650000000000001   0.06990000000000002   0.09230000000000001   0.15360000000000001   0.23490000000000003   0.15060000000000004   0.10700000000000001   0.053500000000000006;
    0.022097790220977905  0.04609539046095391   0.0524947505249475    0.05649435056494351   0.07939206079392061   0.09309069093090691   0.1498850114988501    0.25407459254074594   0.15638436156384364   0.08999100089991001;
    0.0176017601760176    0.0328032803280328    0.040604060406040594  0.0477047704770477    0.059005900590058995  0.0742074207420742    0.0912091209120912    0.1549154915491549    0.3058305830583058    0.1761176117611761;
    0.007301460292058411  0.020404080816163232  0.02090418083616723   0.027505501100220042  0.03250650130026005   0.031906381276255245  0.05251050210042008   0.08681736347269454   0.1807361472294459    0.5394078815763153]
n_grid = [0.17; 0.42; 0.70; 0.76; 0.83; 1.0; 1.22; 1.55; 2.20; 5.53]

inv_P = P^1000
inv_P = inv_P[1,:]

# Compute the average productivity
N = sum(n_grid.*inv_P)
# N = 1

n_grid .= n_grid ./ N
N = sum(n_grid.*inv_P)

### para Q fixo ###
a_min  = 1e-10
a_size = 200
a_vals = range(a_min, a_max; length = a_size)

z_chain = MarkovChain(P, n_grid)
z_size  = length(z_chain.state_values)

n = a_size * z_size
s_i_vals  = gridmake(1:a_size, 1:z_size)
### ### ### ### ###
const Q_FIX = build_Q(a_size, z_size, s_i_vals, z_chain)


### Teste de gráfico do equilíbrio
# Create a grid of r values at which to compute demand and supply of capital
r_vals = range(0.005, 0.05, length = 20)

# Compute supply of capital
k_vals = prices_to_capital_stock.(r_vals; A, N, alpha, beta = 0.98, delta, a_max,
                                 cmin = 1.5, Nivel_IR = 0.525)
r_inverse_demand_vals = r_inverse_demand.(k_vals; A, N, alpha, delta)

# Ploting against demand for capital by firms
labels = ["demand for capital" "supply of capital"]
plot(k_vals, [r_inverse_demand_vals r_vals], label = labels, lw = 2,
     alpha = 0.6, title = "novo com gov: beta = 0.98, cmin = 1.5, \n λ = .525, z_n = 10", dpi=1000)
plot!(xlabel = "capital", ylabel = "interest rate")
png("C:\\Users\\lucinaldosouza\\Pictures\\novo_gov_n7_distr")

# Variáveis de equilíbrio
r_eq, K_eq, L_eq, w_eq, excDem, G_eq, C_eq, Y_eq, c_vec,
 y_vec, y_bruto, ratio, converged, iter, model_moments, resultados, 
 stationary_probs, results, am, T_model, G_model, Cgov, IRgov, dist = find_equilibrium(r0 = 0.005, beta = 0.875480813067229, 
 cmin = 1.3751372838375977, A=1, N = N, alpha = 0.45,
 delta = 0.06, a_max = 40, Nivel_IR = 0.9999)


# Para o equilíbrio
am = Household(; beta= 0.875480813067229, a_max = 40, w = w_eq, r = r_eq, cmin = 1.3751372838375977, τ_k = 0.15, τ_c = 0.29, Nivel_IR = 0.9999)
aiyagari_ddp = DiscreteDP(am.R, am.Q, am.beta)
results = QuantEcon.solve(aiyagari_ddp, VFI)
stationary_probs = stationary_distributions(results.mc)[:, 1][1]


frac_c = sum( (c_vec .< C_eq) .* stationary_probs )

# Compute model moments
model_moments = [
    K_eq / Y_eq,  # Capital-output ratio
    frac_c,  # Consumption ratio # mudar p/ target cmin ou outro
    G_eq / Y_eq   # Debt ratio
]


#
z_vals = am.z_chain.state_values

(; z_size, a_size, n, a_vals) = am
a_star_steady = reshape([am.a_vals[results.sigma[s_i]] for s_i in 1:am.n], am.a_size, am.z_size)

# Plot capital
plot(a_vals, a_star_steady[:, 1], label = "z = $(z_vals[1])", xlabel = "Current assets (a_vals)",
 ylabel = "Next period assets (a_star)",
  title = "Next Period Capital vs. Current Capital \n Novo com beta = 0.96 e cmin = 3.5", 
  lw = 2, dpi = 1000)
plot!(a_vals, a_star_steady[:, 2], label = "z = $(z_vals[2])", lw = 2)
plot!(a_vals, a_star_steady[:, 3], label = "z = $(z_vals[3])", lw = 2)
plot!(a_vals, a_star_steady[:, 4], label = "z = $(z_vals[4])", lw = 2)
plot!(a_vals, a_star_steady[:, 5], label = "z = $(z_vals[5])", lw = 2)
plot!(a_vals, a_star_steady[:, 6], label = "z = $(z_vals[6])", lw = 2)
plot!(a_vals, a_star_steady[:, 7], label = "z = $(z_vals[7])", lw = 2)
plot!(a_vals, a_star_steady[:, 8], label = "z = $(z_vals[8])", lw = 2)
plot!(a_vals, a_star_steady[:, 9], label = "z = $(z_vals[9])", lw = 2)
plot!(a_vals, a_star_steady[:, 10], label = "z = $(z_vals[10])", lw = 2)
plot!(a_vals, a_vals, label = "", color = :black, linestyle = :dash)
plot!(xlabel = "current assets", ylabel = "next period assets", grid = true)

png("C:\\Users\\lucinaldosouza\\Pictures\\ativos_logkappa2")

# Para o gráfico da distribuição estacionária
asset_levels = am.s_vals[:, 1]
K = dot(am.s_vals[:, 1], stationary_probs)

# Plotando como gráfico de linha
plot(
    asset_levels,
    # c_vec,
    stationary_probs,
    seriestype = :line,               # Line plot to visualize distribution
    xlabel = "Nível de ativos",   # Label for x-axis
    ylabel = "Probabilidade",           # Label for y-axis
    title = "Distribuição estacionária de ativos de capital", # Plot title
    #ylim = (0.0, 0.02),
    legend = false,                    # No legend needed for a single series
    dpi=1000
)
png("C:\\Users\\lucinaldosouza\\Pictures\\grafico_dist_ativos_gamma1")

# Para a renda
gini = gini_weighted(y_vec, stationary_probs)
gini_weighted(c_vec, stationary_probs)
top10, bottom50 = top_share(y_bruto, stationary_probs; p = 0.10, pr = 0.50)
top10, bottom50 = top_share(c_vec, stationary_probs; p = 0.10, pr = 0.50)
top10, bottom50 = top_share(asset_levels, stationary_probs; p = 0.10, pr = 0.50)

# Para o consumo
top_share(c_vec, stationary_probs; p = 0.10, pr = 0.10)
top_share(c_vec, stationary_probs; p = 0.01, pr = 0.50)

# top
function share_top(x::AbstractVector, w::AbstractVector; p = 0.10)
    idx     = sortperm(x)              # ordem crescente
    x, w    = x[idx], copy(w[idx]) ./ sum(w)
    cdf_w   = cumsum(w)

    # posição onde começa o top-p
    thr     = 1 - p
    i       = searchsortedfirst(cdf_w, thr)

    # fração do peso nessa observação que realmente fica no top
    frac_i  = (cdf_w[i] - thr) / w[i]

    share   = sum(x[i+1:end] .* w[i+1:end]) + x[i] * w[i] * frac_i
    return share / sum(x .* w)
end

# bottom
function share_bottom(x::AbstractVector, w::AbstractVector; q = 0.50)
    idx     = sortperm(x)              # ordem crescente
    x, w    = x[idx], copy(w[idx]) ./ sum(w)
    cdf_w   = cumsum(w)

    # posição onde termina o bottom-q
    i       = searchsortedfirst(cdf_w, q)

    frac_i  = (q - (i > 1 ? cdf_w[i-1] : 0.0)) / w[i]

    share   = sum(x[1:i-1] .* w[1:i-1]) + x[i] * w[i] * frac_i
    return share / sum(x .* w)
end

share_top(y_vec, stationary_probs; p = 0.10)
share_bottom(y_vec, stationary_probs; q = 0.50)

share_top(c_vec, stationary_probs; p = 0.10)
share_bottom(c_vec, stationary_probs; q = 0.50)

sum( (c_vec .< 1.3751372838375977) .* stationary_probs )

# wealth-income ratio
K_eq / (Y_eq - delta * K_eq)

# Taxa real líquida de imposto
r_eq * (1 - 0.15)


# Welfare
dot(results.v, stationary_probs)

########### Calibração ##############
using Optim
# beta 0.88 ate 0.98 com 11 pontos

lower_bounds = [0.83, 0.2, 0.3]  # beta, cs e lambda_IR mínimo
upper_bounds = [0.99, 2.0, 1.0]  # beta, cs e lambda_IR máximo


data_moments = [2.55, 0.80, 0.196]
W = LinearAlgebra.I(length(data_moments))

params_initial = [0.88806, 1.44798, 0.768577]

function gmm_objective(params)
    
    best_obj = 100
    beta, cs, lambda = params
    best_params = []
    best_moments = []
    # @show params
    r0 = 0.005  # Initial guess for interest rate

    A = 1
    N = 1.0  # Valor calculado de Nbar anteriormente
    alpha = 0.45
    delta = 0.06
    a_max = 40

    # Resolvendo pro eq dado os parâmetros atuais
    r_eq, K_eq, L_eq, w_eq, excDem, G_eq, C_eq, Y_eq, c_vec, 
    y_vec, y_bruto, ratio, converged, iter, model_moments, resultados, 
    stationary_probs, results, am, T_model, Gmodel, Cgov, IRgov, dist = find_equilibrium(; r0, beta, 
    A, N, alpha, delta, a_max, cmin=cs, Nivel_IR=lambda)
    
    # Penalidade por não convergência
    penalty = converged ? 0.0 : 5.0 + ℯ^(dist)

    # Fração de indivíduos abaixo da média (< C_eq)
    frac_c = sum( (c_vec .< C_eq) .* stationary_probs )

    @show model_moments
    
    # fução objetivo do GMM
    diff = (model_moments .- data_moments)./data_moments
    # return diff' * W * diff
    Q = dot(diff, W * diff) + penalty^2
    @show obj = log(Q) # pegar log(obj)

    if obj < best_obj
        best_obj = obj
        best_params = params
        best_moments = model_moments
    end
    

    return obj
end

gmm_objective(params_initial)


# Otimizando por Fminbox
result = optimize(gmm_objective, lower_bounds, upper_bounds, params_initial, Fminbox())
println("Calibrated parameters: ", result.minimizer)

## Testando LBFGS()
using ForwardDiff, Optim

inner = LBFGS();          # algoritmo rápido em alta dimensão
outer = Fminbox(inner)
#

opt = Optim.Options(g_tol = 1e-5, f_tol = 1e-6, iterations = 120)

result = optimize(gmm_objective, lower_bounds, upper_bounds,
                  params_initial, outer, opt; autodiff = :finite)
#
# :finite ou forward

############### Gráfico de perdas ###############
using Plots

# ────────────────────────────────────────────────────────────────────────────────
# Parâmetros best
params_best  = [0.875480813067229, 1.3751372838375977, 0.9999]
best_loss = -2.318672633822192
# Gerando grid em torno de cada parâmetro (11 pontos, passo h)
function profile_axis(idx; h=0.02, n=11)
    p0         = copy(params_best)
    half       = (n-1) ÷ 2
    grid       = p0[idx] .+ (-half:half) .* h
    losses     = Float64[]
    for x in grid
        p0[idx] = x
        push!(losses, gmm_objective(p0))
    end
    return grid, losses
end

# Teste p/ β
β_grid, β_loss = profile_axis(1; h = 0.020, n = 11)
plot(β_grid, β_loss, marker = :o, lw = 2,
     xlabel = "β", ylabel = "log(Q)",
     title  = "Perfil da função-objetivo em torno de β",
     legend = false,
     dpi = 1000)
#png("C:\\Users\\lucinaldosouza\\Pictures\\graf_beta_perda1")

# Para cmin e λ (nível_IR)
cmin_grid, cmin_loss = profile_axis(2; h = 0.05, n = 11)
λ_grid,    λ_loss    = profile_axis(3; h = 0.05, n = 11)

plot(cmin_grid, cmin_loss, marker = :square, lw = 2,
     xlabel = "cₘᵢₙ", ylabel = "log(Q)",
     title  = "Perfil da função-objetivo em torno de cₘᵢₙ",
     legend = false,
     dpi = 1000)
# png("C:\\Users\\lucinaldosouza\\Pictures\\graf_cmin_perda1")

plot(λ_grid, λ_loss, marker = :diamond, lw = 2,
     xlabel = "λ (nível_IR)", ylabel = "log(Q)",
     title  = "Perfil da função-objetivo em torno de λ",
     legend = false,
     #ylim = (0.0, 1.0),
     dpi = 1000)
png("C:\\Users\\lucinaldosouza\\Pictures\\graf_3d")

# ---------------------------------------------------------------------------------
# --------------------- Utilizando um Otimizador Global ---------------------------
using BlackBoxOptim

search_range = [(0.79, 0.99), (1.10, 1.50), (0.700, 0.9999)]


# Otimização
result_bb = bboptimize(gmm_objective,
                       SearchRange = search_range,
                       NumDimensions = 3,
                       # Method = :de_rand_1_bin, # algorítimo
                       Method = :generating_set_search,
                       MaxTime = 25200.0,      # tempo
                       TraceMode = :verbose)

println("Melhores parâmetros (BlackBoxOptim): ", best_candidate(result_bb))
println("Melhor perda (BlackBoxOptim): ", best_fitness(result_bb))


# ---------- Visualizando a Função Objetivo ----------------
using Plots

# Fixe lambda no limite que ele atingiu
lambda_fixed = 0.9999

# grid para beta e cs
beta_grid = 0.85:0.01:0.95
cs_grid = 1.3:0.02:1.5

# armazenando os resultados em uma matriz
loss_surface = zeros(length(beta_grid), length(cs_grid))

# perda para cada ponto do grid
for (i, beta) in enumerate(beta_grid)
    for (j, cs) in enumerate(cs_grid)
        params = [beta, cs, lambda_fixed]
        loss_surface[i, j] = gmm_objective(params)
        println("Calculando para beta=$beta, cs=$cs...")
    end
end

# plot da superfície
plot(beta_grid, cs_grid, loss_surface', st=:surface, camera=(-10, 10),

     xlabel="Beta", ylabel="cs", zlabel="Função Objetivo (log(Q))", dpi = 1000)
