# UBI
# Valor flat do RBU = θ_rbu * (Y_eq - delta*K_eq) com θ_rbu = 0.10

using LinearAlgebra, Statistics, StatsBase
using LaTeXStrings, Plots, QuantEcon
using BenchmarkTools
using SparseArrays
using NLsolve

#= 
=#

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

    #= for s in 1:n
        cur_z = s_i_vals[s, 2]                # índice de produtividade atual

        for a_idx in 1:a_size                 # cada ação = ativo escolhido
            # estados de destino têm SEMPRE esse ativo
            for z′ in 1:z_size
                next_s = (a_idx - 1) * z_size + z′   # (ativo = a_idx , z = z′)
                Q[s, a_idx, next_s] = z_chain.p[cur_z, z′]
            end
        end
    end =#
    return Q
end

# Construção das matrizes de transição e utilidade
function Household(; r = 0.01,
                   w = 1.0,
                   sigma = 1.0,
                   beta = 0.96,
                   τ_k = 0, 
                   τ_c = 0.29,
                   cmin = 0,
                   Nivel_IR,
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
    T(y_e, λ, τ = 0.087) = max(y_e - λ * y_e^(1 - τ), 0) - ω()

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
            Tr = T(y_e, Nivel_IR)
            

            # c = w * z + (1 + r) * a - a_new
            c = (w * z + (1 + (1 - τ_k)*r)*a - Tr - a_new) / (1 + τ_c)
            # c = w * z + (1 + (1 - τ_k)*r)*a - a_new
            # tw = w * z + (1 + r) * a
            tw = w * z - Tr + (1 + (1 - τ_k) * r) * a

            if c < cmin
                tw_rel = (cmin - c) / cmin
            else
                tw_rel = 1
            end

            minimum = min(cmin, tw * tw_rel)

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

# Cálculo do equilíbrio
function find_equilibrium(; r0, t0, beta, A, N, alpha, delta, a_max, cmin, tolEq = 1e-4, maxIterEq = 50, peso = 0.1)
    
    # Iniciando as variáveis
    K_supply = 0.0
    L = N
    w = 1.0
    dist = 1.0
    iter = 0
    τ_c = 0.29
    τ_k = 0.15
    Nivel_IR = t0
    #t0 = 0.825
    sumw = 0

    G_new, C, Y, G_model, T_model = 0, 0, 0, 0, 0
    c_vec, y_vec, ratio = [], [], []
    y_bruto = Float64[]
    converged = false
    model_moments = [0,0,0]
    dist_prev = 0
    r_prev = 0
    t_prev = 0
    resultados = Float64[]
    stationary_probs = Float64[]
    am = Float64[]
    results = Float64[]
    contagem = Float64[]
    cont = 0
    V0 = nothing
    S = 0.0
    cond = 0.0

    @show beta, cmin, t0

    excDem = 0.0
    
    while dist > tolEq && iter < maxIterEq
        iter += 1

        # Calculando Demanda de Capital e salário inicial
        Kd = (alpha / (r0 + delta))^(1 / (1 - alpha)) * L # do matlab
        w = (1 - alpha) * (Kd / L)^alpha # do matlab

        @show t0
        # Create an instance of Household given the parameters
        am = Household(; beta, a_max, w, r = r0, cmin = cmin, τ_k = τ_k, Nivel_IR = t0)
        
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
        total_consumption, c_vec, y_vec, y_bruto, ratio, resultados = compute_total_consumption(results, am, stationary_probs, cmin, t0)
        
        #total_consumption, consumption, y_vec, ratio
        C = total_consumption

        # Compute total consumption tax revenue
        total_consumption_tax = τ_c * total_consumption

        # Compute total capital income tax revenue
        total_capital_income_tax = τ_k * r0 * K_supply

        # Compute total labor income tax revenue
        total_labor_income_tax, S, cond, sumw, cont = compute_total_labor_income_tax(stationary_probs, am.s_vals, w, t0)
        push!(contagem, cont)
        @show cond
        # Total tax revenue
        total_tax_revenue = total_consumption_tax + total_capital_income_tax + total_labor_income_tax

        # Compute total transfers
        total_transfers = compute_total_transfers(stationary_probs, am.s_vals, w, am.r, am.y_med)
        # total_transfers = compute_total_transfers(stationary_probs, am.s_vals, w, am.r, 1)

        G_model = total_transfers # despesa bruta do governo
        T_model = total_tax_revenue # carga tributária do governo
        
        G_new = total_tax_revenue - total_transfers
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
        frac_c,  # % de consumidores
        G_new / Y   # Debt ratio
        ]

        # Valores anteriores salvos
        dist_prev = dist
        r_prev    = r0
        t_prev    = t0
        
        # Novo r
        r_new = A * alpha * (L / K_supply)^(1 - alpha) - delta

        # G + ubi  = K * r * τ_k + C * τ_c + ∑T(y)*Π
        # ∑(∑yπ = ̄y = 1 - λy^(1-τ))*Π = G - total_consumption_tax - total_capital_income_tax
        # 
        # dif = G_model - (total_tax_revenue)
        # dif = G_model - T_model
        # ta = (1 - (G_model - total_consumption_tax - total_capital_income_tax))
        # G = 0.5452724712772316 sendo o gasto do governo
        # ta = min((1 - (0.5452724712772316 + G_model - total_consumption_tax - total_capital_income_tax))/S, cond)
        ta = min((sumw - 0.5452724712772316 - G_model + total_consumption_tax + total_capital_income_tax)/S, cond)
        # ta = min((1 - 0.5452724712772316 - G_model + total_consumption_tax + total_capital_income_tax)/S, cond) 
        

        println("r_dif = $(r_new - r0) e ta_dif = $(ta - t0)")
        
        dist = maximum(abs.([r_new - r0, ta - t0]))
        

        r0 = peso*r0 + (1-peso)*r_new
        # t0 = min(peso*t0 + (1-peso)*ta, 0.9999)
        # t0 = peso*t0 + (1-peso)*ta
        t0 = min(peso*t0 + (1-peso)*ta, cond)
        
        
        # Flexibilidade do peso 
        if iter > 1 && abs(dist) > abs(dist_prev)  # piorou
            peso  = min(0.99, peso + 0.05)         # aumenta amortecimento
            r0    = (r_prev + r0)/2                # retrocede meio passo
            t0    = (t_prev + t0)/2 
        else
            r0 = peso*r_prev + (1-peso)*r_new       # passo normal
            # t0 = peso*t_prev + (1-peso)*ta;
            t0 = min(peso*t0 + (1-peso)*ta, cond)
        end
        
    
        println("Iteration $iter: r = $r0, t0 = $t0, dist = $dist")
        
    end

    converged = dist ≤ tolEq 

    return r0, t0, K_supply, L, w, excDem, G_new, C, Y, c_vec, y_vec, y_bruto, 
    ratio, converged, iter, model_moments, resultados, stationary_probs, 
    results, am, T_model, G_model, sumw, contagem
end

# min 0.001 e max 0.11
# Transferências e T(y) (θ_rbu = 0.10 e Y_eq - delta*K_eq = 1.86) ----
ω() = 0.12 * 2.21 # 0.10 * 2.21 = 0.22 com Y_eq = 2.21

function Taxa_y(y_e, λ, τ = 0.087)
    return max(y_e - λ * y_e^(1 - τ), 0)
end
#-------------------------------------------------------

# Government Functions
function compute_total_consumption(results, am, stationary_probs, cmin, Nivel_IR)
    policy_indices = results.sigma  # Optimal policy indices
    n_states       = length(am.s_vals[:, 1])
    consumption = zeros(n_states)
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
        Tr_value = ω()

        net_transfer[s_i] = Tr_value - T_value

        # Gross_income ≡ tw
        y_i = am.w * z + (1 + (1 - am.τ_k) * am.r) * a - T_value + Tr_value
        y_b = am.w * z + (1 + am.r) * a  # sem T e Tr
        c = (y_i - a_new) / (1 + am.τ_c)

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
    S = 0.0
    sumw = 0.0
    cont = 0.0
    for s_i in 1:length(stationary_probs)
        z = s_vals[s_i, 2]
        y_e = w * z
        sumw += y_e * stationary_probs[s_i]
        y_s = y_e^(1 - 0.087)
        # c_s = y_e^(0.087)
        labor_tax = Taxa_y(y_e, Nivel_IR)
        total_tax += labor_tax * stationary_probs[s_i]
        S += y_s * stationary_probs[s_i]
        # cond += c_s
        if (y_e + Nivel_IR * y_s) < 0
           cont += 1
        end
    end
    #return total_tax, S, cond
    return total_tax, S, sumw/S, sumw, cont
end

function compute_total_transfers(stationary_probs, s_vals, w, r, y_med)
    total_tr = 0.0
    for s_i in 1:length(stationary_probs)
        a = s_vals[s_i, 1]
        z = s_vals[s_i, 2]
        y_e = w * z
        y_k = r * a
        y_i = y_e + y_k
        transfer = ω()
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

### Firms' parameters
A = 1
alpha = 0.45
delta = 0.06
a_max = 40

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


r_eq, t_eq, K_eq, L_eq, w_eq, excDem, G_eq, C_eq, Y_eq, c_vec,
 y_vec, y_bruto, ratio, converged, iter, model_moments, resultados, 
 stationary_probs, results, am, T_model, G_model, sumw, contagem = find_equilibrium(r0 = 0.005, t0 = 0.800, 
 beta = 0.875480813067229, cmin = 1.3751372838375977, A=1, N = N, 
 alpha = 0.45, delta = 0.06, a_max = 40)

# Com w = 0.11
# r_eq = 0.14
# K_eq = 4.54
# t_eq = 0.71
# G_eq = 0.63

# Com w = 0.22

dot(results.v, stationary_probs)

gini = gini_weighted(y_vec, stationary_probs)
gini = gini_weighted(c_vec, stationary_probs)


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


top10, bottom50 = top_share(y_bruto, stationary_probs; p = 0.10, pr = 0.50)
# para 0% (0.2334260299594888, 0.20413303603261612)
# para 6% (0.2538652750848566, 0.15847636460279593)
# para 12% (0.2673369283643604, 0.13133936966296786)

top_share(c_vec, stationary_probs; p = 0.10, pr = 0.50)
# (0.20137401165542898, 0.33306851329319315)
# (0.19316556004959962, 0.3226568466506214)
#