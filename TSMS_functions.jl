function gen_stmat(TSMS)
    NO_ST, DMNSION = TSMS.NO_ST, TSMS.DMNSION
    
    st_mat = zeros(DMNSION, NO_ST)
    j = 1
    st2 = 1
    while st2 <=3
        st1 = 1
        while st1 <= 3
            st = 1
            while st <= 3
                st_mat[j,1] = st2
                st_mat[j,2] = st1
                st_mat[j,3] = st

                j = j +1 
                st = st+1
            end
            st1 = st1+1
        end
        st2 = st2 +1 
    end

    s1t_mat = zeros(DMNSION, NO_ST)
    s2t_mat = zeros(DMNSION, NO_ST)
    s3t_mat = zeros(DMNSION, NO_ST)

    j = 1 
    while j <= DMNSION
        i = 1
        while i <= NO_ST
            if st_mat[j,i]==1
                s1t_mat[j,i] = 1
            end
            i = i+1
        end
        j = j+1
    end

    j = 1 
    while j <= DMNSION
        i = 1 
        while i <= NO_ST
            if  st_mat[j,i]==2
                s2t_mat[j,i] = 1
            end
            i = i+1 
        end
        j = j+1
    end

    j = 1 
    while j <= DMNSION
        i = 1 
        while i <= NO_ST
            if  st_mat[j,i]==3
                s3t_mat[j,i] = 1
            end
            i = i+1 
        end
        j = j+1
    end
    
    return s1t_mat, s2t_mat, s3t_mat
    
end

###
function bingen(p0::Float64, p1::Float64, m::Int64)
    pr0=p0/(p0+p1) # prob(s=0)
    u= rand(Uniform(), (m, 1))
    s= u .>= pr0
    return s
end

###
function rg2(a::Float64)
    b=a-1.0
    c=3*a - 0.75
    accept=0
    x=0
    while accept == 0
        u=rand(Uniform())
        v=rand(Uniform())
        w=u * (1-u)
        y=sqrt(c / w) * (u - 0.5)
        x=b+y
        
        if x >= 0
            z=64 * (w^3) * (v^2)
            accept= z <= ( 1-(2*y^2)/x )
            
            if accept == 0
                accept = log(z) <= 2*(b*log(x/b) - y)
            end
        end
    end
    return x
end

function rg1(a::Float64)
    x=0
    if a > 1
        x=rg2(a)
    elseif a < 1
        a= a+1.0
        u= rand(Uniform())
        x = rg2(a)*u^(1/a)
        
    elseif a == 1
        x = -log(rand(Uniform()))
    end
    
    return x
end
        

function rndc(a::Int)
    a=a/2
    w=rg1(a)
    return w * 2
end


###
function rndb(a::Float64,b::Float64)
    a1n = rg1(a)
    a1d = rg1(b)
    return a1n / (a1n + a1d)
end


###
function gen_st(TSMS, P_TT::Array{Float64,1}, ΣTT::Array{Float64,1}, ΦTT::Array{Float64,1}, μTT::Array{Float64,1})
    t0, yy, x_mat, DMNSION  = TSMS.t0, TSMS.yy, TSMS.x_mat, TSMS.DMNSION
    s1t_mat, s2t_mat, s3t_mat = gen_stmat(TSMS)
    
    STT = Any[]
    
    TRYIT=1    
    while TRYIT == 1    
    
        tstar=t0-2
        flt_pr=zeros(Float64, 3, tstar)  

        prmtr= [P_TT[1]; P_TT[2]; P_TT[4]; P_TT[5]; P_TT[7]; P_TT[8];  
                ΦTT[1]; ΦTT[2]; 
                ΣTT[1]; ΣTT[2]; ΣTT[3]
                μTT[1]; μTT[2]; μTT[3]]         
                    

        #INITIAL PROB. Pr[S0/Y0]
        pr_tr=zeros(Float64,3,3)
        pr_tr[1:2,1]=prmtr[1:2]
        pr_tr[1:2,2]=prmtr[3:4]
        pr_tr[1:2,3]=prmtr[5:6]

        pr_tr[3,1]=1-sum(prmtr[1:2])
        pr_tr[3,2]=1-sum(prmtr[3:4])
        pr_tr[3,3]=1-sum(prmtr[5:6])    

        
        ϕ = prmtr[7:8,1]
        var = prmtr[9:11,1]
        μ = prmtr[12:14,1]

        var_l = var[1,1]*s1t_mat[:,size(s1t_mat,2)] + var[2,1]*s2t_mat[:,size(s2t_mat,2)] + var[3,1]*s3t_mat[:,size(s3t_mat,2)]
        μ_mat = μ[1,1]*s1t_mat + μ[2,1]*s2t_mat + μ[3,1]*s3t_mat

        A=ones(Float64, 4, 3)
        A[1:3,1:3]=Matrix{Float64}(I, 3, 3) - pr_tr


        EN=[0,0,0,1] 

        PROB__T=inv(A'A)*A'EN  #PR[S_t=0]|PR[S_t=1], 3x1 UNCONDITIONAL PROBABILITIES



        FLTPR0=PROB__T

        PR_TRF0=reshape(pr_tr,size(pr_tr,1)*size(pr_tr,2),1)
        PR_TRF =[PR_TRF0; PR_TRF0; PR_TRF0]


        PROB__T=[transpose(PROB__T); transpose(PROB__T); transpose(PROB__T)]      #@3^2 x 1@

        PROB__T = reshape(PROB__T,size(PROB__T,1)*size(PROB__T,2),1) .* PR_TRF0 

        FLTPR1=PROB__T[1:3]+PROB__T[4:6]+PROB__T[7:9]

        PROB__=[transpose(PROB__T); transpose(PROB__T); transpose(PROB__T)]                #@3^3 x 1@
        PROB__ = reshape(PROB__,size(PROB__,1)*size(PROB__,2),1)


        j=1
        while j<= tstar

            F_CAST1=(yy[j,1] - transpose(x_mat[j,:])*ϕ) * ones(Float64,DMNSION,1) - (μ_mat[:,3]  - μ_mat[:,2] * ϕ[1,1] - μ_mat[:,1]*ϕ[2,1])

            PROB_DD=PR_TRF .* PROB__   #@Pr[S_t,....S_{t-4}|Y_{t-1}]

            PR_VL = (1 ./ sqrt.(2 .* π .* var_l)) .* exp.(-0.5 * F_CAST1 .* F_CAST1 ./ var_l) .* PROB_DD;

            PR_VAL=sum(PR_VL)
            PRO_=PR_VL/PR_VAL #Pr[S1t,S2t,S1tl,S2tl|Y_t-1]


            PROB__T=PRO_[1:Int(DMNSION/3), 1] + PRO_[Int(DMNSION/3+1):Int(DMNSION*2/3), 1] + PRO_[Int(DMNSION*2/3+1):DMNSION, 1] #@Pr[S_t, S_{t-1}| Y_t]@    

            flt_pr[:,j] = PROB__T[1:3] + PROB__T[4:6] + PROB__T[7:9]   


            PROB__=[transpose(PROB__T); transpose(PROB__T); transpose(PROB__T)]                #@3^3 x 1@
            PROB__ = reshape(PROB__,size(PROB__,1)*size(PROB__,2),1)       

            j=j+1

        end


        flt_pr= [FLTPR0 FLTPR1 flt_pr]
        tstar=tstar+2


        #GENERATE S_T BASED ON CARTER & KOHN (1994)

        s_t=zeros(Float64, tstar,1)
        s_t___=zeros(Float64, tstar,3)

        s_tmp= bingen(flt_pr[1, tstar], flt_pr[2, tstar]+ flt_pr[3, tstar],1) 
#         print("\n s_tmp:", s_tmp)
        if s_tmp[1]==0 
            s_t[tstar,1]=1
        else 
            s_tmp1= bingen(flt_pr[2, tstar], flt_pr[3, tstar],1) 
#             print("\n s_tmp1:", s_tmp1)
            if s_tmp1[1] == 0 
                s_t[tstar,1]=2
            else
                s_t[tstar,1]=3
            end
        end

        s_t___[tstar,Int(s_t[tstar,1])]=1


        P1_1TT, P1_2TT, P1_3TT, P2_1TT, P2_2TT, P2_3TT, P3_1TT, P3_2TT, P3_3TT = P_TT[1], P_TT[2], P_TT[3], P_TT[4], P_TT[5], P_TT[6], P_TT[7], P_TT[8], P_TT[9]
        j=tstar-1    
        while j>=1

            if s_t[j+1,1]==1
                p1=P1_1TT*flt_pr[1,j]
                p2=P2_1TT*flt_pr[2,j]
                p3=P3_1TT*flt_pr[3,j]     

            elseif s_t[j+1,1]==2
                p1=P1_2TT*flt_pr[1,j]
                p2=P2_2TT*flt_pr[2,j]
                p3=P3_2TT*flt_pr[3,j]

            elseif s_t[j+1,1]==3
                p1=P1_3TT*flt_pr[1,j]
                p2=P2_3TT*flt_pr[2,j]
                p3=P3_3TT*flt_pr[3,j]  
            end

            s_tmp= bingen(p1,p2+p3,1) 
            if s_tmp[1]==0
                s_t[j,1]=1
            else 
                s_tmp1= bingen(p2,p3,1)
                if s_tmp1[1]==0
                    s_t[j,1]=2
                else 
                    s_t[j,1]=3
                end
            end

            s_t___[j, Int(s_t[j,1])]=1

            j=j-1

        end
      
        STT = [s_t, s_t___[:,1], s_t___[:,2], s_t___[:,3]]
    
        if sum(STT[2][3:t0])<4 || sum(STT[3][3:t0])<4 || sum(STT[4][3:t0])<4
            TRYIT=1

        else
            TRYIT=0
        end   
    end
    
    return STT    
end


###
function gen_phi(TSMS, MUTT, SIGMATT::Array{Float64,1})
    yy0, t0 = TSMS.yy0, TSMS.t0
    R0_, t0_ = TSMS.R0_, TSMS.t0_
    
    YSTARR=yy0- MUTT
    YSTAR= YSTARR[3:t0,1] ./ sqrt.(SIGMATT[3:t0])
    XSTAR=[YSTARR[2:t0-1,1]./sqrt.(SIGMATT[3:t0])  YSTARR[1:t0-2,1]./sqrt.(SIGMATT[3:t0])]
    
    V= inv(R0_ + XSTAR'XSTAR)
    PHI =  V * (R0_ * t0_ + XSTAR'YSTAR)
    C = cholesky(Hermitian(V))
    accept = 0
    counter = 0
    PHI_G = zeros(size(PHI,1))
    while accept == 0
        PHI_G = PHI + C.L*rand(Normal(), 2) 
        coef = [1;-(PHI_G)]
        root = roots(Polynomial(coef))
        rootmod = abs.(root)
        if minimum(rootmod) >= 1.0001
            accept = 1
        else 
            counter += 1
            accept = 0
        end
    end
    return [PHI_G[1],PHI_G[2]]    
end


###
function gen_mu(TSMS, ΦTT::Array{Float64,1}, STT, SIGMATT::Array{Float64,1})
    yy0, t0 = TSMS.yy0, TSMS.t0
    r0_m, t0_m = TSMS.r0_m, TSMS.t0_m
    
    S2TT, S3TT = STT[3], STT[4]
    
    YSTARR = yy0[3:t0]- ΦTT[1]*yy0[2:t0-1]-ΦTT[2]*yy0[1:t0-2]
    X1STARR = ones(t0-2,1)*(1-ΦTT[1]-ΦTT[2])
    X2STARR = S2TT[3:t0]- ΦTT[1]*S2TT[2:t0-1]-ΦTT[2]*S2TT[1:t0-2]
    X3STARR = S3TT[3:t0]- ΦTT[1]*S3TT[2:t0-1]-ΦTT[2]*S3TT[1:t0-2]
    
    YSTAR = YSTARR ./ sqrt.(SIGMATT[3:t0])
    XSTAR = [X1STARR./sqrt.(SIGMATT[3:t0]) X2STARR./sqrt.(SIGMATT[3:t0]) X3STARR./sqrt.(SIGMATT[3:t0])]
    
    V = inv(r0_m + XSTAR'XSTAR)
    MU =  V*(r0_m*t0_m + XSTAR'YSTAR)
    C = cholesky(Hermitian(V))
    
    MU_G = zeros(size(MU,1))
    
    accept = 0
    while accept == 0
        MU_G = MU + C.L * rand(Normal(), 3)    
        if MU_G[2] > MU_G[1] && MU_G[3] > MU_G[2]
            accept = 1
        end
    end

    MU2TT=MU_G[1]+MU_G[2]
    MU3TT=MU_G[1]+MU_G[3]    
    
    return [MU_G[1],MU2TT,MU3TT]
end


###
function gen_ss(TSMS, MUTT, STT, HTT::Array{Float64,1}, ΦTT::Array{Float64,1})
    t0, yy0 = TSMS.t0, TSMS.yy0
    v0_, d0_ = TSMS.v0_, TSMS.d0_ 
    
    YSTARR=yy0- MUTT
    YSTAR=YSTARR[3:t0]
    XSTAR=[YSTARR[2:t0-1] YSTARR[1:t0-2]]
    FORCST=XSTAR*[ΦTT[1]; ΦTT[2]]
    
    ERR=YSTAR-FORCST  
    
    S2TT, S3TT = STT[3], STT[4]
    H2TT, H3TT = HTT[1], HTT[2] 
    e_mat=ERR ./ sqrt.(1 .+ S2TT[3:t0]*H2TT) ./ sqrt.((1 .+ S3TT[3:t0]*H2TT) .* (1 .+ S3TT[3:t0]*H3TT));

    TSTAR=size(e_mat,1)
    nn = TSTAR + v0_
    d = d0_ + e_mat'e_mat
    
    c = rndc(nn) 
    t2=c/d
    σ2= 1/t2
    return σ2, FORCST+MUTT[3:t0]    
end 


###
function gen_h(TSMS, MUTT, ΣTT::Array{Float64,1}, STT, ΦTT::Array{Float64,1}, HTT::Array{Float64,1})
    yy0, t0 = TSMS.yy0, TSMS.t0
    v0_, d0_ = TSMS.v0_, TSMS.d0_
    
    #H2TT
    YSTARR = yy0- MUTT

    YSTAR = YSTARR[3:t0]
    XSTAR = [YSTARR[2:t0-1] YSTARR[1:t0-2]]

    FORCST=XSTAR*[ΦTT[1]; ΦTT[2]]

    ERR=YSTAR-FORCST   
    
    e_mat_= ERR ./ sqrt.(ΣTT[1]) ./ sqrt.(1 .+ STT[4][3:t0]*HTT[2])
    
    e_mat=Any[]
    for i=3:t0
        if STT[3][i]+STT[4][i] == 1
            push!(e_mat, e_mat_[i-2])
        end
    end
    
    tstar=size(e_mat,1)
    nn = tstar + v0_
    d = d0_ + e_mat'e_mat
    
    accept =0
    σ2=0
    while accept == 0
        c= rndc(nn)
        t2=c/d
        σ2=1/t2
        
        if σ2 > 0 
            accept = 1
        end
        
    end    
    
    HTT[1] = σ2-1
    
    #H3TT    
    YSTARR = yy0- MUTT #Ver si se pueden borrar, es posible que no sean necesarios
    YSTAR = YSTARR[3:t0] ##Ver si se pueden borrar, es posible que no sean necesarios
    XSTAR = [YSTARR[2:t0-1] YSTARR[1:t0-2]] ##Ver si se pueden borrar, es posible que no sean necesarios
    
    FORCST = XSTAR * [ΦTT[1]; ΦTT[2]] #Ver si se pueden borrar, es posible que no sean necesarios
    ERR = YSTAR-FORCST #Ver si se pueden borrar, es posible que no sean necesarios
    
    e_mat_ = ERR ./ sqrt.(ΣTT[1] .* (1 .+ STT[4][3:t0] .* HTT[1]))

    e_mat = Any[]
    for i=3:t0
        if STT[4][i] == 1
            push!(e_mat, e_mat_[i-2])
        end
    end

    tstar = size(e_mat,1) 
    nn = tstar + v0_  
    d = d0_ + e_mat'e_mat 
        
    accept =0
    σ2=0
    while accept == 0
        c=  rndc(nn) 
        t2=c/d
        σ2 = 1/t2
        
        if σ2 > 0
            accept = 1
        end
    end
    
    HTT[2] = σ2-1 
    
    return   HTT
    
end


###
function switchg(s,g::Array{Int64,1})
    n=size(s,1)
    m=size(g,1)
    switch = zeros(Float64,m,m)
    t=2
    
    while t <= n
        st1 = s[t-1]
        st = s[t]
        switch[Int(st1),Int(st)] = switch[Int(st1),Int(st)] +1
        t = t+1
    end
    return switch
end


###
function gen_PTT(TSMS, TRANMAT::Array{Float64,2})
    U1_1, U1_23, U1_2_23, U1_3_23 = TSMS.U1_1, TSMS.U1_23, TSMS.U1_2_23, TSMS.U1_3_23
    U2_2, U2_13, U2_1_13, U2_3_13 = TSMS.U2_2, TSMS.U2_13, TSMS.U2_1_13, TSMS.U2_3_13
    U3_3, U3_12, U3_2_12, U3_1_12 = TSMS.U3_3, TSMS.U3_12, TSMS.U3_2_12, TSMS.U3_1_12
    
    P1_1TT=rndb(TRANMAT[1,1]+U1_1, TRANMAT[1,2]+TRANMAT[1,3]+U1_23) #0.2
    P1_2TT=rndb(TRANMAT[1,2]+U1_2_23,TRANMAT[1,3]+U1_3_23)*(1-P1_1TT) #0.3
    P1_3TT=1-P1_1TT-P1_2TT        

    P2_2TT=rndb(TRANMAT[2,2]+U2_2,TRANMAT[2,1]+TRANMAT[2,3]+U2_13); #0.2
    P2_1TT=rndb(TRANMAT[2,1]+U2_1_13,TRANMAT[2,3]+U2_3_13)*(1-P2_2TT); #0.3
    P2_3TT=1-P2_1TT-P2_2TT;

    P3_3TT=rndb(TRANMAT[3,3]+U3_3, TRANMAT[3,2]+TRANMAT[3,1]+U3_12); #0.2
    P3_2TT=rndb(TRANMAT[3,2]+U3_2_12,TRANMAT[3,1]+U3_1_12)*(1-P3_3TT); #0.3
    P3_1TT=1-P3_2TT-P3_3TT
    
    P_TT = [P1_1TT, P1_2TT, P1_3TT, P2_1TT, P2_2TT, P2_3TT, P3_1TT, P3_2TT, P3_3TT]
    
    return P_TT
end


###########
function store_iterations!(i::Int64, ΣTT::Array{Float64,1}, μTT::Array{Float64,1}, P_TT::Array{Float64,1}, Φ_TT::Array{Float64,1}, ΣMM::Array{Float64,2}, μMM::Array{Float64,2}, PMM::Array{Float64,2}, ΦMM::Array{Float64,2})
    #STORAGE SPACES
    ΣMM[i, 1] = ΣTT[1]
    ΣMM[i, 2] = ΣTT[2]    
    ΣMM[i, 3] = ΣTT[3]
    
    μMM[i, 1] = μTT[1]
    μMM[i, 2] = μTT[2]    
    μMM[i, 3] = μTT[3]   
    
    PMM[i, 1] = P_TT[1]
    PMM[i, 2] = P_TT[2]
    PMM[i, 3] = P_TT[3]
    PMM[i, 4] = P_TT[4]
    PMM[i, 5] = P_TT[5]
    PMM[i, 6] = P_TT[6]
    PMM[i, 7] = P_TT[7]
    PMM[i, 8] = P_TT[8]
    PMM[i, 9] = P_TT[9]    
    
    ΦMM[i, 1] = Φ_TT[1]
    ΦMM[i, 2] = Φ_TT[2] 
end

###


###
function gibs_s3(TSMS)
    t0 = TSMS.t0
    P_TT, ΣTT, ΦTT, μTT, CAPN, HTT, N0, MM, MUTT  = TSMS.P_TT, TSMS.ΣTT, TSMS.ΦTT, TSMS.μTT, TSMS.CAPN, TSMS.HTT, TSMS.N0, TSMS.MM, TSMS.MUTT
    ΣMM, μMM, PMM, ΦMM, S1TTMM, S2TTMM, S3TTMM = TSMS.ΣMM, TSMS.μMM, TSMS.PMM, TSMS.ΦMM, TSMS.S1TTMM, TSMS.S2TTMM, TSMS.S3TTMM
    SIGMATT, sigmamm, MUMM, EX_AMM, EX_AMM2, stetmm = TSMS.SIGMATT, TSMS.sigmamm, TSMS.MUMM, TSMS.EX_AMM, TSMS.EX_AMM2, TSMS.stetmm
    yy = TSMS.yy
    i=1
    while i <= CAPN 
        #Start sampling
        
        STT = gen_st(TSMS, P_TT, ΣTT, ΦTT, μTT)

        SIGMATT = ΣTT[1]*STT[2] + ΣTT[2]*STT[3] + ΣTT[3]*STT[4] 
        
        ΦTT = gen_phi(TSMS, MUTT, SIGMATT)
        
        μTT = gen_mu(TSMS, ΦTT, STT, SIGMATT)
        MUTT=μTT[1]*STT[2] + μTT[2]*STT[3] + μTT[3]*STT[4]
        
        
        ΣTT[1], EX_ATT = gen_ss(TSMS, MUTT, STT, HTT, ΦTT)
        
        HTT = gen_h(TSMS, MUTT, ΣTT, STT, ΦTT, HTT)

        
        ΣTT[2] = ΣTT[1]*(1+HTT[1])
        ΣTT[3]=ΣTT[1]*(1+HTT[1])*(1+HTT[2])        
        
        TRANMAT = switchg(STT[1],[1,2,3])
     
        P_TT = gen_PTT(TSMS, TRANMAT)

#         #END OF ONE ITERATION
        store_iterations!(i, ΣTT, μTT, P_TT, ΦTT, ΣMM, μMM, PMM, ΦMM)
        
        if i>N0
            EX_AMM = EX_AMM + [0; 0; EX_ATT]
            EX_AMM2 = EX_AMM2 + [0; 0 ; EX_ATT] .* [0; 0; EX_ATT]
            
            S1TTMM = S1TTMM+STT[2]
            S2TTMM = S2TTMM+STT[3]
            S3TTMM = S3TTMM+STT[4]
            sigmamm=sigmamm+SIGMATT
            MUMM=MUMM+MUTT
#             print("ACA:",(size(yy), size(MUTT)))
#             print("ACA:", (size(stetmm), size((yy .- MUTT)), size(SIGMATT)))
#             stetmm=stetmm + (yy .- MUTT) ./ sqrt.(SIGMATT)
        end
        
        i=i+1
    end
    
    S1TTMM = S1TTMM/MM; S2TTMM = S2TTMM/MM; S3TTMM = S3TTMM/MM
    sigmamm=sigmamm/MM ; EX_AMM = EX_AMM/MM
    MUMM = MUMM/MM
#     stetmm=stetmm/MM
    
    xt_mn = EX_AMM
    xt_sd = sqrt.((EX_AMM2-MM*xt_mn.*xt_mn) ./ MM)
    low_b = xt_mn - 1.96*xt_sd
    up_b = xt_mn + 1.96*xt_sd 
  
    return ΣMM, μMM, PMM, ΦMM, sigmamm, S1TTMM, S2TTMM, S3TTMM#, stetmm

end













###############TABLES AND PLOTS #########################
function table_output(Prior, ΣMM, μMM, PMM, ΦMM; latex_format=true)
    table_out=Array{Any,2}(undef, 17, 6)
    if latex_format==true
        table_out[:,1]=[L"p_{11}", L"p_{12}", L"p_{13}",L"p_{21}", L"p_{22}", L"p_{23}", L"p_{31}", L"p_{32}", L"p_{33}", L"\phi_1", L"\phi_2" ,L"\sigma_1^2", L"\sigma_2^2", L"\sigma_3^2",L"\mu_1", L"\mu_2", L"\mu_3" ]
    else
        table_out[:,1]=["p11", "p12", "p13", "p21", "p22", "p23", "p31", "p32", "p33", "ϕ1", "ϕ2" ,"σ1", "σ2", "σ3","μ1", "μ2", "μ3" ]
    end
    
    #Domain
#     table_out[:,2]=["Beta", "Beta", "Beta", "Beta", "Beta", "Beta", "Beta", "Beta", "Beta", "ϕ1", "ϕ2" ,"σ1", "σ2", "σ3","μ1", "μ2", "μ3"]
    
    #Distribution
    #p_ii(p222)
    table_out[:,2]=["Beta", "Beta", "Beta", "Beta", "Beta", "Beta", "Beta", "Beta", "Beta", "Normal", "Normal" ,"InvGamma", "InvGamma", "InvGamma","Normal", "Normal", "Normal" ]

    @unpack ΣTT_prior, μTT_prior, ΦTT_prior, HTT_prior, P1_1, P1_2, P1_3, P2_1, P2_2, P2_3, P3_1, P3_2, P3_3 = Prior
    
    #Prior
    t1 = 3
    PP = [P1_1, P1_2, P1_3, P2_1, P2_2, P2_3, P3_1, P3_2, P3_3]
    table_out[1:9, t1] = PP
    table_out[10:11, t1] = ΦTT_prior
    table_out[12:14, t1] = ΣTT_prior
    table_out[15:17, t1] = μTT_prior
    
    #Mean
    t2 = 4
    table_out[1,t2]=mean(PMM[:,1]);table_out[2,t2]=mean(PMM[:,2]); table_out[3,t2]=mean(1 .- PMM[:,1] .- PMM[:,2]); 
    
    table_out[4,t2]=mean(PMM[:,4]); table_out[5,t2]=mean(PMM[:,5]); table_out[6,t2]=mean(1 .- PMM[:,4] .- PMM[:,5]); 
    
    table_out[7,t2]=mean(PMM[:,7]);table_out[8,t2]=mean(PMM[:,8]); table_out[9,t2]=mean(1 .- PMM[:,7] .- PMM[:,8]); 
    
    table_out[10,t2]=mean(ΦMM[:,1]);table_out[11,t2]=mean(ΦMM[:,2]);table_out[12,t2]=mean(ΣMM[:,1]);
    table_out[13,t2]=mean(ΣMM[:,2]);table_out[14,t2]=mean(ΣMM[:,3])
    table_out[15,t2]=mean(μMM[:,1]); table_out[16,t2]=mean(μMM[:,2]); table_out[17,t2]=mean(μMM[:,3]); 
    
#     #Mode
#     t2 = 5
#     table_out[1,t2]=mode(PMM[:,1]);table_out[2,t2]=mode(PMM[:,2]); table_out[3,t2]=mode(1 .- PMM[:,1] .- PMM[:,2]); 
    
#     table_out[4,t2]=mode(PMM[:,4]); table_out[5,t2]=mode(PMM[:,5]); table_out[6,t2]=mode(1 .- PMM[:,4] .- PMM[:,5]); 
    
#     table_out[7,t2]=mode(PMM[:,7]);table_out[8,t2]=mode(PMM[:,8]); table_out[9,t2]=mode(1 .- PMM[:,7] .- PMM[:,8]); 
    
#     table_out[10,t2]=mode(ΦMM[:,1]);table_out[11,t2]=mode(ΦMM[:,2]);table_out[12,t2]=mode(ΣMM[:,1]);
#     table_out[13,t2]=mode(ΣMM[:,2]);table_out[14,t2]=mode(ΣMM[:,3])
#     table_out[15,t2]=mode(μMM[:,1]); table_out[16,t2]=mode(μMM[:,2]); table_out[17,t2]=mode(μMM[:,3]);     
    
#     #median
#     t2 = 6
#     table_out[1,t2]=median(PMM[:,1]);table_out[2,t2]=median(PMM[:,2]); table_out[3,t2]=median(1 .- PMM[:,1] .- PMM[:,2]); 
    
#     table_out[4,t2]=median(PMM[:,4]); table_out[5,t2]=median(PMM[:,5]); table_out[6,t2]=median(1 .- PMM[:,4] .- PMM[:,5]); 
    
#     table_out[7,t2]=median(PMM[:,7]);table_out[8,t2]=median(PMM[:,8]); table_out[9,t2]=median(1 .- PMM[:,7] .- PMM[:,8]); 
    
#     table_out[10,t2]=median(ΦMM[:,1]);table_out[11,t2]=median(ΦMM[:,2]);table_out[12,t2]=median(ΣMM[:,1]);
#     table_out[13,t2]=median(ΣMM[:,2]);table_out[14,t2]=median(ΣMM[:,3])
#     table_out[15,t2]=median(μMM[:,1]); table_out[16,t2]=median(μMM[:,2]); table_out[17,t2]=median(μMM[:,3]);     
        

#     #Standar Dev
#     ts = 7   
#     table_out[1,ts]=std(PMM[:,1]);table_out[2,ts]=std(PMM[:,2]); table_out[3,ts]=std(1 .- PMM[:,1] .- PMM[:,2]); 
#     table_out[4,ts]=std(PMM[:,4]);table_out[5,ts]=std(PMM[:,5]); table_out[6,ts]=std(1 .- PMM[:,4] .- PMM[:,5]); 
#     table_out[7,ts]=std(PMM[:,7]);table_out[8,ts]=std(PMM[:,8]); table_out[9,ts]=std(1 .- PMM[:,7] .- PMM[:,8]); 
#     table_out[10,ts]=std(ΦMM[:,1]);table_out[11,ts]=std(ΦMM[:,2]);table_out[12,ts]=std(ΣMM[:,1]);table_out[13,ts]=std(ΣMM[:,2]);table_out[14,ts]=std(ΣMM[:,3])
#     table_out[15,ts]=std(μMM[:,1]); table_out[16,ts]=std(μMM[:,2]); table_out[17,ts]=std(μMM[:,3]);
        
    #Median    
#     table_out[1,5]=median(PMM[:,1]);table_out[2,5]=median(PMM[:,2]); table_out[3,5]=median(1 .- PMM[:,1] .- PMM[:,2]); 
#     table_out[4,5]=median(PMM[:,4]);table_out[5,5]=median(PMM[:,5]); table_out[6,5]=median(1 .- PMM[:,4] .- PMM[:,5]); 
#     table_out[7,5]=median(PMM[:,7]);table_out[8,5]=median(PMM[:,8]); table_out[9,5]=median(1 .- PMM[:,7] .- PMM[:,8]); 
#     table_out[10,5]=median(ΦMM[:,1]);table_out[11,5]=median(ΦMM[:,2]);table_out[12,5]=median(ΣMM[:,1]);
#     table_out[13,5]=median(ΣMM[:,2]);table_out[14,5]=median(ΣMM[:,3])
     #table_out[15,5]=median(μMM[:,1]); #table_out[16,5]=median(μMM[:,2]);      #table_out[17,5]=median(μMM[:,3]); 

    
#     #95%l 95%h 
    t3 = 5; t4 = 6
    b_p11=sort(PMM[:,1]); b_p12=sort(PMM[:,2]); b_p13=sort(1 .- PMM[:,1] .- PMM[:,2]); 
    b_p21=sort(PMM[:,4]); b_p22=sort(PMM[:,5]); b_p23=sort(1 .- PMM[:,4] .- PMM[:,5]); 
    b_p31=sort(PMM[:,7]); b_p32=sort(PMM[:,8]); b_p33=sort(1 .- PMM[:,7] .- PMM[:,8]); 
    b_s1=sort(ΣMM[:,1]); b_s2=sort(ΣMM[:,2]); b_s3=sort(ΣMM[:,3]); b_ϕ1=sort(ΦMM[:,1]); b_ϕ2=sort(ΦMM[:,2]); 
    b_mu1=sort(μMM[:,1]); b_mu2=sort(μMM[:,2]); b_mu3=sort(μMM[:,3]); 

    table_out[1,t3]=b_p11[Int(size(b_p11,1)*0.05)]; table_out[1,t4]=b_p11[Int(size(b_p11,1)*0.95)]; 
    table_out[2,t3]=b_p12[Int(size(b_p12,1)*0.05)]; table_out[2,t4]=b_p12[Int(size(b_p12,1)*0.95)]; 
    table_out[3,t3]=b_p13[Int(size(b_p13,1)*0.05)]; table_out[3,t4]=b_p13[Int(size(b_p13,1)*0.95)]; 
    
    table_out[4,t3]=b_p21[Int(size(b_p21,1)*0.05)]; table_out[4,t4]=b_p21[Int(size(b_p21,1)*0.95)]; 
    table_out[5,t3]=b_p22[Int(size(b_p22,1)*0.05)]; table_out[5,t4]=b_p22[Int(size(b_p22,1)*0.95)]; 
     table_out[6,t3]=b_p23[Int(size(b_p23,1)*0.05)]; table_out[6,t4]=b_p23[Int(size(b_p23,1)*0.95)]; 
    
    table_out[7,t3]=b_p31[Int(size(b_p31,1)*0.05)]; table_out[7,t4]=b_p31[Int(size(b_p31,1)*0.95)]; 
    table_out[8,t3]=b_p32[Int(size(b_p32,1)*0.05)]; table_out[8,t4]=b_p32[Int(size(b_p32,1)*0.95)]; 
    table_out[9,t3]=b_p33[Int(size(b_p33,1)*0.05)]; table_out[9,t4]=b_p33[Int(size(b_p33,1)*0.95)]; 
    
    table_out[10,t3]=b_ϕ1[Int(size(b_ϕ1,1)*0.05)]; table_out[10,t4]=b_ϕ1[Int(size(b_ϕ1,1)*0.95)]; 
    table_out[11,t3]=b_ϕ2[Int(size(b_ϕ2,1)*0.05)]; table_out[11,t4]=b_ϕ2[Int(size(b_ϕ2,1)*0.95)]; 
    table_out[12,t3]=b_s1[Int(size(b_s1,1)*0.05)]; table_out[12,t4]=b_s1[Int(size(b_s1,1)*0.95)]; 
    table_out[13,t3]=b_s2[Int(size(b_s2,1)*0.05)]; table_out[13,t4]=b_s2[Int(size(b_s2,1)*0.95)]; 
    table_out[14,t3]=b_s3[Int(size(b_s3,1)*0.05)]; table_out[14,t4]=b_s3[Int(size(b_s3,1)*0.95)]; 
    table_out[15,t3]=b_mu1[Int(size(b_mu1,1)*0.05)]; table_out[15,t4]=b_mu1[Int(size(b_mu1,1)*0.95)]; 
    table_out[16,t3]=b_mu2[Int(size(b_mu2,1)*0.05)]; table_out[16,t4]=b_mu2[Int(size(b_mu2,1)*0.95)]; 
    table_out[17,t3]=b_mu3[Int(size(b_mu3,1)*0.05)]; table_out[17,t4]=b_mu3[Int(size(b_mu3,1)*0.95)]; 
    
    table_out[:,3:end] = round.(table_out[:,3:end], digits=3)
    return table_out
#     return pretty_table(table_out, ["", "Mean", "SD", "MD", "95%l", "95%h"];formatters = ft_printf("%5.4f"))    
end

function table_latex(table_out, path; OIL=true)
    if OIL == true
        table_name = "table_OIL.tex"
    else
        table_name = "table_CUP.tex"
    end
    
    latex_tabular(path*table_name,
                  Tabular("lccccc"),
                  [Rule(:top),
                   ["", "", "", MultiColumn(4, :l, "{Posterior")],
                   Rule(:mid),
                   ["Parameter", "Density","Prior", "Mean",MultiColumn(2, :c, L"95\% \text{posterior bands}")],
                   Rule(:mid),
                   table_out,
                   Rule(:bottom)])

end

###
function plot_output_oil(S1TTMM, S2TTMM, S3TTMM, sigmamm, path, dates)
    
#     dates = Date(1948,2,1):Month(1):Date(2021,3,1); dates = collect(dates)
    tick_years = Date.(range(1950, 2020, length=15))
    DateTick = Dates.format.(tick_years, "yyyy")    
    
    linealpha = 0.5; linecolor=:blue; gridalpha=0.8
    a=plot(dates, S1TTMM, label="", title="Probability of a low-variance", linealpha = linealpha, linecolor=linecolor, foreground_color_grid=:grey, gridalpha=gridalpha)
    b=plot(dates, S2TTMM, label="", title="Probability of a medium-variance", linealpha = linealpha, linecolor=linecolor, foreground_color_grid=:grey, gridalpha=gridalpha)
    c=plot(dates, S3TTMM, label="", title="Probability of a high-variance", linealpha = linealpha, linecolor=linecolor, foreground_color_grid=:grey, gridalpha=gridalpha)

    plot1 = plot(a,b,c, layout=(3,1))    
    plot1 = plot!(xticks=(tick_years,DateTick), xtickfontsize=6, ytickfontsize=6, titlefontsize=8)
    
    plot!(size=(700,600))
    savefig(path*"OIL_prob.pdf")
    
    plot2 = plot(dates, sigmamm, label="", title="Estimated variance of historical returns", linealpha = linealpha, linecolor=linecolor)    
    plot2 = plot!(xticks=(tick_years,DateTick), xtickfontsize=7, ytickfontsize=7, titlefontsize=9)
    plot!(size=(500,300))
    savefig(path*"OIL_var.pdf")

    return (plot1=plot1, plot2=plot2)
end



function plot_output_cu(S1TTMM, S2TTMM, S3TTMM, sigmamm, path, dates)
    
#     dates = Date(1960,2,1):Month(1):Date(2021,3,1); dates = collect(dates)
    tick_years = Date.(range(1960, 2020, length=13))
    DateTick = Dates.format.(tick_years, "yyyy")
    
    linealpha = 0.5; linecolor=:blue; gridalpha=0.8
    a=plot(dates, S1TTMM, label="", title="Probability of a low-variance", linealpha = linealpha, linecolor=linecolor, foreground_color_grid=:grey, gridalpha=gridalpha)
    b=plot(dates, S2TTMM, label="", title="Probability of a medium-variance", linealpha = linealpha, linecolor=linecolor, foreground_color_grid=:grey, gridalpha=gridalpha)
    c=plot(dates, S3TTMM, label="", title="Probability of a high-variance", linealpha = linealpha, linecolor=linecolor, foreground_color_grid=:grey, gridalpha=gridalpha)

    plot1 = plot(a,b,c, layout=(3,1))
    plot1 = plot!(xticks=(tick_years,DateTick), xtickfontsize=7, ytickfontsize=7, titlefontsize=7)
    plot!(size=(700,600))
    savefig(path*"CUP_prob.pdf")    
    
    plot2 = plot(dates, sigmamm, label="", title="Estimated variance of historical returns", linealpha = linealpha, linecolor=linecolor)   
    plot2 = plot!(xticks=(tick_years,DateTick), xtickfontsize=6, ytickfontsize=6, titlefontsize=8)
    plot!(size=(500,300))
    savefig(path*"CUP_var.pdf")    
    return (plot1=plot1, plot2=plot2)
end



function plot_output_oil_cu(S1_O, S2_O, S3_O, S1_C, S2_C, S3_C , path, dates)
    
#     dates = Date(1960,2,1):Month(1):Date(2021,3,1); dates = collect(dates)
    tick_years = Date.(range(1959, 2023, length=17))
    DateTick = Dates.format.(tick_years, "yyyy")    
    
    linealpha = 0.5; linecolor=:blue; gridalpha=0.4
    dates_aux = length(dates)-1
    
    
    a=plot(dates, S1_O[end-dates_aux:end], label="Oil", title="Probability of a low-variance", linealpha = linealpha, linecolor=linecolor, foreground_color_grid=:grey, gridalpha=gridalpha)
    a=plot!(dates, S1_C, label="Copper", line = (:steppre, :dot, 1, 1, :red))      
    
    b=plot(dates, S2_O[end-dates_aux:end], label="", title="Probability of a medium-variance", linealpha = linealpha, linecolor=linecolor, foreground_color_grid=:grey, gridalpha=gridalpha)   
    b=plot!(dates, S2_C, label="", line = (:steppre, :dot, 1, 1, :red))   
    
    c=plot(dates, S3_O[end-dates_aux:end], label="", title="Probability of a high-variance", linealpha = linealpha, linecolor=linecolor, foreground_color_grid=:grey, gridalpha=gridalpha)    
    c=plot!(dates, S3_C, label="", line = (:steppre, :dot, 1, 1, :red))  
    
    plot1 = plot(a,b,c, layout=(3,1))    
    plot1 = plot!(xticks=(tick_years,DateTick), xtickfontsize=7, ytickfontsize=7, titlefontsize=9, foreground_color_legend = nothing, legend=(0.9, 0.90))
    
    plot!(size=(700,600))
    savefig(path*"OIL_CUP_prob.pdf")
    
    return plot1
end



function plot_prices(POIL, PCU, path, dates; stand=true)
    
#     dates = Date(1960,1,1):Month(1):Date(2021,3,1); dates = collect(dates)
    tick_years = Date.(range(1959, 2023, length=17))
    DateTick = Dates.format.(tick_years, "yyyy")   
    linealpha = 0.5; linecolor=:blue; gridalpha=0.6  
    dates_aux = length(dates)-1
    
    if stand==true
#         stetmm_POIL = (POIL.POIL) ./ std(POIL.POIL)
#         stetmm_CUP = (PCU.PCU) ./ std(PCU.PCU) 
        stetmm_POIL = (POIL.POIL .- mean(POIL.POIL)) ./ std(POIL.POIL)
        stetmm_CUP = (PCU.PCU .- mean(PCU.PCU)) ./ std(PCU.PCU) 
        
        p1=plot(dates, stetmm_POIL[end-dates_aux:end], label="Oil", linealpha = linealpha, linecolor=linecolor, foreground_color_grid=:grey, gridalpha=gridalpha)
        p1=plot!(dates, stetmm_CUP, label="Copper", line = (:steppre, :dot, 1, 1, :red))
        p1=plot!(xticks=(tick_years,DateTick), xtickfontsize=7, ytickfontsize=7, titlefontsize=9)

    else

        p1=plot(dates, log.(POIL.POIL[end-dates_aux:end]), label="Oil", linealpha = linealpha, linecolor=linecolor, foreground_color_grid=:grey, gridalpha=gridalpha)        
        p1=plot!(dates, log.(PCU.PCU), label="Copper", line = (:steppre, :dot, 1, 1, :red))
        p1=plot!(xticks=(tick_years,DateTick), xtickfontsize=7, ytickfontsize=7, titlefontsize=9)        
    
    end
    plot!(legend=(0.9, 0.95), size=(800,300), foreground_color_legend = nothing)
    savefig(path*"OIL_CUP_prices.pdf")    
    
    return p1 
end



##########################################################3
function gibs_s3(TSMS)
    t0 = TSMS.t0
    P_TT, ΣTT, ΦTT, μTT, CAPN, HTT, N0, MM, MUTT  = TSMS.P_TT, TSMS.ΣTT, TSMS.ΦTT, TSMS.μTT, TSMS.CAPN, TSMS.HTT, TSMS.N0, TSMS.MM, TSMS.MUTT
    ΣMM, μMM, PMM, ΦMM, S1TTMM, S2TTMM, S3TTMM = TSMS.ΣMM, TSMS.μMM, TSMS.PMM, TSMS.ΦMM, TSMS.S1TTMM, TSMS.S2TTMM, TSMS.S3TTMM
    SIGMATT, sigmamm, MUMM, EX_AMM, EX_AMM2, stetmm = TSMS.SIGMATT, TSMS.sigmamm, TSMS.MUMM, TSMS.EX_AMM, TSMS.EX_AMM2, TSMS.stetmm
    yy = TSMS.yy
    i=1
    while i <= CAPN 
        #Start sampling
        
        STT = gen_st(TSMS, P_TT, ΣTT, ΦTT, μTT)

        SIGMATT = ΣTT[1]*STT[2] + ΣTT[2]*STT[3] + ΣTT[3]*STT[4] 
        
        ΦTT = gen_phi(TSMS, MUTT, SIGMATT)
        
        μTT = gen_mu(TSMS, ΦTT, STT, SIGMATT)
        MUTT=μTT[1]*STT[2] + μTT[2]*STT[3] + μTT[3]*STT[4]
        
        
        ΣTT[1], EX_ATT = gen_ss(TSMS, MUTT, STT, HTT, ΦTT)
        
        HTT = gen_h(TSMS, MUTT, ΣTT, STT, ΦTT, HTT)

        
        ΣTT[2] = ΣTT[1]*(1+HTT[1])
        ΣTT[3]=ΣTT[1]*(1+HTT[1])*(1+HTT[2])        
        
        TRANMAT = switchg(STT[1],[1,2,3])
     
        P_TT = gen_PTT(TSMS, TRANMAT)

#         #END OF ONE ITERATION
        store_iterations!(i, ΣTT, μTT, P_TT, ΦTT, ΣMM, μMM, PMM, ΦMM)
        
        if i>N0
            EX_AMM = EX_AMM + [0; 0; EX_ATT]
            EX_AMM2 = EX_AMM2 + [0; 0 ; EX_ATT] .* [0; 0; EX_ATT]
            
            S1TTMM = S1TTMM+STT[2]
            S2TTMM = S2TTMM+STT[3]
            S3TTMM = S3TTMM+STT[4]
            sigmamm=sigmamm+SIGMATT
            MUMM=MUMM+MUTT
#             print("ACA:",(size(yy), size(MUTT)))
#             print("ACA:", (size(stetmm), size((yy .- MUTT)), size(SIGMATT)))
#             stetmm=stetmm + (yy .- MUTT) ./ sqrt.(SIGMATT)
        end
        
        i=i+1
    end
    
    S1TTMM = S1TTMM/MM; S2TTMM = S2TTMM/MM; S3TTMM = S3TTMM/MM
    sigmamm=sigmamm/MM ; EX_AMM = EX_AMM/MM
    MUMM = MUMM/MM
#     stetmm=stetmm/MM
    
    xt_mn = EX_AMM
    xt_sd = sqrt.((EX_AMM2-MM*xt_mn.*xt_mn) ./ MM)
    low_b = xt_mn - 1.96*xt_sd
    up_b = xt_mn + 1.96*xt_sd 
  
    return ΣMM, μMM, PMM, ΦMM, sigmamm, S1TTMM, S2TTMM, S3TTMM, xt_mn, xt_sd

end



function res_sum(EX_AMM, yy0)
    return sum(abs.(EX_AMM .- yy0))
    
end


function time_function(f, Σ_values, μ_values; N0=500, MM=1000)
    time_in = now()
    TSMS = gen_TSMS(POIL, N0=N0, MM=MM, ΣTT_prior=ΣTT, μTT_prior=μTT, stand=false);
    ΣMM, μMM, PMM, ΦMM, sigmamm, S1TTMM, S2TTMM, S3TTMM, xt_mn, xt_sd = f(TSMS);
end


function gibs_sampling_with_prior(Data, Prior; N0=100, MM0=500, N1=100, MM1=500, stand=false, prueba=true)
    
    if prueba==true
        #Prior original
        TSMS = gen_TSMS(Data, Prior, N0=N0, MM=MM0, stand=stand);
        ΣMM, μMM, PMM, ΦMM, sigmamm, S1TTMM, S2TTMM, S3TTMM, xt_mn, xt_sd = gibs_s3(TSMS);
        print("First step ok \n")   
        #TABLA
        table_out=Array{Any,2}(undef, 17, 2)
        table_out[:,1]=["p11", "p12", "p13", "p21", "p22", "p23", "p31", "p32", "p33", "ϕ1", "ϕ2" ,"σ1", "σ2", "σ3","μ1", "μ2", "μ3" ]
        table_out[1,2]=mean(PMM[:,1]);table_out[2,2]=mean(PMM[:,2]); table_out[3,2]=mean(1 .- PMM[:,1] .- PMM[:,2]); 
        table_out[4,2]=mean(PMM[:,4]); table_out[5,2]=mean(PMM[:,5]); table_out[6,2]=mean(1 .- PMM[:,4] .- PMM[:,5]); 
        table_out[7,2]=mean(PMM[:,7]);table_out[8,2]=mean(PMM[:,8]); table_out[9,2]=mean(1 .- PMM[:,7] .- PMM[:,8]); 
        table_out[10,2]=mean(ΦMM[:,1]);table_out[11,2]=mean(ΦMM[:,2]);table_out[12,2]=mean(ΣMM[:,1]);
        table_out[13,2]=mean(ΣMM[:,2]);table_out[14,2]=mean(ΣMM[:,3])    
        table_out[15,2]=mean(μMM[:,1]); table_out[16,2]=mean(μMM[:,2]); table_out[17,2]=mean(μMM[:,3]); 
        
        return table_out, ΣMM, μMM, PMM, ΦMM, sigmamm, S1TTMM, S2TTMM, S3TTMM, xt_mn, xt_sd, Prior
        
        
    else
        #Prior original
        TSMS = gen_TSMS(Data, Prior, N0=N0, MM=MM0, stand=stand);
        ΣMM, μMM, PMM, ΦMM, sigmamm, S1TTMM, S2TTMM, S3TTMM, xt_mn, xt_sd = gibs_s3(TSMS);
        print("First step ok \n")

        #Calcula nuevos prior
        p_gen = @with_kw (ΣTT_prior = [mean(ΣMM[:,1]), mean(ΣMM[:,2]), mean(ΣMM[:,3])], 
                          μTT_prior = [mean(μMM[:,1]), mean(μMM[:,2]), mean(μMM[:,3])], 
                          ΦTT_prior = [mean(ΦMM[:,1]), mean(ΦMM[:,2])],
                          HTT_prior = [0.0, 0.0],
                          P1_1=mean(PMM[:,1]), P1_2=mean(PMM[:,2]), P1_3=1-P1_1-P1_2, 
                          P2_1=mean(PMM[:,4]), P2_2=mean(PMM[:,5]), P2_3=1-P2_1-P2_2,    
                          P3_1=mean(PMM[:,7]), P3_2=mean(PMM[:,8]), P3_3=1-P3_1-P3_2);
        Prior2 = p_gen()
        TSMS = gen_TSMS(Data, Prior2, N0=N1, MM=MM1, stand=stand);
        #Nuevo cálculo
        ΣMM, μMM, PMM, ΦMM, sigmamm, S1TTMM, S2TTMM, S3TTMM, xt_mn, xt_sd = gibs_s3(TSMS);
        print("Second step ok")
        return ΣMM, μMM, PMM, ΦMM, sigmamm, S1TTMM, S2TTMM, S3TTMM, xt_mn, xt_sd, Prior2
    end
    
end

function gibs_sampling(Data, Prior; N0=100, MM0=500, stand=false)
    
    #Prior original
    TSMS = gen_TSMS(Data, Prior, N0=N0, MM=MM0, stand=stand);
    ΣMM, μMM, PMM, ΦMM, sigmamm, S1TTMM, S2TTMM, S3TTMM, xt_mn, xt_sd = gibs_s3(TSMS);
    return ΣMM, μMM, PMM, ΦMM, sigmamm, S1TTMM, S2TTMM, S3TTMM, xt_mn, xt_sd
end
