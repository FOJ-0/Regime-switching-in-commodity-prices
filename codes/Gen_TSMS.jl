mutable struct gen_TSMS{}
    #Data
    yy::Array{Float64,1}
    x_mat::Array{Float64,2}
    T::Int64
    t0::Int64
    NO_ST::Int64
    DMNSION::Int64
    yy0::Array{Float64,1}
    
    #Prior
    #FOR TRAN. PROB.    
    U1_1::Float64
    U1_23::Float64
    U1_2_23::Float64
    U1_3_23::Float64
    U2_13::Float64
    U2_2::Float64
    U2_1_13::Float64
    U2_3_13::Float64
    U3_3::Float64
    U3_12::Float64
    U3_2_12::Float64
    U3_1_12::Float64
    
    #FOR SIGMA'S
    v0_::Int64
    d0_::Int64
    R0_::Array{Float64,2}
    
    #FOR PHI1,PHI2
    t0_::Array{Int64,1}
    
    #For MU1, tmp2, tmp3
    r0_m::Array{Float64,2}
    t0_m::Array{Float64,1}
    
    #INITIAL VALUES FOR THE GIBBS SAMPLER
    P1_1TT::Float64
    P1_2TT::Float64
    P1_3TT::Float64
    P2_1TT::Float64
    P2_2TT::Float64
    P2_3TT::Float64
    P3_1TT::Float64
    P3_2TT::Float64
    P3_3TT::Float64
    P_TT::Array{Float64,1}
    ΦTT::Array{Float64,1}
    MU1TT::Float64
    tmp2tt::Float64
    tmp3tt::Float64
    MU2TT::Float64
    MU3TT::Float64
    MUTT::Array{Float64,2}
    μTT::Array{Float64,1}
    ΣTT::Array{Float64,1}
    HTT::Array{Float64,1}
    
    #DESIGN  FOR THE GIBBS SAMPLER    
    N0::Int64
    MM::Int64
    LL::Int64
    CAPN::Int64
    
    #Store space
    ΣMM::Array{Float64,2}
    μMM::Array{Float64,2}
    PMM::Array{Float64,2}
    ΦMM::Array{Float64,2}
    S1TTMM::Array{Float64,2}
    S2TTMM::Array{Float64,2}
    S3TTMM::Array{Float64,2}
    SIGMATT::Array{Float64,1}
    sigmamm::Array{Float64,2}
    MUMM::Array{Float64,2}
    EX_AMM::Array{Float64,2}
    EX_AMM2::Array{Float64,2}
    stetmm::Array{Float64,2}
    
    
    
    
    
    function gen_TSMS(data, Prior; N0=1000, MM=5000, LL=5, stand=false ) #data
        
        #1. Gen data
    #    data=CSV.read("C:\\Users\\ordon\\Dropbox\\Markov_Python_2020\\Kim-Nelson\\Ch9\\Data\\SUMMARY\\POIL_PCU.csv", header=true, decimal='.') #@Quarter, U.S. real GNP@ 1947.1 -- 1995.3
#          data=CSV.read("/home/felix/Dropbox/Markov_Python_2020/Kim-Nelson/Ch9/Data/SUMMARY/POIL_PCU.csv", DataFrame, header=true)

        if stand==false
            yy0 = log.(data[2:end,2] ./ data[1:end-1,2])*100   
#             yy0 = (data[2:end,2] ./ data[1:end-1,2] .-1)*100   
        else
            data[:,2] = (data[:,2]) ./std(data[:,2]) #.- mean(data[:,2])
            yy0 = log.(data[2:end,2] ./ data[1:end-1,2])*100   
#             yy0 = (data[2:end,2] ./ data[1:end-1,2] .-1)*100   
        end
        
        @unpack ΣTT_prior, μTT_prior, ΦTT_prior, HTT_prior, P1_1, P1_2, P1_3, P2_1, P2_2, P2_3, P3_1, P3_2, P3_3 = Prior
            
        t0=size(yy0,1)

        yy_d=yy0[2:t0]-yy0[1:t0-1]

        lag_ar=2
        NO_ST=lag_ar+1 #@ NUMBER OF STATES TO BE CONSIDERED@
        DMNSION=3^NO_ST

        yy = yy0[lag_ar+1:t0,1]

        x_mat = zeros(t0-lag_ar,2)
        x_mat[:,1] = yy0[lag_ar:t0-1,1] 
        x_mat[:,2] = yy0[lag_ar-1:t0-2,1]

        T=size(yy,1)
        
        
        #2. Gen prior
        #FOR TRAN. PROB.
        U1_1 = .1; U1_23 = .1;  U1_2_23 = .1;   U1_3_23 = .1;
        U2_13 = .1; U2_2 = .1;  U2_1_13 = .1;   U2_3_13 = .1; 
        U3_3 = .1; U3_12 = .1;  U3_2_12 = .1;   U3_1_12 = .1;

        #FOR SIGMA'S
        v0_=0; d0_=0;

        R0_ = Matrix{Float64}(I, 2, 2) / 25

        #FOR PHI1,PHI2
        t0_ = [0;0]

        #For MU1, tmp2, tmp3
        r0_m = Matrix{Float64}(I, 3, 3) / 25
        t0_m = [0; 0.1; 0.2]

        #INITIAL VALUES FOR THE GIBBS SAMPLER
        
#         P1_1TT=0.9  ;  P1_2TT=0.05 ;  P1_3TT=1-P1_1TT-P1_2TT; 
#         P2_1TT=0.05 ; P2_2TT=0.9; P2_3TT=1-P2_1TT-P2_2TT;
#         P3_1TT=0.1  ;  P3_2TT=0.1; P3_3TT=1-P3_1TT-P3_2TT;
        
        P_TT = [P1_1, P1_2, P1_3, P2_1, P2_2, P2_3, P3_1, P3_2, P3_3]
        ΦTT = similar(ΦTT_prior)
        ΦTT .= ΦTT_prior #[0.0, 0.0]

#         MU1TT = 0.0; tmp2tt = 0.0; tmp3tt = 0.0; #    MU1TT = -1; tmp2tt = 1.5; tmp3tt = 3;
#         MU2TT = MU1TT+tmp2tt; MU3TT = MU1TT+tmp3tt;


#         μTT = [MU1TT, MU2TT, MU3TT]
        μTT = similar(μTT_prior)
        μTT .= μTT_prior
        MUTT= zeros(Float64, t0,1) #MU1TT*S1TT+MU2TT*S2TT+MU3TT*S3TT
        
        ΣTT  = similar(ΣTT_prior)
        ΣTT .= ΣTT_prior #[0.5, 50.0, 300.0]
        #ΣTT = [0.5, 40, 200]
        #ΣTT = [1.0, 1.0, 1.0] 
        HTT = similar(HTT_prior)
        HTT .= HTT_prior #[0.0, 0.0]

        #DESIGN  FOR THE GIBBS SAMPLER
        CAPN = N0 + MM;  #TOTAL NUMBER OF DRAWS   
        
        #3. Gen store space
        ΣMM = zeros(Float64, CAPN, 3)
        μMM = zeros(Float64, CAPN, 3)
        PMM = zeros(Float64, CAPN, 9)
        ΦMM = zeros(Float64, CAPN, 2)

        S1TTMM=zeros(Float64, t0,1)
        S2TTMM=zeros(Float64, t0,1)
        S3TTMM=zeros(Float64, t0,1)

        SIGMATT=zeros(Float64, t0,1)
        sigmamm=zeros(Float64, t0,1)
        stetmm=zeros(Float64, t0,1)

        MUMM = zeros(Float64, t0,1)
        MUTT = zeros(Float64, t0,1)
        EX_AMM = zeros(Float64, t0,1)
        EX_AMM2 = zeros(Float64, t0,1)           
        
        

        (yy=yy, x_mat=x_mat, T=T, t0=t0, NO_ST=NO_ST, DMNSION=DMNSION, yy0=yy0, 
         U1_1=U1_1, U1_23=U1_23, U1_2_23=U1_2_23, U1_3_23=U1_3_23, U2_13=U2_13, U2_2=U2_2, U2_1_13=U2_1_13, U2_3_13=U2_3_13, 
         U3_3=U3_3, U3_12=U3_12, U3_2_12=U3_2_12, U3_1_12=U3_1_12, v0_=v0_ , d0_=d0_, t0_=t0_, R0_=R0_, r0_m=r0_m, t0_m=t0_m, 
         P_TT=P_TT, ΦTT=ΦTT, μTT=μTT, ΣTT=ΣTT, HTT=HTT, N0=N0, MM=MM, LL=LL, CAPN=CAPN, MUTT=MUTT,
         ΣMM=ΣMM, μMM=μMM, PMM=PMM, ΦMM=ΦMM, S1TTMM=S1TTMM, S2TTMM=S2TTMM, S3TTMM=S3TTMM, SIGMATT=SIGMATT, 
         sigmamm=sigmamm, MUMM=MUMM, EX_AMM=EX_AMM, EX_AMM2=EX_AMM, stetmm=stetmm)
#         gen_TSMS(yy, x_mat, T, t0, NO_ST, DMNSION, yy0)        
    end    
end

