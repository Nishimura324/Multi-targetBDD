import scipy.special as sc

f_inf = float('inf')
delta0 = lambda beta: RR(beta * (pi*beta)**(1/beta) / RR(2*pi*e)) ** (1 / (2*(beta-1))) #root Hermite factor of BKZ-blocksize beta
log2 = lambda x: log(x,2)
Vl = lambda l,R: RR(RR(pi**(l/2)) * R**l / RR(gamma(1+l/2))) #Volume of R radius l-dimension ball

def chernoff_bounds(p1,p2): #[BBTV23, H.1]
    alpha = log(p2/(1-p1)) / log(p2*p1/(1-p2)/(1-p1))
    return ceil(-log(2) / (alpha*log(p1/alpha) + (1-alpha)*log((1-p1)/(1-alpha))))

def TOTALTIME(n,k,m,q,sigma,beta_range,itr_range,P,qc,mlwe,qset='s'):
    d = n*m #lattice dimension
    vol = RR(q**(n*m-n*k)) #lattice volume

    RESULT=[[0,0,f_inf]]

    for itr in itr_range:

        BKZ_data=[f_inf for beta in range(d)]
        ENUM_data=[f_inf for beta in range(d)]   #only classical
        DETECT_data=[f_inf for beta in range(d)] #only quantum
        FIND_data=[f_inf for beta in range(d)]   #only quantum
        TOTAL_data=[f_inf for beta in range(d)]

        #Unique setting (For BDD) [ANSS18]
        b_l = [0 for l in range(d+1)]
        lowerVolC_l = [0 for l in range(d+1)]  
        N_beta_l = [[0 for l in range(d+1)] for beta in range(d)] #number of nodes at depth l
        N_beta = [0 for beta in range(d+1)] #total number of nodes

        if mlwe == 'm':
            p = P/(n*itr)
        else:
            p = P/itr

        for l in range(1,d+1):
            b_l[l] = sc.gammaincinv(l/2,p) * 2 * sigma**2
            lowerVolC_l[l] = Vl(l,RR(sqrt(b_l[l])))

        for beta in beta_range:
            for l in range(1,d+1):
                N_beta_l[beta][l] = lowerVolC_l[l] * delta0(beta) ** (l*(d-l)) / vol**(l/d)
            N_beta[beta] = RR(1+sum(N_beta_l[beta]))

            T_BKZ = lambda beta,d,qc:2**(0.292*beta) if qc == 'c' else 2**(0.265*beta) #BKZ cost(Kyber)
            BKZ_data[beta] = log2(T_BKZ(beta,d,qc)*itr)

            if qc == 'c':
                T_node = 200 #clock cycle per node of enumeration tree[CN11]
                T_ENUM = N_beta[beta]*T_node*P/p
                ENUM_data[beta] = log2(T_ENUM)
                TOTAL_data[beta] = log2(T_BKZ(beta,d,qc)*itr + T_ENUM)
                if (beta == beta_range[-1]):
                    RESULT.append([itr, beta, TOTAL_data[beta], BKZ_data[beta], ENUM_data[beta]])
                    break
                elif(TOTAL_data[beta] > TOTAL_data[beta-1]):
                    RESULT.append([itr, beta-1, TOTAL_data[beta-1], BKZ_data[beta-1], ENUM_data[beta-1]])
                    break
            else:
                b = 1/8 #<1/4
                gam = chernoff_bounds(2*b,0.5) #=20
                delta_D = 0.1 #failure rate of DETECTMV
                delta_F = 0.1 #failure rate of FINDMV

                prec = RR(4*0.3*d+4*log(q)+(log(d))**2+7*log(d)) 

                if mlwe == 'm' and qset == 's':
                    delta = delta_D/(n*itr)
                else:
                    delta = delta_D/itr

                if mlwe == 'm' and qset == 'p':
                    Tn = n * N_beta[beta] * d
                    P_depth = prec*(log2(q) + 2*d + (n+2)*prec^(0.158))
                    P_count = prec*d*((q+1)*d/2+98*(n+2)*prec^(0.585))
                    C = 2 * sqrt(n) + 28
                else:
                    Tn = N_beta[beta] * d
                    P_depth = prec*(log2(q) + d + 2*prec^(0.158))
                    P_count = prec*d*((q+1)*d/2+196*prec^(0.585))
                    C = 40
                #U_depth = 16*P_depth + 16*q**2 * log(q*sqrt(Tn)) + 8*q**2*log(q)
                U_count = 32*(q+1)*P_count + 16*q**2 * log(q*sqrt(Tn)) + 32*q**2*log(q)

                T_DETECT = RR(gam * delta_D/delta * log(1/delta) * sqrt(Tn)/b * U_count)
                T_FIND = RR(C * sqrt(Tn) * log(d) / (-log(1-delta_F)) * U_count)

                DETECT_data[beta] = log2(T_DETECT)
                FIND_data[beta] = log2(T_FIND)
                TOTAL_data[beta] = log2(T_BKZ(beta,d,qc)*itr + T_DETECT + T_FIND)
                if (TOTAL_data[beta] > TOTAL_data[beta-1]):
                    RESULT.append([itr,beta-1, TOTAL_data[beta-1], BKZ_data[beta-1], DETECT_data[beta-1], FIND_data[beta-1]])
                    break
                elif (beta == beta_range[-1]):
                    RESULT.append([itr,beta, TOTAL_data[beta], BKZ_data[beta], DETECT_data[beta], FIND_data[beta]])
                    break

        if (RESULT[-1][2] > RESULT[-2][2] ):
            return RESULT[-2]
        elif (RESULT[-1][0] == itr_range[-1]):
            return RESULT[-1]

#param = (n,k,m,q,sigma,beta_range,itr_range)

#param = (256,2,3,3329,sqrt(3/2),range(550,620,1),range(1,10,1)) #Kyber512
#param = (256,3,4,3329,sqrt(2/2),range(990,1024,1),range(1,10,1)) #Kyber768
#param = (256,4,5,3329,sqrt(2/2),[1279],range(1,300,10)) #Kyber1024

#P = 1/2

#TOTALTIME(*param,P,'c','n')
#TOTALTIME(*param,P,'c','m')
#TOTALTIME(*param,P,'q','n')
#TOTALTIME(*param,P,'q','m')
#TOTALTIME(*param,P,'q','m','p')