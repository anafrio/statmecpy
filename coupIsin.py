#Classes and functions needed for the coupled ising models simulation
import scipy
import itertools
#import animated2dLattice
import pylab
import lattice

class ConstantCoupledIsingModel:
    """
    Coupled Ising model class with constant intra and inter coupling and external field.
    Nspin: number of different spin variables
    N: lattice size on each simension
    d: lattice dimension
    T: k*temperature
    vJ0: vector with the inter coupling constant for each spin variable
    vh: vector of external fields on each spin variable
    MK: symmetric matrix with intra coupling between every two different spin variables
    Intrange: 0 for infinite range internal interactions, 1 for nearest neighbors
    Algo: 0 for heat bath
    seed: seed for random number generator
    """

    def __init__(self, Nspin=2, N=10, d=2, T=0.2, vJ0=0.5*scipy.ones(2, float),vh=scipy.zeros(2, float), MK=scipy.array([[0., 1.],[1., 0.]]) ,Intrange=1, Algo = 0, seed=1):
        """
        Check Nspin, n and d are compatible with vJ0, vh dimensions
        Check MK is symmetric
        Sets objects characteristics to Nspim, T, vJ0, vh, MK, Intrange, Algo
        Builds array of possible spin configurations for each agent spinComb
        Call scipy.random.seed with argument seed
	Construct associcated lattice class
        Build self.Spinlattice to be random 50/50 0/1 for each spin
       	Call self.SetProbabilites 
	"""

        # Check that all arrays passed have the correct dimensions
        if Nspin != len(vJ0) or Nspin != len(vh):
            raise ValueError, "Coupling or external field vector don't have appropiate dimensions"
       

        if Nspin != 1 and (Nspin,Nspin) != MK.shape:
            raise ValueError, "Matrix of coupling between spin variables doesn't have appropriate dimensions"

        
        if scipy.transpose(MK).all()  != MK.all():
            raise ValueError, "Interaction between different spin variables must be symmetric"


        #Set the Ising model characteristics which are directly given by parameters
        self.Nspin = Nspin
        self.Intrange = Intrange
        self.Algorithm = Algo
        self.Temperature = T
        self.FieldVec = vh
        self.IntraCoupVec = vJ0
        self.InterCoupMat = MK
        #self.kB = 1.3806503*10**-23
        self.kB = 1

        #Array with possible combinations of the different spin variables for every agent
        self.SpinComb = 2*(scipy.array([i for i in itertools.product(range(2), repeat=self.Nspin)])-0.5)

  

        #print "[ConstantCoupledIsingModel] SpinComb", self.SpinComb

   
        #Initialization of the random seed
        if seed==None:
	    scipy.random.seed()
	else:
            scipy.random.seed(seed)
       

        #Lattice class initialization
        self.Lattice = lattice.Lattice(N,d)

        #List with each agent's possible state depending on how many of its neighbors has each spin variable up
        self.SetNeiUpComb()
       
        #List of arrays (one for each spin variable) initializing each variable for each agent randomly(0/1) 
        spinLattAux = []
        for i in range(Nspin):
            spinLattAux.append(scipy.random.random_integers(0,1,self.Lattice.Nag))

        #List of arrays (one for each agent) with each agent's value for each variable (0/1)
        self.SpinLattice = []
        for i in range(self.Lattice.Nag):
            spv = scipy.zeros(self.Nspin, float)
            for j in range(self.Nspin):
                spv[j] = spinLattAux[j][i] 
            self.SpinLattice.append(spv) 

        
        #Boltzmann probabilities computed depending to each agents own and neighboring configurations
        self.SetProbabilities()


   #INITIALIZATION FUNCTIONS

    def SetNeiUpComb(self):
        """
        Constructs the list of possible states depending on how many of the neighbors has each spin variable UP. Each element is an array:(N1Up, N2Up,...) 
        """
        if self.Intrange == 0:
            neiUpComb = [0]
        elif self.Intrange == 1:
            Nn = 2*self.Lattice.dim
            neiUpComb = [i for i in itertools.product(range(Nn+1), repeat=self.Nspin)]
        else:
            raise ValueError, "[SetNeiUpComb] Not such interaction range defined" 
        self.NeiUpComb = neiUpComb


    def SetProbabilities(self):
        """
        Builds vectors of fraction of agents in each possible state and of average magnetizations
         Calls appropiate function to calculate updating probabilities depending on the interaction range (Intrange is 0 for infinite range and 1 for nearest neighbors) and algorithm (0 for heat bath algorithm)
        """
        T = self.GetTemperature()
        vh = self.GetFieldVector()
        vJ0 = self.GetIntraCoupVector()
        MK = self.GetInterCoupMatrix()

        
        #Fraction per state (spin combination)(magnetization per spin variable is then dalculated from it) (len(SpinComb array))
        self.SetSpinCombFraVec()
        

        if self.Algorithm == 0 : #Heat bath algorithm 
            self.SetProbHB()
        else:
            raise ValueError, "[SetProbabilities] Not defined for such algorithm"
       
    def SetSpinCombFraVec(self):
        """
        Constructs the dictionary with the fraction of agents on each possible spin configuration
        """
        #Array of different possible states(-1/1)
        spinC = self.GetSpinCons()
        #Spin lattice (-1/1)
        lattice = [2*(i - 0.5) for i in self.SpinLattice]

        nag = self.GetNumberOfAgents()
        nspin = self.GetNspin()

        #Initialize to zero
        nCon = scipy.zeros(len(spinC), float)

        #Count number of agents in each state
        for i in range(len(spinC)):
            sc = spinC[i]
            cou = 0
            for j in lattice:
                if tuple(j) == tuple(sc):
                    cou = cou + 1
            nCon[i] = cou

        #Array with the fraction of agents in each state    
        self.SpinCFraVec = nCon/self.Lattice.Nag


        print "[SetSpinCombFraVec] SpinCFraVec :", self.SpinCFraVec
         
        # Average magnetization for each spin variable
        self.SetAvgMagFromSpinCFra()    
 
       

#   def SetAvgMagVec(self):
#         """
#         Constructs vector with the average magnetization for every spin variable
#         """
#         
         # We construct the vector of fraction of upspins
#         NUp = scipy.zeros(self.Nspin)
        
#         for i in range(self.Nspin):
#             for j in range(self.Lattice.Nag):
#                NUp[i] = NUp[i]+self.SpinLattice[j][i]
                
#         upFraVec = NUp/self.Lattice.Nag
        #print "[SetAvgMagVec] upFraVec: ", upFraVec

#        self.AvgMagVec = 2*upFraVec - 1.
         #print "[SetAvgMagVec] AvgMagVec: ", self.AvgMagVec
            


    
    def SetAvgMagFromSpinCFra(self):
        """
        Construncts vector with the average magnetization for every spin variable from the vector of fraction of agents in each possible state
        """
        #Array of different possible states(-1/1)
        spinC = self.GetSpinCons()
        #Array with the fraction of agents in each state 
        spinCFra = self.GetSpinCFraVector()

        nspin = self.GetNspin()

        #Initialization to zero of auxiliary variable
        avgMag = scipy.zeros(nspin,float) 
        #Average magnetization per spin variable
        for i in range(nspin):
            #Count for every possible state
            for j in range(len(spinC)):
                avgMag[i] = avgMag[i]+spinC[j][i]*spinCFra[j]

        
        self.AvgMagVec = avgMag 




#       conlist = []
     
#       for i in range(nspin):
#           conlist.append([]) 
#           for s in spinC:
#               if s[i] == 1:
#                   conlist[i].append(s)
                
#       self.AvgMagVec = scipy.zeros(nspin,float) 
#       for i in range(len(conlist)):
#           spinlis = conlist[i] 
#            for j in range(len(spinlis)):
#                for k in range(len(spinC)):
#                    if  tuple(spinC[k]) == tuple(spinlis[j]):
#                        self.AvgMagVec[i] = self.AvgMagVec[i]+spinCFra[k]

        #print "[SetAvgMagFromSpinCFra] AvgMagVec: ", self.AvgMagVec

    def SetProbHB(self):
        """
        For the heat bath algorith:
        Constructs array with boltzmann probabilities for every neighborhood configuration NeiComUp (only one in infinite range case) and every spin configuration
        Hamiltonian considered (2 spin variables)
        Ei = -J1/2* m1* si -J2/2*m2* ti - ksiti -h1si -h2ti
        where m1 and m2 are the average magnetization in the inifnite range case and the neighborhood absolute magnetization in the short range case
        """


        neiUpComb = self.GetNeiUpComb()
            
        hbpro = []
 
        for i in range(len(neiUpComb)):
            nei = neiUpComb[i]
            ECon = self.GetEnergySpinC(nei)
            HBp = list(self.GetHBProb(ECon))
            hbpro.append(HBp)

        self.HBProbC = scipy.array(hbpro)
        #print "[SetProbHB] HBProbC: ", self.HBProbC 

       
#        self.HBProbC = []
#        self.ProbUp = []
#        neiUpComb = self.GetNeiUpComb()
        
#        for i in range(len(neiUpComb)):
#            nei = neiUpComb[i]
#            ECon = self.GetEnergySpinC(nei)
#            HBp = self.GetHBProb(ECon)
#            self.HBProbC.append(HBp)
#            self.ProbUp.append(self.GetProbUpFromHBProb(HBp))


    def GetEnergySpinC(self, neiUp):
        """
        Fir a given neighbors up configuration it returns an array with the energy of each possible spin configuration
        """
        avgMag = self.GetAverageMagVector()
        sc = self.GetSpinCons()
        N = self.GetNumberOfAgents()
        vJ0 = self.GetIntraCoupVector()
        vh = self.GetFieldVector()
        T = self.GetTemperature()
        MK = self.GetInterCoupMatrix()
        Nn = 2*self.GetDimension() 
        nspin = self.GetNspin()

        if neiUp == 0:
            mag = avgMag
        else:
            mag = 2*scipy.array(neiUp)-Nn

#        ECon = {}
        ECon = scipy.zeros(len(sc))
        for n in range(len(sc)):
            edec = 0
            ecou = 0 
       
            for j in range(nspin):
                #decoupled model energy
                edec = edec - (0.5*vJ0[j]*mag[j]+vh[j])*sc[n][j]
          
                #coupled model energy
                if self.Nspin != 1 and MK.any() != 0.: 
                    for k in range(MK.shape[0]):
                        for l in range(k+1, MK.shape[1]):
                            ecou = ecou - MK[k,l]*sc[n,k]*sc[n,l]

            #pre = tuple(sc[n])
            #ECon[pre] = edec + ecou
            ECon[n] = edec + ecou
     
        #print "[GetEnergySpinC] ECon:", ECon
        return ECon


    def GetHBProb(self, ECon):
        """
        For a given  energy for each possible spin configuration dictionary, it returns an array with the Boltzmann probability for each spinC configuration
        """
       
        sc = self.GetSpinCons()
        T = self.GetTemperature()

        #bolCon = {}
        bolCon = scipy.zeros(len(sc))
        bolTot = 0
        #HBProbC = {}
        HBProbC = scipy.zeros(len(sc))
      
        if T != 0.:
            for j in range(len(sc)):
                #pre = tuple(sc[j])
                #bolCon[pre] = scipy.exp(-scipy.array(ECon[pre])/(self.kB*T))
                bolCon[j] = scipy.exp(-scipy.array(ECon[j])/(self.kB*T))
                #bolTot = bolTot + bolCon[pre]
                bolTot = bolTot + bolCon[j]
                           
            for j in range(len(bolCon)):
                HBProbC[j] = bolCon[j]/bolTot

        else:
            eminInfo = self.GetMinEnergy(ECon)
            emin = eminInfo[0]
            eminconf = eminInfo[1]
            for j in range(len(sc)):
                if tuple(sc[j]) in eminconf:
                    #HBProbC[tuple(sc[j])] = 1./len(eminconf)
                    HBProbC[j] = 1./len(eminconf)
                else:
                    #HBProbC[tuple(sc[j])] = 0.
                    HBProbC[j] = 0.
                        
               
        print "[GetHBProb] HBProbC:",HBProbC
        return HBProbC


    def GetMinEnergy(self,ECon):
        """
        For a given ECon array, returns a list with two elements: the first one is the minimum possible energy and the second one a list with the spin configurations with this energy
        """

             
        emin = 1000. 
        eminconf = []
        sc = self.GetSpinCons()
       
        for j in range(len(sc)):
            #pre = tuple(sc[j])
            #if ECon[pre] < emin:
            if ECon[j] < emin:
                #emin = ECon[pre]
                emin = ECon[j]

        for j in range(len(sc)):
            #pre = tuple(sc[j])
            #if ECon[pre] == emin:
            if ECon[j] == emin:
                eminconf.append(tuple(sc[j]))
       
        #print "[GetMinEnergy] emin:", emin
        #print "[GetMinEnergy] eminconf: ", eminconf
        
  
        return [emin, eminconf]

       

#    def GetProbUpFromHBProb(self, HBProbC):
#        """
#        Returns an array with the up probability of each spin variable as constructed form a given dictionary of spin configuration probabilities
#        """
#        ProbUp = scipy.zeros(self.Nspin, float)

#        cons = HBProbC.keys()
#        nspin = self.GetNspin()

#        for i in cons:
#            pro = HBProbC[i]
#            for k in range(nspin):
#                if  i[k] == 1:
#                    ProbUp[k] = ProbUp[k] + pro
           
                         
        #print "[GetProbUpFromHBProb] ProbUp:", ProbUp   
#        return ProbUp



# FUNCTIONS THAT RETURN THE OBJECT'S CHARACTERISTICS 
  

    #Returns objects properties
    def GetNumberOfAgents(self):
        """
        Returns the number of agents on the lattice
        """
        return self.Lattice.Nag

    def GetNspin(self):
        """
        Returns the number of agents on the lattice
        """
        return self.Nspin

    def GetDimension(self):
        """
        Returns lattice dimension
        """
        return self.Lattice.dim
 
    def GetTemperature(self):
        """
        Returns temperature (in kb)
        """
        return self.Temperature

    def GetFieldVector(self):
        """
        Returns vector of external fields
        """
        return self.FieldVec

    def GetIntraCoupVector(self):
        """
        Returns vector of internal couplings
        """
        return self.IntraCoupVec

    def GetInterCoupMatrix(self):
        """
        Returns matrix of intra interaction couplings
        """
        return self.InterCoupMat

    def GetAverageMagVector(self):
        """
        Returns vector of average magnetizations
        """
        return self.AvgMagVec

    def GetSpinCons(self):
        """
        Returns the different possible spin combinations for the model's Nspin 
        """
        return self.SpinComb

    def GetSpinComFraVector(self):
        """
        Returns fraction of the agents in avery possible states depending on the the sign of the different spin variables 
        """
        return self.SpinCombFraVec

    def GetNeiUpComb(self):
        """
        Returns a list with the possible values of neighborhood values according to the number of neighbors with each spin variable up
        """
        return self.NeiUpComb

    def GetSpinCFraVector(self):
        """
        Returns a vector with the fraction of agents in each state of SpinComb
        """
        return self.SpinCFraVec
 


    #Returns other quantities

    def GetUpFraVector(self):
        """
        Returns vector of fractions of the agents with spin up
        """
        return 0.5*(self.AvgMagVec-1)


    def GetUpNeiVector(self,agent):
        """
        Returns the number of up neighbors vector for the given agent
        """
        lattice = self.Lattice
        d = self.GetDimension() 
        hop = lattice.Hop
        spinLatt = self.SpinLattice

        upNeiVec = scipy.zeros(self.Nspin)

        agentnei = []

        for i in hop.keys():
            if i[0] == agent:
                agentnei.append(hop[i])

        for j in agentnei:
            for i in range(self.Nspin):
                if spinLatt[int(j)][i]==1:
                    upNeiVec[i] = upNeiVec[i]+1
        return upNeiVec    
        
    
        


    #MONTE CARLO SIMULATION FUNCTIONS
   
    #Agent updates

    def UpdateRndAgent(self):
        """
        Updates a randomly picked agent by calling the appropiate method depending on interaction range and algorithm 
        """
        agent = scipy.random.randint(0,self.Lattice.Nag)
        algo = self.Algorithm
        iran = self.Intrange

        if algo == 0:
            self.UpdateAgentHB(agent)
        else:
            raise ErrorValue, "UpdateRndAgent not defined for that algorithm and/or interaction range"
        

    def UpdateAgentHB(self, agent):
        """
        Updates given agent using heat bath algorithm and total probability for each spin up
        """
        
        neiCom = self.GetNeiUpComb()
        nspin = self.GetNspin()
        agentSpin = scipy.zeros(nspin, int)
        for i in range(nspin):
            agentSpin[i] = self.SpinLattice[agent][i] 

        
        rnd = 1. - scipy.random.random(nspin)

        agentNUpNei = self.GetUpNeiVector(agent)
        for i in range(len(neiCom)):
            if neiCom[i] == 0:
                ind = 0
            elif neiCom[i] == tuple(agentNUpNei):
                ind = i
          
        prob = scipy.zeros(self.Nspin)
     
        for i in range(self.Nspin):
            prob[i] = self.ProbUp[ind][i]

        for i in range(self.Nspin):
            if rnd[i] < prob[i]:
                self.SpinLattice[agent][i] = 1
            else:
                self.SpinLattice[agent][i] = 0
           
        self.SetSpinCombFraVec()
        
        if self.Intrange == 0:
            self.SetProbHB()
        
   #Monte Carlo step    

    def MCStep(self):
        """
        Runs a monte carlo step, i.e., updates self.Nag randomly chosen agents
        """
        n = self.Lattice.Nag

        for i in range(n):
            self.UpdateRndAgent()


    def Sweep(self, Ntimes = 1):
        """
        Runs Ntimes MC steps
        """
        for i in range(Ntimes):
            self.MCStep()

   #Sweeps writing evolution to files

    def WriteSweep2File(self, Ntimes = 1, file = "resultsIsingEvo.dat"):
        """
        Runs Ntimes MC steps writing to a file the models average magnetization and fraction of ups vector after every step
        """
        outfile = open(file, 'a')
        outfile.write('Nspin='+str(self.Nspin)+'\n')
 
        for i in range(Ntimes):
            self.MCStep()
            towrite = str(i)+'];'+str(self.GetAverageMagVector())+';'+str(self.GetSpinCFraVector())+';[\n'
            outfile.write(towrite)
     
        outfile.close()

    def WriteFineSweep2File(self, Ntimes = 1, file = "resultsFineIsingEvo.dat"):
            """
            Runs Ntimes MC steps writing to a file the models average magnetization and up fraction vectors after every agent update
            """
            outfile = open(file, 'a')
            n = self.Lattice.Nag
            for i in range(Ntimes):
                for j in range(n):
                    self.UpdateRndAgent()
                    towrite = str(i+j)+'];'+str(self.GetAverageMagVector())+';'+str(self.GetSpinCFraVector())+';[\n'
                    outfile.write(towrite)
            outfile.close()


    def Sweep2Equilibrium(self, Neq = 10, Ntol = 0.0001):
        """
        Runs MC steps until statistictal equilibrium is achieved (no variation in magnetizations) and then for Neq MC steps more returning a list of all average magnetization vectors and a list with all fractions in each state
        """
        avgMag = self.GetAverageMagVector()
        avgMagpre = avgMag
        fraSc = self.GetSpinCFraVector()
        nspin = self.GetNspin()
        
        avgMagList = []
        spinCFraList = []

        while True:
            avgMagprepre = avgMagpre
            avgMagpre = avgMag
            self.MCStep()
            avgMag = self.GetAverageMagVector()
            cond = scipy.zeros(nspin)
            condtot = 1
            avgDif = abs(avgMag-avgMagpre)
            avgDif2 = abs(avgMag-avgMagprepre)
            print "[Sweep2Equilibrium] avgDif:", avgDif
            print "[Sweep2Equilibrium] avgDif2:", avgDif2
            for j in range(len(avgMag)):
                if avgDif[j] < Ntol and avgDif2[j] < Ntol:
                    cond[j] = 1
                else:
                    cond[j] = 0
                condtot = condtot*cond[j]
            if condtot == 1:
                print "Equilibrium achieved"
                break
            
        for i in range(Neq):
            self.MCStep()
            avgMag = self.GetAverageMagVector()
            fraSc = self.GetSpinCFraVector()
            avgMagList.append(tuple(avgMag))
            spinCFraList.append(tuple(fraSc))

        return avgMagList, spinCFraList



    
