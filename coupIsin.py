#Classes and functions needed for the coupled ising models simulation
import scipy
import itertools
#import dynaLatt
import pylab
import defLatt

class ConsCoupIsin:
    """
    Coupled Ising model class with constant intra and inter coupling and external field.
    Initialization for:
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
    vacDen: constant fraction of unoccupied sites 
    """

    def __init__(self, Nspin=2, N=10, d=2, T=0.2, vJ0=0.5*scipy.ones(2, float),vh=scipy.zeros(2, float), MK=scipy.array([[0., 1.],[1., 0.]]) ,Intrange=1, Algo = 0, seed=1, vacDen = 0.0):
        """
        Calls SetChar to set objects characteristics to Nspin, T, vJ0, vh, MK, Intrange, Algo
        Calls SetSpinConf to build array of possible spin configurations for each agent spinConf
        Call scipy.random.seed with argument seed
	Construct associcated lattice object
        Calls SetNeiUpCon to build array of possible neigbor magnetization
        Calls SetSpinLatt to build self.SpinLatt to be random 50/50 0/1 for each spin
       	Call SetProb to calculate Boltzmann probabilities for each spin and neighbor configuration
	"""

        
        #Set the Ising model characteristics which are directly given by parameters
        self.SetChar(Nspin,Intrange,Algo,T,vh,vJ0,MK,vacDen,1)
        
        #Array with possible combinations of the different spin variables for every agent
        self.SetSpinConf()
       
        #Initialization of the random seed
        if seed==None:
	    scipy.random.seed()
	else:
            scipy.random.seed(seed)
  
        #Lattice class initialization
        self.Lattice = defLatt.Lattice(N,d)

        #List with each agent's possible state depending on how many of its neighbors has each spin variable up
        self.SetNeiUpComb()

        #List of arrays (one for each agent) with each agent's value for each variable (0/1)
        self.SetSpinLattRnd()

                
        #Boltzmann probabilities computed depending to each agents own and neighboring configurations
        self.SetProbabilities()


   #INITIALIZATION FUNCTIONS

    def SetChar(self,Nspin,Intrange,Algo,T,vh,vJ0,MK,vacDen,kB):
        """
        Check Nspin, n and d are compatible with vJ0, vh dimensions
        Check MK is symmetric
        Sets the object's characteristics that are directly equal to parameters passed on initialization: number of spin variables, interaction renge, algorith to be used in evolution, temperature, vector of external fields, vector of intra coupling, matrix of intra coupling and denity of vacancies.
        """
        # Check that all arrays passed have the correct dimensions
        if Nspin != len(vJ0) or Nspin != len(vh):
            raise ValueError, "Coupling or external field vector don't have appropiate dimensions"
            
        if Nspin != 1 and (Nspin,Nspin) != MK.shape:
            raise ValueError, "Matrix of coupling between spin variables doesn't have appropriate dimensions"
        
        if scipy.transpose(MK).all()  != MK.all():
            raise ValueError, "Interaction between different spin variables must be symmetric"
            
        self.Nspin = Nspin
        self.Intrange = Intrange
        self.Algorithm = Algo
        self.Temperature = T
        self.FieldVec = vh
        self.IntraCoupVec = vJ0
        self.InterCoupMat = MK
        self.VacDen = vacDen
        #self.kB = 1.3806503*10**-23
        self.kB = 1
        print "[SetCar] Nspin ", self.Nspin, " Intrange ", self.Intrange," Algorithm ", self.Algorithm," kb*Temperature ", self.Temperature*self.kB," FieldVec ", self.FieldVec," IntraCoupVec ", self.IntraCoupVec," InterCoupMat ", self.InterCoupMat," VacDen ", self.VacDen


    def SetSpinConf(self):
        """
        Constructs the array of possible states depending on the value of each spin variable. Array:(spin1, spin2,...) with spini=+/-1 
        """
        nspin = self.GetNspin()
        vacden = self.GetVacDen()

        spinconaux =[i-0.5*scipy.ones(nspin) for i in itertools.product(range(2), repeat=nspin)]
        if vacden != 0.0:
            spinconaux.append(scipy.zeros(nspin,int))
        
        self.SpinConf = 2*scipy.array(spinconaux)
        print "[SetSpinConf] SpinConf", self.SpinConf


    def SetNeiUpComb(self):
        """
        Constructs the list of possible states (which can be degenerate if there are vacancies) depending on how many of the neighbors has each spin variable UP. Each element is an array:(N1Up, N2Up,...). For infinite range it is [0]. 
        """
        intr = self.GetIntRange()
        dim = self.GetDimension()
        nspin = self.GetNspin()
        if intr == 0:
            neiUpComb = [0]
        elif intr == 1:
            Nn = 2*dim
            neiUpComb = [i for i in itertools.product(range(Nn+1), repeat=nspin)]
        else:
            raise ValueError, "[SetNeiUpComb] Not such interaction range defined" 
        self.NeiUpComb = neiUpComb

        print "[SetNeiUpComb] NeiUpComb ", self.NeiUpComb
        

    def SetSpinLattRnd(self):
        """
        Constructs the spin lattice, ie, randomly assigns a value -1/1 to each agent and the includes as many vacancies as needed at random sites
        """
        nspin = self.GetNspin()
        spinlattaux = []
        vacden = self.GetVacDen()
        nag = self.GetNumberOfAgents()
        for i in range(nag):
            spv = scipy.random.random_integers(0,1,nspin)
            spinlattaux.append(2*(spv-0.5))
            
        if vacden != 0.0:
            #Number of vacancies to be introduced
            nholes = int(nag*vacden)
            print "[SetSpinLattRnd] nholes ", nholes
            #Insert holes at random sites. These can be repeated yielding less hole than expected
            posarray = scipy.random.random_integers(0,nag-1,nholes)
            print "[SetSpinLattRnd] posaray ", posarray
            for i in range(len(posarray)):
                print i, posarray[i]
                spinlattaux[posarray[i]] = scipy.zeros(nspin,int)
            #Real vacancy density calculated
            def counthole(x):
                return x.all()==scipy.zeros(len(x),int).all()
            newNholes = float(len(filter(counthole,spinlattaux)))
            if newNholes/nag != vacden:
                self.VacDen = newNholes/nag
                Warning("Number of vacancies introduced is smaller than required. VACANCY DENSITY HAS BEEN REDIFINED to its real vale. If you need to use your previous attempt you may want to try to increase the systems size (or try lower densities)")
                print "WARNING:Number of vacancies introduced is smaller than required. VACANCY DENSITY HAS BEEN REDIFINED to its real vale ",self.VacDen,". If you need to use your previous attempt you may want to try to increase the systems size (or try lower densities)" 
        self.SpinLatt = spinlattaux
        print "[SetSpinLattRnd] SpinLatt ", self.SpinLatt

     
    def SetProbabilities(self):
        """
        Builds vectors of fraction of agents in each possible state and of average magnetizations
        Calls appropiate function to calculate updating probabilities depending on the interaction range (Intrange is 0 for infinite range and 1 for nearest neighbors) and algorithm (0 for heat bath algorithm)
        """
        T = self.GetTemperature()
        vh = self.GetFieldVector()
        vJ0 = self.GetIntraCoupVector()
        MK = self.GetInterCoupMatrix()
        algo = self.GetAlgorithm()
       
        #Fraction per state (spin combination)(magnetization per spin variable is then calculated from it) (len(SpinConf array))
        self.SetSpinConfFraVec()
 
        if algo == 0 : #Heat bath algorithms (probability to be in each state) 
            self.SetProbHB()
        else:
            raise ValueError, "[SetProbabilities] Not defined for such algorithm"
       
    def SetSpinConfFraVec(self):
        """
        Constructs the array with the fraction of agents on each possible spin configuration
        """
        #Array of different possible states(-1/1)
        spinC = self.GetSpinCons()
        #Spin lattice (-1/1)
        lattice = self.GetSpinLatt()
       
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
        self.SpinCFraVec = nCon/nag


        print "[SetSpinConfFraVec] SpinCFraVec :", self.SpinCFraVec
         
        # Average magnetization for each spin variable
        self.SetAvgMagFromSpinCFra()    

    
    def SetAvgMagFromSpinCFra(self):
        """
        Construncts vector with the average magnetization for every spin variable from the vector of fraction of agents in each possible state
        """
        spinC = self.GetSpinCons()
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
        print "[SetAvgMagFromSpinCFra] AvgMagVec: ", self.AvgMagVec


    def SetProbHB(self):
        """
        For the heat bath algorith:
        Constructs dictionary with boltzmann probabilities for every neighborhood magnetization  and every spin configuration (only present one in infinite range case)
        Hamiltonian considered (2 spin variables)
        Ei = -J1/2* m1* si -J2/2*m2* ti - ksiti -h1si -h2ti
        where m1 and m2 are the average magnetization in the inifnite range case and the neighborhood absolute magnetization in the short range case
        """
        neiUpComb = self.GetNeiUpComb()
        avgMag = self.GetAverageMagVector()
        dim = self.GetDimension()
        nn = 2*dim
        self.HBProb = {}

        if neiUpComb == [0]:
            mag = avgMag
            self.WriteHBProb(mag)
        else:
            for i in range(len(neiUpComb)):
                nei =list(neiUpComb[i])
                # Neighboring magnetization will also dependend on neighboring vacancies
                maxvac = nn-max(nei)                
                for j in range(maxvac):
                    mag = 2*scipy.array(nei)-j*scipy.ones(len(nei))
                    self.WriteHBProb(mag)
                   
        print "[SetProbHB] HBProb: ", self.HBProb

       
    def WriteHBProb(self,mag):
        """
        For a given neighboring magnetization, it appends conrresponding boltzmann probabilities to that neighboring magnetization and every spinC configuration
        """
        sc = list(self.GetSpinCons())
        T = self.GetTemperature()
        nspin = self.GetNspin()
                
        #Exclude unoccupied state as final possible sites
        for i in range(len(sc)):
            if max(list(sc[i]))==0:
                sc.pop(i)
        
        ECon = scipy.zeros(len(sc))
        bolCon = scipy.zeros(len(sc))
        bolTot = 0
        HBProbC = scipy.zeros(len(sc))
        
      
        if T != 0.:
            for j in range(len(sc)):
                ECon[j] = self.GetEnergy(mag, sc[j])
                bolCon[j] = scipy.exp(-ECon[j]/(self.kB*T))
                bolTot = bolTot + bolCon[j]
                           
            for j in range(len(sc)):
                HBProbC[j] = bolCon[j]/bolTot
                pre = (tuple(mag),tuple(sc[j]))
                if pre not in self.HBProb.keys():
                    self.HBProb[pre] = []
                self.HBProb[pre].append(HBProbC[j]) 

        else:
            for j in range(len(sc)):
                ECon[j] = self.GetEnergy(mag, sc[j])
                
            eminInfo = self.GetMinEnergy(ECon)
            emin = eminInfo[0]
            eminconf = eminInfo[1]
            for j in range(len(sc)):
                pre = (tuple(mag),tuple(sc[j]))
                if pre not in self.HBProb.keys():
                    self.HBProb[pre] = []
                if tuple(sc[j]) in eminconf:
                    HBProbC[j] = 1./len(eminconf)
                    self.HBProb[pre].append(HBProbC[j])

                else:
                    self.HBProb[pre].append(0.)
                        
               
        print "[WriteHBProb] HBProb:", self.HBProb

    def GetEnergy(self, mag, sc):
        """
        Gets energy value for a neighboring magnetization and  spin configuration
        """
        vJ0 = self.GetIntraCoupVector()
        vh = self.GetFieldVector()
        T = self.GetTemperature()
        MK = self.GetInterCoupMatrix()
        Nn = 2*self.GetDimension() 
        nspin = self.GetNspin()

        edec = 0
        ecou = 0 
       
        for j in range(nspin):
            #decoupled model energy
            edec = edec - (0.5*vJ0[j]*mag[j]+vh[j])*sc[j]
            #coupled model energy
            if nspin != 1 and MK.any() != 0.: 
                for k in range(MK.shape[0]):
                    for l in range(k+1, MK.shape[1]):
                        ecou = ecou - MK[k,l]*sc[k]*sc[l]

        Ene = edec+ecou
      
        print "[GetEnergy] Ene: ",Ene
        return Ene


    def GetMinEnergy(self,ECon):
        """
        For a given ECon array, returns a list with two elements: the first one is the minimum possible energy and the second one a list with the spin configurations with this energy
        """
             
        emin = 1000. 
        eminconf = []
        sc = self.GetSpinCons()
                   
        for j in range(len(sc)):
            if ECon[j] < emin:
                emin = ECon[j]

        for j in range(len(sc)):
            if ECon[j] == emin:
                eminconf.append(tuple(sc[j]))
       
        print "[GetMinEnergy] emin:", emin
        print "[GetMinEnergy] eminconf: ", eminconf
        
  
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
        return self.SpinConf

    def GetSpinComFraVector(self):
        """
        Returns fraction of the agents in avery possible states depending on the the sign of the different spin variables 
        """
        return self.SpinConfFraVec

    def GetNeiUpComb(self):
        """
        Returns a list with the possible values of neighborhood values according to the number of neighbors with each spin variable up
        """
        return self.NeiUpComb

    def GetSpinCFraVector(self):
        """
        Returns a vector with the fraction of agents in each state of SpinConf
        """
        return self.SpinCFraVec

    def GetVacDen(self):
        """
        Returns the real fraction of vacancy densitie
        """
        return self.VacDen

    def GetIntRange(self):
        """
        Returns interaction range for the system: 0 for infinite range, 1 for nearest neighbors
        """
        return self.Intrange

    def GetAlgorithm(self):
        """
        Returns algorithm to used for evolution: 0 for pure heat bath
        """
        return self.Algorithm
 
    def GetSpinLatt(self):
        """
        Returns list with the values of all spin variables for every agent
        """
        return self.SpinLatt
 


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
        spinLatt = self.SpinLatt

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


    
def GetVacs(self,agent):
        """
        Returns the number of vacancies in the  for the given agent
        """
        lattice = self.Lattice
        d = self.GetDimension() 
        hop = lattice.Hop
        spinLatt = self.SpinLatt

        vacs = 0

        agentnei = []

        for i in hop.keys():
            if i[0] == agent:
                agentnei.append(hop[i])

        for j in agentnei:
            if spinLatt[int(j)].all()==0:
                    vacs = vacs + 1
        return vacs
        
    

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
            agentSpin[i] = self.SpinLatt[agent][i] 

        
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
                self.SpinLatt[agent][i] = 1
            else:
                self.SpinLatt[agent][i] = 0
           
        self.SetSpinConfFraVec()
        
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



    
