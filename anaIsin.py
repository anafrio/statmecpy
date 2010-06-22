#Functions for plotting evolution under sweep and distribution in statistical equilibrium

import scipy

def ReadSweepFromFile(file = "resultsIsingEvo.dat"):
        """
        Extracts data from evolution of up fractions and magnetizations from a file cretaed by ConstantCoupledIsingModel.WriteSweep2File or ConstantCoupledIsingModel.WriteFineSweep2File  
        """
        infile = open(file, 'r')
        iniline = infile.readline()
        nspin = int(iniline.split('=')[1])
        lines = infile.readlines()
        linenum = len(lines)

        times = []
        avgMag = {}
        for i in range(nspin):
            avgMag[i] = []
            
        fraSc = {}
        for i in range(2**nspin):
            fraSc[i] = []
        
        for i in range(linenum):
            line = lines[i].split('][\n')[0]
      
            mags = line.split('];[')
    
            times.append(int(mags[0]))
            for j in range(nspin):
                Avgmags = mags[1].split()
                avgMag[j].append(float(Avgmags[j]))
            for j in range(2**nspin):
                Frasc = mags[2].split()
                fraSc[j].append(float(Frasc[j]))
     
        infile.close()
      
        return [nspin, times, avgMag, fraSc]

def GenSweepEvoGraphsFromFile(file = "resultsIsingEvo.dat"):
    """
    Calls ReadSweepFromFile and draws a graph with all magnetizations and a graph with all up fractions
    """
    info = ReadSweepFromFile(file)
    nspin = info[0]
    times = info[1]
    mags = info[2]
    fras = info[3]
    print mags
    print fras
    GenSweepEvoGraph(times, mags, 1, 'AvgMag')
    GenSweepEvoGraph(times, fras, 2, 'FraScon')
    pylab.show()



def GenSweepEvoGraph(times, mags, fignum = 1, filesave = 'AvgMag'):
    """
    Plots and saves graph of all elements of mags vs time vs time
    """
    nspin = len(mags)
    pylab.figure(fignum)
    for i in range(nspin):
        pylab.plot(times, mags[i])
    pylab.savefig(filesave)
  
    


def runHeatBath2D(Nspin=2, N=10,  T=0.5, vJ0=1*scipy.ones(2, float),vh=scipy.zeros(2, float), MK=-0.5*(scipy.ones((2,2))-scipy.eye(2)) ,Intrange=1, nSweepsPerShow=1, nShow=50, Algo = 0):
    ising = ConstantCoupledIsingModel(Nspin, N, 2, T, vJ0,vh, MK ,Intrange, Algo)

    dl1 = animated2dLattice.DynamicLattice((N,N))
    dl2 = animated2dLattice.DynamicLattice((N,N))
    for t in range(nShow):
        ag = 0
        lattice1 = scipy.zeros((N,N))
        lattice2 = scipy.zeros((N,N))
        for i in range(N):

           for j in range(N):
                lattice1[i,j] = ising.SpinLattice[ag][0]
                lattice2[i,j] = ising.SpinLattice[ag][1]
                ag = ag+1

        dl1.display(lattice1)
        dl2.display(lattice2)
        ising.Sweep(nSweepsPerShow)
  
                

