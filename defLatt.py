#Classes and functions needed for the deffinition of a lattice

import scipy

class Lattice:
    """
    Lattice class
    """

    def __init__(self, N=10, dim=2):
        """
        Calls SetHop
        """
        self.SetHop(N,dim)


    def SetHop(self, N, dim):
        """
        Sets self.N to N and self.dim to dim
        Calls SetHop for the appropiate dimension
        """
        self.N = N
        self.dim = dim

        if dim==1:
            self.SetHop1d(N)
        elif dim==2:
            self.SetHop2d(N)
        else:
            print "SetHop not defined for such dimensions"
            return None
        
    def SetHop1d(self, N):
        """
        Builds AgGrid lattice array and hop dictionary for 1d lattices
        """
              
        self.Hop = {}

        self.Nag = N
        self.AgGrid = scipy.array(range(N),int)

        #print "[Lattice.SetHop1d] AgGrid: ", self.AgGrid

        for i in self.AgGrid:
            preM = (i, 1)
            prem = (i, -1)
            self.Hop[preM] = (i+1)%N
            self.Hop[prem] = (i-1)%N

        #print "[Lattice.SetHop1d] Hop: ", self.Hop
        
    def SetHop2d(self, N):
        """
        Builds AgGrid lattice array and hop dictionary for 1d lattices
        """
        self.Hop = {}

        self.Nag = N**2

        self.AgGrid = scipy.zeros((N,N), int)
        
        count = 0
        for i in range(N):
            for j in range (N):
                self.AgGrid[i,j]= count
                count = count + 1
        
        print "[Lattice.SetHop2d] AgGrid: ", self.AgGrid
        
        for i in range(N):
            for j in range(N):
                preM1 = (self.AgGrid[i,j], 1)
                prem1 = (self.AgGrid[i,j], -1)
                preM2 = (self.AgGrid[i,j], 2)
                prem2 = (self.AgGrid[i,j], -2)
           
                self.Hop[preM1] = self.AgGrid[(i+1)%N,j]
                self.Hop[prem1] = self.AgGrid[(i-1)%N,j]

                self.Hop[preM2] = self.AgGrid[i,(j+1)%N]
                self.Hop[prem2] = self.AgGrid[i,(j+1)%N]

        print "[Lattice.SetHop2d] Hop: ", self.Hop
    
                
