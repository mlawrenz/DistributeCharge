import heapq
import optparse
import numpy, pylab
from mpl_toolkits.mplot3d import Axes3D
import itertools
import pylab

class Cell(object):
    def __init__(self, x, y, z, id):
        """
        Initialize new cell

        @param x cell x coordinate
        @param y cell y coordinate
        @param z cell z coordinate
        """
        self.x = x
        self.y = y
        self.z = z
        self.id = id



class AStar(object):
    def __init__(self, basic):
        self.closed = []
        self.opened = []
        self.results=[]
        self.cells=[]
        self.basic=basic

    def init_grid(self, reinitial=False):
        # want to assign the gridpoints
        n=0
        for (n, cell) in enumerate(self.basic):
            x=cell[0]
            y=cell[1]
            z=cell[2]
            index=n
            self.opened.append(Cell(x, y, z, n))
            self.cells.append(Cell(x, y, z, n))
        return

    def get_heuristic(self, cell1, cell2):
        """
        Compute the heuristic value H for a cell: distance between
        two cells being evaluated for charging

        @param cell
        @returns heuristic value H
        """
        #return numpy.sqrt((cell1.x - cell2.x)**2 + (cell1.y - cell2.y)**2)
        r=numpy.sqrt((cell1.x - cell2.x)**2 + (cell1.y - cell2.y)**2+((cell1.z - cell2.z)**2))
        force=1.0/(r**2)
        return force


    def process(self, nresidue):
        # start: looping through distance cells?
        # add tuples with (cost, cell1, cell2) for all unique combos
        #heapq.heappush(local_results, (-1*self.start.cost, (self.start.x, self.start.y)))
        self.closed=[]
        reference=self.opened[0]
        self.opened.pop(0)
        self.results.append(reference)
        follow=0
        while self.opened:
            cum_cost=10000
            for (index, cell) in enumerate(self.opened):
                cost_array=[self.get_heuristic(reference, cell) for reference in self.results]
                if sum(cost_array) == cum_cost:
                    print "equal:", cell, follow
                if sum(cost_array) < cum_cost: #minimie columbic repulsive force
                    cum_cost=sum(cost_array)
                    follow=cell
            if follow not in self.opened:
                print "BREAKING"
                break # no more choices at a lower cost
            self.opened.pop(index)
            self.results.append(follow)
            if len(self.results)==nresidue:
                break
        print [(i.x, i.y, i.z, i.id) for i in self.results]
                
        
    def visualize3d(self):
        #visualize test
        fig=pylab.figure()
        ax=fig.add_subplot(111, projection='3d')
        H=numpy.zeros((6,6,6))
        ax = fig.add_subplot(111, projection='3d')
        x_r=[i.x for i in self.results]
        y_r=[i.y for i in self.results]
        z_r=[i.z for i in self.results]
        x_b=[i.x for i in self.cells]
        y_b=[i.y for i in self.cells]
        z_b=[i.z for i in self.cells]
        ax.scatter(x_b, y_b, z_b, alpha=0.5, c='k', marker='o')
        ax.scatter(x_r, y_r, z_r, alpha=1.0, c='r', marker='o')
        pylab.show()

def parse_cmdln():
    import os
    parser=optparse.OptionParser()
    parser.add_option('-n','--nresidue',dest='nresidue',type='string')
    parser.add_option('-p','--pdb',dest='pdb',type='string')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__=="__main__":
    #(options,args)=parse_cmdln()
    #basic=[]
    #for i in range(0,6):
    #    for j in range(0,6):
    #        for k in range(0,6):
    #            basic.append((i,j,k))
    basic = ((0, 5, 2), (1, 0 , 1), (1, 1, 1), (1, 5, 4), (2, 3, 3), 
             (3, 1, 4), (3, 2, 1), (3, 5, 1), (4, 1, 2), (4, 4, 4), (5, 1, 2))
    nresidue=3
    a = AStar(basic)
    a.init_grid()
    a.process(nresidue)
    a.visualize3d()


