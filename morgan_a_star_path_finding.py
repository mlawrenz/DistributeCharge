import heapq
import optparse
import numpy, pylab
import itertools
import pylab

class Cell(object):
    def __init__(self, x, y, index):
        """
        Initialize new cell

        @param x cell x coordinate
        @param y cell y coordinate
        
        @param distance is distance to reference at a given time
        """
        self.index = index
        self.x = x
        self.y = y



class AStar(object):
    def __init__(self):
        self.opened = []
        #heapq._heapify_max(self.opened)
        heapq.heapify(self.opened)
        self.cells = []
        self.results=0

    def init_grid(self, basic):
        # want to assign 4 charges to the gridpoints
        n=0
        for (n, cell) in enumerate(basic):
            x=cell[0]
            y=cell[1]
            self.cells.append(Cell(x, y, index=n))

    def get_heuristic(self, cell1, cell2):
        """
        Compute the heuristic value H for a cell: distance between
        two cells being evaluated for charging

        @param cell
        @returns heuristic value H
        """
        return numpy.sqrt((cell1.x - cell2.x)**2 + (cell1.y - cell2.y)**2)


    def process(self, nresidue):
        # start: looping through distance cells?
        # add tuples with (cost, cell1, cell2) for all unique combos
        combos = itertools.combinations(self.cells, nresidue)
        for combo in combos:
            cost=0
            subcombos = itertools.combinations(combo, 2)
            for subcombo in subcombos:
                reference = combo[0]
                next = combo[1]
                cost+=self.get_heuristic(reference, next)
            # make negative to get min heap (hack)
            heapq.heappush(self.opened, (-1*cost, [(i.x, i.y) for i in combo]))
        print "n= ", len(self.cells)
        results=[heapq.heappop(self.opened) for _ in range(0, len(self.opened))]
        self.results=results


def visualize(results, nresidue):
    #visualize test
    rank=0
    for result in results:
        rank+=1
        H=numpy.zeros((6,6))
        for coor in basic:
            ind1=coor[0]
            ind2=coor[1]
            H[ind1,ind2]=0.5
        print result
        n=0
        point=result[1]
        while n < nresidue: 
            ind1=point[n][0]
            ind2=point[n][1]
            H[ind1,ind2]=1.0
            n+=1
        pylab.figure()
        pylab.pcolor(H)
        pylab.colorbar()
        pylab.title('config %s' % rank)
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
    basic=[]
    for i in range(0,6):
        for j in range(0,6):
            basic.append((i,j))
    #basic = ((0, 5), (1, 0), (1, 1), (1, 5), (2, 3), 
    #         (3, 1), (3, 2), (3, 5), (4, 1), (4, 4), (5, 1))
    nresidue=5
    a = AStar()
    a.init_grid(basic)
    a.process(nresidue)
    visualize(a.results, nresidue)
