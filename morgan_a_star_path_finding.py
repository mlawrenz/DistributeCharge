import heapq
import optparse
import numpy, pylab
import itertools
import pylab

class Cell(object):
    def __init__(self, x, y, cost, follow):
        """
        Initialize new cell

        @param x cell x coordinate
        @param y cell y coordinate
        
        @param distance is distance to reference at a given time
        """
        self.cost = cost
        self.follow=follow
        self.x = x
        self.y = y



class AStar(object):
    def __init__(self, basic):
        self.closed = []
        #heapq._heapify_max(self.opened)
        self.cells = []
        self.opened = []
        self.results=[]
        self.basic=basic

    def init_grid(self, reinitial=False):
        # want to assign 4 charges to the gridpoints
        n=0
        for (n, cell) in enumerate(self.basic):
            x=cell[0]
            y=cell[1]
            self.cells.append((x,y))
            self.opened.append((x,y))
            #self.cells.append(Cell(x, y, cost=0, follow=0))
            #self.opened.append(Cell(x, y, cost=0, follow=0))
        self.start=self.cells[0]
        self.end=self.cells[-1]

    def get_heuristic(self, cell1, cell2):
        """
        Compute the heuristic value H for a cell: distance between
        two cells being evaluated for charging

        @param cell
        @returns heuristic value H
        """
        #return numpy.sqrt((cell1.x - cell2.x)**2 + (cell1.y - cell2.y)**2)
        return numpy.sqrt((cell1[0] - cell2[0])**2 + (cell1[1] - cell2[1])**2)


    def process(self, nresidue):
        # start: looping through distance cells?
        # add tuples with (cost, cell1, cell2) for all unique combos
        #heapq.heappush(local_results, (-1*self.start.cost, (self.start.x, self.start.y)))
        self.closed=[]
        reference=self.start
        index=self.opened.index(reference)
        self.opened.pop(index)
        self.results.append(reference)
        follow=0
        while self.opened:
            cum_cost=0
            for i in self.opened:
                cost=0
                for reference in self.results:
                    tmp_cost=self.get_heuristic(reference, i)
                    if tmp_cost < 2:
                        break
                    else:
                        cost+=tmp_cost
                if cost > cum_cost:
                    cum_cost=cost
                    follow=i
            if follow not in self.opened:
                break
            else:
                ind=self.opened.index(follow)
                self.opened.pop(ind)
                self.results.append(follow)
        print self.results
                
        
    def visualize(self, nresidue):
        #visualize test
        pylab.figure()
        H=numpy.zeros((6,6))
        for coor in self.basic:
            ind1=coor[0]
            ind2=coor[1]
            H[ind1,ind2]=0.5
        rank=0
        for result in self.results:
            rank+=1
            ind1=result[0]
            ind2=result[1]
            H[ind1,ind2]=1.0
            if rank==nresidue:
                break
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
    nresidue=6
    a = AStar(basic)
    a.init_grid()
    a.process(nresidue)
    #heapq.heapify(a.results)
    #best=[heapq.heappop(a.results) for _ in range(0, len(a.results))]
    #print best
    a.visualize(nresidue)
