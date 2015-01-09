import heapq
import optparse
import numpy, pylab
from mpl_toolkits.mplot3d import Axes3D
import itertools
import pylab

# helper functions





class Cell(object):
    def __init__(self, x, y, z, id, name):
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
        self.name=name



class AStar(object):
    def __init__(self, basic):
        self.opened = []
        self.results=[]
        self.basic=basic
        self.opened=basic

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
        # compute cumulative Coulombic force cost for each point to reference
        # points, deplete opened list of availables cells
        reference=self.opened[0]
        self.opened.pop(0)
        self.results.append(reference)
        follow=0
        while self.opened:
            cum_cost=10000
            for (index, cell) in enumerate(self.opened):
                cost_array=[self.get_heuristic(reference, cell) for reference in self.results]
                print (cell.x, cell.y, cell.z), cost_array
                if sum(cost_array) == cum_cost:
                    print "equal:", cell, follow
                total=sum(cost_array)
                if total < cum_cost: #minimize total columbic repulsive force
                    cum_cost=total
                    follow=cell
            if follow not in self.opened:
                print "BREAKING"
                break # no more choices at a lower cost
            print "FOLLOW:", (follow.x, follow.y, follow.z)
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
        x_b=[i.x for i in self.basic]
        y_b=[i.y for i in self.basic]
        z_b=[i.z for i in self.basic]
        ax.scatter(x_b, y_b, z_b, alpha=0.5, c='k', marker='o')
        ax.scatter(x_r, y_r, z_r, alpha=1.0, c='r', marker='o')
        pylab.show()

def find_basic(pdbfile):
    basic=[]
    fhandle=open(pdbfile)
    resid=0
    basic_resnames=['ARG', 'LYS', 'HIS', 'HIE', 'HIP', 'HID']
    for line in fhandle.readlines():
        if len(line.split()) > 10:
            type=line.split()[0]
            atom=line.split()[2]
            resname=line.split()[3]
            x=float(line.split()[6])
            y=float(line.split()[7])
            z=float(line.split()[8])
            if 'HET' in line.split()[0]:
                type=line.split()[0]
                resname=line.split()[3]
                resid+=1
                print "skipping unknown residue %s %s" % (type, resname)
            if type=='ATOM':
                resid+=1
                if atom == 'CA':
                    if resname in basic_resnames:
                        basic.append(Cell(x, y, z, resid, resname))
    return basic 

def parse_cmdln():
    import argparse
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-n','--nresidue',dest='nresidue', help='num. basic residues to be charges on the surface')
    parser.add_argument('-p','--pdb',dest='pdb', help='PDB file to be parsed for basic residues')
    args = parser.parse_args()
    return args

if __name__=="__main__":
    args=parse_cmdln()
    basic=find_basic(args.pdb)
    print [(i.id, i.name) for i in basic]
    nresidue=int(args.nresidue)
    a = AStar(basic)
    a.process(nresidue)
    a.visualize3d()


