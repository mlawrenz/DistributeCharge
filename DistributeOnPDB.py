import sys
import numpy, pylab
from mpl_toolkits.mplot3d import Axes3D
import itertools
import pylab

# helper functions





class Cell(object):
    def __init__(self, x, y, z, id, resname, pdbnum, sasa=None):
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
        self.pdbnum=pdbnum
        self.resname=resname
        self.sasa=sasa



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
                total=sum(cost_array)
                if cell.sasa < 100:
                    continue
                if total < cum_cost: #minimize total columbic repulsive force
                    cum_cost=total
                    follow=cell
            if follow not in self.opened:
                print "NO MORE CHOICES AT LOW COST"
                break # no more choices at a lower cost
            self.opened.pop(index)
            self.results.append(follow)
            if len(self.results)==nresidue:
                break
        final=[(i.pdbnum, i.sasa) for i in self.results]
        return final
                
        
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
        if 'HET' in line.split()[0]:
            resid+=1
            continue
            print "skipping unknown residue %s %s" % (type, resname)
        type=line.split()[0]
        if type=='ATOM': # parse normal PDB line
            atom=line.split()[2]
            resname=line.split()[3]
            x=float(line.split()[6])
            y=float(line.split()[7])
            z=float(line.split()[8])
            try:
                int(line.split()[4])
                pdbnum=int(line.split()[4]) # a check for columns
            except ValueError:
                print "PDB contains chain column"
                pdbnum=int(line.split()[5]) # a check for columns
            if atom == 'CA': #only count alpha carbond
                resid+=1
                if resname in basic_resnames:
                    basic.append(Cell(x, y, z, resid, resname, pdbnum))
    return basic 


def order_by_sasa(basic, sasa_list, res):
    import heapq
    basic_reordered=[]
    heapq.heapify(basic_reordered) # use heap to keep ordering of sasa 
    for cell in basic:
        location=numpy.where(res==cell.pdbnum)[0]
        cell.sasa=sasa_list[location][0]
        heapq.heappush(basic_reordered, (-1*sasa_list[location][0], cell)) #use (-1)*sasa bc its a min_heap
    surface_basic=[]
    while basic_reordered:
        sasa, cell=heapq.heappop(basic_reordered)
        surface_basic.append(cell)
    return surface_basic
    


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
    res=numpy.loadtxt('residue_sasa.dat', usecols=(0,), dtype=int)
    sasa_list=numpy.loadtxt('residue_sasa.dat', usecols=(1,))
    basic=order_by_sasa(basic, sasa_list, res)
    #print "-----------------------------------"
    #print [(i.pdbnum, i.resname, i.sasa) for i in basic]
    nresidue=int(args.nresidue)
    a = AStar(basic)
    final=a.process(nresidue)
    numpy.savetxt('%s_charge_basic_sasa.dat' % args.pdb.split('.pdb')[0], final, fmt='%i')
    #a.visualize3d()


