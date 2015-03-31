import sys
import numpy, pylab
from mpl_toolkits.mplot3d import Axes3D
import itertools
import pylab

# helper functions





class Cell(object):
    def __init__(self, x, y, z, id, resname, pdbnum, pka, sasa=None):
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
        self.pka=pka



class AStar(object):
    def __init__(self, residues):
        self.opened = []
        self.results=[]
        self.residues=residues
        self.opened=residues

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
            if len(self.results)==nresidue:
                break
            cum_cost=10000
            for (index, cell) in enumerate(self.opened):
                cost_array=[self.get_heuristic(reference, cell) for reference in self.results]
                total=sum(cost_array)
                if total < cum_cost: #minimize total columbic repulsive force
                    cum_cost=total
                    follow=cell
            if follow not in self.opened:
                print "NO MORE CHOICES AT LOW COST"
                break # no more choices at a lower cost
            self.opened.pop(index)
            self.results.append(follow)
        return self.results
                
        
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
    acidic=[]
    intermed=[]
    fhandle=open(pdbfile)
    resid=0
    charge=0
    basic_charge=0
    basic_resnames=['ARG', 'LYS', 'HIS', 'HIE', 'HIP', 'HID']
    acidic_resnames=['ASP', 'GLU']
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
                    if resname=='ARG' or resname=='LYS':
                        charge+=1
                        basic.append(Cell(x, y, z, resid, resname, pdbnum, 'high'))
                    else:
                        intermed.append(Cell(x, y, z, resid, resname, pdbnum, 'low'))
                elif resname in acidic_resnames:
                    if resname=='GLU' or resname=='ASP':
                        charge+=-1
                    acidic.append(Cell(x, y, z, resid, resname, pdbnum, 'low'))
                else:
                    pass
            else:
                pass
    return basic, acidic, intermed, charge,


def order_by_sasa(list, sasa_list, res):
    import heapq
    reordered=[]
    heapq.heapify(reordered) # use heap to keep ordering of sasa 
    for cell in list:
        location=numpy.where(res==cell.pdbnum)[0]
        cell.sasa=sasa_list[location][0]
        if cell.sasa==0:
            continue
        heapq.heappush(reordered, (-1*sasa_list[location][0], cell)) #use (-1)*sasa bc its a min_heap
    surface=[]
    while reordered:
        sasa, cell=heapq.heappop(reordered)
        surface.append(cell)
    return surface
    
def parse_cmdln():
    import argparse
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-n','--input_charge',dest='input_charge', help='charge on the surface')
    parser.add_argument('-p','--pdb',dest='pdb', help='PDB file to be parsed for charge assignment')
    parser.add_argument('-s','--sasa',dest='sasa', help='per residue SASA for PDB file')
    args = parser.parse_args()
    return args

def print_charge_change(filename, array1, array2=0):
    file=open('charge.sh', 'w')
    if array2!=0:
        loops=[array1, array2]
    else:
        loops=[array1,]
    for loop in loops:
        for res in loop:
            if 'HI' in res.resname:
                target='HIP'
            if res.resname=='ASP':
                target='ASH'
            if res.resname=='GLU':
                target='GLH'
            if res.resname=='LYS':
                target='LYN'
            if res.pdbnum > 9:
                if res.pdbnum > 99:
                    if res.pdbnum > 999:
                        file.write('sed -i "s/%s X%s/%s X%s/g" noh-%s\n' % (res.resname, res.pdbnum, target, res.pdbnum, filename))
                    else:
                        file.write('sed -i "s/%s X %s/%s X %s/g" noh-%s\n' % (res.resname, res.pdbnum, target, res.pdbnum, filename))
                else:
                    file.write('sed -i "s/%s X  %s/%s X  %s/g" noh-%s\n' % (res.resname, res.pdbnum, target, res.pdbnum, filename))
            else:
                file.write('sed -i "s/%s X   %s/%s X   %s/g" noh-%s\n' % (res.resname, res.pdbnum, target, res.pdbnum, filename))
    file.close()
    print "wrote sed file for %s" % filename


if __name__=="__main__":
    args=parse_cmdln()
    pdbfile=args.pdb
    orig_basic, orig_acidic, orig_intermed, charge,=find_basic(args.pdb)
    input_charge=int(args.input_charge)
    maxcharge=charge+len(orig_intermed)
    if 'CytC' in pdbfile: # add for HEME
        "ADDING -2 for CytC"
        charge=charge-2
    print "target charge: %s" % input_charge
    print "solution charge: %s" % charge
    print "all basic sites charge: %s" % maxcharge
    res=numpy.loadtxt(args.sasa, usecols=(0,), dtype=int)
    sasa_list=numpy.loadtxt(args.sasa, usecols=(1,))
    basic=order_by_sasa(orig_basic, sasa_list, res)
    acidic=order_by_sasa(orig_acidic, sasa_list, res)
    intermed=order_by_sasa(orig_intermed, sasa_list, res)
    if input_charge==charge:
        print "SAME CHARGE AS SOLUTION"
        sys.exit()
    if input_charge > maxcharge: # not enough basic charges to add
        print "need to neutralize acidic sites, reduce negative"
        a = AStar(acidic)
        target=(input_charge-maxcharge)
        prot_acidic=a.process(target)
        print "protonate all %s basic residues: " % len(basic)
        #print [(i.pdbnum, i.resname, i.sasa) for i in basic]
        print "protonate all %s histidine residues: " % len(intermed)
        #print [(i.pdbnum, i.resname, i.sasa) for i in intermed]
        print "protonate these %s acidic residues out of %s: " % (len(prot_acidic), len(orig_acidic))
        #print [(i.pdbnum, i.resname, i.sasa) for i in prot_acidic ]
        print_charge_change(args.pdb, intermed, prot_acidic)
        sys.exit()
    if input_charge < maxcharge and input_charge > charge:
        print "add protonated histidines"
        a = AStar(intermed)
        target=(input_charge-charge)
        prot_hist=a.process(target)
        print "protonate all %s basic residues: " % len(basic)
        #print [(i.pdbnum, i.resname, i.sasa) for i in basic]
        print "protonate these %s histidine residues: " % len(prot_hist)
        #print [(i.pdbnum, i.resname, i.sasa) for i in prot_hist ]
        print_charge_change(args.pdb, prot_hist)
        sys.exit()
    if input_charge < charge: #too many highly basic charges
        print "need to neutralize basic sites, reduce positive"
        basic_lys=[i for i in basic if i.resname=='LYS' ]
        target=charge-input_charge
        a = AStar(basic_lys)
        only_lys=a.process(target)
        print "only protonate these %s basic residues: " % len(only_lys)
        print_charge_change(args.pdb, only_lys)
        sys.exit()
    #a.visualize3d()


