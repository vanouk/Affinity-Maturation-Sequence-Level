import compute_descriptors_kayla
import csv
import sys
import modeller
import os
import importlib
import numpy as np                          # numerical tools
from copy import deepcopy                   # deepcopy copies a data structure without any implicit references
from timeit import default_timer as timer   # timer for performance
import importlib
import random
from Bio.Seq import Seq
import random
from collections import Counter
from collections import defaultdict
import warnings
warnings.filterwarnings("ignore")

###### Global parameters ######
    
p_mut        = 0.14                             # probability of mutation per division round
p_CDR        = 0.85                             # probability of mutation in the CDR region
p_CDR_lethal = 0.30                             # probability that a CDR mutation is lethal
p_CDR_silent = 0.00                             # probability that a CDR mutation is silent
p_CDR_affect = 1. - p_CDR_lethal - p_CDR_silent # probability that a CDR mutation affects affinity
p_recycle    = 0.70            			# probability that a B cell is recycled
p_exit       = 1. - p_recycle  			# probability that a B cell exits the GC

nb_GC_founders_I1 = 1
nb_GC_founders_I2or3 = 3
epsilon           = 1e-16
activation_energy = -9.011959719155849
nb_trial          = 1
energy_scale = 3.5/0.6           			# inverse temperature : 1/(8.3145/1000*310)
help_cutoff  = 0.70            			# only B cells in the top (help_cutoff) fraction of binders receive T cell help

# Upload dictionary of antigens
import dictionary_little_code
importlib.reload(dictionary_little_code)
from dictionary_little_code import dicAgs
from dictionary_little_code import dicconc
from dictionary_little_code import dicGCDur
from dictionary_little_code import flag
print('flag = ',flag)

#Upload seeding cells
import seedingBcell
importlib.reload(seedingBcell)
from seedingBcell import seedingCells, VRC01GL_nt, VRC01GL_nt_str
 
#Upload founder cells
if flag == 2 or flag == 3: 
    import founders
    importlib.reload(founders)
    from founders import founderseq

lenFRWH1=len(VRC01GL_nt[0][0])
lenCDRH1=len(VRC01GL_nt[0][1])
lenFRWH2=len(VRC01GL_nt[0][2])
lenCDRH2=len(VRC01GL_nt[0][3])
lenFRWH3=len(VRC01GL_nt[0][4])
lenCDRH3=len(VRC01GL_nt[0][5])
lenFRWL1=len(VRC01GL_nt[1][0])
lenCDRL1=len(VRC01GL_nt[1][1])
lenFRWL2=len(VRC01GL_nt[1][2])
lenCDRL2=len(VRC01GL_nt[1][3])
lenFRWL3=len(VRC01GL_nt[1][4])
lenCDRL3=len(VRC01GL_nt[1][5])   

ftrack = open('output-track.txt', 'w')
ftrack2 = open('output-track2.txt', 'w')
ftrack3 = open('output-track3.txt', 'w')
ftrack4 = open('output-track4.txt', 'w')
ftotTEST = open('output-totalTEST.csv', 'w')
fseed = open('seed.csv', 'w')

###### B cell clone class ######

class BCell:

    def __init__(self, nb = 128, **kwargs): ###FIX###
        """ Initialize clone-specific variables. 
            nb          - population size
            res         - array of the values for each residue of the B cell 
            E           - total binding energy
            nb_CDR_mut  - number of accumulated CDR mutations
            antigens    - list of antigens
            breadth     - breadth against test panel
            nb_Ag       - number of antigens available
            last_bound  - number of individuals that last bound each Ag
            mut_res_id  - index of mutated residue
            delta_res   - incremental change at residue
            generation  - generation in the GC reaction
            cycle_number - cycle number
            history     - history of mutations, generation occurred, and effect on Q/Ec and breadth against test panel """
        
        self.nb = nb    
	
        if 'nb_CDR_mut' in kwargs: self.nb_CDR_mut = kwargs['nb_CDR_mut']      ####NOTE: this is the number of mutations in the GL CODON sequence, not AA sequence (# aa is < or = # nt mutations)
        else:                      self.nb_CDR_mut = 0

        if 'res' in kwargs: self.res  = np.array(kwargs['res']) 
        else:
            print('res not recognized as an input argument')    

        if 'cycle_number' in kwargs: self.cycle_number = kwargs['cycle_number']
        else:                        self.cycle_number = 2

        if 'antigens' in kwargs: self.antigens = np.array(kwargs['antigens'])
        else: 			 self.antigens = np.array(dicAgs[self.cycle_number])

        if 'E' in kwargs: self.E = kwargs['E']  
        else:             
            self.res = str(self.res)
            self.E = self.energy(self.antigens[0],self.cycle_number) #assuming single Ag case!!!!### ###FIX###

        if 'breadth' in kwargs: self.breadth = np.array(kwargs['breadth'])
        else:                   self.breadth = 0
            
        if 'nb_Ag' in kwargs:              self.nb_Ag = kwargs['nb_Ag']
        else:                              self.nb_Ag = np.shape(self.antigens)[0]
        
        if 'last_bound' in kwargs: self.last_bound = kwargs['last_bound']
        else:                      self.last_bound = np.random.multinomial(self.nb, pvals = [1/float(self.nb_Ag)] * self.nb_Ag)
        
        if 'mut_res_id' in kwargs: self.mut_res_id = kwargs['mut_res_id']
        else:                      self.mut_res_id = 0
        
        if 'delta_res' in kwargs: self.delta_res = kwargs['delta_res']
        else:                     self.delta_res = 0

        if 'generation' in kwargs: self.generation = kwargs['generation']
        else:                      self.generation = 0

        if 'history' in kwargs: self.history = kwargs['history']
        else:                   self.history = {'generation' : [self.generation], 'cycle_number' : [self.cycle_number], 'res' : [self.res], 'nb_CDR_mut' : [self.nb_CDR_mut], 'mut_res_id' : [self.mut_res_id], 'E' : [self.E], 'delta_res' : [self.delta_res]}

    """ Return a new copy of the input BCell"""
    @classmethod
    def clone(cls, b):
        return cls(1, res = deepcopy(b.res), antigens = deepcopy(b.antigens), cycle_number = b.cycle_number, nb_Ag = b.nb_Ag, E = b.E, breadth = b.breadth, generation = b.generation, mut_res_id = b.mut_res_id, nb_CDR_mut = b.nb_CDR_mut, delta_res = b.delta_res, last_bound = deepcopy(b.last_bound), history = deepcopy(b.history))
                   
    def update_history(self):
        """ Add current parameters to the history list. """
        self.history['generation'].append(self.generation)
        self.history['res'].append(self.res)
        self.history['nb_CDR_mut'].append(self.nb_CDR_mut)
        self.history['mut_res_id'].append(self.mut_res_id)
        self.history['E'].append(self.E)      
        self.history['delta_res'].append(self.delta_res)
        self.history['cycle_number'].append(self.cycle_number)
    
    def energy(self, Ag, cycle_number):
        """ Return binding energy with input antigen. """ 
        Ab_nt = Seq(self.res)
        Ab_aa = str(Ab_nt.translate())
        Ab = Ab_aa[0:118]+'/'+Ab_aa[118:]
        print(cycle_number,self.nb_CDR_mut,Ab_aa)
        if "*" not in Ab_aa:
                print('Ag=',Ag)
                print('Ab=',Ab)
                dg = compute_descriptors_kayla.main(Ag, Ab)
        else: dg = 100
        return dg

    def divide(self, cycle_number):
        """ Run one round of division. """
        self.nb *= 2
        self.generation += 1
        self.cycle_number = cycle_number

    def pick_Antigen(self):
        """ Assuming nb_Ag > 1, return one antigen randomly chosen. """
        return self.antigens[np.random.randint(self.nb_Ag)]

    def calculate_breadth(self, testpanel, threshold, panelSize):   
        test_energies = [self.energy(testpanel[j]) for j in range(testpanelSize)]
        return float(sum(x > threshold for x in test_energies))/panelSize
        
    def update_antigens(self, newAntigens):
        self.antigens = deepcopy(newAntigens)
        self.nb_Ag = np.shape(newAntigens)[0]   

    def mutate_CDR(self, Ag, cycle_number):
        """ Change in energy due to affinity-affecting CDR mutation. Only one residue mutates."""
        """ Targeting model: choose nt to mutate based on 5mer mutability and sequence occurrences """
        """ Substitution model: choose nt to mutate to based on probability of middle nt in 5mer mutating to any other nt """

	#VRC01GL_str = []
	#VRC01GL_str.append(self.res)
        Fivemer = []
        Fivemer_index = []
        
        for i in range(len(str(self.res))):
            if i<(len(str(self.res))-4):
                if lenFRWH1<=(i+2)<lenFRWH1+lenCDRH1:
                    Fivemer_index.append(i+2)
                    Fivemer.append(str(self.res)[i+0:i+5])
                elif lenFRWH1+lenCDRH1+lenFRWH2<=(i+2)<lenFRWH1+lenCDRH1+lenFRWH2+lenCDRH2: 
                    Fivemer_index.append(i+2)
                    Fivemer.append(str(self.res)[i+0:i+5])
                elif lenFRWH1+lenCDRH1+lenFRWH2+lenCDRH2+lenFRWH3<=(i+2)<lenFRWH1+lenCDRH1+lenFRWH2+lenCDRH2+lenFRWH3+lenCDRH3: 
                    Fivemer_index.append(i+2)
                    Fivemer.append(str(self.res)[i+0:i+5])
                elif lenFRWH1+lenCDRH1+lenFRWH2+lenCDRH2+lenFRWH3+lenCDRH3+lenFRWL1<=(i+2)<lenFRWH1+lenCDRH1+lenFRWH2+lenCDRH2+lenFRWH3+lenCDRH3+lenFRWL1+lenCDRL1: 
                    Fivemer_index.append(i+2)
                    Fivemer.append(str(self.res)[i+0:i+5])
                elif lenFRWH1+lenCDRH1+lenFRWH2+lenCDRH2+lenFRWH3+lenCDRH3+lenFRWL1+lenCDRL1+lenFRWL2<=(i+2)<lenFRWH1+lenCDRH1+lenFRWH2+lenCDRH2+lenFRWH3+lenCDRH3+lenFRWL1+lenCDRL1+lenFRWL2+lenCDRL2: 
                    Fivemer_index.append(i+2)
                    Fivemer.append(str(self.res)[i+0:i+5])	
                elif lenFRWH1+lenCDRH1+lenFRWH2+lenCDRH2+lenFRWH3+lenCDRH3+lenFRWL1+lenCDRL1+lenFRWL2+lenCDRL2+lenFRWL3<=(i+2)<lenFRWH1+lenCDRH1+lenFRWH2+lenCDRH2+lenFRWH3+lenCDRH3+lenFRWL1+lenCDRL1+lenFRWL2+lenCDRL2+lenFRWL3+lenCDRL3: 
                    Fivemer_index.append(i+2)
                    Fivemer.append(str(self.res)[i+0:i+5])	

        reader    = csv.DictReader(open('Mutability.csv'))
        mutTable  = [obs for obs in reader]
        a = [float(o['Mutability']) for o in mutTable]
        mutSum = np.sum(a)
	
        probs = []
        for i in Fivemer:
            FivemerMutTable = selectObservations(mutTable, ['Fivemer'], [str(i)])
            a      = [float(o['Mutability'])           for o in FivemerMutTable]
            probs.append(a[0])
        probsNorm = probs
        probsNorm /= np.sum(probsNorm)
        picked_int = np.random.choice(range(len(Fivemer)), 1, p=probsNorm)
        picked_5mer = Fivemer[picked_int[0]]
        picked_5mer_index = Fivemer_index[picked_int[0]]
        
        reader2    = csv.DictReader(open('Substitution.csv'))
        mutTable2  = [obs for obs in reader2]
        bases = ['A', 'C', 'G', 'T']
        probs2 = []
        temp=[]
        temp.append(picked_5mer)

        for i in temp:
            FivemerMutTable2 = selectObservations(mutTable2, ['Fivemer'], [str(i)])
            A      = [float(o['A'])           for o in FivemerMutTable2]
            C      = [float(o['C'])           for o in FivemerMutTable2]
            G      = [float(o['G'])           for o in FivemerMutTable2]
            T      = [float(o['T'])           for o in FivemerMutTable2]
            probs2.extend([A, C, G, T])
            probs3=np.concatenate(probs2, axis=0)
            picked_list2=np.random.choice(bases, 1, p=probs3)
        picked_ntmut = picked_list2[0]
	
        temp_res1 = list(deepcopy(str(self.res)))
        index = picked_5mer_index
        self.mut_res_id = index
        delta = picked_ntmut
        self.delta_res = delta
        temp_res1[index] =  delta
        self.nb_CDR_mut += 1
        self.breadth = 0
        #self.breadth = self.calculate_breadth(testpanel, breadth_threshold,testpanelSize)
        temp_res1 = "".join(map(str, temp_res1))
        self.res = deepcopy(temp_res1) 
        self.E = self.energy(Ag,cycle_number)
        
        Ab_nt = Seq(str(self.res))
        Ab_aa = str(Ab_nt.translate())
        ftrack.write('%d,%lf,%d,%s' % (self.cycle_number, self.E, self.nb, str(Ab_aa)))
        ftrack.write('\n')
        ftrack.flush()
        self.update_history()

    def shm(self,cycle_number):
        """ Run somatic hypermutation and return self + new B cell clones. """
        
        # get number of cells that mutate, remove from the clone
        new_clones = []
        n_mut      = np.random.binomial(self.nb, p_mut)
        self.nb   -= n_mut
            
        # get number of CDR vs framework (FR) mutations
        n_CDR = np.random.binomial(n_mut, p_CDR)
        n_FR  = n_mut - n_CDR
            
        # process CDR mutations
        n_die, n_silent, n_affect  = np.random.multinomial(n_CDR, pvals = [p_CDR_lethal, p_CDR_silent, p_CDR_affect])
        self.nb                   += n_silent #add silent mutations to the parent clone
        print("naffect= %d" %(n_affect))
        ftrack2.write('n_affect = %d,%d,%d' %(n_affect, cycle_number, self.nb))
        ftrack2.write('\n')
        ftrack2.flush()
        for i in range(n_affect):
            b = BCell.clone(self)
            if b.nb_Ag > 1: Ag = b.pick_Antigen()
            else: Ag = b.antigens[0]
            b.mutate_CDR(Ag, cycle_number)
            new_clones.append(b)
        
        # return the result
        if (self.nb>0): new_clones.append(self)
        return new_clones


###### Main functions ######


def usage():
    print("")


def main(verbose=False):
    """ Simulate the affinity maturation process in a single germinal center (GC) and save the results to a CSV file. """
    
    # Run multiple trials and save all data to file
        
    start    = timer()
    
    fend  = open('output-end.csv', 'w')
    ftot  = open('output-total.csv',  'w')
    fbig  = open('output-largest-clone.csv', 'w')
    fend.write('trial,exit_cycle,number,generation,CDR_mutations,E,breadth,sequence\n')
    ftot.write('trial,cycle,number recycled,number exit,mean E,mean breadth,mean nb CDR mut\n')
    fbig.write('trial,cycle,index,CDR_mutations,E,delta_res,mut_res_id,sequence\n')
    
    # Events of a trial
    for t in range(nb_trial):
    
        print_update(t, nb_trial)   # status check

        # INITIALIZATION - DEFINE DATA STRUCTURES

        recycled_cells   = []
        exit_cells       = [] # cells at the end of the simulation
        memory_cells     = [] # exit cells from previous cycles
        nb_recycled      = []
        nb_exit          = []
        memory_founders  = []

        # CYCLES 1 + 2 - CREATE FOUNDERS AND REPLICATE WITHOUT MUTATION
        
        if flag == 1: 
            id_seeding_cells = np.random.choice(len(seedingCells), nb_GC_founders_I1, replace = False)
            B_cells = [BCell(res = seedingCells[id_seeding_cells[i]]) for i in range(nb_GC_founders_I1)]
            print(list(b.res for b in B_cells))
        elif flag == 2 or flag == 3: 
            id_founding_cells= np.random.choice(len(founderseq), nb_GC_founders_I2or3, replace = False)
            cdr_mut = []
            for seq in founderseq:
                mutations_aa=[i for i in range(len(VRC01GL_nt_str)) if VRC01GL_nt_str[i] != seq[i]]
                cdr_mut.append(len(mutations_aa))
            B_cells = [BCell(res = founderseq[id_founding_cells[i]], nb_CDR_mut = cdr_mut[id_founding_cells[i]]) for i in range(nb_GC_founders_I2or3)]
            print(list(b.res for b in B_cells))
        # Update data
        #cycle 0
        if flag == 1: nb_recycled.append(nb_GC_founders_I1) # all founders are recycled
        elif flag == 2 or flag == 3: nb_recycled.append(nb_GC_founders_I2or3)
        nb_exit.append(0)                                   # no founders exit the GC
        recycled_cells.append([deepcopy(b) for b in B_cells]) # add all cells of all 3 clones       
        
        #cycle 1
        nb_recycled.append(np.sum([b.nb for b in B_cells])) # all founders replicate and are recycled
        recycled_cells.append([deepcopy(b) for b in B_cells]) # add all cells of all 3 clones
        nb_exit.append(0)                                   # no founders exit
        
        # AFFINITY MATURATION
        
        GC_size_max  = nb_recycled[-1]  # maximum number of cells in the GC (= initial population size)
        cycle_number = 2
        nb_cycle_max = len(dicAgs)+ cycle_number -1     # maximum number of GC cycles
        
        while (cycle_number < nb_cycle_max): 
             
            cycleAntigens = np.array(dicAgs[cycle_number])
            nb_Ag = find_nb_Ag(cycleAntigens)           
            cycleconc = dicconc[cycle_number]
            cycledur  = dicGCDur[cycle_number]

            if cycle_number < cycledur:
                # keep same GC
                B_cells, out_cells = run_GC_cycle(B_cells, cycleAntigens, cycleconc, nb_Ag, cycle_number)
#            elif cycle_number == cycledur:
#                # start new GC
#                print('starting new GC at cycle number %d' % (cycle_number))
#
#                recycled_cells   = []
#                exit_cells       = []
#                nb_recycled      = []
#                nb_exit          = []
#                
#                memory_founders = pick_memCells_for_new_GC(memory_cells, nb_GC_founders_I2or3) 
#                for b in memory_founders: b.nb = 128
#                nb_recycled.append(nb_GC_founders_I2or3)                  
#                nb_exit.append(0)                                  
#                recycled_cells.append([deepcopy(b) for b in memory_founders]) 
#                nb_recycled.append(np.sum([b.nb for b in memory_founders])) 
#                recycled_cells.append([deepcopy(b) for b in memory_founders]) 
#                nb_exit.append(0)                                  
#                 
#                B_cells, out_cells = run_GC_cycle(memory_founders, cycleAntigens, cycleconc, nb_Ag, cycle_number)
#                memory_cells     = []
#                memory_founders  = []
#            else: 
#                print('error in starting a GC')
#                print(cycle_number)                
            
            GC_size = np.sum([b.nb for b in B_cells])       # total number of cells in the GC
            
            if (cycle_number == nb_cycle_max-1) or (GC_size >= GC_size_max): # at the end, all B cells exit the GC
                out_cells += B_cells
                nb_exit.append(np.sum([b.nb for b in out_cells]))
                memory_founders = pick_memCells_for_new_GC(memory_cells, nb_GC_founders_I2or3)
                print(len(memory_cells),len(memory_founders),nb_GC_founders_I2or3)
                for b in memory_founders:
                    fseed.write('%s\n' %(b.res))
                    fseed.flush()
            else:
                memory_cells += out_cells
                nb_exit.append(np.sum([b.nb for b in out_cells]))
                out_cells = []
            
            recycled_cells.append([deepcopy(b) for b in B_cells])
            exit_cells.append(out_cells)
            nb_recycled.append(GC_size)
           
            if (nb_recycled[-1] == 0):break
            elif (GC_size>GC_size_max): cycle_number = cycledur
            else: cycle_number += 1

            ftotTEST.write('trial,cycle,number recycled,number exit,mean E,mean breadth,mean nb CDR mut\n')
            for i in range(len(recycled_cells)):
                meanE = 0
                meanBreadth = 0
                meanCDRMutations = 0
                count_clones = 0
                num_cells = []
                if nb_recycled[i] > 0:
                    for b in recycled_cells[i]:
                        count_clones += 1
                        meanE += b.E
                        meanBreadth += b.breadth
                        meanCDRMutations += b.nb_CDR_mut
                        num_cells.append(b.nb)
                    indmax_num_cells = np.argmax(num_cells)
		    #bigClone = recycled_cells[indmax_num_cells][-1]
		    #Ab_nt = Seq(str(bigClone.res))
		    #Ab_aa = str(Ab_nt.translate())
                    meanE /= count_clones
                    meanBreadth /= count_clones
                    meanCDRMutations /= count_clones
                    cycle = recycled_cells[i][0].cycle_number
                ftotTEST.write('%d,%d,%d,%d,%lf,%lf,%lf\n' % (t, cycle, nb_recycled[i],nb_exit[i], meanE, meanBreadth, meanCDRMutations))
            ftotTEST.flush()
            
	#print(exit_cells)
        for i in range(len(exit_cells)):
	    #print(i)
            for b in exit_cells[i]:
		#print(i,b)
                Ab_nt = Seq(str(b.res))
                Ab_aa = str(Ab_nt.translate())
                
                fend.write('%d,%d,%d,%d,%d,%lf,%lf,%s' % (t, b.cycle_number, b.nb, b.generation, b.nb_CDR_mut, b.E, b.breadth, Ab_aa))
                fend.write('\n')
                for j in range(len(b.history['E'])):
                    Ab_nt = Seq(str(b.history['res'][j]))
                    Ab_aa = str(Ab_nt.translate())
                    fbig.write('%d,%d,%d,%d,%1f,%s,%s,%s,%s' % (t, b.history['cycle_number'][j], j, b.history['nb_CDR_mut'][j], b.history['E'][j], b.history['delta_res'][j], b.history['mut_res_id'][j], Ab_nt, Ab_aa))
                    fbig.write('\n')
                fbig.write('%d,%d,%d,%d,%1f,%s,%s,%s,%s' % (t, b.cycle_number, len(b.history['E']), b.nb_CDR_mut, b.E, b.delta_res, b.mut_res_id, Ab_nt, Ab_aa))
                fbig.write('\n')
        fend.flush()
        fbig.flush()
        
        for i in range(len(recycled_cells)):    
            meanE = 0
            #meanEc = 0
            meanBreadth = 0
            meanCDRMutations = 0
            count_clones = 0
            if nb_recycled[i] > 0:
                for b in recycled_cells[i]:
                    count_clones += 1
                    meanE += b.E
                    #meanEc += b.Ec
                    meanBreadth += b.breadth
                    meanCDRMutations += b.nb_CDR_mut          
                meanE /= count_clones
                #meanEc /= count_clones
                meanBreadth /= count_clones
                meanCDRMutations /= count_clones
                cycle = recycled_cells[i][0].cycle_number
                
            ftot.write('%d,%d,%d,%d,%lf,%lf,%lf\n' % (t, cycle, nb_recycled[i],nb_exit[i], meanE, meanBreadth, meanCDRMutations))
        ftot.flush()
            
    # End and output total time
    fend.close()
    ftot.close()
    ftotTEST.close()  
    ftrack.close()  
    ftrack2.close()
    ftrack3.close()
    ftrack4.close()
    fbig.close()
    fseed.close()
    end = timer()
    print('\nTotal time: %lfs, average per cycle %lfs' % ((end - start),(end - start)/float(nb_trial)))

def selectObservations(observations, keys, values):
    assert len(keys)==len(values), 'Mistmatched lengths for selectObservations!'
    selected = []
    for obs in observations:
        allowed = True
        for i in range(len(keys)):
            if obs[keys[i]]!=values[i]:
                allowed = False
                break
        if allowed: selected.append(obs)
    return selected

def find_nb_Ag(antigens):
    nb_Ag=np.shape(antigens)[0]
    return nb_Ag
       
def print_update(current, end, bar_length=20):
    """ Print an update of the simulation status. h/t Aravind Voggu on StackOverflow. """
    
    percent = float(current) / end
    dash    = ''.join(['-' for k in range(int(round(percent * bar_length)-1))]) + '>'
    space   = ''.join([' ' for k in range(bar_length - len(dash))])

    sys.stdout.write("\rSimulating: [{0}] {1}%\n".format(dash + space, int(round(percent * 100))))
    sys.stdout.flush()

def updating_antigens(B_cells, cycleAntigens):
    """ The antigens for all B cells are updated with the beginning of a new cycle. """
    for b in B_cells:
        b.update_antigens(cycleAntigens)    
    return B_cells
    
def run_dark_zone(B_cells, cycle_number, nb_rounds = 2):
    """ B cells proliferate and undergo SHM in the dark zone. """
    
    for i in range(nb_rounds):
        ftrack2.write('division/mutation round = %d' %(i+1))
        ftrack2.write('\n')
        ftrack2.flush()
        new_cells = []
        for b in B_cells:
            b.divide(cycle_number)
            new_cells += b.shm(cycle_number)
        B_cells = new_cells
    return B_cells

def run_binding_selection(B_cells, cycleconc, nb_Ag, cycle_number):
    """ Select B cells for binding to antigen. """
    
    new_cells=[]
    #print(B_cells)
    for b in B_cells:
	#print(b)
	#print("in run_bind_sel")
	#print(b.nb_CDR_mut, b.res)
        b.res=str(b.res)
        b.last_bound = np.random.multinomial(b.nb, pvals = [1./float(nb_Ag)] * nb_Ag)
        ftrack3.write('%d,%lf,%d' % (b.cycle_number, b.E, b.nb))       
        ftrack3.write('\n')
        for i in range(nb_Ag):
            # compute binding energy and chance of death ( = 1 - chance of survival )
            Ag_bound      = np.exp(-1 * energy_scale * (b.energy(b.antigens[i], cycle_number)-activation_energy)) ##changed from +1
            factor        = cycleconc * Ag_bound
            langmuir_conj = 1. / (1. + factor)
            
            # remove dead cells and update binding details
            n_die            = np.random.binomial(b.last_bound[i], langmuir_conj)
            b.nb            -= n_die
            b.last_bound[i] -= n_die
        ftrack3.write('%d,%lf,%d' % (b.cycle_number, b.E, b.nb))
        ftrack3.write('\n')
        if b.nb>0:new_cells.append(b)
    ftrack3.flush()
    return new_cells

def run_help_selection(B_cells, nb_Ag, cycle_number):
    """ Select B cells to receive T cell help. """
    #nb_Ag = B_cells[0].nb_Ag
    
    # get binding energies
    binding_energy     = [[b.energy(b.antigens[i], cycle_number) for i in range(nb_Ag)] for b in B_cells]
    binding_energy_tot = []
    for i in range(len(B_cells)):
        for j in range(nb_Ag): binding_energy_tot += [binding_energy[i][j]] * B_cells[i].last_bound[j]
    
    # cells in the top (help_cutoff) fraction of binders survive
    if len(binding_energy_tot)>0:
        cut_idx       = np.max([0, int(np.floor(help_cutoff * len(binding_energy_tot)))-1])
        energy_cutoff = np.array(binding_energy_tot)[np.argsort(binding_energy_tot)][::1][cut_idx]  ##changed from -1
        n_die_tie     = len(binding_energy_tot) - cut_idx - np.sum(binding_energy_tot > energy_cutoff) ##changed from <
        ftrack4.write('%d,%d,%lf,%d' % (cycle_number, cut_idx, energy_cutoff, n_die_tie))       
        ftrack4.write('\n')
        ftrack4.flush()
        # kill all B cells above threshold ##changed from "below"
        for i in np.random.permutation(len(B_cells)):
            for j in np.random.permutation(nb_Ag):
                energy = binding_energy[i][j]
                if energy > energy_cutoff:  ##changed from <
                    B_cells[i].nb            -= B_cells[i].last_bound[j]
                    B_cells[i].last_bound[j]  = 0
                elif (energy == energy_cutoff) and (n_die_tie > 0):
                    if B_cells[i].last_bound[j] < n_die_tie:
                        B_cells[i].nb            -= B_cells[i].last_bound[j]
                        n_die_tie                -= B_cells[i].last_bound[j]
                        B_cells[i].last_bound[j]  = 0
                    else:
                        B_cells[i].nb            -= n_die_tie
                        B_cells[i].last_bound[j] -= n_die_tie
                        n_die_tie                 = 0
    cells_surv = np.sum([b.nb for b in B_cells])   
    return B_cells, cells_surv
    

def run_recycle(B_cells):
    """ Randomly select B cells to be recycled back into the GC or to exit. """

    new_cells  = []                                 # cells that will remain in the GC
    exit_cells = []                                 # cells that will exit the GC
    n_tot      = np.sum([b.nb for b in B_cells])    # total number of cells currently in the GC
    n_exit     = int(np.floor(p_exit * n_tot))      # number of cells that will exit the GC
    b_exit     = np.array([])                       # index of cells that will exit the GC

    if (n_tot > 0) and (n_exit > 0):
        b_exit = np.random.choice(n_tot, n_exit, replace=False)

    idx = 0
    for b in B_cells:
    
        # find which cells exit the GC
        n_exit  = np.sum((idx <= b_exit) * (b_exit < idx + b.nb))
        idx    += b.nb
        b.nb   -= n_exit
        
        # add remainder to recycled cells
        if (b.nb > 0):
            new_cells.append(b)
    
        # record exit cells
        if (n_exit > 0):
            exit_cells.append(deepcopy(b))
            exit_cells[-1].nb = n_exit

    return new_cells, exit_cells

def pick_memCells_for_new_GC(memory_cells, nb_GC_founders):
    n_mem_cells = len(memory_cells)
    id_new_founders = np.random.choice(n_mem_cells, nb_GC_founders, replace=False)
    new_founders = [memory_cells[id_new_founders[i]] for i in range(nb_GC_founders)]
    return new_founders

def run_breadth_calculation(panel_energies, threshold, panelSize):
    average  = np.mean(panel_energies)
    variance = np.var(panel_energies)
    breadth  = float(sum(x > threshold for x in panel_energies))/panelSize 
    return average, variance, breadth

def run_GC_cycle(B_cells, cycleAntigens, cycleconc, nb_Ag, cycle_number):
    """ Run one cycle of the GC reaction. """
    B_cells = updating_antigens(B_cells, cycleAntigens)         # UPDATE antigens
    B_cells = run_dark_zone(B_cells, cycle_number)              # DARK  ZONE - two rounds of division + SHM + updates cycle_number
    total_cells = np.sum([b.nb for b in B_cells])
    
    if total_cells == 0: 
       print('GC extinct at cycle ', cycle_number)
       exit_cells = [] 
    else: 
        B_cells = run_binding_selection(B_cells, cycleconc, nb_Ag, cycle_number)  # LIGHT ZONE - selection for binding to Ag
        B_cells, cells_surv = run_help_selection(B_cells, nb_Ag, cycle_number)    # LIGHT ZONE - selection to receive T cell help
        B_cells, exit_cells = run_recycle(B_cells)
    
    return B_cells, exit_cells               # RECYCLE    - randomly pick exiting cells from the surviving B cells

if __name__ == '__main__': main()
