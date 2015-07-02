import random
import simplesbml as ss
import tellurium as te
import os

numSpecies = 10
maxReactions = 30 # Total starting number of reactions (some will be deleted)
maxTries = 200 # If no model is found by then, execution halts
S0Concentration = 10.0 # Concentration of starting species
SxConcentration = 5.0 # Concentration of ending species
k0 = 0.3 # Rate constant for S0 --> S1 
kx = 0.3 # Rate constant for S10 --> Sx

printReactions = True # Prints a summary of the reactions in the finished model
printDecayReactions = True # Prints a summary of the decay reactions

forbidZeros = False # Tosses models where 1+ species fall to zero concentration

limitConnections = True # Removes reactions that make nodes with too many edges
maxConnections = 5

class Reactions:
    """ This class handles protein decay, as well as printing out all of the 
    reactions at the end of a successful attempt. (Protein decay is needed 
    because the species aren't destroyed in the random reactions.)
    """    
    def proteindecay(self):
        for x in range(1,numSpecies+1):
            rate = round(random.random(), 2)
            k = "".join(["kd", str(x)])
            m.addParameter(k, rate)
            sp1 = "".join(["S",  str(x)])
            m.addReaction([sp1], [], "*".join([k, sp1]))
            if printDecayReactions:
                decayData.append([sp1, rate])          
    def printreactions(self):
        print("$S0 --> S1, rate constant = " + str(k0))
        for i in data:           
            sp1 = "S" + str(i[0])
            sp2 = "S" + str(i[1])
            sp3 = "S" + str(i[2]) 
            if i[6] == 1:
                print(sp1 + " activates " + sp2 + ", k = " + str(i[3]))
            elif i[6] == 2:
                print(sp1 + " represses " + sp2 + ", k = " + str(i[3]))
            elif i[6] == 3:
                print(sp1 + " AND " + sp2 + " activate " + sp3 + ", k1 = " 
                    + str(i[3]) + ", k2 = " + str(i[4]))
            elif i[6] == 4:
                print(sp1 + " OR " + sp2 + " activate " + sp3 + ", k1 = " + 
                    str(i[3]) + ", k2 = " + str(i[4]))
            elif i[6] == 5:
                print(sp1 + " XOR " + sp2 + " activate " + sp3 + ", k1 = " +
                    str(i[3]) + ", k2 = " + str(i[5]) + ", k3 = " + str(i[6]))
            elif i[6] == 6:
                print(sp1 + " NAND " + sp2 + " activate " + sp3 + ", k1 = " +
                    str(i[3]) + ", k2 = " + str(i[5]) + ", k3 = " + str(i[6]))
            elif i[6] == 7:
                print(sp1 + " NOR " + sp2 + " activate " + sp3 + ", k1 = " +
                    str(i[3]) + ", k2 = " + str(i[5]) + ", k3 = " + str(i[6]))
            elif i[6] == 8:
                print(sp1 + " EQ " + sp2 + " activate " + sp3 + ", k1 = " +
                    str(i[3]) + ", k2 = " + str(i[5]) + ", k3 = " + str(i[6]))
            elif i[6] == 9:
                print(sp1 + " counter " + sp2 + " activate " + sp3 + ", k1 = "+
                    str(i[3]) + ", k2 = " + str(i[5]) + ", k3 = " + str(i[6]))
        print("S10 --> $Sx, rate constant = " + str(kx))           
    def printdecayreactions(self):
        for d in decayData:
            print(d[0] + " decays with k = " + str(d[1]))
                  
re = Reactions()
zcount = 0 # The total number of tries

while True: # Sets up a loop so that multiple models can be generated 
    data = [] # Holds the master list of reactions
    zcount += 1    
    for x in range(1, maxReactions + 1):
        reactionType = random.randrange(1,10) # Select which type reaction next
        speciesChoice = range(1, numSpecies+1) 
        species1 = random.choice(speciesChoice) # Grab random species
        speciesChoice.remove(species1)
        species2 = random.choice(speciesChoice)
        rate1 = round(random.random(), 2)
        if reactionType == 1: # activation
            data.append([species1, species2, 0, rate1, 0, 0, 1])
        elif reactionType == 2: # repression
            data.append([species1, species2, 0, rate1, 0, 0, 2])
        else:
            speciesChoice.remove(species2)
            species3 = random.choice(speciesChoice)
            rate2 = round(random.random(), 2)
            rate3 = round(random.random(), 2)
            if reactionType == 3:                   # AND
                data.append([species1, species2, species3, rate1, rate2, 0, 3])
            elif reactionType == 4:                 # OR
                data.append([species1, species2, species3, rate1, rate2, 0, 4])
            elif reactionType == 5:                 # XOR
                data.append([species1, species2, species3, rate1, rate2, rate3, 5])
            elif reactionType == 6:                 # NAND
                data.append([species1, species2, species3, rate1, rate2, rate3, 6])
            elif reactionType == 7:                 # NOR
                data.append([species1, species2, species3, rate1, rate2, rate3, 7])
            elif reactionType == 8:                 # EQ
                data.append([species1, species2, species3, rate1, rate2, rate3, 8])
            elif reactionType == 9:                 # counter
                data.append([species1, species2, species3, rate1, rate2, rate3, 9])
              
    for x in range(1, numSpecies + 1):
        count = 0
        for i in data:
            if x in i:
                if count >= maxConnections:
                    data.remove(i)
                else:
                    count += 1  
                                                                        
    m = ss.sbmlModel()
    m.addSpecies('$S0', S0Concentration)
    m.addSpecies('$Sx', SxConcentration)
    m.addSpecies('S1', 5.0)
    m.addSpecies("S" + str(numSpecies), 5.0)
    m.addParameter('k0', k0)
    m.addReaction(['S0'], ['S1'], 'k0*S0')
    m.addParameter('kx', kx)
    m.addReaction(["S" + str(numSpecies)], ['Sx'], "kx*S" + str(numSpecies))
    
    for x in range(2, numSpecies):
        m.addSpecies("S" + str(x), 5.0)
        
    count = 0
    
    for i in data:
        count += 1
        k1 = "k" + str(count) + "a"
        m.addParameter(k1, i[3])
        sp1 = "S" + str(i[0])
        sp2 = "S" + str(i[1])
        sp3 = "S" + str(i[2])
        s = "*"
        if i[6] == 1:
            m.addReaction([], [sp2], "(" + k1 + "*" + sp1 + ")/(1+" +
                k1 + "*" + sp1 + ")")                        
        elif i[6] == 2:
            m.addReaction([], [sp2], "(1)/(1+" + k1 + "*" + sp1 + ")")
        elif i[6] in range(3,5):
            k2 = "k" + str(count) + "b"
            m.addParameter(k2, i[4])
            if i[6] == 3:
                m.addReaction([], [sp3], "(" + s.join([k1, k2, sp1, sp2]) 
                    + ")/(1+" + s.join([k1, sp1]) + "+" + s.join([k2, sp2]) + 
                    "+" + s.join([k1, k2, sp1, sp2]) + ")")                         
            elif i[6] == 4:
                m.addReaction([], [sp3], "(" + s.join([k1, sp1]) + "+" +
                    s.join([k2, sp2]) + ")/(1+" + s.join([k1, sp1]) + "+" +
                    s.join([k2, sp2]) + ")")
        else:
            k2 = "k" + str(count) + "b"
            m.addParameter(k2, i[4])
            k3 = "k" + str(count) + "c"
            m.addParameter(k3, i[5])
            if i[6] == 5:
                m.addReaction([], [sp3], "(1)/(1+" + s.join([k1, sp1]) + "+" +
                    s.join([k2, sp2]) + "+" + s.join([k3, sp1, sp2]) + ")")
            elif i[6] == 6:
                m.addReaction([], [sp3], "(1+" + s.join([k1, sp1]) + "+" +
                    s.join([k2, sp2]) + ")/(1+" + s.join([k1, sp1]) + "+" +
                    s.join([k2, sp2]) + "+" + s.join([k3, sp1, sp2]) + ")")
            elif i[6] == 7:
                m.addReaction([], [sp3], "(" + s.join([k1, sp1]) + "+" + 
                    s.join([k2, sp2]) + ")/(1+" + s.join([k1, sp1]) + "+" +
                    s.join([k2, sp2]) + "+" + s.join([k3, sp1, sp2]) + ")")
            elif i[6] == 8:
                m.addReaction([], [sp3], "(1+" + s.join([k1, sp1, sp2]) +  
                    ")/(1+" + s.join([k1, sp1]) + "+" + s.join([k2, sp2]) + "+"
                    + s.join([k3, sp1, sp2]) + ")")                     
            elif i[6] == 9:
                m.addReaction([], [sp3], "(" + s.join([k1, sp1]) + ")/(1+" +
                    s.join([k1, sp1]) + "+" + s.join([k2, sp2]) + "+" +
                    s.join([k3, sp1, sp2]) + ")")
                    
    decayData = [] # Holds protein decay information
    re.proteindecay()        
                
    te.saveToFile('GeneticModel.xml', str(m))

    r = te.loadSBMLModel('GeneticModel.xml')
    try:
        if forbidZeros:
            r.simulate()
            for i in r.getSteadyStateValues():
                if i < 1e-5:
                    raise ValueError
        steady = r.steadyState() # Determine if final concentrations steady out
        if steady < 10e-6: # If they do, then it's a successful model           
            print "Success!"
            print ("Number of tries = " + str(zcount))
            zcount = 0 
            for i in data:
                zcount += 1
            print ("Number of reactions = " + str(zcount))
            if printReactions:
                print("\n====Printing reactions====\n")
                re.printreactions()    
            if printDecayReactions:
                print("\n====Printing decay reactions====\n")
                re.printdecayreactions()
            break
        else:
            raise ValueError           
    except BaseException as be:
        if zcount >= maxTries:
            print ("Exceeded maximum number of attempts. Try increasing " 
                    "maxConnections or numReactions.")
            os.remove('Model2.xml')
            break
        else:
            pass # If maxTries hasn't been reached, ignore error and try again
