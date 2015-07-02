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

limitConnections = True # Removes reactions that make nodes with too many edges
maxConnections = 5 # Set the maximum number of edges for each node

forbidZeros = True # Tosses models where 1+ species fall to zero concentration

limitBoundarySpecies = True # Tosses models that have too many boundary species
maxBoundarySpecies = 4

class Reactions:
    """ This class adds each reaction as it is created to the data array. 
        It also handles changing species to boundary species and printing out
        all of the reactions at the end of execution.
    """
    def reaction1(self, species1, species2):
        rate = round(random.random(), 1)
        data.append([species1, 0, species2, 0, rate])
    def reaction2(self, species1, species2, species3):
        rate = round(random.random(), 1)
        data.append([species1, species2, species3, 0, rate])
    def reaction3(self, species1, species2, species3):
        rate = round(random.random(), 1)
        data.append([species1, 0, species2, species3, rate])
    def reaction4(self, species1, species2, species3, species4):
        rate = round(random.random(), 1)
        data.append([species1, species2, species3, species4, rate])
    def changetoboundary(self, species):
        speciesType[species-1] = 1
    def printreactions(self):
        if speciesType[0] == 1:
            print("$S0 --> $S1, rate constant = " + str(k0))
        else:
            print("$S0 --> S1, rate constant = " + str(k0))
        for i in data:           
            sp1 = "S" + str(i[0])
            sp2 = "S" + str(i[1])
            sp3 = "S" + str(i[2])
            sp4 = "S" + str(i[3]) 
            if speciesType[i[0]-1] == 1: sp1 = "$" + sp1
            if speciesType[i[1]-1] == 1: sp2 = "$" + sp2    
            if speciesType[i[2]-1] == 1: sp3 = "$" + sp3
            if speciesType[i[3]-1] == 1: sp4 = "$" + sp4  
            if i[1] == 0 and i[3] == 0:
                print(sp1 + " --> " + sp3 + ", rate constant = " + str(i[4]))
            elif i[1] == 0:
                print(sp1 + " --> " + sp3 + " + " + sp4 + 
                        ", rate constant = " + str(i[4]))
            elif i[3] == 0:
                print(sp1 + " + " + sp2 + " --> " + sp3 + ", rate constant = " 
                        + str(i[4]))
            else:
                print(sp1 + " + " + sp2 + " --> " + sp3 + " + " + sp4 + 
                        ", rate constant = " + str(i[4]))
        if speciesType[numSpecies-1] == 1:
            print("$S10 --> $Sx, rate constant = " + str(kx))
        else:
            print("S10 --> $Sx, rate constant = " + str(kx))
                          
re = Reactions()
zcount = 0 # The total number of tries

while True: # Sets up a loop so that multiple models can be generated 
    data = [] # Holds the master list of reactions
    speciesType = [0]*numSpecies # Holds info on which species are boundary
    zcount += 1     
    for x in range(1, maxReactions + 1):
        reactionType = random.randrange(1,5) # Select which type reaction next
        speciesChoice = range(1, numSpecies+1) 
        species1 = random.choice(speciesChoice) # Grab random species
        if reactionType == 1:
            speciesChoice.remove(species1)
            species2 = random.choice(speciesChoice) 
            re.reaction1(species1, species2)
        elif reactionType in range(2,4):
            species2 = random.choice(speciesChoice)
            speciesChoice.remove(species1)
            if species2 in speciesChoice:
                speciesChoice.remove(species2)
            species3 = random.choice(speciesChoice)
            if reactionType == 2:
                re.reaction2(species1, species2, species3)
            elif reactionType == 3:
                re.reaction3(species3, species2, species1)
        else:
            species2 = random.choice(speciesChoice)
            speciesChoice.remove(species1)
            if species2 in speciesChoice:
                speciesChoice.remove(species2) 
            species3 = random.choice(speciesChoice)            
            speciesChoice.remove(species3) 
            species4 = random.choice(speciesChoice)
            re.reaction4(species1, species2, species3, species4)
   
    if limitConnections:       
        for x in range(1, numSpecies + 1):
            count = 0
            for i in data:
                if x in i:
                    if count >= maxConnections:
                        data.remove(i) # Remove if involves overcrowded nodes
                    else:
                        count += 1
    
    for x in range(1, numSpecies + 1):
        count = 0
        for i in data:
            k = i[:2]
            count += k.count(x)
        if count <= 1:
            re.changetoboundary(x)
        else:
            count = 0
            for i in data:
                k = i[2:4]
                count += k.count(x)
            if count <= 1:
                re.changetoboundary(x)    
                                                                        
    m = ss.sbmlModel()
    m.addSpecies('$S0', S0Concentration)
    m.addSpecies('$Sx', SxConcentration)
    m.addSpecies('S1', 5.0)
    m.addSpecies("S" + str(numSpecies), 5.0)
    m.addParameter('k0', k0)
    m.addReaction(['$S0'], ['S1'], 'k0*S0')
    m.addParameter('kx', kx)
    m.addReaction(["S" + str(numSpecies)], ['$Sx'], "kx*S" + str(numSpecies))
    
    for x in range(2, numSpecies):
        if speciesType[x-1] == 1:
            m.addSpecies("$S" + str(x), 5.0)
        else:
            m.addSpecies("S" + str(x), 5.0)
        
    count = 0
    
    for i in data:
        count += 1
        m.addParameter("k" + str(count), i[4])
        sp1 = "S" + str(i[0])
        sp2 = "S" + str(i[1])
        sp3 = "S" + str(i[2])
        sp4 = "S" + str(i[3])
        if i[1] == 0 and i[3] == 0:
            m.addReaction([sp1], [sp2], "k" + str(count) + "*" + sp1)
        elif i[3] == 0:
            m.addReaction([sp1, sp2], [sp3], "k" + str(count) + "*" + sp1 
                            + "*" + sp2)
        elif i[1] == 0:
            m.addReaction([sp1], [sp2, sp3], "k" + str(count) + "*" + sp1)
                          
        else:
            m.addReaction([sp1, sp2], [sp3, sp4], "k" + str(count) + "*" + sp1 
                            + "*" + sp2)
                          
                
    te.saveToFile('Model.xml', str(m))

    r = te.loadSBMLModel('Model.xml')
    try:
        if forbidZeros:
            r.simulate()
            for i in r.getSteadyStateValues():
                if i < 1e-5:
                    raise ValueError
        if limitBoundarySpecies:
            if speciesType.count(1) > maxBoundarySpecies:
                raise ValueError
        steady = r.steadyState() #Determine if concentrations converge
        if steady < 10e-6: # If they do, it's a successful model
            print "Success!"
            print ("Number of tries = " + str(zcount))
            zcount = 0 
            for i in data:
                zcount += 1
            print ("Number of reactions = " + str(zcount))
            if printReactions:
                re.printreactions()                   
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
            pass #If maxTries hasn't been reached, ignore error and try again




