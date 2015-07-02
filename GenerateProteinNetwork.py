import random
import simplesbml as ss
import tellurium as te

numSpecies = 10
maxReactions = 30 # Total starting number of reactions (some will be deleted)
maxTries = 200 # If no model is found by then, execution halts
S0Concentration = 10.0 # Concentration of starting species
SxConcentration = 5.0 # Concentration of ending species
k0 = 0.3 # Rate constant for S0 --> S1 
kx = 0.3 # Rate constant for S10 --> Sx
kRatio = 1000 # Ratio of kcat to kM, in general

forbidZeros = False # Tosses models where 1+ species fall to zero concentration
limitConnections = True # Removes reactions that make nodes with too many edges
maxConnections = 7 

printReactions = True # Prints a summary of the reactions


class Reactions:
    def reaction1(self, species1, species2, species3, data):
        rate1 = round(kRatio*random.random(), 1)
        rate2 = round(random.random(), 1)
        data.append([species1, species3, species2, 0, rate1, rate2])
    def reaction2(self, species1, species2, species3, data):
        rate = round(random.random(), 1)
        data.append([species1, species2, species3, 0, rate])
    def reaction3(self, species1, species2, species3, data):
        rate = round(random.random(), 1)
        data.append([species1, 0, species2, species3, rate])
    def DeleteRow(self, species):
        for i in data:
            if species in i:
                data.remove(i)
    def printreactions(self):
        print("$S0 --> S1, rate constant = " + str(k0))
        for i in data:           
            sp1 = "S" + str(i[0])
            sp2 = "S" + str(i[1])
            sp3 = "S" + str(i[2])
            sp4 = "S" + str(i[3]) 
            if len(i) == 6:
                print(sp1 + " --> " + sp3 + ", " + sp2 + ", kcat = "
                        + str(i[4]) + ", KM = " + str(i[5]))
            elif i[1] == 0:
                print(sp1 + " --> " + sp3 + " + " + sp4 + ", rate constant = " 
                        + str(i[4]))
            else:
                print(sp1 + " + " + sp2 + " --> " + sp3 + ", rate constant = " 
                        + str(i[4]))
        print("S" + str(numSpecies) + " --> $Sx, rate constant = " + str(kx))
        
       
re = Reactions()
zcount = 0 # The total number of tries

while True: # Sets up a loop so that multiple models can be generated
    data = [] # Holds the master list of reactions
    zcount += 1    
    for x in range(1, maxReactions + 1):
        reactionType = random.randrange(1,4) # Select which type reaction next
        numbers = range(1, numSpecies+1)
        species1 = random.choice(numbers) # Grab random species
        numbers.remove(species1)
        species2 = random.choice(numbers)
        numbers.remove(species2)
        species3 = random.choice(numbers)
        if reactionType == 1:
            re.reaction1(species1, species2, species3, data)
        elif reactionType == 2:
            re.reaction2(species1, species2, species3, data)
        elif reactionType == 3:
            re.reaction3(species1, species2, species3, data)
     
    if limitConnections:          
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
    m.addParameter('k0', 0.3)
    m.addReaction(['S0'], ['S1'], 'k0*S0')
    m.addParameter('kx', 0.3)
    m.addReaction(["S" + str(numSpecies)], ['Sx'], "kx*S" + str(numSpecies))
    
    for x in range(2, numSpecies):
        m.addSpecies("S" + str(x), 5.0)
        
    count = 0
    
    for i in data:
        count += 1
        m.addParameter("k" + str(count), i[4])
        if len(i) == 6:
            m.addParameter("KM" + str(count), i[5])
            m.addReaction(["S" + str(i[0])], ["S" + str(i[2])], 
                          "(k" + str(count) + "*S" + str(i[0]) + 
                          "*S" + str(i[1]) + ")/(KM" + str(count) + "+S" + 
                          str(i[0]) + ")\n")                         
        elif i[3] == 0:
            m.addReaction(["S" + str(i[0]), "S" + str(i[1])], 
                          ["S" + str(i[2])], "k" + str(count) + "*S" + 
                          str(i[0]) + "*S" + str(i[1]))
        elif i[1] == 0:
            m.addReaction(["S" + str(i[0])], ["S" + str(i[1]), 
                          "S" + str(i[2])], "k" + str(count) + "*S" + 
                          str(i[0]))
        else:
            m.addReaction(["S" + str(i[0]), "S" + str(i[1])], ["S" + str(i[2]),
                          "S" + str(i[3])], "k" + str(count) + "*S" + 
                          str(i[0]) + "*S" + str(i[1]))
            
    
    te.saveToFile('ProteinModel.xml', str(m))
    r = te.loadSBMLModel('ProteinModel.xml')
    
    try:
        if forbidZeros:
            r.simulate()
            for i in r.getSteadyStateValues():
                if i < 1e-5:
                    raise ValueError
        r.simulate(end=1000)
        steady = r.steadyState()
        if steady < 10e-6:            
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
    except:
        if zcount >= maxTries:
            print ("Exceeded maximum number of attempts. Try increasing " 
                    "maxConnections or numReactions.")
            break
        else:
            pass
        




