import os
import json
import glob


class ParseOutput(object):
    '''
    This class contains all the methods required to extract the data 
    specified from "aims.out" and write the data to json files. See 
    documentation for a description of each method and its parameters.
    '''

    def __init__(self, list_of_Properties_to_get, working_dir):
        '''
        This method retrieves the kinds of values the user wanted to save
        and initializes variables.
        '''
        self.working_dir = working_dir
        #The possible flags in the config file that would be included if the
        # user wanted to keep the data for that property
        list_of_Possible_Properties = ['cartesian_atomic_coordinates',\
        'Total_energy','vdW_energy_correction','Total_atomic_forces',\
        'Hirshfeld_volume','lattice_vector','MBD_at_rsSCS_energy',\
        'Maximum_force_component','after_each_convergence_cycle']
        # Make directory to save relaxation data
        out_dir = os.path.join(self.working_dir, "relaxation_data")
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        #store True for properties the user wants to save and False otherwise
        propertyBools = []

        for prop in list_of_Possible_Properties:
            if prop in list_of_Properties_to_get:
                propertyBools.append(True)
            else:
                propertyBools.append(False)

        [self.bcartesian_atomic_coordinates, self.bTotal_energy,\
         self.bvdW_energy_correction, self.bTotal_atomic_forces, \
         self.bHirshfeld_volume, self.blattice_vector,\
         self.bMBD_at_rsSCS_energy, self.bMaximum_force_component,\
         self.bafter_each_convergence_cycle\
         ] = propertyBools

        self.convergenceCycleNum = self.getInitCycleNum()
        self.lastOne = False
        self.resetVars()

    def resetVars(self):
        '''
        This method resets many of the data variables to their initial values
        for use in the next self-consistent convergence cycle.
        data is a dict to easily keep track of all data being parsed and to
        easily write the data to the json file.
        '''
        self.data = {}
        self.data['cartesian_atomic_coordinates'] = []
        self.data['Total_atomic_forces'] = []
        self.data['Hirshfeld_volume'] = []
        self.data['lattice_vector'] = []
        self.data['properties'] = {}
        self.Hirshfeld_volume = []
        self.converged = False

    def getInitCycleNum(self):
        dataDir = (os.path.join(self.working_dir, "relaxation_data", 
                  "data_for_convergence_cycle_*"))
        dataFileNames = glob.glob(dataDir)
        return len(dataFileNames)
 
    def parseFile(self, sFilePath):
        '''
        This func iterates through each line in the file and passes each line
        to switchCase
        '''
        searchStrings = ['Parsing geometry.in',\
                         'Self-consistency cycle converged',\
                         '|   Hirshfeld volume        :',\
                         '| vdW energy correction         :',\
                         '| MBD@rsSCS energy              :',\
                         '| Total energy                  :',\
                         'Total atomic forces (unitary forces cleaned) [eV/Ang]:',\
                         'Maximum force component is',\
                         'Updated atomic structure:',\
                         'Final atomic structure:'\
                         ]

        with open(sFilePath, 'r') as f:
            self.lines = f.readlines()

            self.findLastEntry(searchStrings)

            self.lineNum = 0
            while self.lineNum < len(self.lines):
                line = self.lines[self.lineNum]
                self.switchCase(line, searchStrings)
                self.lineNum += 1

    def incrementLineNumAndReturnSplitLine(self):
        '''
        This method splits the next line in the file into strings in a list
        e.g.
        atom        -4.21034890        5.31230613       -3.50454216  S
        becomes
        ['atom','-4.21034890','5.31230613','-3.50454216','S']
        '''
        self.lineNum += 1
        line = self.lines[self.lineNum]
        return line.split(), line
        
    def getToLineOfDesiredValue(self, subStrInDesiredLine, indexOfSubStrInSplitLine):
        '''
        The target value may be a few lines down from the search string used 
        to find the value. This method searches from the line of the search
        string line-by-line until the lines containing the target value is
        reached.
        '''
        splitLine, line = self.incrementLineNumAndReturnSplitLine()
    
        while self.lineNum < len(self.lines):
            try:
                bLineFound = subStrInDesiredLine == splitLine[indexOfSubStrInSplitLine]
            except IndexError:
                #IndexError here means line is an empty line ('').
                splitLine, line = self.incrementLineNumAndReturnSplitLine()
                continue
            if bLineFound:
                break
            else:
                splitLine, line = self.incrementLineNumAndReturnSplitLine()
                
        #exit if reached EOF
        if self.lineNum == len(self.lines):
            exit
        
    def getAtomList(self):
        '''
        This method populates atomList with a list of atomic symbols in the 
        molecule being evaluated.
        '''
        lineNum = self.lineNum
        try:
            #return if atomList is already populated
            return self.atomList
        except:
            #atomList isn't defined so this is the first find of 'geometry.in'
            pass
        #first set of coordinates has 'atom' as a substring in position 0 in 
        #the split line
        self.getToLineOfDesiredValue('atom', 0)
            
        self.atomList = []
        line = self.lines[self.lineNum]
        splitLine = line.split()
        
        #Because the lines with atomic forces will contain the atom #, this 
        #goes for all atoms dynamically
        while splitLine != [] and splitLine[0] == 'atom':
            #write one entry of x,y,z coordinate data per atom
            self.atomList.append(splitLine[4])
            #write one entry of x,y,z coordinate data per atom
            if self.bcartesian_atomic_coordinates:
                self.data['cartesian_atomic_coordinates'].append(self.getXYZ(line, \
                'cartesian_atomic_coordinates') + [self.atomList[-1]])
            
            splitLine, line = self.incrementLineNumAndReturnSplitLine()
        self.lineNum = lineNum
        if self.blattice_vector:
            self.latticeVectors() 
    def getXYZ(self, line, valueToGet):
        '''
        This func optionally returns the x,y,z coordinates of an atom or
        the forces acting in the x, y, and z directions on an atom
        '''
        if valueToGet == 'cartesian_atomic_coordinates':
            #_ is a place holder. See template.
            _,x,y,z,_ = line.split()
 
        elif valueToGet == 'Total_atomic_forces':
            _,_,x,y,z = line.split()
            x = self.convertSciNotationToNumber(x)
            y = self.convertSciNotationToNumber(y)
            z = self.convertSciNotationToNumber(z)
            
        elif valueToGet == 'lattice_vector':
            _,x,y,z = line.split()
        else:
            return []
                      
        return list(map(float,[x, y, z]))
    
    def HirshfeldVol(self, line):
        '''This func adds a hirshfeld volume of an atom to a list'''
        
        lineSplit = line.split()
        self.Hirshfeld_volume.append(float(lineSplit[4]))
        
        #only write hirshfeld volume data when data for all atoms has been
        #included
        if len(self.Hirshfeld_volume) != len(self.atomList):
            return
            
        for atomNum, volSet in enumerate(self.Hirshfeld_volume):
            self.data['Hirshfeld_volume'].append([volSet] + [self.atomList[atomNum]])
    
    def vdWCorrection(self, line):
        '''This func gets the vdWaals correction value'''
        #save the vdW energy correction value
        splitLine = line.split()
        self.vdWCorrectionValue = float(splitLine[7])
        self.data['vdW_energy_correction'] = self.vdWCorrectionValue

    def MBDEnergy(self, line):
        '''This func gets the MBD energy'''
        splitLine = line.split()
        self.MBDenergy = float(splitLine[6])
        self.data['MBD_at_rsSCS_energy'] = self.MBDenergy
    
    def totEnergy(self, line):
        '''This func gets the total energy''' 
        #save the Total energy value
        splitLine = line.split()
        self.totEnergyValue = float(splitLine[6])
        #include the Total energy
        self.data['Total_energy'] = self.totEnergyValue
    
    def convertSciNotationToNumber(self, numberInSciNotationAsString):
        num = numberInSciNotationAsString
        e='E'
        try:
            if num.find(e) == -1:
                #E not found
                if num.find('e') != -1:
                    e='e'
                else:
                    return str(num)
        except:
            #not a string
            return str(num)
        
        numberPortion = num[:num.find(e)]
        sign = num[num.find(e)+1:num.find(e)+2]
        exp = num[num.find(e)+2:]
                  
        if sign == '-':
            return float(numberPortion)*(10.0**(-1*float(exp)))
        elif sign == '+':
            return float(numberPortion)*(10.0**float(exp))
        else:
            exp = num[num.find(e)+1:]
            return float(numberPortion)*(10.0**float(exp))
        
    def getMaxForceComponent(self, line):
        '''This func gets the max force component'''
        splitLine = line.split()
        #convert scientific notation to number
        self.Maximum_force_component = self.convertSciNotationToNumber(splitLine[4])
        #include the max force component
        self.data['Maximum_force_component'] = self.Maximum_force_component
    
    def latticeVectors(self):
        '''This function saves the cartesian x,y,z of the lattice vectors'''
        self.getToLineOfDesiredValue('lattice_vector', 0)
        line = self.lines[self.lineNum]
        splitLine = line.split()
        while splitLine[0] == 'lattice_vector':
            self.data['lattice_vector'].append(self.getXYZ(line, 'lattice_vector'))
            splitLine, line = self.incrementLineNumAndReturnSplitLine()
            if len(splitLine) == 0:
                break
    
    def newAtomicStructure(self):
        '''This writes the cartesian atomic coordinates of all atoms to data'''
        self.getToLineOfDesiredValue('atom', 0)
        line = self.lines[self.lineNum]
        splitLine = line.split()
        
        atomNum = 0
        while splitLine[0] == 'atom':
            #write one entry of x,y,z coordinate data per atom
            self.data['cartesian_atomic_coordinates'].append(self.getXYZ(line, \
            'cartesian_atomic_coordinates') + [self.atomList[atomNum]])
            atomNum+=1
            splitLine, line = self.incrementLineNumAndReturnSplitLine()
            if len(splitLine) == 0:
                break
    
    def totForces(self):
        '''
        This writes the forces in the x,y,and z directions on all 
        atoms to data
        '''
        splitLine, line = self.incrementLineNumAndReturnSplitLine()
        atomNum = 0
        
        #Because the lines with atomic forces will contain the atom #, this 
        # goes for all atoms dynamically
        while int(splitLine[1]) == atomNum+1:
            self.data['Total_atomic_forces'].append(self.getXYZ(line, \
            'Total_atomic_forces') + [self.atomList[atomNum]])
            atomNum+=1
            splitLine, line = self.incrementLineNumAndReturnSplitLine()
    
            # a blank line will have len = 0 (which occurs in the line after
            # the forces list)
            if len(splitLine) == 0:
                break
    
    def goToLastCycleIfDesired(self):
        if self.bafter_each_convergence_cycle == False and self.lastOne == False:
            self.lineNum = self.lastConvergenceLineNum - 1
            self.lastOne = True
    
    def writeData(self):
        '''write data to JSON file'''
        fileName = os.path.join(self.working_dir, "relaxation_data", "data_for_convergence_cycle_"+\
        str(self.convergenceCycleNum)+".json")

        # Make compatible with GA structure format
        output_data = self.data
        try:
            output_data['geometry'] = self.data['cartesian_atomic_coordinates']
            output_data['properties']['lattice_vector_a'] = self.data['lattice_vector'][0]
            output_data['properties']['lattice_vector_b'] = self.data['lattice_vector'][1]
            output_data['properties']['lattice_vector_c'] = self.data['lattice_vector'][2]
            output_data['properties']['energy'] = self.data['Total_energy']
            del output_data['cartesian_atomic_coordinates']
            del output_data['lattice_vector']
        except: pass
             
        with open(fileName, 'w') as outfile:
            #json.dump(self.data, outfile, indent=4)
            json.dump(output_data, outfile, indent=4)
        #A new self-consistency cycle begins...
        self.convergenceCycleNum += 1
        self.resetVars()    
    
    def getGeometry(self):
        if self.blattice_vector:
            self.latticeVectors()
        if self.bcartesian_atomic_coordinates:
            self.newAtomicStructure()
    
    def findLastEntry(self, searchStrings):
        '''Find last occurrence of "Self-consistency cycle converged"'''
        self.foundFinal = False
        foundNextToLast = False
        self.lineNum = len(self.lines)-1
        while self.lineNum > 0:
            line = self.lines[self.lineNum]

            if searchStrings[9] in line:
                self.foundFinal = True
            
            if searchStrings[8] in line and (self.foundFinal or foundNextToLast):
                self.lastConvergenceLineNum = self.lineNum
                break
            elif searchStrings[8] in line and self.foundFinal == False:
                foundNextToLast = True
                self.lineNum -= 1
            else:
                self.lineNum -= 1
    
    def switchCase(self, line, searchStrings):
        '''
        This func searches the given line for the presence a search string and 
        executes the corresponding function when the substr is found in the
        converged data and if the user wanted to save the corresponding value.
        '''
        for strNum, searchStr in enumerate(searchStrings):
            if searchStr in line:
                if strNum == 0:
                    self.getAtomList()
                    #self.getGeometry()
               
                if strNum == 1:
                    self.converged = True
              
                if strNum == 2 and self.bHirshfeld_volume:
                    self.HirshfeldVol(line)
                
                elif strNum == 3 and self.converged \
                and self.bvdW_energy_correction:
                    self.vdWCorrection(line)
                    
                elif strNum == 4 and self.converged \
                and self.bMBD_at_rsSCS_energy:
                    self.MBDEnergy(line)
                    
                elif strNum == 5 and self.converged \
                and self.bTotal_energy:
                    self.totEnergy(line)
                    
                elif strNum == 6 and self.bTotal_atomic_forces:
                    self.totForces()
                    
                elif strNum == 7 and self.bMaximum_force_component:
                    self.getMaxForceComponent(line)
                    
                elif strNum == 8 and (self.lineNum != self.lastConvergenceLineNum \
                or self.bafter_each_convergence_cycle):
                    self.writeData()
                    self.goToLastCycleIfDesired()
                    self.getGeometry()
                    
        if self.lineNum == len(self.lines)-1 and self.foundFinal:
            self.writeData()
