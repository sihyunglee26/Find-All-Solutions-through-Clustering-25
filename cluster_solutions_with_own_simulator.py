'''
Authored by Sihyung Lee
'''
import math
import random


def generateAnswersWithOneCluster(searchSpaceSize, numTargets):
    '''
    Generate solutions for Scenario-1
        cluster 1: 111..111xxx..xxx
    '''
    assert math.log2(searchSpaceSize) == int(math.log2(searchSpaceSize))
    assert math.log2(numTargets) == int(math.log2(numTargets))

    answerSet = set()
    for i in range(searchSpaceSize-1, searchSpaceSize-numTargets-1, -1): answerSet.add(i) # cluster 1: 111..111xxx..xxx
    return answerSet


def generateAnswersWithTwoEqualClusters(searchSpaceSize, numTargets):
    '''
    Generate solutions for Scenario-2
        cluster 1: 111..111xxx..xxx
        cluster 2: 000..000xxx..000
    '''
    assert math.log2(searchSpaceSize) == int(math.log2(searchSpaceSize))
    assert math.log2(numTargets) == int(math.log2(numTargets))

    answerSet = set()
    numTargetsInCluster = int(numTargets/2)    
    for i in range(searchSpaceSize-1, searchSpaceSize-numTargetsInCluster-1, -1): answerSet.add(i) # cluster 1: 111..111xxx..xxx
    for i in range(numTargetsInCluster): answerSet.add(i) # cluster 2: 000..000xxx..000
    return answerSet


def generateAnswersWithFourEqualClusters(searchSpaceSize, numTargets):
    '''
    Generate solutions for Scenario-4
        cluster 1: 111..111xxx..xxx
        cluster 2: 000..000xxx..000
        cluster 3: 000..111xx..x (equal numbers of 0's and 1's exist in the prefix)    
        cluster 4: 111..000xx..x (equal numbers of 0's and 1's exist in the prefix)    
    '''
    assert math.log2(searchSpaceSize) == int(math.log2(searchSpaceSize))
    assert math.log2(numTargets) == int(math.log2(numTargets))

    answerSet = set()
    numTargetsInCluster = int(numTargets/4)    
    for i in range(searchSpaceSize-1, searchSpaceSize-numTargetsInCluster-1, -1): answerSet.add(i) # cluster 1: 111..111xxx..xxx   
    for i in range(numTargetsInCluster): answerSet.add(i) # cluster 2: 000..000xxx..000
    
    postfixLength = int(math.log2(numTargetsInCluster))
    prefixLength = int(math.log2(searchSpaceSize)) - postfixLength

    startIndex = 0
    for i in range(postfixLength, postfixLength + prefixLength//2): startIndex += 2**i
    for i in range(startIndex, startIndex + numTargetsInCluster): answerSet.add(i) # cluster 3: 000..111xx..x (equal numbers of 0's and 1's exist in the prefix)    

    startIndex = 0
    for i in range(postfixLength + prefixLength//2, postfixLength + prefixLength): startIndex += 2**i
    for i in range(startIndex, startIndex + numTargetsInCluster): answerSet.add(i) # cluster 4: 111..000xx..x (equal numbers of 0's and 1's exist in the prefix)    
    
    return answerSet


def generateAnswersWithEvenDistribution(searchSpaceSize, numTargets):
    '''
    Generate solutions for Scenario-even
        randomly generate solutions, so that no explicit cluter exists
    '''
    assert math.log2(searchSpaceSize) == int(math.log2(searchSpaceSize))
    assert math.log2(numTargets) == int(math.log2(numTargets))

    answerSet = set() 
    for _ in range(numTargets):
        r = random.randint(0, searchSpaceSize - 1)
        while r in answerSet:
            r = random.randint(0, searchSpaceSize - 1)
        answerSet.add(r)

    return answerSet


def clusterSolutions(numQubits, numTargets, answersSetOriginal, numSimulations=1, debug=False):
    '''
    This function simulates the discovery of all solutions considering clustered solutions
        and it stores a limited number of clusters, while flushing out those with low priority
        'x' is considered a match
        when a new solution s is found, choose a cluster that s matches the most, among those that do not grow bigger than the max cluster size
        An entire cluster is searched for 
            (1) when the number of solutions found in the cluster is greater than or equal to a variable threhold, and
                (this threshold is a function of the cluster's size)
            (2) when size of the cluster > threshold, which we use sqrt(M)
                (this prevents searching a small cluster of size <= 4, as doing it does not improve performance significantly)
        
    Input:
        numQubits - the number of qubits that represents the search space (n)
        numTargets - the number of targets within the search space (M)
        answersSetOriginal - the set of M solutions
        numSimulations - the number of simulations to get numShots statistics        

    Output:
        Each element in the following lists represent one simulation result        
        numOracleQueriesList - # of oracle queries    
    '''
    # 
    # Initialize variables for the entire simulations
    #    
    searchSpaceSize = 2 ** numQubits
    numOracleQueriesList = []
    simulationCount = 0    
    numIterations = math.floor(math.sqrt(searchSpaceSize / numTargets) * math.pi / 4)
    successProbability = math.sin((2 * numIterations + 1) * math.sqrt(numTargets/searchSpaceSize)) ** 2    
    numTargetsSqrt = math.sqrt(numTargets)
    minMatchingBits = math.floor(numQubits/2)
    maxNumClusters = math.floor(math.log2(searchSpaceSize)) / 2

    while simulationCount < numSimulations:
        #
        # Initialize variables for one simulation run        
        #  
        answersSet = answersSetOriginal.copy()    
        answersSortedList = sorted(list(answersSet))
        if debug: print(f"{numTargets} answers = {answersSortedList}")        
        numMeasurements, numTotalIterations, numOracleQueries = 0, 0, 0
        answersSetFound = set()        
        numRemainingSolutions = numTargets    # the number of remaining solutions yet to be found        

        #
        # Initialize data structures used for clustering
        #
        bitAndCluster = [] # stores all patterns and the corresponding cluster IDs
        bitAndCluster.append([[] for _ in range(numQubits)]) # add lists for bit 0
        bitAndCluster.append([[] for _ in range(numQubits)]) # add lists for bit 1
        clusters = {} # hashmap that stores all clusters. key: cluster id, value: Cluster instance
        nextClusterID = 0        

        #
        # Find all of the remaining targets
        #                   
        while numRemainingSolutions > 0:                
            # Simulate the application and measurement on the quantum circuit            
            if random.random() < successProbability:
                # one of the solutions is measured
                measuredValueInt = random.choice(answersSortedList)
                measuredValueBinary = f"{measuredValueInt:0{numQubits}b}"                    
            else:                
                # one of the non-solutions is measured
                measuredValueInt = None
                measuredValueBinary = "one of non-solutions"                
            numMeasurements += 1
            numTotalIterations += numIterations
            numOracleQueries += numIterations        
            if debug: print(f"measured value: {measuredValueBinary}({measuredValueInt}) with success probability {successProbability:.3f} and {numIterations} Grover iterations")
                        
            if measuredValueInt is not None and measuredValueInt not in answersSetFound:
                # If the measured value is an undiscovered solution
                if debug: print(f"a new solution found and thus remove from answers")                
                answersSetFound.add(measuredValueInt)                
                numRemainingSolutions -= 1                
                numOracleQueries += 1

                #
                # Compare the newly-found answer with the cluster table
                #
                for _, cluster in clusters.items(): cluster.count = 0 # initialize all counts to 0
                for i in range(numQubits):
                    bit = int(measuredValueBinary[i])
                    relevantClusterList = bitAndCluster[bit][i]
                    for clusterID in relevantClusterList:
                        cluster = clusters.get(clusterID)
                        cluster.count += 1
                
                #
                # Find the cluster that matches the most bits with the newly found solution
                #   among the clusters that do not grow larger than the maximum allowed cluster size, when the new solution is included
                #
                clusterWithMaxCount = None
                for _, cluster in clusters.items():
                    if cluster.count - (len(cluster.pattern) - cluster.patternLength) < minMatchingBits: continue
                        # len(cluster.pattern) - cluster.patternLength: number of 'None' bits
                        # (cluster.count - (len(cluster.pattern) - cluster.patternLength)): number of matching bits that are not None    
                        # including the newly found sonlution into this cluster can grow this cluster larger than the maximum allowed cluster size, and thus continue
                    if clusterWithMaxCount == None or cluster.count > clusterWithMaxCount.count:
                        clusterWithMaxCount = cluster
                
                if clusterWithMaxCount != None:
                    #
                    # include the measured value into the cluster with a max count
                    #                        
                    if debug: print(f"{measuredValueBinary} matches cluster {clusterWithMaxCount}")
                    clusterWithMaxCount.numFoundSolutions += 1
                    clusterWithMaxCount.lastSolutionFoundTime = numMeasurements

                    # adjust the pattern according to the measured value
                    for i in range(numQubits):
                        bit = int(measuredValueBinary[i])
                        if clusterWithMaxCount.pattern[i] == bit or clusterWithMaxCount.pattern[i] == None: pass
                        else:
                            clusterWithMaxCount.pattern[i] = None
                            clusterWithMaxCount.patternLength -= 1
                            clusterWithMaxCount.size *= 2                            
                            bitAndCluster[bit][i].append(clusterWithMaxCount.id)
                    if debug: print(f"  this cluster is adjusted to {clusterWithMaxCount}")

                    if clusterWithMaxCount.numFoundSolutions >= math.log2(clusterWithMaxCount.size)\
                    and clusterWithMaxCount.size >= numTargetsSqrt:                                                               
                        # if (1) the # of discovered solutions >= log_2(cluster size)
                        # and (2) cluster size > threshold
                        #
                        # go through all potential solutions in the cluster and 
                        # see if they pass the oracle query
                        #
                        allElements = clusterWithMaxCount.generateAllelements()
                        if debug: print(f"  search all elements in the cluster: {allElements}")
                        for e in allElements:
                            if e not in answersSetFound:
                                numOracleQueries += 1
                                if e in answersSet:    
                                    if debug: print(f"      a new solution {e} found and thus remove from answers")                                    
                                    answersSetFound.add(e)                                
                                    numRemainingSolutions -= 1

                        #
                        # remove this cluster
                        #
                        for i in range(numQubits):
                            if clusterWithMaxCount.pattern[i] == None: 
                                bitAndCluster[0][i].remove(clusterWithMaxCount.id)
                                bitAndCluster[1][i].remove(clusterWithMaxCount.id)
                            else:
                                bitAndCluster[clusterWithMaxCount.pattern[i]][i].remove(clusterWithMaxCount.id)
                        clusters.pop(clusterWithMaxCount.id)

                else:
                    #
                    # store the measured value as a new cluster
                    #
                    clusters[nextClusterID] = Cluster(nextClusterID, measuredValueBinary, numMeasurements)
                    for i in range(numQubits):
                        bit = int(measuredValueBinary[i])
                        bitAndCluster[bit][i].append(nextClusterID)                    
                    if debug: print(f"{measuredValueBinary} does not match any cluster, so create a new cluster {clusters[nextClusterID]}")
                    nextClusterID += 1  

                    # 
                    # if the number of clusters exceed the specified maximum
                    #   then flush out the cluster with the minimum priority
                    #
                    if len(clusters) > maxNumClusters:
                        #
                        # select the cluster with the minimum priority to flush out
                        #                        
                        clusterWithMinPriority = None
                        for _, cluster in clusters.items():
                            if cluster.id == nextClusterID - 1: continue  # do not remove the newly-created cluster
                            if clusterWithMinPriority == None or cluster < clusterWithMinPriority:
                                clusterWithMinPriority = cluster
                        
                        if debug: print(f"remove cluster {clusterWithMinPriority}")

                        #
                        # remove the selected cluster
                        #
                        for i in range(numQubits):
                            if clusterWithMinPriority.pattern[i] == None:
                                bitAndCluster[0][i].remove(clusterWithMinPriority.id)
                                bitAndCluster[1][i].remove(clusterWithMinPriority.id)
                            else:
                                bitAndCluster[clusterWithMinPriority.pattern[i]][i].remove(clusterWithMinPriority.id)
                        clusters.pop(clusterWithMinPriority.id)

            else: # If the measured value is NOT an undiscovered solution
                if measuredValueInt not in answersSetFound: numOracleQueries += 1

        if debug:
            print(f"Algorithm terminated with {numTargets - len(answersSet)}/{numTargets} solutions found")
            print(f"    performed {numMeasurements} measurements, {numTotalIterations} Grover iterations, and {numOracleQueries} oracle queries")
            print()
        
        numOracleQueriesList.append(numOracleQueries)
        simulationCount += 1

    return numOracleQueriesList


class Cluster:
    def __init__(self, id, binaryNumber, solutionFoundTime=None):
        self.id = id
        self.pattern = []  # each element is 0, 1, or None (None means that the corresponding bit is not part of the pattern)
        self.patternLength = len(binaryNumber) # number of bits that are not "None" in the pattern
        for i in range(self.patternLength):
            bit = int(binaryNumber[i])
            if bit == 0: self.pattern.append(0)
            elif bit == 1: self.pattern.append(1)

        self.count = 0      # variable used to determine the most-matching cluster among many
        self.numFoundSolutions = 1      # number of discovered solutions that belong to this cluster
        self.size = 1   # total number of elements that this cluster represents, including both confirmed solutions and those not confirmed yet
        self.lastSolutionFoundTime = solutionFoundTime # store the value of numMeasurements when the last solution was added

    def __str__(self) -> str:
        result = []
        for e in self.pattern:
            if e == 0 or e == 1: result.append(str(e))
            else: result.append('x') # None means "do not care"
        return ''.join(result)
    
    def __lt__(self, other):
            if self.numFoundSolutions != other.numFoundSolutions: return self.numFoundSolutions < other.numFoundSolutions
            #elif self.size != other.size: return self.size < other.size
            else: return self.lastSolutionFoundTime < other.lastSolutionFoundTime  # these two values are always distinct

    def generateAllelements(self):
        # 
        # generate all elements that belong to this cluster as integers and
        # return them as a list
        #
        result = []
        xPositions = []  # list of "don't care" positions
        baseSum, baseMult = 0, 1
        for i in range(len(self.pattern)-1, -1, -1):
            if self.pattern[i] == 0: pass
            elif self.pattern[i] == 1: baseSum += baseMult
            else: xPositions.append([baseMult, 0])
            baseMult *= 2
        result.append(baseSum)
        for _ in range(self.size - 1):
            nextSum = baseSum
            carry = True
            for xPosition in xPositions:
                if carry:
                    if xPosition[1] == 0: carry = False
                    xPosition[1] ^= 1
                if xPosition[1] == 1: nextSum += xPosition[0]
            result.append(nextSum)
        return result


def printSimulationResults(message, results):
    def findMinAvgMax(elements):
        return min(elements), sum(elements)/len(elements), max(elements)
    
    numOracleQueriesList = results    
    print(f"{message} Algorithm performed {findMinAvgMax(numOracleQueriesList)} oracle queries")    


if __name__ == "__main__":
    '''
    Unit test for clusterSolutions
    '''
    for numQubits in range(6, 28, 2):
        searchSpaceSize = 2 ** numQubits
        numTargets = 2 ** int(numQubits / 2)    # max # of targets available for the given searchSpaceSize        
        numSimulations = 100

        print(f"search space size = {searchSpaceSize}:")
        answersSet = generateAnswersWithOneCluster(searchSpaceSize, numTargets)
        printSimulationResults(f"scenario-1:", clusterSolutions(numQubits, numTargets, answersSet, numSimulations=numSimulations, debug=False))
        answersSet = generateAnswersWithTwoEqualClusters(searchSpaceSize, numTargets)
        printSimulationResults(f"scenario-2:", clusterSolutions(numQubits, numTargets, answersSet, numSimulations=numSimulations, debug=False))
        answersSet = generateAnswersWithFourEqualClusters(searchSpaceSize, numTargets)
        printSimulationResults(f"scenario-4:", clusterSolutions(numQubits, numTargets, answersSet, numSimulations=numSimulations, debug=False))
        answersSet = generateAnswersWithEvenDistribution(searchSpaceSize, numTargets)
        printSimulationResults(f"scenario-even:", clusterSolutions(numQubits, numTargets, answersSet, numSimulations=numSimulations, debug=False))

        numOracleQueriesWithoutClustering = 0
        for numRemainingSolutions in range(numTargets, 0, -1):
            numOracleQueriesWithoutClustering += math.floor(math.sqrt(searchSpaceSize / numRemainingSolutions) * math.pi / 4)
            numOracleQueriesWithoutClustering += 1 # additional query required to confirm that the measured value is a solution
        print(f"previous work (clustering not used): it takes {numOracleQueriesWithoutClustering} oracle queries\n")

    