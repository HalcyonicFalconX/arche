import random
import sys
import copy
import timeit

import os
sys.path.append(os.path.abspath(".."))

from .lib.anneal import Annealer
from sortedcontainers import SortedList
from igraph import *

def colour(chromosome, l, total_len):           #BottleNeck

    k = len(l)  #Number of graphs
    clist = []  #Colour list
    new = SortedList()  #Temporary list to create 2d clist
    for i in range(k):
        for j in range(l[i]):
            new.add(j+1)
        clist.append(new)   #Creates list of colours for each graph
        new = SortedList()

    for i in range(total_len):

        arr = []        
        Max = 0
        for j in range(k):
            if chromosome[j][i] != 0:
                arr.append(j)
                if clist[j][0] > Max:
                    Max = clist[j][0]
        chromosome[k+1][i] = Max
        for z in arr:
            clist[z].discard(Max)

    return chromosome

class Movement:

    def __init__(self, d, t):
        self.divider = d
        self.total_len = t

    def mutate2(self, args):

        """
        Chooses a random row (graph) and inserts a random number of 0's (blank spaces)
        at a random position, while preserving solution integrity.
        """

        inp = args[0]
        l = args[1]
        chromosome = inp[0]
        k = len(l)
        c = random.randint(0, k-1)

        newlist = []
        for i in range(self.divider, self.total_len):
            if chromosome[c][i] != 0:
                newlist.append(i)

        l2 = len(newlist)
        idx = newlist[random.randint(0, l2-1)]
        i = idx
        rightzeros = []
        rightpos = []
        zero = 0
        while (i < self.total_len):
            if chromosome[c][i] == 0:
                if zero == 0:
                    rightpos.append(i)
                zero += 1
            elif (chromosome[c][i] != 0 or i == total_len-1) and zero > 0:
                rightzeros.append(zero)
                zero = 0
            if i == self.total_len-1 and zero>0:
                rightzeros.append(zero)
            i+=1

        i = idx
        leftzeros = []
        leftpos = []
        zero = 0
        while (i >= self.divider):
            if chromosome[c][i] == 0:
                if zero == 0:
                    leftpos.append(i)
                zero += 1
            elif (chromosome[c][i] != 0 or i == 0) and zero > 0:
                leftzeros.append(zero)
                zero = 0
            if i == self.divider and zero>0:
                leftzeros.append(zero)
            i-=1
    
        zeros = leftzeros + rightzeros
        pos = leftpos + rightpos
        mutationpt = random.randint(0, len(pos)-1)
        mutationamt = random.randint(1, zeros[mutationpt])

        if pos[mutationpt] > idx:
            for i in range(pos[mutationpt]+mutationamt-1, idx+mutationamt-1, -1):
                chromosome[c][i] = chromosome[c][i-mutationamt]
            for i in range(idx, idx+mutationamt):
                chromosome[c][i] = 0
        else:
            for i in range(pos[mutationpt]-mutationamt+1, idx-mutationamt+1):
                chromosome[c][i] = chromosome[c][i+mutationamt]
            for i in range(idx-mutationamt+1, idx+1):
                chromosome[c][i] = 0

        args[0][0] = chromosome
        return args

    def crossover(self, args):

        inp = args[0]
        l = args[1]
        g = args[2]

        chromosome = inp[0]

        k = len(l)
        c = random.randint(0, k-1)
        total_len = 0
        for i in l:
            self.total_len += int(i)

        idxset = set(range(self.divider+1, self.total_len))
        idx = random.choice(tuple(idxset))
        while chromosome[c][idx] != 0 and chromosome[c][idx-1] != 0 and g[c].are_connected(chromosome[c][idx-1]-1, chromosome[c][idx]-1):
            idxset.remove(idx)
            idx = random.choice(tuple(idxset))

        temp = chromosome[c][idx]
        chromosome[c][idx] = chromosome[c][idx-1]
        chromosome[c][idx-1] = temp

        args[0][0] = chromosome
        return args

    def crossover2(self, args):

        """
        Chooses a random row (graph) and a random valid node (vertex) in that row.
        Moves the chosen node a random number of places to the left.
        Rightshifts all nodes between initial and final positions of the chosen node.
        """

        inp = args[0]
        l = args[1]
        g = args[2]
        chromosome = inp[0]
        k = len(l)
        c = random.randint(0, k-1)

        newlist = []
        for i in range(self.divider+1, self.total_len):
            if chromosome[c][i] != 0:
                newlist.append(i)

        idx = newlist[random.randint(0, len(newlist)-1)]
        leftmost = idx
        for i in range(idx-1, self.divider-1, -1):
            if chromosome[c][i] == 0:
                leftmost = i
            elif not g[c].are_connected(int(chromosome[c][i]-1), int(chromosome[c][idx]-1)):
                leftmost = i
            else:
                break

        if leftmost == idx:
            args[0][0] = chromosome
            return args

        idx2 = random.randint(leftmost, idx-1)  
        temp = chromosome[c][idx]
        for i in range(idx, idx2, -1):
            chromosome[c][i] = chromosome[c][i-1]
        chromosome[c][idx2] = temp

        args[0][0] = chromosome
        return args

    def condense(self, args):

        """
        Chooses a random row (graph) and two random indices.
        Condenses all nodes between the indices and places them at a random starting index.
        Preserves mutual ordering.
        """

        inp = args[0]
        l = args[1]
        g = args[2]
        chromosome = inp[0]
        k = len(l)
        c = random.randint(0, k-1)

        choices = set(range(self.divider, self.total_len))
        leftpos = random.choice(tuple(choices))
        choices.remove(leftpos)
        rightpos = random.choice(tuple(choices))

        if leftpos > rightpos:
            leftpos, rightpos = rightpos, leftpos

        templist = []
        for i in range(leftpos, rightpos+1):
            if chromosome[c][i] != 0:
                templist.append(int(chromosome[c][i]))

        if len(templist) != rightpos-leftpos+1:
            idx = random.randint(leftpos, rightpos-len(templist)+1)
            for i in range(leftpos, idx):
                chromosome[c][i] = 0
            for num in templist:
                chromosome[c][idx] = num
                idx += 1
            for i in range(idx, rightpos+1):
                chromosome[c][i] = 0

        args[0][0] = chromosome
        return args

    def expand(self, args):

        """
        Chooses a random row (graph) and two random indices.
        Distributes all nodes between the indices randomly within the limits of the said indices.
        Preserves mutual ordering.
        """
        inp = args[0]
        l = args[1]
        g = args[2]
        chromosome = inp[0]
        k = len(l)
        c = random.randint(0, k-1)

        choices = set(range(self.divider, self.total_len))
        leftpos = random.choice(tuple(choices))
        choices.remove(leftpos)
        rightpos = random.choice(tuple(choices))

        if leftpos > rightpos:
            leftpos, rightpos = rightpos, leftpos
        
        templist = []
        universe = []
        for i in range(leftpos, rightpos+1):
            universe.append(int(i))
            if chromosome[c][i] != 0:
                templist.append(int(chromosome[c][i]))

        if len(templist) != rightpos-leftpos+1:
            positions = random.sample(universe, len(templist))
            positions.sort()
            universe = list(set(universe) - set(positions))
            for position, item in zip(positions, templist):
                chromosome[c][position] = item
            for position in universe:
                chromosome[c][position] = 0

        args[0][0] = chromosome
        return args

    def shufflecolor(self, args):

        inp = args[0]
        l = args[1]
        g = args[2]
        chromosome = inp[0]
        k = len(l)

        #clist = list(range(1, max(l)+1))
        choicelist = []
        for i in range(sum(l)):
            if chromosome[k+1][i] != 0:
                choicelist.append(i)

        idx1 = random.choice(choicelist)

        graphs = set([])
        for i in range(k):
            if chromosome[i][idx1] != 0:
                graphs.add(i)

        choicelist = []
        for i in range(sum(l)):
            tempset = set([])
            for j in range(k):
                if chromosome[j][i] != 0:
                    tempset.add(j)
            if tempset == graphs and i != idx1:
                choicelist.append(i)

        if len(choicelist) > 0:
            idx2 = random.choice(choicelist)
            t = chromosome[k+1][idx1]
            chromosome[k+1][idx1] = chromosome[k+1][idx2]
            chromosome[k+1][idx2] = t
#

        args[0][0] = chromosome
        return args

    def magiccolor(self, args):

        inp = args[0]
        l = args[1]
        g = args[2]
        chromosome = inp[0]
        k = len(l)

        choicelist = []
        for i in range(sum(l)):
            if chromosome[k+1][i] != 0:
                choicelist.append(i)

        idx = random.choice(choicelist)

        graphs = []
        for i in range(k):
            if chromosome[i][idx] != 0:
                graphs.append(i)

        clist = set(range(1, max(l)+1))
        usedclist = set([])
        for i in range(sum(l)):
            include = 0
            for graph in graphs:
                if chromosome[graph][i] != 0:
                    include = 1
                    break
            if include == 1:
                usedclist.add(chromosome[k+1][i])

        choicelist = set(clist - usedclist)
        if len(choicelist) > 0:
            newcolor = random.sample(choicelist, 1)[0]
            chromosome[k+1][idx] = newcolor

        args[0][0] = chromosome
        return args

    def hardshuffle(self, args):

        inp = args[0]
        l = args[1]
        g = args[2]
        chromosome = inp[0]
        k = len(l)

        #clist = list(range(1, max(l)+1))
        choicelist = []
        for i in range(sum(l)):
            if chromosome[k+1][i] != 0:
                choicelist.append(i)

        idx1 = random.choice(choicelist)
        c1 = chromosome[k+1][idx1]
        choicelist.remove(idx1)
        
        if len(choicelist) == 0:
            sys.exit()

        idx2 = random.choice(choicelist)
        c2 = chromosome[k+1][idx2]
        graph1 = []
        graph2 = []

        for i in range(k):
            if chromosome[i][idx1] != 0:
                graph1.append(i)
            if chromosome[i][idx2] != 0:
                graph2.append(i)

        for i in range(sum(l)):
            if chromosome[k+1][i] == c1:
                tempset = set([])
    
        args[0][0] = chromosome
        return args

class Solution(Annealer, Movement):

    def __init__(self, x):
        self.state = x
        self.breakeven = -1
        self.total_len = self.state[3]
        self.Mover = Movement(0, self.total_len)
        self.CXPB = 0.5
        self.MUTPB = 0.5
        self.EXPB = 0.4
        self.CDPB = 0.4

    def assign(self, chromosome, l):
        assign = 1  #Issues assignments based on chromosome structure
        k = len(l)
        for i in range(self.total_len):
            Bool = 0
            for j in range(k):
                if chromosome[j][i] != 0:
                    Bool = 1
                    break
            if Bool == 1:
                chromosome[k][i] = assign
                assign += 1
            else:
                chromosome[k][i] = 0
        return chromosome

    def colour(self, chromosome, l):            #BottleNeck

        k = len(l)  #Number of graphs
        clist = []  #Colour list
        new = SortedList()  #Temporary list to create 2d clist
        for i in range(k):
            for j in range(l[i]):
                new.add(j+1)
            clist.append(new)   #Creates list of colours for each graph
            new = SortedList()

        for i in range(self.total_len):

            arr = []        
            Max = 0
            for j in range(k):
                if chromosome[j][i] != 0:
                    arr.append(j)
                    if clist[j][0] > Max:
                        Max = clist[j][0]
            chromosome[k+1][i] = Max
            for z in arr:
                clist[z].discard(Max)

        return chromosome

    def move(self):

        l = self.state[1]

        choice = random.randint(1,5)
        if choice == 1:
            self.state = self.Mover.crossover2(self.state)
            chromosome = self.state[0][0]
            chromosome = self.colour(chromosome, l)
            self.state[0][0] = chromosome
        if choice == 2:
            self.state = self.Mover.mutate2(self.state)
            chromosome = self.state[0][0]
            chromosome = self.colour(chromosome, l)
            self.state[0][0] = chromosome
        if choice == 3:
            self.state = self.Mover.condense(self.state)
            chromosome = self.state[0][0]
            chromosome = self.colour(chromosome, l)
            self.state[0][0] = chromosome
        if choice == 4:
            self.state = self.Mover.expand(self.state)
            chromosome = self.state[0][0]
            chromosome = self.colour(chromosome, l)
            self.state[0][0] = chromosome

        energymin = self.energy()
        beststate = self.state

        k = random.randint(0, self.total_len)
        for i in range(k):
            if random.random() > 0.5:
                self.state = self.Mover.magiccolor(self.state)
            else:
                self.state = self.Mover.shufflecolor(self.state)
            energynew = self.energy()
            if energynew < energymin:
                energymin = energynew
                beststate = self.state

        self.state = beststate

    def energy(self):

        inp = self.state[0]
        l = self.state[1]
        g = self.state[2]
        chromosome = inp[0]

        #chromosome = self.colour(chromosome, l)
        chromosome = self.assign(chromosome, l)

        k = len(l)
        c = random.randint(0, k-1)

        x = 0
        color = []
        for i in range(k): 
            newlist = [0]*l[i]
            for j in range(self.total_len):
                if chromosome[i][j] != 0:
                    newlist[chromosome[i][j]-1] = chromosome[k+1][j]
            color.append(newlist)
            newlist = []

        for i in range(self.total_len):
            for j in range(k):
                if chromosome[j][i] != 0:
                    for z in range(j+1, k):
                        if chromosome[z][i] != 0 and g[j].indegree(chromosome[j][i]-1) != g[z].indegree(chromosome[z][i]-1):
                            x += 1
                        elif chromosome[z][i] != 0: 
                            list1 = g[j].predecessors(chromosome[j][i]-1)
                            list2 = g[z].predecessors(chromosome[z][i]-1)
                            clist1 = []
                            clist2 = []
                            for p in range(len(list1)):
                                clist1.append(color[j][list1[p]])
                            for p in range(len(list2)):
                                clist2.append(color[z][list2[p]])
                            if set(clist1) != set(clist2):
                                x+=1
        y = self.total_len - max(chromosome[k])
        C0 = self.total_len - max(l)
        if self.breakeven == -1:
            self.breakeven = C0

        ret = (y - C0*x)/1.0
        #if ret < 0:
        #   ret /= C0

        return ret*-1.0+C0

class Create:

    def __init__(self, x, y):
        self.createcycle = 0 
        self.toporesort = x
        self.divider = 0
        self.total_len = y

    def DFS(graph, visited, current_vertex, t1):

        visited[current_vertex] = 1
        seq = graph.vs.select(lambda vertex: graph.are_connected(current_vertex, vertex.index) and visited[vertex.index] == 0)
        indices = SortedList()
        for vs in seq:
            indices.add(vs.index)
        while len(indices) > 0:
            next = random.choice(indices)
            if visited[next] == 1:
                indices.discard(next)
                continue
            indices.discard(next)
            t1, visited = DFS(graph, visited, next, t1)

        t1.append(current_vertex)

        return t1, visited
    

    def topological_sort(g, l):

        t = []
        i = 0

        for graph in g:
            t1 = []
            visited = [0 for j in range(l[i])]
            seq = graph.vs.select(lambda vertex: graph.indegree(vertex) == 0 and visited[vertex.index] == 0)
            indices = SortedList()
            for vs in seq:
                indices.add(vs.index)
            while len(indices) > 0:
                start = random.choice(indices)
                if visited[start] == 1:
                    indices.discard(start)
                    continue
                indices.discard(start)
                visited[start] = 1
                t1, visited = DFS(graph, visited, start, t1)
            t1.reverse()
            t.append(t1)
            del t1
            visited = []
            i += 1

        return t

    def createseed3(self, T, l, g):

        if self.createcycle >= self.toporesort:
            temp = topological_sort(g, l)
            T[:] = temp
            self.createcycle = 0
        else:
            self.createcycle += 1

        t = copy.deepcopy(T)

        seed1 = []
        for i in range(len(l)+2):
            temp = [0]*sum(l)
            seed1.append(temp)

        indeg0 = 0

        for i in range(len(l)):
            k = 0
            for j in range(l[i]):
                if g[i].indegree(t[i][j]) == 0:
                    seed1[i][k] = t[i][j]+1
                    t[i][j] = -1
                    indeg0 += 1
                    k+=1

        for i in range(len(l)):
            j = 0
            while j < len(t[i]):
                if t[i][j] == -1:
                    t[i].pop(j)
                else:
                    j += 1

        if self.divider == 0:
            self.divider = indeg0

        include = set(range(indeg0, self.total_len))

        for i in range(len(l)):
            j = 0
            positions = random.sample(include, len(t[i]))
            include = include - set(positions)
            positions.sort()
            for j in range(len(t[i])):
                seed1[i][positions[j]] = t[i][j]+1

        #devicelist = list(range(1, max(l)+1))

        #for i in range(len(l)+2, 2*len(l)+2):
        #   random.shuffle(devicelist)
        #   for j in range(max(l)):
        #       seed1[i][j].add(devicelist[j])

        return seed1;

class MIMDSA:

    def __init__(self,graphs,graphcount,outname):
        
        self.__graphs = graphs 
        self.__graphCount = graphcount
        self.__outname = outname

        self.__g = []
        self.__t = []
        self.__l = []
        self.__Tl = 0

    def readGraph(self):

        for i in range(self.__graphCount):
            temp = Graph(directed = True)
            self.__g.append(temp)

        for i in range(self.__graphCount):
            f = open(self.__graphs[i], "r")
            f1 = f.readlines()
            f.close()

            numbers = []
            temp2 = []

            for elem in f1:
                temp = list(map(int, elem.split()))
                numbers.append(temp)

            self.__g[i].add_vertices(max(list(map(max, numbers)))+1)
            for j in range(max(list(map(max, numbers)))+1):
                temp2.append(j)

            temp = []
            for elem in numbers:
                temp.append(tuple(elem))

            self.__g[i].add_edges(temp)
            t1 = self.__g[i].topological_sorting(mode=OUT)
            self.__t.append(t1)
            len1 = len(t1)
            self.__l.append(len1)

        for graph in self.__l:
            self.__Tl += graph

    def genSolution(self):

        creator = Create(1000, self.__Tl)
        chromosome = creator.createseed3(self.__t, self.__l, self.__g)
        chromosome = colour(chromosome, self.__l, self.__Tl)

        state = [[chromosome], self.__l, self.__g, self.__Tl]

        sol = Solution(state)
        sol.Tmax = 5.0
        sol.Tmin = 0.01
        sol.steps = 50000
        sol.updates = 100
        sol.Timeout = 21600.0
        itinerary, fitness = sol.anneal()

        apex = itinerary[0][0]
        C0 = sol.breakeven
        fitness = (fitness-C0)*-1

        print ('\n')
        print ('Final Fitness = ',fitness)
        #print apex

        f = open(self.__outname, "w+")

        for i in range(self.__Tl):
            if apex[len(self.__l)][i] != 0:
                f.write("%d c%d " % (apex[len(self.__l)][i], apex[len(self.__l)+1][i]))
                for j in range(len(self.__l)):
                    if apex[j][i] != 0:
                        f.write("%d " % int(apex[j][i]-1))
                    else:
                        f.write("- ")
                f.write("\n")
        f.close()

        print("Solution: ", end='\n')

        for i in range(self.__Tl):
            if apex[len(self.__l)][i] != 0:
                print(str(apex[len(self.__l)][i])+" c"+str(apex[len(self.__l)+1][i]), end=' ')
                for j in range(len(self.__l)):
                    if apex[j][i] != 0:
                        print(str(apex[j][i]-1), end=' ')
                    else:
                        print("- ", end='')
                print("\n", end='')
    

if __name__ == '__main__':

    if len(sys.argv) < 3:
        print("Usage: python mimdSA.py <File name of graph1> <File name of graph2> ... (at least 2 files)\n")
        sys.exit()

    techMapper = MIMDSA(sys.argv[1:len(sys.argv)-1], len(sys.argv)-1, sys.argv[len(sys.argv)-1])
    techMapper.readGraph()
    techMapper.genSolution()

