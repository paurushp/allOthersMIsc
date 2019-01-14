# Python Codes

### text pattern search

import re 

def matcher(word, pat, what):
    if what == 'first':
        match = re.search(pat, word)
        if match:
            print match.group()
        else:
            print 'Not found'
    if what == 'all':
        matches = re.findall(pat, word)
        if matches:
	    for match in matches:
                print match
        else:
            print 'Not found'      

  
word= '''The Agilent Biopharma Applications Compendium contains dynamic 
       links to more than 125 biopharma applications for discovery, 
       development and QA/QC. The guide covers workflows for intact 
       protein analysis, aggregation analysis, charge variant analysis, 
       peptide mapping and glycosylation, using LC, LC/MS, CE/MS, and on 10/11/2012'''

pat_repeat_word= ''
pat_date='\d\d.\d\d.\d\d\d\d'


from itertools import permutations
perms = [''.join(p) for p in permutations('stack')]
perms

matcher(word, '\d\d.\d\d.\d\d\d\d', what='all') # 
matcher(word, 'pharma', what='all')
### LCS algorithm
def LCS(string_one, string_two, x, y): 
    alignment_matrix = [[0 for j in range(y+1)] for i in range(x+1)] # Computing DP matrix
    for i in range(1, x+1):
        for j in range(1, y+1):
            if string_one[i-1] == string_two[j-1]:
                alignment_matrix[i][j] = alignment_matrix[i-1][j-1] + 1
            else:
	        alignment_matrix[i][j] = max(alignment_matrix[i][j-1], alignment_matrix[i-1][j])
    return alignment_matrix


def back_tracker(string_one, string_two, alignment_matrix, x, y): # back-track function
    if x == 0 or y == 0: # Empty string inputs
        return set([""])
    elif string_one[x-1] == string_two[y-1]:
        return set([sequence + string_one[x-1] for sequence in back_tracker(string_one, string_two, alignment_matrix, x-1, y-1)])
    else:
        lcs_seq = set()
        if alignment_matrix[x][y-1] >= alignment_matrix[x-1][y]:
	    lcs_seq.update(back_tracker(string_one, string_two, alignment_matrix, x, y-1)) # recursive calling
	if alignment_matrix[x-1][y] >= alignment_matrix[x][y-1]:
	    lcs_seq.update(back_tracker(string_one, string_two, alignment_matrix, x-1, y)) # recursive calling
	return lcs_seq

def find_lcs(string_one, string_two):
    m = len(string_one)
    n = len(string_two)
    dp_mat = LCS(string_one, string_two, m, n)
    seq_len = max(max(dp_mat))
    lcs_seq = back_tracker(string_one, string_two, dp_mat, m, n)
    return [seq_len, lcs_seq]
    


s1 = 'AATCC'
s2 = 'ACACG'
my_lcs = find_lcs(s1, s2)
my_lcs

### Read a FASTA FILE
import re

def read_fasta(file_name):
    read_file = open(file_name)
    sequence = read_file.read()
    if sequence[0] == '>':
        anotation = re.findall(r'>(.*?)\n', sequence)
        sequence = a = ''.join(sequence.splitlines(True)[1:])
        sequence = sequence.upper()
        sequence = sequence.rstrip("\n")
        sequence = sequence.replace('\n','')
        print ("sequence length: %d" % len(sequence))
        return sequence, anotation
    else:
        sequence = sequence.upper()
        sequence = sequence.rstrip("\n")
        sequence = sequence.replace('\n','')
        print ("sequence length: %d" % len(sequence))
        return sequence
        
    
    
read_fasta('/home/praveen/seq1.txt')



### Needleman Wumsch

import numpy as np
arr = np.array(aln)

def compute_matrix(string_one, string_two, x, y, match, mis_match, gap): 
    alignment_matrix = [[0 for j in range(y+1)] for i in range(x+1)] # Computing DP matrix alignment_matrix = np.zeros((x+1, y+1))
    for i in range(1, x+1):
        for j in range(1, y+1):
            if string_one[i-1] == string_two[j-1]:
                alignment_matrix[i][j] = alignment_matrix[i-1][j-1] + match
            else:
	        alignment_matrix[i][j] = max(alignment_matrix[i][j-1]-gap, alignment_matrix[i-1][j-1]-mis_match, alignment_matrix[i-1][j]-gap)
    return alignment_matrix


def compute_back_tracker(string_one, string_two, alignment_matrix, x, y, gap): # back-track function
    S = ''
    T = ''
    U = ''
    i = x
    j = y
    while i > 0 or j > 0 :
        if alignment_matrix[i][j] == alignment_matrix[i-1][j] - gap:
            S = string_one[i-1] + S  # Please remeber we craeted matrix of sixe x+1 times y+1 to include gaps
            T = '-' + T
            U = ' ' + U
            i = i-1
        elif alignment_matrix[i][j] == alignment_matrix[i][j-1] - gap:
            S = '-' + S 
            T = string_two[j-1] + T
            U = ' ' + U
            j = j-1
        else:
	    S = string_one[i-1] + S 
	    T = string_two[j-1] + T
	    U = '|' + U
            i = i-1
            j = j-1
    print S
    print U
    print T
    return [S, T]


 
def align_sequence(string_one, string_two, match, mis_match, gap):
    m = len(string_one)
    n = len(string_two)
    aln = compute_matrix(s1, s2, m, n, match, mis_match, gap)
    ss = compute_back_tracker(s1,s2, aln, m, n, gap)
    score = max(map(max,aln))
    return [seq_len, lcs_seq]

    
    
compute_back_tracker(s1,s2, aln, 8, 8,1)
'''
Complexity: Matrix filling O(m+n) + O(mn)
            Traceback O(m+n) + O(n) :: n > m
            Overall O(mn)
            Same for smith waterman: matrix filling same, but to find max in matrix O(mn)
Possible improvements:
Use of numpy array in place of list of list
saves on space
Another- Use of tuples words of size k ## Four Russian speed up
Another- Inclusion of gap extension penalty

Doing local alignment
Start with i,j = max(alignment_matrix)
and trace back (Smith- Waterman)

To use for protein sequence the gap, match and mismatches has to be called via by a matrix called substitution matrix
'''
def find_max_index(mat):
    arr = np.array(mat)
    temp = 0
    for i in range(0, arr.shape[0]):
        for j in range(0, arr.shape[1]):
	    if arr[i,j] > temp:
	        temp = arr[i,j]
	        index_i = i
	        index_j = j
    return ([i,j])
    
    

    
    
### ORF Finder 

start_codons = ['ATG', 'ACG'] 
stop_codons = ['TAA', 'TGT']

def find_ORF(sequence, treshold, start_codons, stop_codons):
    start_codon_index = 0
    end_codon_index = 0
    start_codon_found = False
    orfs = []
    for j in range(0, 3):
        for indx in range(j, len(sequence), 3):
            current_codon = sequence[indx:indx+3]
            if current_codon in start_codons and not start_codon_found:
                start_codon_found = True
                start_codon_index = indx
            if current_codon in stop_codons and start_codon_found:
                end_codon_index = indx
                print 'ORF is '
                print seq[start_codon_index:end_codon_index]
                length = end_codon_index - start_codon_index + 1
                if length >= treshold * 3:
                    orfs.append(start_codon_index)
                start_codon_found = False
        start_codon_index = 0
        end_codon_index = 0
        start_codon_found = False
    return len(orfs), orfs

    
   
seq= 'ATGCTAGTCGATGCTACTGTGCTGATCGTAGCTATGTAA'
start_codons = 'ATG'
stop_codons = 'TAA'
find_ORF(seq, 3,start_codons, stop_codons) 


#  Complexity: O(n)
# Improvement: Search for multiple codons (1) Iteration (2) searching all the codons

## ## Coin change dynamic programing

coins = [1,2,5,10,20,50,100,200]
change = 243

# Recursive approach
def recursive_change(coins,change):
    minCoins = change
    if change in coins:
        return 1
    else:
        for i in [c for c in coins if c <= change]:
            numCoins = 1 + recursive_change(coins,change-i)
            if numCoins < minCoins:
                minCoins = numCoins
    return minCoins

print(recursive_change([1,5,10,25],63))
      
# Dynamic programing approach

def dp_change(coins, change):
    result = []
    for i in range(0, change+1):
        result.append(999)      
    result[0]=0
    coin_list = []
    for amount in range(1,change+1):
        for j in range(0,len(coins)):
            ref = result[amount - coins[j]] + 1
            if coins[j] <= amount and ref < result[amount]:
                tempA = result[amount - coins[j]]
                tempB = tempA+1
                result[amount] = tempB
    return result[change]

value = dp_change([1,2,5, 10, 20, 50], 243)


import timeit

start = timeit.default_timer()
print(recursive_change([1,5,10,25],63))
stop = timeit.default_timer()
t1 = stop - start 

start = timeit.default_timer()
print(dp_change([1,5,10,25],63))
stop = timeit.default_timer()
t2 = stop - start 



### DFS 



class Graph(object):
    def __init__(self, graph_dict={}):
        """ initializes a graph object """
        self.__graph_dict = graph_dict
    def vertices(self):
        """ returns the vertices of a graph """
        return list(self.__graph_dict.keys())
    def edges(self):
        """ returns the edges of a graph """
        return self.__generate_edges()
    def add_vertex(self, vertex):
        """ If the vertex "vertex" is not in 
            self.__graph_dict, a key "vertex" with an empty
            list as a value is added to the dictionary. 
            Otherwise nothing has to be done. 
        """
        if vertex not in self.__graph_dict:
            self.__graph_dict[vertex] = []
    def add_edge(self, edge):
        """ assumes that edge is of type set, tuple or list; 
            between two vertices can be multiple edges! 
        """
        edge = set(edge)
        (vertex1, vertex2) = tuple(edge)
        if vertex1 in self.__graph_dict:
            self.__graph_dict[vertex1].append(vertex2)
        else:
            self.__graph_dict[vertex1] = [vertex2]
    def __generate_edges(self):
        """ A static method generating the edges of the 
            graph "graph". Edges are represented as sets 
            with one (a loop back to the vertex) or two 
            vertices 
        """
        edges = []
        for vertex in self.__graph_dict:
            for neighbour in self.__graph_dict[vertex]:
                if {neighbour, vertex} not in edges:
                    edges.append({vertex, neighbour})
        return edges
    def __str__(self):
        res = "vertices: "
        for k in self.__graph_dict:
            res += str(k) + " "
        res += "\nedges: "
        for edge in self.__generate_edges():
            res += str(edge) + " "
        return res


graph = {'A' : set(['B','C']), 'B' : set(['A','D', 'E']), 'C' : set(['A','F']), 'D' : set(['B']), 'E' : set(['B','F']), 'F' : set(['C','E'])}
graph = Graph(graph)
print(graph.vertices())
print(graph.edges())


def find_path(self, start_vertex, end_vertex, path=[]):
        """ find a path from start_vertex to end_vertex 
            in graph """
        graph = self.__graph_dict
        path = path + [start_vertex]
        if start_vertex == end_vertex:
            return path
        if start_vertex not in graph:
            return None
        for vertex in graph[start_vertex]:
            if vertex not in path:
                extended_path = self.find_path(vertex, end_vertex, path)
                if extended_path: 
                    return extended_path
        return None



def generate_edges(graph):
    edges = []
    for node in graph:
        for neighbour in graph[node]:
            edges.append((node, neighbour))
    return edges
print(generate_edges(graph))



## BFS algorithm
## Graph as a dictionary

def back_track(parent, start, end):
    path = [end]
    while path[-1] != start:
        path.append(parent[path[-1]])
    path.reverse()
    return path
    
def breadth_first_search(graph, root, end):
    parent = {}
    queue = []
    queue.append(root)
    while queue:
        node = queue.pop(0)
        if node == end:
            return back_track(parent, root, end)
        for adjacent in graph.get(node, []):
            parent[adjacent] = node # <<<<< record its parent 
            queue.append(adjacent)

            
            
depth_first_search(graph, 'A')



### Sorting algorithms

## Bubble sort Complexity O(n^2) , O(1)
## Swaps adjacent elements till sorted
## 

def b_sort(n):
    for i in range(len(n)):
    for j in range(len(n) - 1, i, -1):
        if (n[j] < n[j - 1]):
        swap(n, j, j - 1 )
        print n
 
def swap(n, x, y):
    tmp = n[x]
    n[x] = n[y]
    n[y] = tmp

## Insertation sort Complexity time complexity O(n^2) , space complexity O(1)
## inserts a number to right place
## 

def i_sort(n):
    for i in range(len(n)):
        j = i
        while j > 0 and n[j-1] > n[j]:
	    n[j-1], n[j] = n[j], n[j-1]
	    print n
	    j -= 1
	    

## Quick sort Complexity time complexity O(nlog(n)) , space complexity O(log(n))
## Divide and conquer 
## 

def q_sort(n):
    if len(n) <= 1:
        return n
    else:
        sorted_n = q_sort([x for x in n[1:] if x<n[0]]) + [n[0]] + q_sort([x for x in n[1:] if x>=n[0]])
        print sorted_n
        return sorted_n

## Merge sort Complexity time complexity O(nlog(n)) , space complexity O(log(n))
## Divide and conquer select mid point recursively and sort 
## 
   
   
def m_sort(n):
    if len(n) <= 1:
        return n
    sorted_n = []
    split_point = int(len(n)/2)
    x = m_sort(n[:split_point])
    y = m_sort(n[split_point:])
    print x, y
    while (len(x)>0) and (len(y) > 0):
        if x[0] > y[0]:
	    sorted_n.append(y.pop(0))
	else:
	    sorted_n.append(x.pop(0))
    print sorted_n
    sorted_n.extend(x+y)
    return sorted_n
    

   m_sort(n)     
## Search Algorithms

# Binary search
# Split array into half and search each half
# For sorted array
# An iterative approach is also possible
# Divide and conquer approach
# time O(log(n)) space O(1)
def binary_search(l, value):
    if len(l) == 0:
        return -1
    elif len(l) == 1:
        if l[0] == value:
	    return value
	else:
	    return -1
    else:
        mid = len(l)//2
        if l[mid] == value:
	    return value
	elif l[mid] > value:
	    return binary_search(l[0:mid], value)
	elif l[mid] < value:
	    return binary_search(l[mid:len(l)], value)
    
binary_search(l, 8)
# To return the index as well :) create an initial set of positions and remove half every time you dont find :)
def binary_search_index(l, value, res):
    if len(l) == 0:
        return -1, -1
    elif len(l) == 1:
        if l[0] == value:
	    return value, res
	else:
	    return -1, -1
    else:
        mid = len(l)//2
        if l[mid] == value:
	    res = res[mid]
	    return value, res
	elif l[mid] > value:
	    res = res[0:mid]
	    return binary_search_index(l[0:mid], value, res)
	elif l[mid] < value:
	    res = res[mid:len(l)]
	    return binary_search_index(l[mid:len(l)], value, res)

        
binary_search_index(l,102, range(0,len(l)-1))

# Iterative search

def iterative_search(l, value):
    i = 0
    for values in l:
        #print values
        i += 1
        if values == value:
	    return value, i
	    break

	    
start = timeit.default_timer()
print(iterative_search(l,100))
stop = timeit.default_timer()
t1 = stop - start 

start = timeit.default_timer()
print(binary_search_index(l,100, range(0,len(l))))
stop = timeit.default_timer()
t2 = stop - start 

### Check Palindrome number
##
##
def palindrome(n):
    if str(n) == str(n)[::-1]:
        print '%d is a palindrome' %n

### Check Palindrome sequence
## 
##

def palindrome_checker(seq):
    seq = list(seq)
    rev_seq = reversed(seq)
    if list(seq) == list(rev_seq):
        print("It is palindrome")
    else:
        print("It is not palindrome")


### MCMC

## LDA

## Decision tree



## Implementing data structures in Python

# Network alteration and MCMC
random.seed(123)

def find_possible_actions(mat, x, y, time):
    actions = [0,0,0,0] # remove, increase, decrease, swap (direction change)
    if x != y: # avoiding self loops
        if mat[x,y] == 0: # !remove, increase, !decrease, swap (direction change)
	    actions[1] = 1
	    if mat[x,y] != mat[y,x]:
	        actions[3] = 1
	if mat[x,y] > 0 and mat[x,y] < time: # remove, increase, decrease, swap (direction change)
	    actions[0] = 1
	    actions[1] = 1
	    actions[2] = 1
	    if mat[x,y] != mat[y,x]:
	        actions[3] = 1
	if mat[x,y] == time: # remove, !increase, decrease, swap (direction change)
	    actions[0] = 1
	    actions[2] = 1
	    if mat[x,y] != mat[y,x]:
	        actions[3] = 1
    return actions


def alter_net(mat, time):
    shape_mat = mat.shape
    x = random.sample(range(0,shape_mat[0]), 1) # samples an edge
    y = random.sample(range(0,shape_mat[1]), 1) # samples an edge
    operations = find_possible_actions(mat, x, y, time) # finds all possible operation on the edge
    op = random.sample(range(0,3),1)[0] # random sampling of operations
    new_mat = mat
    if operations[op] == 1 and op == 0:
        new_mat[x,y] = 0
    if operations[op] == 1 and op == 1:
        new_mat[x,y] =  new_mat[x,y] + 1
    if operations[op] == 1 and op == 2:
        new_mat[x,y] =  new_mat[x,y] - 1
    if operations[op] == 1 and op == 3:
        temp_mat = new_mat
        temp_mat[x,y] = new_mat[y,x]
        temp_mat[y,x] = new_mat[x,y]
        new_mat = temp_mat
    return new_mat

def calc_likelihood(dat,net, prior):
    likelihood = np.sum(net)
    likelihood = likelihood + prior
    return likelihood
    
    
    
def mcmc_mh(sample, burnins, time, dat, vertices, prior):
    scores = []
    net = np.zeros((vertices,vertices))
    final_net = np.zeros((vertices,vertices))
    networks = []
    score_zero = calc_likelihood(dat, net, prior) 
    burns = 0
    converged = 0
    for burn in range(0, sample+burnins):
        net = alter_net(net, time)
        #print net
        score = calc_likelihood(dat, net, prior)
        if score >= score_zero:
	    score_zero = score
	    scores.append(score)
            burns += 1
        if burn < burnins and not converged:
            print burns
        if burns >= burnins: # may include another converging condition here otherwise mininmum burnin is assumed
	    burns += 1
	    #print burns
	    converged = 1
	    final_net = np.add(final_net, net)
    return scores, final_net/sample

mcmc_mh(100, 50, 5, 0, 4, 0)


import matplotlib.pyplot as plt
plt.plot(mcmc_mh(10000, 5000, 5, 0, 4, 0)[0])
plt.ylabel('scores')
plt.show()

    
# Poisson process

import math
import random

def next_poisson(rateParameter):
    return -math.log(1.0 - random.random()) / rateParameter

n = []
for i in range(0,100):
    next_poisson(nextTime(0.5))

   
## generate normal distribution
  
def generate_normal_population(mu, N, sigma): #Extract samples from a normal distribution
    return np.random.normal(mu, scale=sigma, size=N)
    
def pdf_model(x, p):
    mu1, sig1, mu2, sig2, pi_1 = p
    return pi_1*py.normpdf(x, mu1, sig1) + (1-pi_1)*py.normpdf(x, mu2, sig2)

plt.hist(s, 10, normed=True)
plt.show()

### Expectation maximization Algorithms
# To guess the distribution parameters that eventually generated mixed observations.
# Assume gaussian distributions ~N(mu1, sigma1) and ~N(mu2, sigma2)
# Genearate data 'normally distributed' 

N = 100
a = 0.40
m1 = 0.5         # true mean 1 this is what we want to guess
m2 = 0.8       # true mean 2 this is what we want to guess
s1 = generate_normal_population(m1, N*a, 0.5)
s2 = generate_normal_population(m2, N*(1-a), 0.5)
sig1 = s1
sig2 = s2
s = np.concatenate([s1, s2])   # put all together
sigma_tot = np.concatenate([sig1, sig2])
 
#plt.hist(s, bins=np.r_[-1:2:0.025], alpha=0.3, color='g', histtype='stepfilled');
#ax = py.twinx(); ax.grid(False)
#ax.plot(s, 0.1/sigma_tot, 'o', mew=0, ms=6, alpha=0.4, color='b')
#plt.xlim(-0.5, 1.5)
#plt.title('Sample to be fitted')
#plt.show()


def expectation_maximization(s):
    p_start = np.array([0.4, 0.5, 0.7, 0.5, 0.4]) # Parameter guess mu:=mean, sig:= SD, pi:= fraction (m1 and m2) 
    mu1, sig1, mu2, sig2, pi1 = p_start
    mu = np.array([sig1, sig2])
    sig = np.array([sig1, sig2])
    pi = np.array([pi1, 1-pi1])
    p_new = p_start
    delta = 0.000001
    improvement = float('inf')
    gamma = np.zeros((2, s.size))
    N_ = np.zeros(2)
    counter = 0
    while improvement > delta:
# Compute the responsibility func. and new parameters
        for k in [0,1]:
            gamma[k,:] = pi[k]*py.normpdf(s, mu[k], sig[k])/pdf_model(s, p_new)   # responsibility
            N_[k] = 1.*gamma[k].sum() # effective number of objects to k category
            mu[k] = sum(gamma[k]*s)/N_[k] # new sample mean of k category
            sig[k] = np.sqrt( sum(gamma[k]*(s-mu[k])**2)/N_[k] ) # new sample var of k category
            pi[k] = N_[k]/s.size # new mixture param of k category
            # updated parameters will be passed at next iter
            p_old = p_new
            p_new = [mu[0], sig[0], mu[1], sig[1], pi[0]]
            # check convergence
            improvement = max(abs(p_old[0] - p_new[0]), abs(p_old[1] - p_new[1]) )
            counter += 1
            print "Means: %6.3f %6.3f" % (p_new[0], p_new[2])
            print "Std dev: %6.3f %6.3f" % (p_new[1], p_new[3])
            print "Mix (1): %6.3f " % p_new[4]
            print "Total iterations %d" % counter
            print pi.sum(), N_.sum()

            
            
p_start = np.array([0.4, 0.5, 0.7, 0.5, 0.4]) # Parameter guess mu:=mean, sig:= SD, pi:= fraction (m1 and m2) 
mu1, sig1, mu2, sig2, pi1 = p_start
mu = np.array([sig1, sig2])
sig = np.array([sig1, sig2])
pi = np.array([pi1, 1-pi1])



### Distance measures
## Let a and b be two vectors

import math*
import operator
import numpy as np

# Euclidean
def euclidean_dist(a,b):
    d=0
    for i in range(len(a)):
        d += (a[i]-b[i])**2
        return math.sqrt(d)

# Hamming

def hamming_dist(a, b):
    if len(a) == len(b):
        d = 0
        for i in range(len(a)):
            a[i] == b[i]
            d += 1
        return d
    else:
        return None

# Cosine

def cosine_dist(a,b):
    a = np.array(a)
    b = np.array(b)
    dp = (a*b).sum()
    mod_a = (a*a).sum()
    mod_b = (b*b).sum()
    d = 1 - (dp/sqrt(mod_a*mod_b))
    return d

# Edit

def edit_dist(a, b):
    len_a = len(a)
    len_b = len(b)
    dp_mat = LCS(a, b, len_a, len_b)
    seq_len = max(max(dp_mat))
    d = (len_a + len_b) - (2*seq_len)
    return d

# Jaccard # In scipy module := scipy.spatial.distance.jaccard(u, v)

def jaccard_dist(a,b):
    uni = list(set(a)|set(b))
    inter = list (set(a) & set(b))
    d = 1 - (len(inter)/len(uni))
    return d


# Page rank

from PageRank import PageRanker

web = ((0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (1, 0, 0, 0))
pr = PageRanker(0.85, web)
pr.improve_guess(100)
print pr.getPageRank()


def page_rank(G, damping):
    n = G.shape[0]
    mat = csc_matrix(G, dtype=np.float)
    sums = np.array(mat.sum(1))[:0]
    r = mat.nonzero()
    c = mat.nonzero()
    mat.data /= sums[r]
    sink = 0
    sums = 0

def pagerank(graph, damping=0.85, epsilon=1.0e-8):
    edge_in = []
    edge_out = []
    for i in range(graph.shape[0]):
        for j in range(graph.shape[1]):
    
    
    def new_node(node):
        if node not in inlink_map: inlink_map[node] = set()
        if node not in outlink_counts: outlink_counts[node] = 0
    
    for tail_node, head_node in graph:
        new_node(tail_node)
        new_node(head_node)
        if tail_node == head_node: continue
        
        if tail_node not in inlink_map[head_node]:
            inlink_map[head_node].add(tail_node)
            outlink_counts[tail_node] += 1
    
    all_nodes = set(inlink_map.keys())
    for node, outlink_count in outlink_counts.items():
        if outlink_count == 0:
            outlink_counts[node] = len(all_nodes)
            for l_node in all_nodes: inlink_map[l_node].add(node)
    
    initial_value = 1 / len(all_nodes)
    ranks = {}
    for node in inlink_map.keys(): ranks[node] = initial_value
    
    new_ranks = {}
    delta = 1.0
    n_iterations = 0
    while delta > epsilon:
        for node, inlinks in inlink_map.items():
            new_ranks[node] = ((1 - damping) / len(all_nodes)) + (damping * sum(ranks[inlink] / outlink_counts[inlink] for inlink in inlinks))
        delta = abs(sum(new_ranks.values()) - sum(ranks.values()))
        ranks = new_ranks
        n_iterations += 1
    
    return ranks, n_iterations


## Fastq files
def parse_file(self):
    with open(self.filename, 'r') as f:
        content = f.readlines()
        
        # Recreate content without lines that start with @ and +
        content = [line for line in content if not line[0] in '@+']
        
        # Now the lines you want are alternating, so you can make a dict
        # from key/value pairs of lists content[0::2] and content[1::2]
        data = dict(zip(content[0::2], content[1::2]))
    
    return data


# fastq
import random

with open("test.fastq") as input:
    with open("sample.fastq", "w") as output:
        for line1 in input:
            line2 = input.next()
            line3 = input.next()
            line4 = input.next()
            if random.randrange(0,10) == 0:
                output.write(line1)
                    output.write(line2)
                    output.write(line3)
                    output.write(line4)

