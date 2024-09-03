
import numpy as np

def viral_centrality(inList, inWeight, outList, Niter = 5, beta = 1.0, tol = 0.0001):  
    ''' User has a choice to either run each simulation until the probabilities have converged within
    a specified relative tolerance, or just iterate for 'Niter' iterations for each seed node. If Niter is less than 1, the former option is selected, and
    if Niter is 1 or greater, the latter option is selected.
    inList[i] is list of all the nodes sending connections to i; inWeight[i] is list of corresponding weights (ie transmission probabilities)
    outList[i] is  a list of all the nodes i sends connections to
    All transmission probabilities are universally multiplied by 'beta' '''

    N=len(inList)
    
    avg_infections = np.zeros(N)
    
    if Niter < 1: #if Niter is less than 1, this means we want to iterate until the uninfected array has converged to within the prescribed relative tolerance
        
        for seed in range(N):
            prev_uninfected = np.ones(N) 
            uninfected = np.ones(N) #probabilty that node hasn't been infected yet. starts at one for each node
            last_infected = np.zeros(N) #probability that node was infected on last timestep
            cur_infected = np.zeros(N) #probability that node was infected on current timestep
            last_infected[seed] = 1
            uninfected[seed] = 0
            
            t = 0
            #breadth-first search (BFS)
            bfs_queue = -1 * np.ones(N, dtype=int) #first-in/first-out buffer to perform breadth-first search (see section 10.3.3 of Mark Newman's textbook)
            bfs_queue[0] = seed
            seed_distance = -1 * np.ones(N, dtype=int) #array of distance from seed node
            seed_distance[seed] = 0
            read=0
            write=1
            while np.nanmax( (prev_uninfected[bfs_queue[0:write]] - uninfected[bfs_queue[0:write]]) / prev_uninfected[bfs_queue[0:write]]) > tol: #iterate until relative tolerance is met; do not need absolute value bc. prev_uninfected will always be geq to uninfected; ; note there will always be an 'nan' in entry corresponding to seed node, hence the nanmax function
                prev_uninfected = np.copy(uninfected)
                
                #expand BFS to find next ring of nodes that are within seed node's reach (see section 10.3.3 of Mark Newman's text)
                if read != write: #if read==write, then breadth-first search is exhausted
                    write_start = write
                    while read < write_start: #when read is equal to write_start, that signifies reaching the end of the "t+1st ring"
                        for neighb in outList[bfs_queue[read]]:
                            if seed_distance[neighb] < 0: #if distance from seed node has not yet been determined, then record it
                                seed_distance[neighb] = t+1
                                bfs_queue[write] = neighb 
                                write += 1
                        read += 1
                    
                    
                for node in bfs_queue[0:write]: #for all nodes within reach of the seed node at this point in time
                    prob_uninfected = 1 #set probability of being uninfected to one
                    for con in range(len(inList[node])): #look at each incoming connection to node
                        prob_uninfected = prob_uninfected*(1-(last_infected[inList[node][con]]*(beta*inWeight[node][con])))
                        #prob that node remains uninfected decreases as we go through each connection. beta again multiplies weights to reset max
                    cur_infected[node] = (1-prob_uninfected)*uninfected[node]
                    #prob of current infection is prob of being infected times the prob that node hasn't been infected yet
                
                for node in bfs_queue[0:write]:
                    last_infected[node] = cur_infected[node] 
                    #update our last infected list for each timestep
                    uninfected[node] = uninfected[node] - cur_infected[node]
                    #new prob of being uninfected is old prob - the prob that node is currently infected
                
                t = t+1
                    
            avg_infections[seed] = sum(1 - uninfected) - 1 #don't want to include seed node infection in total
        
    else: #if Niter is a positive integer, then just iterate for that number of time steps
    
        for seed in range(N):
            prev_uninfected = np.ones(N) 
            uninfected = np.ones(N) #probabilty that node hasn't been infected yet. starts at one for each node
            last_infected = np.zeros(N) #probability that node was infected on last timestep
            cur_infected = np.zeros(N) #probability that node was infected on current timestep
            last_infected[seed] = 1
            uninfected[seed] = 0
            
            t = 0
            #breadth-first search (BFS)
            bfs_queue = -1 * np.ones(N, dtype=int) #first-in/first-out buffer to perform breadth-first search (see section 10.3.3 of Newman's textbook)
            bfs_queue[0] = seed
            seed_distance = -1 * np.ones(N, dtype=int) #array of distance from seed node
            seed_distance[seed] = 0
            read=0
            write=1
            while t < Niter: #in contrast to original algorithm, just iterate for fixed number of time steps
                prev_uninfected = np.copy(uninfected)
                
                #expand BFS to find next ring of nodes that are within seed node's reach (see section 10.3.3 of Newman's text)
                if read != write: #if read==write, then breadth-first search is exhausted
                    write_start = write
                    while read < write_start: #when read is equal to write_start, that signifies reaching the end of the "t+1st ring"
                        for neighb in outList[bfs_queue[read]]:
                            if seed_distance[neighb] < 0: #if distance from seed node has not yet been determined, then record it
                                seed_distance[neighb] = t+1
                                bfs_queue[write] = neighb 
                                write += 1
                        read += 1
                    
                    
                for node in bfs_queue[0:write]: #for all nodes within reach of the seed node at this point in time
                    prob_uninfected = 1 #set probability of being uninfected to one
                    for con in range(len(inList[node])): #look at each incoming connection to node
                        prob_uninfected = prob_uninfected*(1-(last_infected[inList[node][con]]*(beta*inWeight[node][con])))
                        #prob that node remains uninfected decreases as we go through each connection. beta again multiplies weights to reset max
                    cur_infected[node] = (1-prob_uninfected)*uninfected[node]
                    #prob of current infection is prob of being infected times the prob that node hasn't been infected yet
                
                for node in bfs_queue[0:write]:
                    last_infected[node] = cur_infected[node] 
                    #update our last infected list for each timestep
                    uninfected[node] = uninfected[node] - cur_infected[node]
                    #new prob of being uninfected is old prob - the prob that node is currently infected
                
                t = t+1
                    
            avg_infections[seed] = sum(1 - uninfected) - 1 #dont want to include seed node infection in total
    
    return avg_infections

