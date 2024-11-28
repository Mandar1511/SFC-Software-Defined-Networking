from gurobipy import *
import random
import sys
import networkx as nx
model = Model("sfc")


#########################################################################INPUT#######################################################
num_nodes = 20
num_vnf = 15
num_sfc = 5
num_edges = 40

nodes = [i for i in range(num_nodes)]
W = ["vnf"+str(i) for i in range(num_vnf)]
Y = ["sfc"+str(i) for i in range(num_sfc)]

edges = []

while len(edges) < num_edges:
	first = str(random.randint(0,num_nodes))
	second = str(random.randint(0,num_nodes))
	if first != second:
		if (first, second) not in edges and (second, first) not in edges:
			edges.append((first, second))


latency = {}
for x in edges:
	if x not in latency:
		latency[x] = 0
	latency[x] = random.randint(100,200)

c1 = {} # compute available at v
for x in nodes:
	if x not in c1:
		c1[x] = 0
	c1[x] = random.randint(100,200)

m1 = {} # memory available at v
for x in nodes:
	if x not in m1:
		m1[x] = 0
	m1[x] = random.randint(1000,2000)
	
b1 = {} # bandwidth of edges
for x in edges:
	if x not in b1:
		b1[x] = 0
	b1[x] = random.randint(20,40)

SFC = {} # contains list of VNFs for each SFC
for sfc in Y:
	if sfc not in SFC:
		SFC[sfc] = []
	# assign 3-6 vnfs randomly to each sfc
	SFC[sfc] = random.sample(W, random.randint(3,6))

b2 = {} # store bandwidth required for sfc
for x in Y:
	if x not in b2:
		b2[x] = 0
	b2[x] = random.randint(10,20)
	
c2 = {} # compute required for vnfs
for x in W:
	if x not in c2:
		c2[x] = 0
	c2[x] = random.randint(20,80)

m2 = {} # compute required for vnfs
for x in W:
	if x not in m2:
		m2[x] = 0
	m2[x] = random.randint(500,1000)

delay = {} # delay threshold for SFCs
for x in Y:
	if x not in delay:
		delay[x] = 0
	delay[x] = random.randint(100,200)

originate = {} # originating node for SFCs
for x in Y:
	if x not in originate:
		originate[x] = 0
	originate[x] = nodes[random.randint(0, len(nodes)-1)]

############################################################## networkx ################################################
graph = nx.Graph()
for e in edges:
	graph.add_edge(e[0], e[1], weight=latency[e])

print(graph.nodes)

def check(e, v1, v2):
	# this function checks if edge e is on shortest path from v1 to v2
	try:
		path = nx.shortest_path(graph, source=v1, target=v2, weight="weight", method="dijkstra")
		edge_in_path = 0
		for i in range(len(path)-1):
			if (path[i] == v1 and path[i + 1] == v2) or (path[i] == v2 and path[i + 1] == v1):
				edge_in_path = 1
				break
		return edge_in_path
	except nx.NetworkXNoPath:
		print(f"No path exists between {v1} and {v2}.")

def len_shortest_path(v1, v2):
	try:
		path = nx.shortest_path(graph, source=v1, target=v2, weight="weight", method="dijkstra")
		path_len = 0
		for i in range(len(path) - 2):
			path_len += latency[(path[i], path[i+1])]
		return path_len
	except nx.NetworkXNoPath:
		print(f"No path exists between {v1} and {v2}.")
		
############################################################### VARIABLES ###############################################
X = model.addVars(W,nodes,vtype=GRB.BINARY, name="X")
Z = model.addVars(Y, vtype=GRB.BINARY, name="Z")

###################################################### VNF constraints######################################################

# C0 - one VNF can be provisioned on atmost 1 compute node
for vnf in W:
	model.addConstr(quicksum(X[vnf,node] for node in nodes) <=1)
	
# C1 - Compute capacity constraint for each node
for node in nodes:
	model.addConstr(quicksum(X[vnf, node]*c2[vnf] for vnf in W) <= c1[node])

# C2 - Memory constraint for each node
for node in nodes:
	model.addConstr(quicksum(X[vnf, node]*m2[vnf] for vnf in W) <= m1[node])

################################################### SFC constraints #######################################################

# C3 - If SFC is provisioned then all of its VNFs are provisioned
for sfc in Y:
	model.addGenConstrIndicator(Z[sfc], True, quicksum(X[vnf,node] for vnf in SFC[sfc] for node in nodes), GRB.EQUAL, len(SFC[sfc]))

# C4 - if sfc is provisioned, then last and first vnf of that sfc chain will be placed on originate node
for sfc in Y:
	chain = SFC[sfc]
	model.addGenConstrIndicator(Z[sfc], True, X[chain[0], originate[sfc]], GRB.EQUAL, 1)
	model.addGenConstrIndicator(Z[sfc], True, X[chain[-1], originate[sfc]], GRB.EQUAL, 1)

# C5 - Link capacity constraint
# Bandwidth of each link must be greater than or equal to sum of bandwidths of all SFCs using that link
# vnf1-get_next(sfc,vnf1) is link in sfc

def get_next(sfc, vnf1):
	# returns the vnf next to vnf1 in chain
	for i in range(len(SFC[sfc])):
		if SFC[sfc][i] == vnf1:
			if i < len(SFC[sfc]):
				return SFC[sfc][i+1]

for e in edges:
	model.addConstr(b1[e] >= quicksum(Z[sfc] * quicksum(check(e,v1,v2)*X[vnf1,v1]*X[get_next(sfc,vnf1),v2]*b2[sfc] for v1 in nodes for v2 in nodes for vnf1 in SFC[sfc][:-1]) for sfc in Y))


# c6 - latency threshold of each SFC is met
for sfc in Y:
	model.addConstr(delay[sfc] <= quicksum(Z[sfc] * quicksum(len_shortest_path(v1,v2) *X[vnf1,v1]*X[get_next(sfc,vnf1),v2] for v1 in nodes for v2 in nodes for vnf1 in SFC[sfc][:-1])))
	


#################################################### Objective #################################################################
model.setObjective(quicksum(Z[sfc] for sfc in Y), GRB.MAXIMIZE)


status = model.Status
if status == GRB.INFEASIBLE:
	model.computeIIS()
	model.write("infeasible_model.ilp")
if status == GRB.UNBOUNDED:
	print("The model cannot be solved because it is unbounded")
	sys.exit(0)
if status == GRB.OPTIMAL or status == GRB.TIME_LIMIT:
	for vnf in W:
		for node in nodes:
			if X[vnf, node].X == 1:
				print(f"Place {vnf} at node {node}")
	for sfc in Y:
		if Y[sfc].X == 1:
			print(f"{sfc} is fully provisioned")
