###############################
# Load and store drug data in a graph.
# Lowell Milliken
###############################
import networkx as nx
import pickle

default_gene_file = 'drug_data/genes/genes.tsv'
default_relationships_file = 'drug_data/relationships/relationships.tsv'
auto_relationships_file = 'drug_data/genedrug_relationship_100417_sfsu.tsv'

# graph node
class Entity:
    def __init__(self, eid, etype):
        self.eid = eid
        self.etype = etype

    def __hash__(self):
        return hash((self.eid, self.etype))

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) != type(other):
            return False
        else:
            return self.etype == other.etype and self.eid == other.eid

# node for genes
class Gene(Entity):
    def __init__(self, eid, symbols):
        Entity.__init__(self, eid, 'Gene')
        self.symbols = symbols

# node for non-genes
class Other(Entity):
    def __init__(self, eid, etype, name):
        Entity.__init__(self, eid, etype)
        self.name = name


# Gene-Drug interaction graph class
class DrugGraph:
    # load genes and relationships from the given csv files
    # source is "original" for PharmGKB or "auto" for Mallory et al.
    def __init__(self, gene_file=default_gene_file, relationship_file=default_relationships_file, source='original'):
        self.genes, self.symbol_gene = self.load_genes(gene_file)
        self.graph, self.diseases = self.load_relationships(relationship_file, source)
        self.source = source

    # load genes from the csv file
    @staticmethod
    def load_genes(gene_file):
        genes = {}
        symbol_genes = {}
        first = True
        with open(gene_file, 'r') as infile:
            for line in infile:
                if first:
                    first = False
                    continue
                tokens = line.split('\t')
                symbols = [tokens[5]]
                alt_symbols = tokens[7]
                if alt_symbols:
                    alt_symbols = alt_symbols.replace('\"', '')
                    alt_symbols = alt_symbols.split(',')
                else:
                    alt_symbols = []

                symbols = symbols + alt_symbols
                genes[tokens[0]] = Gene(tokens[0], symbols)

                for symbol in symbols:
                    symbol_genes[symbol] = genes[tokens[0]]

        return genes, symbol_genes

    # load relationships from the relationships file
    # source determines how the csv file in read
    def load_relationships(self, relationship_file, source, threshold=0.9):
        first = True

        graph = nx.Graph()
        graph.add_nodes_from(self.genes.values())
        diseases = {}

        with open(relationship_file, 'r') as infile:
            for line in infile:
                if first:
                    first = False
                    continue

                tokens = line.split('\t')

                if source == 'original':
                    e1_id = tokens[0]
                    e1_type = tokens[2]
                    e2_id = tokens[3]
                    e2_type = tokens[5]

                    if e1_type == 'Gene':
                        e1 = Gene(e1_id, [tokens[1]])
                    else:
                        e1 = Other(e1_id, e1_type, tokens[1])
                        if e1.etype == 'Disease':
                            diseases[e1.name] = e1

                    if e2_type == 'Gene':
                        e2 = Gene(e2_id, [tokens[4]])
                    else:
                        e2 = Other(e2_id, e2_type, tokens[4])
                        if e2.etype == 'Disease':
                            diseases[e2.name] = e2

                    graph.add_nodes_from([e1, e2])
                    association = tokens[7]

                    if association == 'associated':
                        graph.add_edge(e1, e2)
                else:
                    if tokens[9] not in self.symbol_gene:
                        continue
                    e1 = self.symbol_gene[tokens[9]]
                    e2_id = tokens[8]
                    e2_type = 'Chemical'

                    e2 = Other(e2_id, e2_type, e2_id)

                    association = tokens[12]
                    if (e1, e2) in graph.edges:
                        graph.edges[e1, e2]['weight'] = graph.edges[e1, e2]['weight'] + float(association)
                        graph.edges[e1, e2]['count'] = graph.edges[e1, e2]['count'] + 1
                    graph.add_edge(e1, e2, weight=float(association), count=1)

        return graph, diseases

    # Get node with given name
    def get_entity(self, name):
        if name in self.diseases:
            return self.diseases[name]
        elif name in self.symbol_gene:
            return self.symbol_gene[name]
        else:
            return None

    # Find neighboring drugs from node with given name with confidence over the threshold
    def find_drugs(self, name, threshold=0.9):
        drugs = []
        entity = self.get_entity(name)
        if entity:
            neighbors = self.graph[entity]
            for neighbor in list(neighbors):
                if neighbor.etype == 'Chemical':
                    if self.source == 'auto':
                        weight = self.graph.edges[entity, neighbor]['weight']
                        count = self.graph.edges[entity, neighbor]['count']
                        if weight/count < threshold:
                            continue

                    drugs.append(neighbor.name)

        return drugs


# create and save drug graph to pickle file
def save_drug_graph(relationship_file, source, outfilename):
    druggraph = DrugGraph(relationship_file=relationship_file, source=source)
    print("writing to file")
    with open(outfilename, 'wb') as outfile:
        pickle.dump(druggraph, outfile)


# load drug graph from pickle file
def load_drug_graphs():
    with open('pharmgkbDG.pickle', 'rb') as infile:
        dgp = pickle.load(infile)
    with open('malloryDG.pickle', 'rb') as infile:
        dgm = pickle.load(infile)

    return dgp, dgm


# create and save drug graphs for both data sources.
def save_all():
    print("In file")
    save_drug_graph(default_relationships_file, 'original', 'pharmgkbDG.pickle')
    save_drug_graph(auto_relationships_file, 'auto', 'malloryDG.pickle')
