import sys
from queue import Queue


class Node:
    def __init__(self, nid, parent_id, level):
        self.id = nid
        self.parent_id = parent_id
        self.children = set()
        self.level = level

    def AddChild(self, nid):
        if nid != self.id:
            self.children.add(nid)

    def IsLeaf(self):
        return len(self.children) == 0

class GTDB:
    def __init__(self, taxonomy_path: str):
        self.path = taxonomy_path
        self.nodes = []
        self.names2id = dict()
        self.ids2name = dict()
        self.current_id = 0
        self.root = 0

        self.AddTaxon("r__Root", "r__Root", 0)

        self.Build()

    def AddTaxon(self, name: str, parent_name: str, level: int):
        if name not in self.names2id:
            self.names2id[name] = self.current_id
            self.ids2name[self.current_id] = name
            parent_id = self.names2id[parent_name] if name != parent_name else self.current_id

            node = Node(self.current_id, parent_id, level)
            self.nodes.append(node)
            self.nodes[parent_id].AddChild(self.current_id)

            self.current_id += 1

    def CountLeavesWorker(self, nid: int):
        node_count = 0

        if len(self.nodes[nid].children) == 0:
            return 1;

        for childid in self.nodes[nid].children:
            node_count += self.CountLeavesWorker(childid)

        return node_count

    def CountLeaves(self, name):
        visited = set()

        q = Queue()
        count = 0

        q.put(self.names2id[name])

        while not q.empty():
            temp = q.queue[0]
            q.get()

            visited.add(temp)

            node = self.nodes[temp]

            for child in node.children:
                if len(self.nodes[child].children) == 0:
                    count += 1
                else:
                    q.put(child)
        
        return count

    def CountLeavesRec(self, name):
        return self.CountLeavesWorker(self.names2id[name])

    def Build(self):
        with open(self.path, 'r') as file:
            for line in file:
                line = line.rstrip()

                tokens = line.split('\t')

                genome = tokens[0]

                lineage = tokens[1].split(';')

                parent_taxon = self.ids2name[self.root]
                idx = 0
                for idx in range(0, len(lineage)):
                    taxon = lineage[idx]
                    if idx > 0:
                        parent_taxon = lineage[idx - 1]

                    self.AddTaxon(taxon, parent_taxon, idx + 1)

                idx += 1
                self.AddTaxon(genome, taxon, idx + 1)

    def Get(self, name):
        nid = self.names2id[name]
        return self.nodes[nid]

    def GetParent(self, name):
        nid = self.names2id[name]
        return self.nodes[self.nodes[nid].parent_id]

    def LCA(self, node1, node2):
        nid1 = 0
        nid2 = 0
        if isinstance(node1, int) and isinstance(node2, int):
            nid1 = node1
            nid2 = node2
        else:
            nid1 = self.names2id[node1]
            nid2 = self.names2id[node2]

        node_a = self.nodes[nid1]
        node_b = self.nodes[nid2]

        while node_a.level > node_b.level:
            node_a = self.nodes[node_a.parent_id]

        while node_b.level > node_a.level:
            node_b = self.nodes[node_b.parent_id]

        while node_b.id != node_a.id:
            node_a = self.nodes[node_a.parent_id]
            node_b = self.nodes[node_b.parent_id]

        return(self.ids2name[node_a.id])

    def IsAncestorOf(self, node1, node2):
        nid1 = self.names2id[node1]
        nid2 = self.names2id[node2]

        node_a = self.nodes[nid1]
        node_b = self.nodes[nid2]

        if node_a.level > node_b.level:
            return False

        while node_b.level > node_a.level:
            node_b = self.nodes[node_b.parent_id]

        return node_a.id == node_b.id

    def GetAt(self, name, rank):
        nid = self.names2id[name]

        node = self.nodes[nid]
        while node.level > 0 and self.ids2name[node.id][0] != rank:
            node = self.nodes[node.parent_id]

        return self.ids2name[node.id]

    def RootName(self):
        return self.ids2name[self.root]

    def LeafsFromNode(self, name: str):
        if name not in self.names2id.keys():
            print("Unknown taxon {}".format(name))
            exit()

        leafs = []

        queue = []

        nid = self.names2id[name]
        node = self.nodes[nid]

        queue.append(node)

        while queue:
            cnode = queue.pop(0)

            for cid in cnode.children:
                child = self.nodes[cid]
                if child.IsLeaf():
                    leafs.append(self.ids2name[cid])
                else:
                    queue.append(self.nodes[cid])

        return leafs


    def NodesFromNode(self, name: str, level: int):
        if name not in self.names2id.keys():
            print("Unknown taxon {}".format(name))
            exit()

        targets = []

        queue = []

        nid = self.names2id[name]
        node = self.nodes[nid]

        if node.level > level: return None;

        queue.append(node)

        while queue:
            cnode = queue.pop(0)

            for cid in cnode.children:
                child = self.nodes[cid]
                if child.level == level:
                    targets.append(self.ids2name[cid])
                else:
                    queue.append(self.nodes[cid])

        return targets


    def NodesAtLevel(self, level: int):
        return self.NodesFromNode(self.ids2name[self.root], level)


    def PrintLineage(self, name):
        nid = self.names2id[name]

        node = self.nodes[self.nodes[nid].parent_id]
        lineage_str = name
        while node.level != 0:
            lineage_str = self.ids2name[node.id] + ";" + lineage_str
            node = self.nodes[node.parent_id]
        return lineage_str


# if __name__ == "__main__":
#     gtdb_path = sys.argv[1]

#     taxonomy = GTDB(gtdb_path)

#     print(taxonomy.ids2name[taxonomy.Get("o__Oscillospirales").id])
#     for child_id in taxonomy.Get("o__Oscillospirales").children:
#         print(taxonomy.ids2name[child_id])

#     taxonomy.LCA(3,5)
#     taxonomy.LCA("f__CAG-382","f__Oscillospiraceae")

#     print(taxonomy.IsAncestorOf("p__Firmicutes", "f__GUT_GENOME256086"))

