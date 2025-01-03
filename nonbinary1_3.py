from ete3 import Tree
from LIS import LIS_len, LIS_seq
from random import randint
from LCS1 import process_trees

# Separator between a label and a node name in a tree representation
SEP_NODE = ":"
# Separators between tree representation elements
SEP_VEC = ","

# Main Class
class TreeVec:
    """
    Vector representation of a tree.

    The topology of a tree with n leaves, ordered from 1 to n is encoded by a
    list of 2n integer labels
    - start with 1,
    - ends with n,
    - contains 2 occurrences of every integer in {1,...,n}
    - the first occurrence of i>1 appears before the second copy of i-1
    - the second occurrence of i>1 appears after the second occurrence of i-1
    - the first occurrence of i encodes an internal node
    - the second occurrence of i encodes a leaf
    The tree is augmented by a root labeld 1 and with a single child called the
    dummy root.
    
    Data structure: list([int,str,float,bool])
    - field 0 (int): label
    - field 1 (str): name of the node in the tree
    - field 2 (float): length of the branch to the parent;
      the root and the dummy root have a branch length equal to 0.0
    - field 3 (bool): True if second occurrence (leaf)
                      False if first occurrence (internal node)

    String encoding
    A tree vector representation can be written in format 1 or 2 and in compact or
    non-compact writing:
    - nodes are separated by SEP_VEC
    - format 1.non-compact: each node is written as label:name:dist
    - format 2.non-compact: each node is written as label:name
    - format 1.compact:
      - each internal node is written as label:name:dist
      - each leaf is written as dist
        to be decoded this requires a mapping idx2leaf (dict int -> str) that
        defines a total order on leaves and allows to recover the leaf name and label
        associated to positions in the vector encoding leaves
    - format 2.compact:
      - each internal node is written as label:name
      - each leaf is written as an empty string whose name and labels can be recovered
        from the mapping idx2leaf as described above
    """

    def __init__(
            self,
            treevec_vec=None,
            tree=None,
            newick_str=None,
            treevec_str=None,
            leaf2idx=None,
            idx2leaf=None,
            format=None,
            compact=None
    ):
        """
        Instantiate a vector representation for a tree on n leaves
        - If treevec_vec is not None, the vector is created using it as vector
        - If tree is not None, tree is a Tree object and the vector is created from it
          using leaf2idx
        - If newick_str is not None it is created from newick_str using idx2leaf and
          expected in Newick format=1   
        - If treevec_str is not None it is created from treevec_str using idx2leaf and
          expected in format defined by format and compact
        - Otherwise an empty vector is created
        - leaf2idx (dict str -> int): leaf name to index in a total order on leaves
          (1-base)
        - idx2leaf (dict int -> str): reverse dictionary
        - format (int in [1,2])
        - compact (bool)
        """
        self.vector = []
        self.simvec = []
        if treevec_vec is not None:
            self.vector = treevec_vec
        elif tree is not None:            
            self.vector = self.tree2treevec(tree, leaf2idx=leaf2idx)
        # elif newick_str is not None:
            # self.vector = self.newick2treevec(newick_str, leaf2idx=leaf2idx)
        # elif treevec_str is not None:
            # self.vector = self.str2treevec(
                # treevec_str, idx2leaf,
                # format=format, compact=compact
            # )    
        self.simvec = [sublist[1] for sublist in self.vector]
    
    def find(self, label: dict, k: int):
        # 从键为1开始遍历字典
        for key in range(1, len(label)):
            value = label.get(key)
            # 如果值是集合并且集合中含有k
            if isinstance(value, set) and k in value:
                return key
        # 如果没有找到，返回一个标识，比如-1
        return -1

    
    def treevec2tree(self):
        """
        Given a tree vector representation, compute a Tree object
        Ouput:
        - (Tree) Tree object
        """
        v = self.vector
        # n = int(((a-1)*len(v)+2-a) / a)
        # Computing the label and name of each node
        __label, __name, __dist, edges = {}, {}, {}, {}
        for i in range(0,len(v)):
            [label,name,dist,leaf] = v[i]
            __label[i] = label
            edges[label] = []
            __name[label] = name
            __dist[label] = dist
        # Decoding into edges
        # edges = {i: [] for i in range(1,len(v)+1)}
        o = {j: (2 if v[j][3] else 1) for j in range(0,len(v))}
        for j in range(0,len(v)-1):
            if o[j] == 1 and o[j+1] == 1:
                edges[__label[j]].append(__label[j+1])
            elif o[j] == 1 and o[j+1] == 2:
                edges[__label[j]].append(__label[j+1])
            elif o[j] == 2 and o[j+1] == 1:
                k = j+2
                while o[k] == 1: k += 1
                edges[__label[self.find(__label, __label[k])]].append(__label[j+1])
            elif o[j] == 2 and o[j+1] == 2:
                edges[__label[self.find(__label, __label[j+1])]].append(__label[j+1])
        # Creating Tree structure
        # nodes = {i: Tree(name=__name[i],dist=__dist[i]) for i in range(1,2*n+1)}
        nodes = {}
        for i in range(0,len(v)):
            nodes[__label[i]] = Tree(name=__name[i],dist=__dist[i])
        for node_from_idx,nodes_to_idx in edges.items():
            node_from = nodes[node_from_idx]
            for node_to_idx in nodes_to_idx:
                node_from.add_child(nodes[node_to_idx])
        root = nodes[__label[0]].children[0]
        return root
    
    def tree2treevec(self, tree, leaf2idx=None):
        """
        Compute the vector representation of a tree with n leaves rom a Tree objet
        Input:
        - t: Tree object with features "name" and "dist" (branch length)
        - leaf2idx: dict(str -> int) leaf name to leaf label
        if None: leaf labels added during a postorder traversal in order of visit.
        """
        # Adding a root labeled 1 and named ""
        T = Tree(name="")
        if leaf2idx is not None:
            for key, value in leaf2idx.items():
                    if value == 1:
                        T.name = {key}
        T.add_feature("label", {1})
        T.add_child(tree)
        # Labeling nodes
        label = 1
        for node in tree.traverse("postorder"):
            if node.is_leaf() and leaf2idx is not None:
                node.add_feature("label", leaf2idx[node.name])
                node.add_feature("min_label", node.label)
            elif node.is_leaf():
                node.add_feature("label", label)
                node.add_feature("min_label", node.label)
                label += 1
            else:
                children_min_label = [child.min_label for child in node.children]
                node.add_feature("min_label", min(children_min_label))
                # node.add_feature("label", max(children_min_label))
                node.add_feature("label", set(children_min_label).difference({node.min_label}))
                new_set = set()
                for key, value in leaf2idx.items():
                    if value in node.label:
                        new_set.add(key)
                node.name = new_set
        # Computing a dictionary from label to leaf
        label2leaf = {}
        for node in T.traverse("postorder"):
            if node.is_leaf():
                label2leaf[node.label] = node
                n = len(label2leaf.keys())
                # Computing paths from leaves to the internal node of same label
        paths = {}
        for i in range(1,n+1):
            paths[i] = []
            node = label2leaf[i].up
            # while node.label != i:
            while i not in node.label:
                paths[i].append([node.label, node.name, node.dist, False])
                node = node.up
        # Concatenating reversed paths
        v = [[{1}, T.name, 0.0, False]]
        for i in range(1,n+1):
            leaf = label2leaf[i]
            v += paths[i][::-1] + [[i,leaf.name,leaf.dist,True]]
        return v
    
    def hop_similarity(self, t2, compute_seq=False):
        """
        Compute the hop smilarity to another tree representations
        Input:
        - t2 (TreeVec)
        assumption: both are on the same leaves order (asserted)
        - compute_seq (bool): if True, returns an actual LCS, if False, returns
          the similarity value
        Output:
        - compute_seq=False: (int) in [0,n]
        - compute_seq=True: list((int,bool)) list of (integers,True if leaf)
          encoding the LCS between v1 and v2
        """
        
        def __relabel_segment(segment1, segment2):
            """
            Relabel the labels of segment1 increasingly from 0 
            and the labels of segment2 according to the relabeling of 
            segment1, excluding labels not present in segment1
            """
            for j in range(0,len(segment1)):
                if isinstance(segment1[j], set):
                    segment1[j] = frozenset(segment1[j])
            for j in range(0,len(segment2)):
                if isinstance(segment2[j], set):
                    segment2[j] = frozenset(segment2[j])
            # map1[x] = position of label x in segment1
            map1 = {segment1[i1]: i1 for i1 in range(0,len(segment1))}
            # Relabeling segment2 according to __map1,
            # excluding labels not in segment1
            relabeled_segment2 = []
            for i2 in range(0,len(segment2)):
                if segment2[i2] in map1.keys():
                    relabeled_segment2.append(map1[segment2[i2]])
            return relabeled_segment2
        
        def __Partition(segment1, segment2):
            result = []
            for A in segment1:
                remaining_elements = A.copy()  
                for B in segment2:
                    intersection = A & B  
                    if intersection:
                        result.append(intersection)  
                    remaining_elements -= intersection  
                if remaining_elements:
                    result.append(remaining_elements)  
            return result
    
        v1,v2 = self.vector,t2.vector
        n=0
        for i in range(0,len(v1)):
            if v1[i][3]:
                n += 1
        # n = int(len(v1)/2)
        second_occ_order = [
            v1[i][0]
            for i in range(0,len(v1)) if v1[i][3]
        ]
        # Compute a list of pairs of subsequences to compare pairwise
        # boundaries1[i] = [j,k]: boundaries of the segment of internal nodes
        # in v1 before leaf i+1 similar for boundaries2 and v2
        # if j>k: empty segment
        boundaries = {1: [], 2: []}
        i1,i2 = 0,0
        for i in range(0,len(v1)):
            leaf1,leaf2 = v1[i][3],v2[i][3]
            if leaf1:
                boundaries[1].append([i1,i-1])
                i1 = i+1
            if leaf2:
                boundaries[2].append([i2,i-1])
                i2 = i+1
        
        def nonbin_sim():
            lcs_len,lcs_seq = 0,[]
            for j in range(0,n):
                b1_start,b1_end = boundaries[1][j][0], boundaries[1][j][1]
                b2_start,b2_end = boundaries[2][j][0], boundaries[2][j][1]
                # Checking that both segments are non-empty (otherwise, no LCS)
                if (b1_end>=b1_start) and (b2_end>=b2_start):
                    # Segments of v1 and v2 to consider
                    __segment1 = [v1[k][0] for k in range(b1_start, b1_end+1)] 
                    __segment2 = [v2[k][0] for k in range(b2_start, b2_end+1)] 
                    # Relabeling __segment2 according to __map1,
                    # excluding labels not in __segment1            
                    segment2 = __Partition(__segment1, __segment2)
                    
        # return (lcs_seq if compute_seq else lcs_len)
        
        # Computes an LCS for each pair of segments using an LIS algorithm
        lcs_len,lcs_seq = 0,[]
        for j in range(0,n):
            b1_start,b1_end = boundaries[1][j][0], boundaries[1][j][1]
            b2_start,b2_end = boundaries[2][j][0], boundaries[2][j][1]
            # Checking that both segments are non-empty (otherwise, no LCS)
            if (b1_end>=b1_start) and (b2_end>=b2_start):
                # Segments of v1 and v2 to consider
                __segment1 = [v1[k][0] for k in range(b1_start, b1_end+1)] 
                __segment2 = [v2[k][0] for k in range(b2_start, b2_end+1)] 
                # Relabeling __segment2 according to __map1,
                # excluding labels not in __segment1            
                segment2 = __relabel_segment(__segment1, __segment2)
                # Computing an LIS in segment2
                a, b = process_trees(segment2)
                if compute_seq: lcs_seq += [
                        (__segment1[i2],False)
                        for i2 in b
                ]
                else: lcs_len += a
            lcs_seq += [(second_occ_order[j],True)]
        return (lcs_seq if compute_seq else lcs_len)
    
    def hop_distance(self, t2):
        v1,v2 = self.vector,t2.vector
        n=0
        for i in range(0,len(v1)):
            if v1[i][3]:
                n += 1
        return (n - self.hop_similarity(self, t2))