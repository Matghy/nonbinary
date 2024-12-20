def process_trees(n, trees):
    """
    处理一个包含 n 个列表的列表。

    参数：
    n (int): 列表 trees 的长度。
    trees (list of list): 一个包含 n 个子列表的列表。

    返回：
    tuple: 包含 LCS 的长度和 LCS 本身的元组。
    """
    # 检查输入是否合法
    if not isinstance(n, int) or n <= 0:
        raise ValueError("参数 n 必须是一个正整数。")

    if not isinstance(trees, list) or len(trees) != n:
        raise ValueError("参数 trees 必须是长度为 n 的列表。")

    for tree in trees:
        if not isinstance(tree, list) or not all(isinstance(x, int) for x in tree):
            raise ValueError("trees 中的每个元素必须是仅包含整数的列表。")

    # 找出所有列表的公共元素
    common_elements = set(trees[0])
    for tree in trees[1:]:
        common_elements &= set(tree)

    # 初始化结果字典
    result = {element: {"predecessors": [], "successors": []} for element in common_elements}
    # print(result.keys())

    # 对公共元素中的所有元素对进行前驱和后继关系处理
    for a in common_elements:
        for b in common_elements:
            if a != b:
                is_predecessor = all(tree.index(a) < tree.index(b) for tree in trees if a in tree and b in tree)
                if is_predecessor:
                    result[a]["successors"].append(b)
                    result[b]["predecessors"].append(a)
    # print(result.items())

    # 定义递归函数来处理末尾元素及其恰好前驱元素
    def process_terminal(element, path, visited):
        """
        递归处理末尾元素，生成最长路径。

        参数：
        element (int): 当前处理的元素。
        path (list): 当前路径。
        visited (set): 已访问的节点。

        返回：
        tuple: 包含路径长度和路径本身的元组。
        """
        # 找到所有的恰好前驱元素
        immediate_predecessors = set()
        for tree in trees:
            if element in tree:
                index = tree.index(element)
                closest_predecessor = None
                for predecessor in result[element]["predecessors"]:
                    if predecessor in tree and tree.index(predecessor) < index:
                        if closest_predecessor is None or tree.index(predecessor) > tree.index(closest_predecessor):
                            closest_predecessor = predecessor
                # print(closest_predecessor)
                if closest_predecessor:
                    immediate_predecessors.add(closest_predecessor)

        # 去重并更新结果字典
        result[element]["immediate_predecessors"] = list(immediate_predecessors)

        # 更新路径
        current_path = path + [element]

        # 递归处理恰好前驱元素
        max_length = len(current_path)
        max_path = current_path

        for predecessor in immediate_predecessors:
            if predecessor not in visited:
                visited.add(predecessor)
                length, pred_path = process_terminal(predecessor, current_path, visited)
                if length > max_length:
                    max_length = length
                    max_path = pred_path

        return max_length, max_path

    # 找到所有末尾元素
    terminal_elements = [key for key, value in result.items() if not value["successors"]]

    # 初始化全局最大长度和路径
    global_max_length = 0
    global_max_path = []

    # 对每个末尾元素计算最长路径
    for terminal in terminal_elements:
        visited = set()
        length, path = process_terminal(terminal, [], visited)
        if length > global_max_length:
            global_max_length = length
            global_max_path = path
    # print(result.items())

    # 返回 LCS 的长度和路径
    return global_max_length, list(reversed(global_max_path))

