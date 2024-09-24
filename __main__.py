import nonbinary
from nonbinary import TreeVec
import ete3
from ete3 import Tree
import argparse
import json

def read_file(file_path):
    T1, T2 = "", "" 
    id1, id2 = [], []
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()  # 读取文件的所有行
            if len(lines) >= 4:
                T1 = lines[0].strip()
                id1 = json.loads(lines[1].strip())
                T2 = lines[2].strip() 
                id2 = json.loads(lines[3].strip())
            elif len(lines) == 2:
                T1 = lines[0].strip()  # 如果文件只有一行，只赋值 T1
            return T1, id1, T2, id2  # 返回 T1 和 T2
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None, None

if __name__ == "__main__":
    # 初始化解析器
    parser = argparse.ArgumentParser(description="Read a file from the command line")

    # 添加文件路径参数
    parser.add_argument("file", type=str, help="The path to the file to read")

    # 解析命令行参数
    args = parser.parse_args()

    # 调用读取文件的函数并将内容保存到 T1 和 T2
    T1, id1, T2, id2 = read_file(args.file)

t1 = Tree(T1)
t2 = Tree(T2)

tree1 = TreeVec(tree = t1, leaf2idx=id1)
tree2 = TreeVec(tree = t2, leaf2idx=id2)

print(tree1.simvec)
print(tree2.simvec)

print(tree1.hop_similarity(tree2))