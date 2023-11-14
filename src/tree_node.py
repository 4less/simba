
class TreeNode:
    def __init__(self, name):
        self.name = name
        self.children = set()

    def AddChild(self, node):
        self.children.add(node)
        