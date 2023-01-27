class DoF():
    """
    class: representing a single Degree of Freedom

        self.type: type of DoF (str)
        self.disp: current displacement (float)
        self.nodeID:: node the DoF belongs to (int)
        self.isFixed: fixity status (bool)
        self.record: indicator whether to save displacement history (bool)
        self.history: displacement history (list of floats)
    """

    def __init__(self, type):
        self.type = type
        self.disp = 0
        self.nodeID = None
        self.isFixed = False
        self.record = False
        self.history = []

    def __str__(self):
        return '{} DoF of node {}, isFixed: {}, recording: {}'.format(self.type, self.nodeID, self.isFixed, self.record)

    def __repr__(self):
        pass

    def setFixity(self, fixity):
        if fixity:
            self.isFixed = True
        else:
            self.isFixed = False

    def getFixity(self):
        return self.isFixed

    def setDisp(self, disp):
        self.disp = disp
        if self.record:
            self.history.append(disp)

    def getDisp(self):
        return self.disp

    def setNode(self, nodeID):
        self.nodeID = nodeID


if __name__ == "__main__":
    dof1 = DoF('ux')
    print(dof1)
    dof1.setFixity(1)
    print(dof1)
    dof1.setFixity(0)
    print(dof1)

