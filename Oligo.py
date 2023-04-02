class Oligo:
    def __init__(self,name):
        self.name = name
        self.neighbours_1=[]
        self.neighbours_2=[]
        self.neighbours_3=[]

    def create_neighbours (self,neighbour,weight):
        if (weight == 1):
            self.neighbours_1.append(neighbour)
        elif (weight == 2):
            self.neighbours_2.append(neighbour)
        elif (weight == 3):
            self.neighbours_3.append(neighbour)
        

    
